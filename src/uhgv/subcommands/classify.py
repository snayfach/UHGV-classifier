#!/usr/bin/env python

import csv
import math
import os
import shutil
import subprocess as sp
import sys
import time
from collections import OrderedDict
from typing import TypedDict

from uhgv import prodigal, utility
from uhgv.utility import DiamondRow, get_n_available_cpus


class SelfAlignmentStats(TypedDict):
    genes: int
    selfscore: float


def blast_to_tophits(
    inpath, querykey="query", refkey="reference", scorekey="norm_score", refs=None
):
    """
    takes
        inpath: tsv file with header and at least 3 fields: query, reference, score
        querykey, refkey, scorekey: names of respective fields
        keep: optional list of references to keep, excluding the rest
    returns
        tophits: dict mapping each query to highest scoring record
    """
    tophits = {}
    for r in csv.DictReader(open(inpath), delimiter="\t"):
        r[scorekey] = float(r[scorekey])
        if refs and r[refkey] not in refs:
            continue
        elif r[querykey] not in tophits:
            tophits[r[querykey]] = r
        elif r[scorekey] > tophits[r[querykey]][scorekey]:
            tophits[r[querykey]] = r
    return tophits


def assign_aai_taxonomy(taxonomy, score):
    cutoffs = {
        "species": 95.0,
        "subgenus": 80.0,
        "genus": 65.0,
        "subfamily": 32.0,
        "family": 5.5,
    }
    prefix2rank = {
        "vOTU": "species",
        "vSUBGEN": "subgenus",
        "vGENUS": "genus",
        "vSUBFAM": "subfamily",
        "vFAM": "family",
    }
    taxa = taxonomy.split(";")
    for i in range(2, len(taxa) + 1):
        lineage = taxa[: -i + 1]
        prefix = taxa[-i].split("-")[0]
        if prefix != "Unclassified":
            rank = prefix2rank[prefix]
            if score >= cutoffs[rank]:
                return ";".join(lineage)
    return None


class ViralClassifier:
    def __init__(self, input, outdir, dbdir, sens, threads, splits, continue_, quiet):
        self.input = input
        self.outdir = outdir
        self.dbdir = dbdir.rstrip("/")
        self.sens = sens
        self.threads = threads if threads is not None else get_n_available_cpus()
        self.splits = splits if splits is not None else self.threads
        self.continue_ = continue_
        self.quiet = quiet
        self.blastcpus = math.ceil(1.0 * self.threads / self.splits)
        self.set_paths()
        self.perform_checks()

    def set_paths(self):
        self.paths = {}
        self.paths["input"] = self.input
        self.paths["outdir"] = self.outdir
        self.paths["dbdir"] = self.dbdir
        self.paths["tmpdir"] = os.path.join(self.paths["outdir"], "tmp")
        self.paths["blastn"] = os.path.join(self.paths["tmpdir"], "blastn.tsv")
        self.paths["blastani"] = os.path.join(self.paths["tmpdir"], "blastani.tsv")
        self.paths["prodigal"] = os.path.join(self.paths["tmpdir"], "prodigal.faa")
        self.paths["selfaai"] = os.path.join(self.paths["tmpdir"], "selfscores.tsv")
        self.paths["diamond"] = os.path.join(self.paths["tmpdir"], "diamond.tsv")
        self.paths["blastaai"] = os.path.join(self.paths["tmpdir"], "blastaai.tsv")
        self.paths["blastdir"] = os.path.join(self.paths["tmpdir"], "blastn")
        self.paths["dmnddir"] = os.path.join(self.paths["tmpdir"], "diamond")
        self.paths["aaidir"] = os.path.join(self.paths["tmpdir"], "aai")
        self.paths["classify_summary"] = os.path.join(
            self.paths["outdir"], "classify_summary.tsv"
        )
        self.paths["taxon_info"] = os.path.join(self.paths["outdir"], "taxon_info.tsv")

    def perform_checks(self):

        # check executables
        utility.check_executables(["diamond", "blastn"])

        # check database files
        if not os.path.exists(self.dbdir):
            sys.exit("\nError: database directory not found: '%s'" % self.dbdir)
        files = [
            "genomes.fna",
            "proteins.faa",
            "proteins.phr",
            "genome_taxonomy.tsv",
        ]
        for file in files:
            if not os.path.exists(os.path.join(self.paths["dbdir"], file)):
                sys.exit("\nError: database file not found: '%s'" % file)

        # check input seqs
        if not os.path.exists(self.input):
            sys.exit("\nError: input file does not exist")

        # check output directory does not exist
        # unless using --continue flag
        if os.path.exists(self.paths["tmpdir"]):
            if not self.continue_:
                sys.exit(
                    "\nError: output directory already exists. Remove directory or use --continue"
                )
        else:
            os.makedirs(self.paths["tmpdir"])

    def load_queries(self):
        self.queries = OrderedDict()
        for r in utility.read_fasta(self.paths["input"]):
            self.queries[r.accession] = {}
            self.queries[r.accession]["length"] = len(r)

    def load_refdb(self):
        self.ref_genomes = {}
        path = os.path.join(self.paths["dbdir"], "genome_taxonomy.tsv")
        for r in csv.DictReader(open(path), delimiter="\t"):
            r["taxonomy"] = ";".join(
                [
                    r["family_vc"],
                    r["subfamily_vc"],
                    r["genus_vc"],
                    r["subgenus_vc"],
                    r["species_vc"],
                ]
            )
            r["taxonomy"] = r["taxonomy"].replace("NULL", "Unclassified")
            self.ref_genomes[r["genome_id"]] = r
        self.ref_clusters = {}
        path = os.path.join(self.paths["dbdir"], "viral_cluster_info.tsv")
        for r in csv.DictReader(open(path), delimiter="\t"):
            self.ref_clusters[r["taxon_id"]] = r
            del self.ref_clusters[r["taxon_id"]]["taxon_id"]

    def blastani(self):
        if os.path.exists(self.paths["blastani"]):
            return

        utility.split_fasta(
            self.paths["input"], self.paths["blastdir"], self.splits, ".fna"
        )

        commands = []
        for file in os.listdir(self.paths["blastdir"]):
            if file.endswith(".fna"):
                inpath = os.path.join(self.paths["blastdir"], file)
                outpath = os.path.join(self.paths["blastdir"], f"{file}.tsv")
                logpath = os.path.join(self.paths["blastdir"], f"{file}.log")
                cmd = [
                    "blastn",
                    "-query",
                    inpath,
                    "-db",
                    f"{self.paths['dbdir']}/genomes",
                    "-out",
                    outpath,
                    "-num_threads",
                    str(self.blastcpus),
                    "-outfmt",
                    "6 std qlen slen",
                    "-max_target_seqs",
                    "1000",
                ]
                commands.append([cmd, logpath])

        return_codes = utility.parallel(utility.run_shell, commands, self.splits)
        if sum(return_codes) != 0:
            msg = "\nError: One or more blastn tasks failed to run\n"
            msg += f"See logs for details: {self.paths['blastdir']}/*.log"
            sys.exit(msg)

        with open(self.paths["blastn"], "w") as out:
            for file in os.listdir(self.paths["blastdir"]):
                if file.endswith(".tsv"):
                    inpath = os.path.join(self.paths["blastdir"], file)
                    for line in open(inpath):
                        out.write(line)

        utility.ani_calculator(self.paths["blastn"], self.paths["blastani"])

        shutil.rmtree(self.paths["blastdir"])

    def call_genes(self):
        if os.path.exists(self.paths["prodigal"]):
            return
        prodigal_obj = prodigal.ProdigalGv(self.paths["input"], self.paths["prodigal"])
        prodigal_obj.run_parallel_prodigal(threads=self.threads)

    def self_protein_alignment(self):

        if not os.path.exists(self.paths["selfaai"]):
            data: dict[str, SelfAlignmentStats] = {}

            # make dir
            selfaln_dir = os.path.join(self.paths["tmpdir"], "selfaln")
            if not os.path.exists(selfaln_dir):
                os.makedirs(selfaln_dir)

            # split fasta
            handle, qname, filenum = None, None, 0
            try:
                for r in utility.read_fasta(self.paths["prodigal"]):
                    genome = r.accession.rsplit("_", 1)[0]
                    if genome not in data:
                        data[genome] = {"genes": 0, "selfscore": 0.0}
                    data[genome]["genes"] += 1
                    if qname is None or genome != qname:
                        if handle is not None:
                            handle.close()
                        filenum += 1
                        qname = genome
                        handle = open(
                            os.path.join(selfaln_dir, str(filenum)) + ".faa", "w"
                        )
                    if handle is None:
                        sys.exit(
                            "\nError: no protein sequences found in Prodigal output"
                        )
                    handle.write(">" + r.accession + "\n" + str(r.seq) + "\n")
            finally:
                if handle is not None:
                    handle.close()
            if filenum == 0:
                sys.exit("\nError: no protein sequences found in Prodigal output")

            # run diamond
            commands = []
            files = [_ for _ in os.listdir(selfaln_dir) if _.endswith(".faa")]
            for file in files:
                path = os.path.join(selfaln_dir, file)
                logpath = f"{path}.log"
                cmd = [
                    "diamond",
                    "blastp",
                    "--masking",
                    "none",
                    "-k",
                    "1000",
                    "-e",
                    "1e-3",
                    f"--{self.sens}",
                    "--query",
                    path,
                    "--db",
                    path,
                    "--out",
                    f"{path}.tsv",
                    "--threads",
                    "1",
                ]
                commands.append([cmd, logpath])

            return_codes = utility.parallel(utility.run_shell, commands, self.splits)
            if sum(return_codes) != 0:
                msg = "\nError: One or more diamond tasks failed to run\n"
                msg += f"See logs for details: {selfaln_dir}/*.log"
                sys.exit(msg)

            # calculate aai
            files = [_ for _ in os.listdir(selfaln_dir) if _.endswith(".tsv")]
            for file in files:
                path = os.path.join(selfaln_dir, file)
                score = 0.0
                genome = None
                with open(path) as fin:
                    for r in csv.reader(fin, delimiter="\t"):
                        row = DiamondRow.from_csv(r)
                        genome = row.qseqid.rsplit("_", 1)[0]
                        if row.qseqid == row.sseqid:
                            score += row.bitscore
                if genome is not None:
                    data[genome]["selfscore"] = score

            # write results
            with open(self.paths["selfaai"], "w") as out:
                header = ["genome_id", "genes", "selfscore"]
                out.write("\t".join(header) + "\n")
                for id in data:
                    rec = [id, data[id]["genes"], data[id]["selfscore"]]
                    out.write("\t".join([str(_) for _ in rec]) + "\n")

            # cleanup
            shutil.rmtree(selfaln_dir)

        # update queries
        for r in csv.DictReader(open(self.paths["selfaai"]), delimiter="\t"):
            self.queries[r["genome_id"]].update(r)

    def db_protein_alignment(self):
        if os.path.exists(self.paths["diamond"]):
            return
        logpath = f"{self.paths['diamond']}.log"
        cmd = [
            "diamond",
            "blastp",
            "-k",
            "1000",
            "-e",
            "1e-3",
            "--masking",
            "none",
            f"--{self.sens}",
            "--query",
            self.paths["prodigal"],
            "--out",
            self.paths["diamond"],
            "--threads",
            str(self.threads),
            "--db",
            f"{self.paths['dbdir']}/proteins",
            "--outfmt",
            "6",
        ]
        with open(logpath, "w") as log:
            result = sp.run(cmd, stdout=log, stderr=sp.STDOUT)
        if result.returncode != 0:
            msg = "\nError: DIAMOND database search failed to run\n"
            msg += f"See log for details: {logpath}"
            sys.exit(msg)

    def blastaai(self):
        if os.path.exists(self.paths["blastaai"]):
            return

        utility.split_dmnd(self.paths["diamond"], self.paths["dmnddir"], self.threads)

        if not os.path.exists(self.paths["aaidir"]):
            os.makedirs(self.paths["aaidir"])

        argument_list = []
        for file in os.listdir(self.paths["dmnddir"]):
            inpath = os.path.join(self.paths["dmnddir"], file)
            outpath = os.path.join(self.paths["aaidir"], file)
            argument_list.append([inpath, outpath, self.paths["selfaai"]])
        utility.parallel(utility.aai_main, argument_list, threads=self.threads)

        with open(self.paths["blastaai"], "w") as out:
            for file in os.listdir(self.paths["aaidir"]):
                handle = open(os.path.join(self.paths["aaidir"], file))
                out.write(next(handle))
                break

            for file in os.listdir(self.paths["aaidir"]):
                handle = open(os.path.join(self.paths["aaidir"], file))
                next(handle)
                for line in handle:
                    out.write(line)
                handle.close()

        shutil.rmtree(self.paths["dmnddir"])
        shutil.rmtree(self.paths["aaidir"])

    def find_top_hits(self):
        for type in ["blastani", "blastaai"]:
            tophits = blast_to_tophits(self.paths[type], refs=self.ref_genomes)
            for qname, hit in tophits.items():
                self.queries[qname][type] = hit

    def assign_taxonomy(self):
        for id in self.queries:
            r = {}
            r["genome_id"] = id
            r["genome_length"] = self.queries[id]["length"]
            r["genome_num_genes"] = self.queries[id].get("genes")
            r["taxon_id"] = None
            r["taxon_lineage"] = None
            r["class_method"] = None
            r["class_rank"] = None
            r["ani_reference"] = None
            r["ani"] = None
            r["ani_query_af"] = None
            r["ani_target_af"] = None
            r["ani_taxonomy"] = None
            r["aai_reference"] = None
            r["shared_genes"] = None
            r["aai"] = None
            r["proteomic_similarity"] = None
            r["aai_taxonomy"] = None

            if "blastani" in self.queries[id]:
                h = self.queries[id]["blastani"]
                r["ani_reference"] = h["reference"]
                r["ani"] = h["ani"]
                r["ani_query_af"] = h["qcov"]
                r["ani_target_af"] = h["tcov"]
                r["ani_taxonomy"] = self.ref_genomes[h["reference"]]["taxonomy"]

            if "blastaai" in self.queries[id]:
                h = self.queries[id]["blastaai"]
                r["aai_reference"] = h["reference"]
                r["shared_genes"] = int(h["hits"])
                r["aai"] = float(h["aai"])
                r["proteomic_similarity"] = float(h["norm_score"])
                r["aai_taxonomy"] = self.ref_genomes[h["reference"]]["taxonomy"]

            classified_by_ani = False
            ani_reference = r["ani_reference"]
            if ani_reference is not None:
                ani = r["ani"]
                ani_query_af = r["ani_query_af"]
                ani_target_af = r["ani_target_af"]
                ani_taxonomy = r["ani_taxonomy"]
                assert ani is not None
                assert ani_query_af is not None
                assert ani_target_af is not None
                assert ani_taxonomy is not None
                if float(ani) >= 95 and (
                    float(ani_query_af) >= 85 or float(ani_target_af) >= 85
                ):
                    while ani_taxonomy.endswith(";Unclassified"):
                        ani_taxonomy = ani_taxonomy.rsplit(";Unclassified", 1)[0]
                    r["ani_taxonomy"] = ani_taxonomy
                    r["taxon_lineage"] = ani_taxonomy
                    r["taxon_id"] = ani_taxonomy.split(";")[-1]
                    r["class_method"] = "nucleotide"
                    r["class_rank"] = "species"
                    classified_by_ani = True

            if not classified_by_ani and r["aai_reference"] is not None:
                r["taxon_lineage"] = assign_aai_taxonomy(
                    r["aai_taxonomy"], r["proteomic_similarity"]
                )
                if r["taxon_lineage"]:
                    r["class_method"] = "protein"
                    rank_dict = {
                        "vFAM": "family",
                        "vSUBFAM": "subfamily",
                        "vGENUS": "genus",
                        "vSUBGEN": "subgenus",
                    }
                    r["taxon_id"] = r["taxon_lineage"].split(";")[-1]
                    r["class_rank"] = rank_dict[r["taxon_id"].split("-")[0]]

            if r["taxon_id"] is not None:
                r.update(self.ref_clusters[r["taxon_id"]])

            self.queries[id]["record"] = r

    @staticmethod
    def _format_row(record, fields):
        vals = (record.get(f) for f in fields)
        return "\t".join("NA" if v is None else str(v) for v in vals)

    def write_results(self):
        fields = [
            "genome_id",
            "genome_length",
            "genome_num_genes",
            "taxon_id",
            "class_method",
            "class_rank",
            "ani_reference",
            "ani",
            "ani_query_af",
            "ani_target_af",
            "ani_taxonomy",
            "aai_reference",
            "shared_genes",
            "aai",
            "proteomic_similarity",
            "aai_taxonomy",
        ]
        with open(self.paths["classify_summary"], "w") as out:
            out.write("\t".join(fields) + "\n")
            for query in self.queries.values():
                out.write(self._format_row(query["record"], fields) + "\n")

        fields = [
            "genome_id",
            "taxon_id",
            "taxon_lineage",
            "host_lineage",
            "ictv_lineage",
            "lifestyle",
            "genome_length_median",
            "genome_length_iqr",
        ]
        with open(self.paths["taxon_info"], "w") as out:
            out.write("\t".join(fields) + "\n")
            for query in self.queries.values():
                out.write(self._format_row(query["record"], fields) + "\n")


def main(
    input,
    outdir,
    dbdir,
    sens="sensitive",
    threads=None,
    splits=None,
    continue_=False,
    quiet=False,
    cleanup=False,
):

    prog_start = time.time()
    vclass = ViralClassifier(
        input, outdir, dbdir, sens, threads, splits, continue_, quiet
    )

    console = utility.ConsoleLogger(quiet)

    console.log("Reading input sequences")
    vclass.load_queries()

    console.log("Reading database sequences")
    vclass.load_refdb()

    with console.status("Computing nucleotide similarity with BLASTN…"):
        vclass.blastani()

    with console.status("Predicting genes with pyrodigal-gv…"):
        vclass.call_genes()

    with console.status("Computing protein self-alignments…"):
        vclass.self_protein_alignment()

    with console.status("Searching protein database with DIAMOND…"):
        vclass.db_protein_alignment()

    with console.status("Computing AAI scores…"):
        vclass.blastaai()

    console.log("Finding top database hits")
    vclass.find_top_hits()

    console.log("Assigning taxonomy")
    vclass.assign_taxonomy()

    console.log("Writing output files")
    vclass.write_results()

    if cleanup and os.path.exists(vclass.paths["tmpdir"]):
        console.log("Removing temporary directory")
        shutil.rmtree(vclass.paths["tmpdir"])

    console.log("Elapsed time (s): %s" % round(time.time() - prog_start, 2))
    console.log("Peak RAM usage (GB): %s" % round(utility.max_mem_usage(), 2))
