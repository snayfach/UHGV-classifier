#!/usr/bin/env python

import csv
import math
import os
import shutil
import subprocess as sp
import sys
import time
from collections import OrderedDict, defaultdict
from typing import TypedDict

import taxopy
from taxopy.exceptions import MajorityVoteError

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
    with open(inpath) as fi:
        for r in csv.DictReader(fi, delimiter="\t"):
            r[scorekey] = float(r[scorekey])
            if refs and r[refkey] not in refs:
                continue
            elif r[querykey] not in tophits:
                tophits[r[querykey]] = r
            elif r[scorekey] > tophits[r[querykey]][scorekey]:
                tophits[r[querykey]] = r
    return tophits


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
            "uhgv_taxdump/nodes.dmp",
            "uhgv_taxdump/names.dmp",
            "genome_metadata.tsv",
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
        self.genome_to_taxid = {}
        map_path = os.path.join(self.paths["dbdir"], "genome_metadata.tsv")
        with open(map_path) as fi:
            for r in csv.DictReader(fi, delimiter="\t"):
                self.genome_to_taxid[r["genome_id"]] = int(r["taxid"])
        self.ref_genomes = {
            gid: {"taxid": tid} for gid, tid in self.genome_to_taxid.items()
        }
        self.taxdb = taxopy.TaxDb(
            nodes_dmp=os.path.join(self.paths["dbdir"], "uhgv_taxdump", "nodes.dmp"),
            names_dmp=os.path.join(self.paths["dbdir"], "uhgv_taxdump", "names.dmp"),
        )
        self.ref_clusters = {}
        path = os.path.join(self.paths["dbdir"], "viral_cluster_info.tsv")
        with open(path) as fi:
            for r in csv.DictReader(fi, delimiter="\t"):
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
                    with open(inpath) as fi:
                        for line in fi:
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
                with open(path) as fi:
                    for r in csv.reader(fi, delimiter="\t"):
                        row = DiamondRow.from_csv(r)
                        genome = row.qseqid.rsplit("_", 1)[0]
                        if row.qseqid == row.sseqid:
                            score += row.bitscore
                if genome is not None:
                    data[genome]["selfscore"] = score

            with open(self.paths["selfaai"], "w") as out:
                header = ["genome_id", "genes", "selfscore"]
                out.write("\t".join(header) + "\n")
                for id in data:
                    rec = [id, data[id]["genes"], data[id]["selfscore"]]
                    out.write("\t".join([str(_) for _ in rec]) + "\n")

            shutil.rmtree(selfaln_dir)

        # update queries
        with open(self.paths["selfaai"]) as fi:
            for r in csv.DictReader(fi, delimiter="\t"):
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
                with open(os.path.join(self.paths["aaidir"], file)) as handle:
                    out.write(next(handle))
                break

            for file in os.listdir(self.paths["aaidir"]):
                with open(os.path.join(self.paths["aaidir"], file)) as handle:
                    next(handle)
                    for line in handle:
                        out.write(line)

        shutil.rmtree(self.paths["dmnddir"])
        shutil.rmtree(self.paths["aaidir"])

    def build_lineage(self, taxon, stop_rank=None):
        rank_order = ["vfam", "vsubfam", "vgenus", "vsubgen", "votu"]
        rank_to_name = dict(zip(taxon.rank_lineage, taxon.name_lineage))
        parts = []
        for rank in rank_order:
            if rank in taxon.rank_taxid_dictionary:
                parts.append(rank_to_name[rank])
            else:
                parts.append("NA")
            if stop_rank is not None and rank == stop_rank:
                break
        return ";".join(parts)

    def get_lineage_string(self, taxid, up_to_rank=None):
        taxon = taxopy.Taxon(taxid, self.taxdb)
        return self.build_lineage(taxon, stop_rank=up_to_rank)

    def get_truncated_taxon(self, taxid, score):
        cutoffs = {
            "vsubgen": 80.0,
            "vgenus": 65.0,
            "vsubfam": 32.0,
            "vfam": 5.5,
        }
        taxon = taxopy.Taxon(taxid, self.taxdb)
        for rank, _ in taxon.ranked_name_lineage:
            if rank in {"no rank", "votu"}:
                continue
            cutoff = cutoffs.get(rank)
            if cutoff is not None and score >= cutoff:
                return taxopy.Taxon(taxon.rank_taxid_dictionary[rank], self.taxdb)
        return None

    def taxon_to_lineage(self, taxon):
        return self.build_lineage(taxon, stop_rank=taxon.rank)

    def find_top_hits(self):
        if os.path.exists(self.paths["blastani"]):
            tophits = blast_to_tophits(self.paths["blastani"], refs=self.ref_genomes)
            for qname, hit in tophits.items():
                self.queries[qname]["blastani"] = hit

        if os.path.exists(self.paths["blastaai"]):
            all_hits = defaultdict(list)
            top_scores = {}
            top_rows = {}
            with open(self.paths["blastaai"]) as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    if r["reference"] not in self.ref_genomes:
                        continue
                    score = float(r["norm_score"])
                    r["norm_score"] = score
                    all_hits[r["query"]].append(r)
                    prev_score = top_scores.get(r["query"])
                    if prev_score is None or score > prev_score:
                        top_scores[r["query"]] = score
                        top_rows[r["query"]] = r
            for qname, hit in top_rows.items():
                self.queries[qname]["blastaai"] = hit
                self.queries[qname]["blastaai_all"] = all_hits[qname]

    def assign_taxonomy(self, aai_similarity_threshold=0.825):
        for id in self.queries:
            r = {}
            r["genome_id"] = id
            r["genome_length"] = self.queries[id]["length"]
            r["genome_n_genes"] = self.queries[id].get("genes")
            r["assigned_taxon"] = None
            r["assignment_method"] = None
            r["references_for_assignment"] = None
            r["top_nucleotide_hit"] = None
            r["top_nucleotide_hit_ani"] = None
            r["top_nucleotide_hit_query_af"] = None
            r["top_nucleotide_hit_target_af"] = None
            r["ani_taxonomy"] = None
            r["top_protein_hit"] = None
            r["top_protein_hit_shared_genes"] = None
            r["top_protein_hit_aai"] = None
            r["top_protein_hit_proteomic_similarity"] = None
            r["assigned_lineage"] = None

            if "blastani" in self.queries[id]:
                h = self.queries[id]["blastani"]
                r["top_nucleotide_hit"] = h["reference"]
                r["top_nucleotide_hit_ani"] = h["ani"]
                r["top_nucleotide_hit_query_af"] = h["qcov"]
                r["top_nucleotide_hit_target_af"] = h["tcov"]
                r["ani_taxonomy"] = self.get_lineage_string(
                    self.ref_genomes[h["reference"]]["taxid"]
                )

            if "blastaai" in self.queries[id]:
                h = self.queries[id]["blastaai"]
                r["top_protein_hit"] = h["reference"]
                r["top_protein_hit_shared_genes"] = int(h["hits"])
                r["top_protein_hit_aai"] = float(h["aai"])
                r["top_protein_hit_proteomic_similarity"] = h["norm_score"]

            classified_by_ani = False
            top_nucleotide_hit = r["top_nucleotide_hit"]
            if top_nucleotide_hit is not None:
                top_nucleotide_hit_ani = r["top_nucleotide_hit_ani"]
                top_nucleotide_hit_query_af = r["top_nucleotide_hit_query_af"]
                top_nucleotide_hit_target_af = r["top_nucleotide_hit_target_af"]
                assert top_nucleotide_hit_ani is not None
                assert top_nucleotide_hit_query_af is not None
                assert top_nucleotide_hit_target_af is not None
                assert r["ani_taxonomy"] is not None
                if float(top_nucleotide_hit_ani) >= 95 and (
                    float(top_nucleotide_hit_query_af) >= 85
                    or float(top_nucleotide_hit_target_af) >= 85
                ):
                    ani_taxid = self.ref_genomes[top_nucleotide_hit]["taxid"]
                    if ani_taxid != 1:
                        taxon = taxopy.Taxon(ani_taxid, self.taxdb)
                        lineage = self.taxon_to_lineage(taxon)
                        r["ani_taxonomy"] = lineage
                        r["assigned_taxon"] = taxon.name
                        r["assigned_lineage"] = lineage
                        r["references_for_assignment"] = top_nucleotide_hit
                        r["assignment_method"] = "nucleotide"
                        classified_by_ani = True

            if not classified_by_ani and "blastaai_all" in self.queries[id]:
                hits = self.queries[id]["blastaai_all"]
                top_score = self.queries[id]["blastaai"]["norm_score"]
                threshold = top_score * aai_similarity_threshold

                taxa, weights, used_refs = [], [], []
                for h in hits:
                    score = h["norm_score"]
                    if score >= threshold:
                        taxid = self.ref_genomes[h["reference"]]["taxid"]
                        truncated = self.get_truncated_taxon(taxid, score)
                        if truncated is not None and truncated.taxid != 1:
                            taxa.append(truncated)
                            weights.append(score)
                            used_refs.append(h["reference"])

                if taxa:
                    if len(taxa) == 1:
                        result = taxa[0]
                    else:
                        try:
                            result = taxopy.find_majority_vote(
                                taxa, self.taxdb, fraction=0.6, weights=weights
                            )
                        except MajorityVoteError:
                            result = None
                    if result is not None and result.taxid != 1:
                        assigned_lineage = self.taxon_to_lineage(result)
                        r["assigned_taxon"] = result.name
                        r["assignment_method"] = "protein"
                        r["assigned_lineage"] = assigned_lineage
                        r["references_for_assignment"] = ",".join(used_refs)

            if r["assigned_taxon"] is not None:
                r.update(self.ref_clusters[r["assigned_taxon"]])

            self.queries[id]["record"] = r

    @staticmethod
    def format_row(record, fields):
        vals = (record.get(f) for f in fields)
        return "\t".join("NA" if v is None else str(v) for v in vals)

    def write_results(self):
        fields = [
            "genome_id",
            "genome_length",
            "genome_n_genes",
            "assigned_taxon",
            "assigned_lineage",
            "assignment_method",
            "references_for_assignment",
            "top_nucleotide_hit",
            "top_nucleotide_hit_ani",
            "top_nucleotide_hit_query_af",
            "top_nucleotide_hit_target_af",
            "top_protein_hit",
            "top_protein_hit_shared_genes",
            "top_protein_hit_aai",
            "top_protein_hit_proteomic_similarity",
        ]
        with open(self.paths["classify_summary"], "w") as out:
            out.write("\t".join(fields) + "\n")
            for query in self.queries.values():
                out.write(self.format_row(query["record"], fields) + "\n")

        fields = [
            "genome_id",
            "assigned_taxon",
            "assigned_lineage",
            "host_lineage",
            "ictv_lineage",
            "lifestyle",
            "genome_length_median",
            "genome_length_iqr",
        ]
        with open(self.paths["taxon_info"], "w") as out:
            out.write("\t".join(fields) + "\n")
            for query in self.queries.values():
                out.write(self.format_row(query["record"], fields) + "\n")


def main(
    input,
    outdir,
    dbdir,
    sens="sensitive",
    threads=None,
    splits=None,
    continue_=False,
    quiet=False,
    keep_tmp=False,
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

    if not keep_tmp and os.path.exists(vclass.paths["tmpdir"]):
        console.log("Removing temporary directory")
        shutil.rmtree(vclass.paths["tmpdir"])

    console.log("Elapsed time (s): %s" % round(time.time() - prog_start, 2))
    console.log("Peak RAM usage (GB): %s" % round(utility.max_mem_usage(), 2))
