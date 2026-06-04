import csv
import multiprocessing as mp
import os
import platform
import resource
import signal
import subprocess as sp
import sys
import textwrap
import time
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from enum import Enum, auto
from pathlib import Path

if sys.version_info >= (3, 14):
    from compression import bz2, gzip, lzma, zstd
else:
    import bz2
    import gzip
    import lzma
import psutil
from rich.console import Console


def check_executables(requirements):
    fails = 0
    for program in requirements:
        found = False
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path.strip('"'), program)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                found = True
                break
        if not found:
            msg = f"\nError: required program '{program}' not executable or not found on $PATH"
            sys.stderr.write(msg)
            fails += 1
    if fails > 0:
        sys.exit("")


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    zstd = auto()
    uncompressed = auto()


def is_compressed(filepath: Path) -> Compression:
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        elif tuple(signature[:4]) == (0x28, 0xB5, 0x2F, 0xFD):
            return Compression.zstd
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath):
    filepath_compression = is_compressed(filepath)
    if filepath_compression is Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression is Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression is Compression.xz:
        fin = lzma.open(filepath, "rt")
    elif filepath_compression is Compression.zstd and sys.version_info >= (3, 14):
        fin = zstd.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    try:
        yield fin
    finally:
        fin.close()


class Sequence:
    def __init__(self, header: str, seq: str, compress: bool = False) -> None:
        self._compress = compress
        self._header = header
        self._seq = seq.encode("ascii")

    @property
    def header(self) -> str:
        return self._header

    @property
    def accession(self) -> str:
        return self._header.split()[0]

    @property
    def seq(self) -> str:
        return self._seq.decode()

    def __str__(self) -> str:
        return (
            f">{self.header}\n{textwrap.fill(self.seq, 60, break_on_hyphens=False)}\n"
        )

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, k: int):
        return Sequence(self.header, self.seq[k], self._compress)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.seq.casefold() == other.seq.casefold()
        if isinstance(other, str):
            return self.seq.casefold() == other.casefold()
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.seq.casefold())


def read_fasta(filepath, uppercase=False, strip_n=False, compress=False):
    with open_file(filepath) as fin:
        last = None
        while True:
            if not last:
                for line in fin:
                    if line[0] == ">":
                        last = line.removesuffix("\n")
                        break
            if not last:
                break
            name, seqs, last = last[1:], [], None
            for line in fin:
                if line[0] == ">":
                    last = line.removesuffix("\n")
                    break
                seqs.append(line.removesuffix("\n"))
            seqs = "".join(seqs)
            if uppercase:
                seqs = seqs.upper()
            if strip_n:
                seqs = seqs.strip("nN")
            if len(seqs):
                yield Sequence(name, seqs, compress)
            if not last:
                break


class ConsoleLogger:
    def __init__(self, quiet=False):
        self.console = Console(quiet=quiet, highlight=False)

    def log(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.console.print(f"[dim][{timestamp}][/dim] {message}")

    def error(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.console.print(f"[dim][{timestamp}][/dim] [red]{message}[/red]")

    def status(self, message):
        return self.console.status(message)


def max_mem_usage():
    """Return max mem usage (GB) of self and child processes"""
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        return (max_mem_self + max_mem_child) / float(1e6)
    else:
        return (max_mem_self + max_mem_child) / float(1e9)


def get_n_available_cpus():
    sched_getaffinity = getattr(os, "sched_getaffinity", None)
    if sched_getaffinity is not None:
        try:
            n_cpus = len(sched_getaffinity(0))
        except OSError:
            n_cpus = os.cpu_count() or 1
    else:
        n_cpus = os.cpu_count() or 1
    return n_cpus


def mean(values):
    return sum(values) / len(values)


def split_dmnd(inpath, outdir, num_splits, ext=""):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    total_size = os.stat(inpath).st_size
    split_size = int(total_size / num_splits)

    last_id = None
    with open(inpath) as infile:
        for line in infile:
            last_id = line.split("\t", 1)[0].rsplit("_", 1)[0]

    split_num = 1
    cursize = 0
    outfile = open(os.path.join(outdir, str(split_num)) + ext, "w")
    with open(inpath) as infile:
        for line in infile:
            cur_id = line.split("\t", 1)[0].rsplit("_", 1)[0]
            if cursize > split_size and cur_id != last_id:
                split_num += 1
                cursize = 0
                outfile.close()
                outfile = open(os.path.join(outdir, str(split_num)) + ext, "w")
            outfile.write(line)
            cursize += len(line)
            last_id = cur_id
    outfile.close()


def split_fasta(inpath, outdir, num_splits, ext):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    total_size = 0
    for r in read_fasta(inpath):
        total_size += len(r)

    split_size = int(total_size / num_splits)

    split_num = 1
    cursize = 0
    out = open(os.path.join(outdir, str(split_num)) + ext, "w")

    for r in read_fasta(inpath):
        if cursize > split_size:
            split_num += 1
            cursize = 0
            out = open(os.path.join(outdir, str(split_num)) + ext, "w")
        out.write(str(r))
        cursize += len(r)

    out.close()


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def terminate_tree(pid, including_parent=True):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.terminate()
    if including_parent:
        parent.terminate()


def run_shell(cmd, log_path=None):
    if log_path:
        with open(log_path, "w") as log:
            result = sp.run(cmd, stdout=sp.DEVNULL, stderr=log)
    else:
        result = sp.run(cmd)
    return result.returncode


def parallel(function, arguments_list, threads):
    pool = mp.Pool(threads, init_worker)
    try:
        results = []
        for arguments in arguments_list:
            p = pool.apply_async(function, args=arguments)
            results.append(p)
        pool.close()
        while True:
            if all(r.ready() for r in results):
                return [r.get() for r in results]
            time.sleep(1)
    except KeyboardInterrupt:
        pid = os.getpid()
        terminate_tree(pid)


@dataclass
class BlastRow:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int

    @property
    def qcoords(self):
        return sorted([self.qstart, self.qend])

    @property
    def tcoords(self):
        return sorted([self.sstart, self.send])

    @classmethod
    def from_csv(cls, row):
        return cls(
            qseqid=row[0],
            sseqid=row[1],
            pident=float(row[2]),
            length=int(row[3]),
            mismatch=int(row[4]),
            gapopen=int(row[5]),
            qstart=int(row[6]),
            qend=int(row[7]),
            sstart=int(row[8]),
            send=int(row[9]),
            evalue=float(row[10]),
            bitscore=float(row[11]),
            qlen=int(row[12]),
            slen=int(row[13]),
        )


@dataclass
class DiamondRow:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float

    @classmethod
    def from_csv(cls, row):
        return cls(
            qseqid=row[0],
            sseqid=row[1],
            pident=float(row[2]),
            length=int(row[3]),
            mismatch=int(row[4]),
            gapopen=int(row[5]),
            qstart=int(row[6]),
            qend=int(row[7]),
            sstart=int(row[8]),
            send=int(row[9]),
            evalue=float(row[10]),
            bitscore=float(row[11]),
        )


def parse_blast(handle):
    for r in csv.reader(handle, delimiter="\t"):
        yield BlastRow.from_csv(r)


def yield_alignment_blocks(handle):
    records = parse_blast(handle)
    first_aln = next(records, None)
    if first_aln is None:
        return

    key = (first_aln.qseqid, first_aln.sseqid)
    alns = [first_aln]

    for aln in records:
        if (aln.qseqid, aln.sseqid) == key:
            alns.append(aln)
        else:
            yield alns
            key = (aln.qseqid, aln.sseqid)
            alns = [aln]
    yield alns


def prune_alns(alns, min_length=0, min_evalue=1e-3):
    keep = []
    cur_aln = 0
    qry_len = alns[0].qlen
    for aln in alns:
        qcoords = aln.qcoords
        aln_len = max(qcoords) - min(qcoords) + 1
        if aln_len < min_length or aln.evalue > min_evalue:
            continue
        if cur_aln >= qry_len or aln_len + cur_aln >= 1.10 * qry_len:
            break
        keep.append(aln)
        cur_aln += aln_len
    return keep


def compute_cov(alns):

    coords = sorted([a.qcoords for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        else:
            nr_coords.append([start, stop])

    alen = sum([stop - start + 1 for start, stop in nr_coords])
    qcov = round(100.0 * alen / alns[0].qlen, 2)

    coords = sorted([a.tcoords for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        else:
            nr_coords.append([start, stop])

    alen = sum([stop - start + 1 for start, stop in nr_coords])
    tcov = round(100.0 * alen / alns[0].slen, 2)

    return qcov, tcov


def ani_calculator(inpath, outpath):
    with open(outpath, "w") as out:
        fields = ["query", "reference", "ani", "qcov", "tcov", "norm_score"]
        out.write("\t".join(fields) + "\n")
        for alns in yield_alignment_blocks(open(inpath)):
            if alns is not None:
                alns = prune_alns(alns)
                if len(alns) > 0:
                    qname, tname = alns[0].qseqid, alns[0].sseqid
                    ani = round(
                        sum(a.length * a.pident for a in alns)
                        / sum(a.length for a in alns),
                        2,
                    )
                    qcov, tcov = compute_cov(alns)
                    norm_score = float(qcov) * ani / 100
                    row = [qname, tname, ani, qcov, tcov, norm_score]
                    out.write("\t".join([str(_) for _ in row]) + "\n")


def yield_diamond_hits(diamond):
    with open(diamond) as f:
        reader = csv.reader(f, delimiter="\t")
        first_row = next(reader, None)
        if first_row is None:
            return
        hits = [DiamondRow.from_csv(first_row)]
        for r in reader:
            row = DiamondRow.from_csv(r)
            query = row.qseqid.rsplit("_", 1)[0]
            last = hits[-1].qseqid.rsplit("_", 1)[0]
            if query != last:
                yield last, hits
                hits = []
            hits.append(row)
        if len(hits) > 0:
            last = hits[-1].qseqid.rsplit("_", 1)[0]
            yield last, hits


def split_hits(hits):
    target_to_hits = defaultdict(list)
    for hit in hits:
        tname = hit.sseqid.rsplit("_", 1)[0]
        target_to_hits[tname].append(hit)
    return target_to_hits


def best_blast_hits(hits):
    bhits = {}
    for hit in hits:
        if hit.qseqid not in bhits:
            bhits[hit.qseqid] = hit
        elif hit.bitscore > bhits[hit.qseqid].bitscore:
            bhits[hit.qseqid] = hit
    return list(bhits.values())


def aai_main(inpath, outpath, selfpath):
    selfaai = {}
    for r in csv.DictReader(open(selfpath), delimiter="\t"):
        selfaai[r["genome_id"]] = float(r["selfscore"])
    with open(outpath, "w") as out:
        header = ["query", "reference", "hits", "aai", "raw_score", "norm_score"]
        out.write("\t".join(header) + "\n")
        for qname, hits in yield_diamond_hits(inpath):
            target_to_hits = split_hits(hits)
            for tname, thits in target_to_hits.items():
                bhits = best_blast_hits(thits)
                aai = mean([_.pident for _ in bhits])
                score = sum([_.bitscore for _ in bhits])
                norm = 100 * score / selfaai[qname]
                row = [
                    qname,
                    tname,
                    len(bhits),
                    round(aai, 2),
                    round(score, 2),
                    round(norm, 2),
                ]
                out.write("\t".join([str(_) for _ in row]) + "\n")
