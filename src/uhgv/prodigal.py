import multiprocessing.pool
from pathlib import Path

import pyrodigal_gv

from uhgv import utility

_orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True, mask=True)


def _predict_genes(seq):
    return (seq.accession, _orf_finder.find_genes(seq.seq))


class ProdigalGv:
    def __init__(self, input_file: Path, prodigal_output: Path) -> None:
        self.input_file = input_file
        self.prodigal_output = prodigal_output

    def run_parallel_prodigal(self, threads: int) -> None:
        input_sequences = utility.read_fasta(self.input_file)
        with (
            multiprocessing.pool.Pool(threads) as pool,
            open(self.prodigal_output, "w") as fout,
        ):
            for seq_i, (seq_acc, predicted_genes) in enumerate(
                pool.imap(_predict_genes, input_sequences), 1
            ):
                for gene_i, gene in enumerate(predicted_genes, 1):
                    header = (
                        f"{seq_acc}_{gene_i} # {gene.begin} # {gene.end} # "
                        f"{gene.strand} # ID={seq_i}_{gene_i};"
                        f"partial={int(gene.partial_begin)}{int(gene.partial_end)};"
                        f"start_type={gene.start_type};rbs_motif={gene.rbs_motif};"
                        f"rbs_spacer={gene.rbs_spacer};"
                        f"genetic_code={gene.translation_table};"
                        f"gc_cont={gene.gc_cont:.3f}"
                    )
                    gene_seq = utility.Sequence(
                        header, gene.translate(include_stop=False)
                    )
                    fout.write(str(gene_seq))
