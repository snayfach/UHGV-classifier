# UHGV-classifier

## Overview

The UHGV-classifier assigns viral genomes to the taxonomy-like clusters of the Unified Human Gut Virome (UHGV), a comprehensive, richly annotated database of viral genomes from the human gut microbiome. UHGV organizes viral genomes into hierarchical taxa (vFAM, vSUBFAM, vGENUS, vSUBGEN, and vOTU), providing a standardized framework for classifying human gut viruses.

Given one or more viral genome sequences, UHGV-classifier will:

- Assign each genome to the most appropriate UHGV genome cluster
- Identify the closest UHGV reference genomes and estimate genomic novelty relative to the database
- Report taxonomic and biological features of the assigned lineage

![](img/classify_workflow.png)

## Installation

UHGV-classifier uses [Pixi](https://pixi.sh/) to manage its software environment and dependencies. Install Pixi first, then clone the repository and run the `--help` command to see usage instructions:

```sh
git clone https://github.com/snayfach/UHGV-classifier.git
cd UHGV-classifier

pixi run uhgv --help
```

Results will be written to the `output/` directory.

## Database setup

UHGV-classifier requires a local copy of a reference database of viral genomes from the human gut. Download it using the `download-database` subcommand:

```sh
pixi run uhgv download-database .
```

## Example workflow

Download a test dataset of 5 phages from [Nishijima _et al._](https://www.nature.com/articles/s41467-022-32832-w) using `curl`:

```sh
curl -O https://raw.githubusercontent.com/snayfach/UHGV-classifier/main/example/viral_sequences.fna
```

Classify the phage genomes in `viral_sequences.fna` using the `classify` subcommand:

```sh
pixi run uhgv classify viral_sequences.fna output uhgv-db-v1.0
```

This command will generate an `output` directory containing two output files:

- `output/classify_summary.tsv`: information related to the classification process.
- `output/taxon_info.tsv`: details about the UHGV clusters that the query genomes were assigned to.

### `classify_summary.tsv`

This file summarizes the classification process and reports the evidence supporting each assignment.

| genome_id        | genome_length | genome_n_genes | assigned_taxon | assigned_lineage                                                | assignment_method | references_for_assignment                | top_nucleotide_hit | top_nucleotide_hit_ani | top_nucleotide_hit_query_af | top_nucleotide_hit_target_af | top_protein_hit | top_protein_hit_shared_genes | top_protein_hit_aai | top_protein_hit_proteomic_similarity |
| ---------------- | ------------- | -------------- | -------------- | --------------------------------------------------------------- | ----------------- | ---------------------------------------- | ------------------ | ---------------------- | --------------------------- | ---------------------------- | --------------- | ---------------------------- | ------------------- | ------------------------------------ |
| 0008_k141_99927  | 96989         | 106            | vGENUS-00180   | vFAM-00050;vSUBFAM-00057;vGENUS-00180                           | protein           | UHGV-0030436,UHGV-0031631,UHGV-0182391,… | UHGV-0030436       | 93.65                  | 86.61                       | 83.58                        | UHGV-0030436    | 93                           | 89.27               | 82.57                                |
| 0009_k141_103019 | 93763         | 97             | vOTU-000013    | vFAM-00033;vSUBFAM-00217;vGENUS-00144;vSUBGEN-00142;vOTU-000013 | nucleotide        | UHGV-2012830                             | UHGV-2012830       | 96.44                  | 94.34                       | 86.39                        | UHGV-2012830    | 92                           | 93.97               | 91.44                                |
| 0021_k141_9110   | 96143         | 90             | vOTU-000001    | vFAM-00055;vSUBFAM-00155;vGENUS-00073;vSUBGEN-00042;vOTU-000001 | nucleotide        | UHGV-0121692                             | UHGV-0121692       | 96.1                   | 97.99                       | 93.93                        | UHGV-0121692    | 87                           | 92.1                | 92.4                                 |
| 0029_k141_61813  | 36943         | 51             | vOTU-000480    | vFAM-00035;vSUBFAM-01016;vGENUS-01771;vSUBGEN-01230;vOTU-000480 | nucleotide        | UHGV-0278246                             | UHGV-0278246       | 98.29                  | 97.42                       | 94.31                        | UHGV-0278246    | 49                           | 99.0                | 96.67                                |
| 0033_k141_114792 | 41425         | 52             | vOTU-000022    | vFAM-00001;vSUBFAM-00035;vGENUS-18848;vSUBGEN-36070;vOTU-000022 | nucleotide        | UHGV-2247228                             | UHGV-2247228       | 97.06                  | 95.75                       | 87.33                        | UHGV-2247228    | 51                           | 95.62               | 92.36                                |

- `genome_id`: query genome identifier
- `genome_length`: length of the query genome in bp
- `genome_n_genes`: number of genes in the query genome
- `assigned_taxon`: identifier of the UHGV genome cluster to which the query was assigned
- `assigned_lineage`: full lineage of the UHGV genome cluster
- `assignment_method`: whether nucleotide- or protein-based similarity was used for assignment
- `references_for_assignment`: reference genomes used for assignment
- `top_nucleotide_hit`: top hit reference based on ANI
- `top_nucleotide_hit_ani`: average nucleotide identity
- `top_nucleotide_hit_query_af`: % of query covered
- `top_nucleotide_hit_target_af`: % of reference covered
- `top_protein_hit`: top hit reference based on AAI
- `top_protein_hit_shared_genes`: number of proteins that are shared between the query and the reference
- `top_protein_hit_aai`: average amino acid identity
- `top_protein_hit_proteomic_similarity`: proteomic similarity

### `taxon_info.tsv`

This file reports metadata associated with the assigned UHGV genome clusters.

| genome_id        | assigned_taxon | assigned_lineage                                                | host_lineage                                                                          | ictv_lineage                                                                                   | lifestyle          | genome_length_median | genome_length_iqr |
| ---------------- | -------------- | --------------------------------------------------------------- | ------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- | ------------------ | -------------------- | ----------------- |
| 0008_k141_99927  | vGENUS-00180   | vFAM-00050;vSUBFAM-00057;vGENUS-00180                           | Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae (88.9%)                | Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Crassvirales;Steigviridae   (100.0%)   | lytic (100.0%)     | 97883.5              | 93684.0-100165.0  |
| 0009_k141_103019 | vOTU-000013    | vFAM-00033;vSUBFAM-00217;vGENUS-00144;vSUBGEN-00142;vOTU-000013 | Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Phocaeicola   (100.0%) | Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Crassvirales;Suoliviridae   (100.0%)   | lytic (100.0%)     | 102504               | 102504.0-102504.0 |
| 0021_k141_9110   | vOTU-000001    | vFAM-00055;vSUBFAM-00155;vGENUS-00073;vSUBGEN-00042;vOTU-000001 | Bacteria;Firmicutes (100.0%)                                                          | Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Crassvirales;Intestiviridae   (100.0%) | lytic (100.0%)     | 100403               | 100403.0-100403.0 |
| 0029_k141_61813  | vOTU-000480    | vFAM-00035;vSUBFAM-01016;vGENUS-01771;vSUBGEN-01230;vOTU-000480 | Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae (100.0%)              | Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes (100.0%)                               | temperate (100.0%) | 38407                | 38407.0-38407.0   |
| 0033_k141_114792 | vOTU-000022    | vFAM-00001;vSUBFAM-00035;vGENUS-18848;vSUBGEN-36070;vOTU-000022 | NA                                                                                    | NA                                                                                             | lytic (100.0%)     | 45698                | 45698.0-45698.0   |

- `genome_id`: user genome identifier
- `assigned_taxon`: UHGV taxon identifier
- `assigned_lineage`: UHGV taxon lineage
- `host_lineage`: Consensus GTDB host lineage
- `ictv_lineage`: Consensus ICTV taxon lineage
- `lifestyle`: Consensus virus lifestyle
- `genome_length_median`: median genome length of viruses in lineage
- `genome_length_iqr`: interquartile range of genome length

## Citation

If you use the UHGV-classifier in your research, please cite both the software and the underlying publication:

**Publication:**

> [**A genomic atlas of the human gut virome elucidates genetic factors shaping host interactions**](https://www.biorxiv.org/content/10.1101/2025.11.01.686033)
>
> Camargo, A. P., Baltoumas, F. A., Ndela, E. O., Fiamenghi, M. B., Merrill, B. D., Carter, M. M., Pinto, Y., Chakraborty, M., Andreeva, A., Ghiotto, G., Shaw, J., Proal, A. D., Sonnenburg, J. L., Bhatt, A. S., Roux, S., Pavlopoulos, G. A., Nayfach, S., & Kyrpides, N. C. — _bioRxiv_ (2025), DOI: 10.1101/2025.11.01.686033

**Software:**
> Nayfach, S. (2025). UHGV classifier (Version 1.0.0) [Software]. Zenodo. <https://doi.org/10.5281/zenodo.17418882>
