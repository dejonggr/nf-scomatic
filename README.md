# nf-scomatic

## Overview

This is a Nextflow implementation of [SComatic](https://github.com/cortes-ciriano-lab/SComatic).
Please refer to the original GitHub for more information about the pipeline.

## Usage

### Samplesheet

First, prepare a comma-delimited samplesheet with your input data. It should 
look like this:

| donor_id       | sample_id      | bam                                                                                                                 |
| -------------- | -------------- | ------------------------------------------------------------------------------------------------------------------- |
| PB_WT_AX001    | PB_WT_AX001    | /lustre/scratch126/casm/team294rr/lm26/whole_transcriptome/mapped/TRAC-2-8004-Cell1/TRAC_2_8004_scisoseq.mapped.bam |
| PB_panel_AX001 | PB_panel_AX001 | /lustre/scratch126/casm/team294rr/lm26/panel/mapped/TRAC-2-8178-Cell1/TRAC_2_8178_scisoseq.mapped.bam               |
| PB_WT_KX004    | PB_WT_KX004    | /lustre/scratch126/casm/team294rr/lm26/whole_transcriptome/mapped/TRAC-2-8005-Cell1/TRAC_2_005_scisoseq.mapped.bam  |
| PB_panel_KX004 | PB_panel_KX004 | /lustre/scratch126/casm/team294rr/lm26/panel/mapped/TRAC-2-8179-Cell1/TRAC_2_8179_scisoseq.mapped.bam               |

The samplesheet must contain the columns `sample_id`, `donor_id`, and `bam`. If
the `bam` value is a local path, instead of an iRODS path, please pass
`--location local`.

### Celltypes

You must also prepare a tab-delimited file containing the celltype annotations, 
with columns `sample_id`, `Index` and `Cell_type`. It should look like this:

| sample_id       | Index           | Cell_type |
|----------------|----------------|-----------|
| PB_panel_KX004 | GAACTAAGAACTCCCT | NK_cell   |
| PB_WT_KX004    | AGAGCTGGACCGTGAC | Monocyte  |
| PB_panel_KX004 | GCCTCTCCTGTTAGGG | B_cell    |
| PB_panel_AX001 | CAAAGCTACCAGACCC | T_cells   |
| PB_WT_AX001    | CCGTACGGATGGGTTT | NK_cell   |
| PB_panel_AX001 | CCAAATTCTATCACCA | NK_cell   |

**N.B.** Make sure the barcodes here match those in your BAM files. For example,
if the BAM files contain the `-1` suffix, include that here as well. You can 
check the barcode formatting in the BAMs by looking at the `CB:Z` tag.

### Reference files

Several reference files must be downloaded from the [SComatic](github.com/cortes-ciriano-lab/SComatic)
GitHub to be used as inputs:

- Required: Panel of normals files (`--pons`) can be found [here](https://github.com/cortes-ciriano-lab/SComatic/blob/main/PoNs).
- Required: RNA editing files (`--editing`) can be found [here](https://github.com/cortes-ciriano-lab/SComatic/tree/main/RNAediting).
- Recommended: High quality regions of the human genome to intersect (`--bed`) can be found [here](https://github.com/cortes-ciriano-lab/SComatic/blob/main/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed).

## Run

Now you can run the pipeline using:

```bash
nextflow run nf-scomatic \
  --samplesheet samplesheet.csv \
  --celltypes celltypes.tsv \
  --genome genome.fa \
  --pons data/PoN.scRNAseq.hg38.tsv \
  --editing data/AllEditingSites.hg38.txt \
  --bed data/UCSC.k100_umap.without.repeatmasker.bed \
  --out_dir out/ \
  -profile pacbio
```

To get a full list of parameters, use the `--help` flag:

```bash
nextflow run nf-scomatic --help
```

## Parameters

### Required options

The following parameters must be passed:

- `--samplesheet`
- `--celltypes`
- `--genome`
- `--pons`
- `--editing`

### Input/output options

```ansi
Script started on 2025-03-27 15:22:15+00:00 [TERM="xterm-256color" TTY="/dev/pts/1196" COLUMNS="161" LINES="31"]
ESC[33mNextflow 24.10.5 is available - Please consider updating your version to itESC(BESC[m

ESC[1mESC[38;5;232mESC[48;5;43m N E X T F L O W ESC[0;2m  ~  ESC[mversion 24.10.4ESC[m
ESC[K
LaunchingESC[35m `./main.nf` ESC[0;2m[ESC[0;1;36mchaotic_kimuraESC[0;2m] DSL2 - ESC[36mrevision: ESC[0;36mbe0fcb1af9ESC[m
ESC[K
--help                    ESC[2m[boolean, string] ESC[0mShow the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed. ESC[2mESC[0m
--helpFull                ESC[2m[boolean]         ESC[0mShow the help message for all non-hidden parameters. ESC[2mESC[0m
--showHidden              ESC[2m[boolean]         ESC[0mShow all hidden parameters in the help message. This needs to be used in combination with `--help` or `--helpFull`. ESC[2mESC[0m

ESC[4mESC[1mInput/output optionsESC[0m
  --samplesheet           ESC[2m[string]  ESC[0mComma-separated samplesheet with columns 'donor_id', 'sample_id', and 'bam'. ESC[2mESC[0m
  --celltypes             ESC[2m[string]  ESC[0mTab-separated file mapping cell barcodes to celltype information. It must contain at least the columns 'Index' and 'Cell_type'. ESC[2mESC[0m
  --genome                ESC[2m[string]  ESC[0mFasta file for genome build. ESC[2mESC[0m
  --bed                   ESC[2m[string]  ESC[0mBED of regions of interest. ESC[2mESC[0m
  --location              ESC[2m[string]  ESC[0mAre the BAMs saved locally or on iRODs? ESC[2m (accepted: irods, local) [default: irods] ESC[0m
  --modality              ESC[2m[string]  ESC[0mSingle cell RNA-seq (GEX) or single cell ATAC-seq (ATAC)? ESC[2m (accepted: GEX, ATAC) [default: GEX] ESC[0m
  --subset_bed            ESC[2m[string]  ESC[0mOptionally subset the BAMs before filtering with a BED file. Make sure the 'chr' conventions match (ie. chr1 or 1)! ESC[2mESC[0m
  --ignore_base_quality   ESC[2m[boolean] ESC[0mOptionally ignore base quality scores in the BAM. The pipeline will artificially set the quality to the maximum value for all bases, to circumvent SComatic's 
filtering. ESC[2mESC[0m 
  --mutations             ESC[2m[string]  ESC[0mA folder with prior scomatic mutations output. ESC[2mESC[0m
  --publish_celltype_bams ESC[2m[boolean] ESC[0mPublish the celltype-split BAMs to the celltype_bams/ subdirectory? ESC[2mESC[0m
  --out_dir               ESC[2m[string]  ESC[0mOutput directory. ESC[2m[default: ./] ESC[0m

ESC[4mESC[1mSComatic options - SplitBamCeltypesESC[0m
  --max_nM                ESC[2m[integer] ESC[0mMaximum number of mismatches permitted to consider reads  for analysis. By default, this filter is switched off, although we recommed using --max_nM 5. If 
applied, this filter requires having the nM tag in the bam file. ESC[2mESC[0m 
  --max_NH                ESC[2m[integer] ESC[0mMaximum number of alignment hits permitted to consider reads for analysis. By default, this filter is switched off, although we recommend using --max_NH 1. This 
filter requires having the NH tag in the bam file. ESC[2mESC[0m 
  --min_MQ                ESC[2m[integer] ESC[0mMinimum mapping quality required to consider reads for analysis. Set this value to 0 to switch this filter off. --min_MQ 255 is recommended for RNA data, and 
--min_MQ 30 for DNA data. ESC[2m[default: 255] ESC[0m 
  --n_trim                ESC[2m[integer] ESC[0mNumber of bases trimmed by setting the base quality to 0 at the beginning and end of each read. ESC[2m[default: 0] ESC[0m

ESC[4mESC[1mSComatic options - BaseCellCounterESC[0m
  --min_ac                ESC[2m[integer] ESC[0mMinimum alt count to consider a genomic site for further analysis. ESC[2m[default: 0] ESC[0m
  --min_af                ESC[2m[number]  ESC[0mMinimum alt allele fraction to consider a genomic site for further analysis. ESC[2m[default: 0] ESC[0m
  --min_dp                ESC[2m[integer] ESC[0mMinimum coverage to consider a genomic site for further analysis. ESC[2m[default: 5] ESC[0m
  --min_cc                ESC[2m[integer] ESC[0mMinimum number of cells required to consider a genomic site for further analysis. ESC[2m[default: 5] ESC[0m
  --min_bq                ESC[2m[integer] ESC[0mMinimum base quality permited for the base counts. ESC[2m[default: 20] ESC[0m
  --max_dp                ESC[2m[integer] ESC[0mMaximum number of reads per genomic site that are read by pysam pileup (to save time and memory. Set this value to 0 to switch this filter off (recommended for 
high-depth sequencing). ESC[2m[default: 8000] ESC[0m 

ESC[4mESC[1mSComatic options - BaseCellCalling.step1ESC[0m
  --max_cell_types        ESC[2m[number] ESC[0mMaximum number of celltypes carrying a mutation to make a somatic call. ESC[2m[default: 1] ESC[0m

ESC[4mESC[1mSComatic options - BaseCellCalling.step2ESC[0m
  --pons                  ESC[2m[string] ESC[0mPanel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts. ESC[2mESC[0m
  --editing               ESC[2m[string] ESC[0mRNA editing file to be used to remove RNA-editing sites. ESC[2mESC[0m

 ESC[2m!! Hiding 3 param(s), use the `--showHidden` parameter to show them !!ESC[0m
-ESC[2m----------------------------------------------------ESC[0m-


Script done on 2025-03-27 15:22:19+00:00 [COMMAND_EXIT_CODE="0"]
```

- `--samplesheet` [string]: Comma-separated samplesheet with columns 'donor_id',
'sampel_id', and 'bam'.
- `--celltypes` [string]: Tab-separated file mapping cell barcodes to celltype
information, with columns 'Index' and 'Cell_type'. 
- `--genome`  [string]: Fasta file for the genome build.
- `--location` [string]: Are the BAMs saved locally or on iRODs? [default:
local] (accepted: local, irods)
- `--out_dir` [string]: Output directory. [default: out/] 
- `--modality` [string]: Single cell RNA-seq (GEX) or single cell ATAC-seq
(ATAC)?  (accepted: GEX, ATAC) 
- `--subset_bed` [string]:  Optionally subset the BAMs before filtering with a 
BED file.
- `--ignore_base_quality` [boolean]: Optionally ignore base quality scores in
the BAM. The pipeline will artificially set the quality to the maximum value for
all bases to circumvent SComatic's filtering.
- `--bed` [string]: BED of regions of interest.
- `--mutations` [string]: A folder with prior SComatic mutations output.
- `--publish_celltype_bams` [boolean]: Publish the celltype-split BAMs to the
`celltype_bams/` subdirectory? [default: false]

### SComatic options - SplitBamCeltypes

- `--max_nM` [integer]: Maximum number of mismatches permitted to consider reads
for analysis. By default, this filter is switched off, although we recommed
using --max_nM 5. If applied, this filter requires having the nM tag in the bam
file.  
- `--max_NH` [integer]: Maximum number of alignment hits permitted to consider
reads for analysis. By default, this filter is switched off, although we
recommend using --max_NH 1. This filter requires having the NH tag in the bam
file.  
- `--min_MQ` [integer]: Minimum mapping quality required to consider reads for
analysis. Set this value to 0 to switch this filter off. `--min_MQ 255` is
recommended for RNA data, and `--min_MQ 30` for DNA data. [default: 255]  
- `--n_trim` [integer]: Number of bases trimmed (by setting the base quality to
0) at the beginning and end of each read. [default: 0]

### SComatic options - BaseCellCounter

- `--min_ac` [integer]: Minimum alt count to consider a genomic site for further
analysis. [default: 0] 
- `--min_af` [number]:  Minimum alt allele fraction to consider a genomic site
for further analysis. [default: 0] 
- `--min_dp` [integer]: Minimum coverage to consider a genomic site for further
analysis. [default: 5] 
- `--min_cc` [integer]: Minimum number of cells required to consider a genomic
site for further analysis. [default: 5] 
- `--min_bq` [integer]: Minimum base quality permited for the base counts.
[default: 20] 
- `--max_dp` [integer]: Maximum number of reads per genomic site that are read
by pysam pileup (to save time and memory). Set this value to 0 to switch this 
filter off (recommended for high-depth sequencing). [default: 8000]

### SComatic options - BaseCellCalling.step1

- `--max_cell_types` [number]: Maximum number of celltypes carrying a mutation
to make a somatic call. [default: 1] 

### SComatic options - BaseCellCalling.step2

- `--pons` [string]: Panel of normals (PoN) file to be used to remove germline
polymorphisms and recurrent artefacts. 
- `--editing` [string]: RNA editing file to be used to remove RNA-editing sites. 

## Profiles

### `pacbio`

PacBio BAMs have different sequencing quality metrics than 10X BAMs. When using
PacBio BAMs, please pass `-profile pacbio` so that the appropriate filters can
be used:

```bash
nextflow run nf-scomatic \
  --samplesheet samplesheet.csv \
  --celltypes celltypes.tsv \
  -profile pacbio
```

This sets the parameters `--ignore_base_quality` and `--min_MQ 60`. These
differences are explained below.

#### Base quality

In long-read sequencing data from PacBio, base-level quality scores are not
meaningful and the may be filled with placeholders (!) or omitted (*).
Because SComatic filters on base quality (BQ), these placeholders can lead all 
reads to be filtered out. In order to prevent this, you must set
`--ignore_base_quality` when running `nf-scomatic` on PacBio samples. This will
replace the placeholder characters in the base quality column with the maximum
value to prevent filtering.

#### Mapping quality

Additionally, the default minimum mapping quality in the SComatic pipeline is 
255. However, the minimum mapping quality for PacBio reads is capped at 60,
so in order to prevent them being filtered out you must lower the `--min_MQ` to 
60.

### `GRCh38`

If you are running this pipeline with GRCh38, you can use `-profile GRCh38` 
and the pipeline will use the appropriate local files for `--genome`, `--pons`,
`--editing`, and `--bed`:

```bash
nextflow run nf-scomatic \
  --samplesheet samplesheet.csv \
  --celltypes celltypes.tsv \
  -profile GRCh38
```

If this does not work for you (e.g. due to permission issues), see the above
section on reference files. You will have to manually download these yourself.

## Output

An example of an output directory, for the sample `TX_WT_KX004`:

```
out/
└── TX_WT_KX004
    ├── cell_callable_sites
    │   ├── TX_WT_KX004.B_cell.SitesPerCell.tsv
    │   ├── TX_WT_KX004.Monocyte.SitesPerCell.tsv
    │   ├── TX_WT_KX004.NK_cell.SitesPerCell.tsv
    │   └── TX_WT_KX004.T_cells.SitesPerCell.tsv
    ├── TX_WT_KX004.calling.step2.intersect.tsv
    ├── TX_WT_KX004.calling.step2.pass.tsv
    ├── TX_WT_KX004.calling.step2.tsv
    ├── TX_WT_KX004.coverage_cell_count.per_chromosome.report.tsv
    └── TX_WT_KX004.coverage_cell_count.report.tsv
```

