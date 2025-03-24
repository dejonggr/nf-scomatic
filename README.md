# nf-scomatic

## Overview

This is a Nextflow implementation of [SComatic](github.com/cortes-ciriano-lab/SComatic).
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
`--location local`. The `sample_id` values cannot contain dashes.

### Celltypes

You must also prepare a tab-delimited file containing the celltype annotations, 
with columns `Index` and `Cell_type`.

```
Index   Cell_type
PB_panel_KX004_GAACTAAGAACTCCCT NK_cell
PB_WT_KX004_AGAGCTGGACCGTGAC    Monocyte
PB_panel_KX004_GCCTCTCCTGTTAGGG B_cell
PB_panel_AX001_CAAAGCTACCAGACCC T_cells
PB_WT_AX001_CCGTACGGATGGGTTT    NK_cell
PB_panel_AX001_CCAAATTCTATCACCA NK_cell
```

The `Index` column must be in the format `{sample_id}_{barcode}` for PacBio BAMs
or `{sample_id}_{barcode}-1` for 10X BAMs. 

### Reference files

Several reference files must be downloaded from the [SComatic](github.com/cortes-ciriano-lab/SComatic)
GitHub to be used as inputs:

- Panel of normals files (`--pons`) can be found [here](https://github.com/cortes-ciriano-lab/SComatic/blob/main/PoNs).
- RNA editing files (`--editing`) can be found [here](https://github.com/cortes-ciriano-lab/SComatic/tree/main/RNAediting).

Additionally, it is recommended to only look in high-quality regions of the
human genome (masking repeats, etc...). This can be set with the `--bed`
parameter. A repeat-masked BED containg high quality regions of the genome
can be found at `/lustre/scratch125/casm/team268im/at31/nextflow/nf-scomatic/data/GRCh38/UCSC.k100_umap.without.repeatmasker.bed`.

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