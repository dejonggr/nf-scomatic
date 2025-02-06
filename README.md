## nf-scomatic

This is a Nextflow implementation of SComatic. All of the SComatic parameters can be passed to `nf-scomatic`. 

### Parameters

```
$ nextflow run nf-scomatic --help
N E X T F L O W   ~  version 24.10.0
Launching nf-scomatic/main.nf [friendly_goldwasser] DSL2 - revision: d744a3d5e6

--help                    [boolean, string] Show the help message for all top level parameters. When a parameter is given to --help, the full help message of that parameter will be printed.
--helpFull                [boolean]         Show the help message for all non-hidden parameters.
--showHidden              [boolean]         Show all hidden parameters in the help message. This needs to be used in combination with --help or --helpFull.

Input/output options
  --samplesheet           [string]  Comma-separated samplesheet with columns 'donor_id', 'sample_id', and 'bam'.
  --location              [string]  Location of the input files.  (accepted: irods, local) [default: irods]
  --modality              [string]  Single cell RNA-seq (GEX) or single cell ATAC-seq (ATAC)?  (accepted: GEX, ATAC) [default: GEX]
  --celltypes             [string]  Tab-separated file mapping cell barcodes to celltype information. It must contain at least the columns 'Index' and 'Cell_type'.
  --genome                [string]  Fasta file for genome build.
  --subset_bed            [string]  Optionally subset the BAMs before filtering with a BED file. Make sure the 'chr' conventions match (ie. chr1 or 1)!
  --ignore_base_quality   [boolean] Optionally ignore base quality scores in the BAM. The pipeline will artificially set the quality to the maximum value for all bases, to circumvent SComatic's filtering.
  --bed                   [string]  BED of regions of interest.
  --mutations             [string]  A folder with prior scomatic mutations output.
  --publish_celltype_bams [boolean] Publish the celltype-split BAMs to the celltype_bams/ subdirectory?
  --out_dir               [string]  Output directory. [default: ./]SComatic options - SplitBamCeltypes
  --max_nM                [integer] Maximum number of mismatches permitted to consider reads  for analysis. By default, this filter is switched off, although we recommed using --max_nM 5. If applied, this filter requires having the nM tag in the bam file.
  --max_NH                [integer] Maximum number of alignment hits permitted to consider reads for analysis. By default, this filter is switched off, although we recommend using --max_NH 1. This filter requires having the NH tag in the bam file.
  --min_MQ                [integer] Minimum mapping quality required to consider reads for analysis. Set this value to 0 to switch this filter off. --min_MQ 255 is recommended for RNA data, and --min_MQ 30 for DNA data. [default: 255]
  --n_trim                [integer] Number of bases trimmed by setting the base quality to 0 at the beginning and end of each read. [default: 0]

SComatic options - BaseCellCounter
  --min_ac                [integer] Minimum alt count to consider a genomic site for further analysis. [default: 0]
  --min_af                [number]  Minimum alt allele fraction to consider a genomic site for further analysis. [default: 0]
  --min_dp                [integer] Minimum coverage to consider a genomic site for further analysis. [default: 5]
  --min_cc                [integer] Minimum number of cells required to consider a genomic site for further analysis. [default: 5]
  --min_bq                [integer] Minimum base quality permited for the base counts. [default: 20]
  --max_dp                [integer] Maximum number of reads per genomic site that are read by pysam pileup (to save time and memory. Set this value to 0 to switch this filter off (recommended for high-depth sequencing). [default: 8000]

SComatic options - BaseCellCalling.step1
  --max_cell_types        [number] Maximum number of celltypes carrying a mutation to make a somatic call. [default: 1]

SComatic options - BaseCellCalling.step2
  --pons                  [string] Panel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts.
  --editing               [string] RNA editing file to be used to remove RNA-editing sites.


```

### Input

The main input files myst be passed to the `--samplesheet` command. The samplesheet must be a comma-separated file (*.csv) containing the columns `sample_id`, `donor_id`, and `bam`. The `bam` should specify full paths.

| donor_id       | sample_id      | bam                                                                                                                 |
| -------------- | -------------- | ------------------------------------------------------------------------------------------------------------------- |
| PB_WT_AX001    | PB_WT_AX001    | /lustre/scratch126/casm/team294rr/lm26/whole_transcriptome/mapped/TRAC-2-8004-Cell1/TRAC_2_8004_scisoseq.mapped.bam |
| PB_panel_AX001 | PB_panel_AX001 | /lustre/scratch126/casm/team294rr/lm26/panel/mapped/TRAC-2-8178-Cell1/TRAC_2_8178_scisoseq.mapped.bam               |
| PB_WT_KX004    | PB_WT_KX004    | /lustre/scratch126/casm/team294rr/lm26/whole_transcriptome/mapped/TRAC-2-8005-Cell1/TRAC_2_005_scisoseq.mapped.bam  |
| PB_panel_KX004 | PB_panel_KX004 | /lustre/scratch126/casm/team294rr/lm26/panel/mapped/TRAC-2-8179-Cell1/TRAC_2_8179_scisoseq.mapped.bam               |
