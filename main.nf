nextflow.enable.dsl=2

// TODO: implement extra params in scomatic commands

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${sample_id}"
  maxForks 10
  label "normal4core"
  input:
    tuple val(sample_id), val(irods), val(bam), val(donor_id)
  output:
    tuple val(sample_id), val(donor_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
  script:
    """
    iget -K ${irods}/${bam} ${sample_id}.bam
    if [[ `ils ${irods} | grep "${bam}.bai" | wc -l` == 1 ]]
    then
        iget -K ${irods}/${bam}.bai ${sample_id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${sample_id}.bam
    fi
    """
}

// The equivalent of an irods download, but for a local copy of samplesheet
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${sample_id}"
  label "normal4core"
  input:
    tuple val(sample_id), path(local), val(bam), val(donor_id)
  output:
    tuple val(sample_id), val(donor_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
  script:
    """
    ln -s ${local}/${bam} ${sample_id}.bam
    if [ -f "${local}/${bam}.bai" ]
    then
        ln -s ${local}/${bam}.bai ${sample_id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${sample_id}.bam
    fi
    """
}

// Optionally subset BAMs based on a BED file of coords
process subset_bams {
  tag "${sample_id}"
  label 'normal4core'
  input:
    tuple val(sample_id), val(donor_id), path(bam), path(bai), path(celltypes)
  output:
    tuple val(sample_id), val(donor_id), path("subset/*.bam"), path("subset/*.bam.bai"), path(celltypes)
  script:
    """
    mkdir subset/
    samtools view -b -h -L ${params.subset_bed} ${bam} > subset/${sample_id}.bam
    (cd subset ; samtools index -@ ${task.cpus} ${sample_id}.bam)
    """
}

// Optionally artificially set base qualities to maximum in the BAM
// This should be used when running on PacBio data, where base quality
// scores are set to `*`, which cannot be handled by various SComatic scripts.
process ignore_base_quality {
  tag "${sample_id}"
  label 'normal4core'
  input:
    tuple val(sample_id), val(donor_id), path(bam), path(bai), path(celltypes)
  output:
    tuple val(sample_id), val(donor_id),
          path("ignore_bq/${sample_id}.bam"),
          path("ignore_bq/${sample_id}.bam.bai"),
          path(celltypes)
  script:
    """
    mkdir ignore_bq/
    samtools view -h ${bam} |
    awk 'BEGIN{OFS="\\t"} /^@/ {print; next} { \$11 = ""; for (i=1; i<=length(\$10); i++) \$11 = \$11 "~"; print }' \
    > ignore_bq/${sample_id}.bam
    (cd ignore_bq ; samtools index -@ ${task.cpus} ${sample_id}.bam)
    """
}

// Extract the cell type definitions for the specified sample
// Returns the sample as part of the tuple for easy joining with sample BAM/BAI lists
// Finds the sample ID at the start of the line as per SAMPLE_BARCODE-1 specification
process getSampleCelltypes {
  tag "${sample_id}"
  label "normal"
  input:
    tuple val(sample_id), path(celltypes)
  output:
    tuple val(sample_id), path("${sample_id}-celltypes.tsv")
  script:
    """
    head -n 1 ${celltypes} > ${sample_id}-celltypes.tsv
    grep "^${sample_id}" ${celltypes} >> ${sample_id}-celltypes.tsv
    """
}

// Step one of scomatic - split the BAM into per cell type BAMs
// However, in the interest of efficiency, we do this per sample
// So we have to merge them on a donor level across samples later
// GEX/ATAC alternate parameterisation (NH/nM vs MQ) prepared in modargs
// The output for this is a little weird, as the process makes a "list of lists"
// Meanwhile I need to create a nice clean tuple for subsequent alternate grouping
// As such, this emits BAM/BAI "list of lists", which have all the necessary info
// (sample.donor.celltype) embedded in the file names
// The typical scomatic convention is sample.celltype, so I have to rename files
process splitBams {
  tag "${sample_id}"
  label "long"
  input:
    tuple val(sample_id), val(donor_id), path(bam), path(bai), path(celltypes)
  output:
    path("output/${sample_id}.${donor_id}*.bam"), emit: "bam"
    path("output/${sample_id}.${donor_id}*.bam.bai"), emit: "bai"
  script:
    def nm_arg = params.max_nM == null ? "" : "--max_nM ${params.max_nM}"
    def nh_arg = params.max_NH == null ? "" : "--max_NH ${params.max_NH}"
    """
    mkdir -p output
    python3 ${params.scomatic}/SplitBam/SplitBamCellTypes.py \\
      --bam ${bam} \\
      --meta ${celltypes} \\
      --id ${sample_id} \\
      --min_MQ ${params.min_MQ} \\
      --n_trim ${params.n_trim} \\
      --outdir output $nm_arg $nh_arg
    for file in output/${sample_id}* ; do
      new_file=\$(echo \$file | sed "s/${sample_id}/${sample_id}.${donor_id}/g")
      mv \$file \$new_file
    done
    """
}

// Merge the per-sample BAMs for each cell type into one per-donor file
// Pass the donor and cell type along for a while for downstream process use
process mergeCelltypeBams {
  tag "${donor_id}"
  label "normal4core"
  input:
    tuple val(donor_id), val(celltype), path(bams), path(bais)
  output:
    tuple val(donor_id), val(celltype), path("${donor_id}.${celltype}.bam")
  script:
    """
    samtools merge -@ {task.cpus} ${donor_id}.${celltype}.bam ${bams}
    """
}

// Index the donor-celltype BAMs, returning both the BAM and BAI
process indexCelltypeBams {
  tag "${donor_id}"
  label "normal4core"
  publishDir "${params.out_dir}/${donor_id}/celltype_bams/", 
    mode: "copy",
    enabled: params.publish_celltype_bams
  input:
    tuple val(donor_id), val(celltype), path(bam)
  output:
    tuple val(donor_id), val(celltype), path(bam), path("*.bam.bai")
  script:
    """
    samtools index -@ {task.cpus} ${bam}
    """
}

// Step two of scomatic - convert the donor-celltype BAMs to scomatic's TSVs
// Up the bin parameter to decrease the number of temporary files made
// Also the LSF config for this uses scratch to not pollute the drive with those files
// (Even if they are just temporary and get removed at the end of the process)
// There can sometimes be no output file, so specify it as optional
// Can't pass out donor info because tuples and optional don't get along
process bamToTsv {
  tag "${donor_id}"
  cpus 16
  label "long16core"
  input:
    tuple val(donor_id), val(celltype), path(bam), path(bai), path(fasta), path(fai)
  output:
    path("output/*.tsv", optional: true)
  script:
    """
    mkdir -p temp
    mkdir -p output
    python3 ${params.scomatic}/BaseCellCounter/BaseCellCounter.py --bam ${bam} \
      --ref ${fasta} \
      --chrom all \
      --out_folder output \
      --min_bq ${params.min_bq} \
      --min_mq ${params.min_MQ} \
      --min_cc ${params.min_cc} \
      --min_dp ${params.min_dp} \
      --tmp_dir temp \
      --bin 1000000 \
      --nprocs ${task.cpus}
    rm -rf temp
    """
}

// Step three of scomatic - merge the celltype TSVs into a single mega TSV
// Important to stage the input in a directory here so scomatic sees it like it wants to
process mergeTsvs {
  tag "${donor_id}"
  label "week"
  input:
    tuple val(donor_id), path(tsvs, stageAs: "input/*")
  output:
    tuple val(donor_id), path("${donor_id}.BaseCellCounts.AllCellTypes.tsv")
  script:
    """
    python3 ${params.scomatic}/MergeCounts/MergeBaseCellCounts.py --tsv_folder input \
        --outfile ${donor_id}.BaseCellCounts.AllCellTypes.tsv
    """
}

// Step 4.1 of scomatic - perform the initial mutation calling
// Set max_cell_types to a stratospheric value to disable that filter
// As we are interested in germline mutations spanning multiple populations
process callMutations {
  tag "${donor_id}"
  label "week"
  input:
    tuple val(donor_id), path(tsv), path(fasta), path(fai)
  output:
    tuple val(donor_id), path("${donor_id}.calling.step1.tsv")
  script:
    """
    python3 ${params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py \
        --infile ${tsv} \
        --outfile ${donor_id} \
        --ref ${fasta} \
        --max_cell_types ${params.max_cell_types}
    """
}

// Step 4.2 of scomatic - perform initial mutation filtering
// Need to do GEX and ATAC separately due to the editing file
// This needs a bunch of memory apparently for reasons that I can't identify from the code
process filterMutationsGex {
  tag "${donor_id}"
  label "long10gb"
  publishDir "${params.out_dir}/${donor_id}", mode:"copy"
  input:
    tuple val(donor_id), path(tsv), path(pons), path(editing)
  output:
    tuple val(donor_id), path("${donor_id}.calling.step2.tsv")
  script:
    """
    python3 ${params.scomatic}/BaseCellCalling/BaseCellCalling.step2.py \
        --infile ${tsv} \
        --outfile ${donor_id} \
        --editing ${editing} \
        --pon ${pons}
    """
}

process filterMutationsAtac {
  tag "${donor_id}"
  label "long10gb"
  publishDir "${params.out_dir}/${donor_id}", mode:"copy"
  input:
    tuple val(donor_id), path(tsv), path(pons)
  output:
    tuple val(donor_id), path("${donor_id}.calling.step2.tsv")
  script:
    """
    python3 ${params.scomatic}/BaseCellCalling/BaseCellCalling.step2.py \
        --infile ${tsv} \
        --outfile ${donor_id} \
        --pon ${pons}
    """
}

// Intersect the filtered mutations with a BED region of interest
process intersectBed {
  tag "${donor_id}"
  label "normal"
  publishDir "${params.out_dir}/${donor_id}", mode:"copy"
  input:
    tuple val(donor_id), path(tsv), path(bed)
  output:
    tuple val(donor_id), path("${donor_id}.calling.step2.intersect.tsv")
  script:
    """
    bedtools intersect -header -a ${tsv} -b ${bed} > ${donor_id}.calling.step2.intersect.tsv
    """
}

// Final mutation filter, just keep the PASS ones
// Standard $ escaping applies
process passMutations {
  tag "${donor_id}"
  label "normal"
  publishDir "${params.out_dir}/${donor_id}", mode:"copy"
  input:
    tuple val(donor_id), path(tsv)
  output:
    tuple val(donor_id), path("${donor_id}.calling.step2.pass.tsv")
  script:
    """
    awk '\$1 ~ /^#/ || \$6 == "PASS"' ${tsv} > ${donor_id}.calling.step2.pass.tsv
    """
}

// Get callable sites on a per cell type level
// This makes a couple files, both end in .report.tsv
process callableSitesCellType {
  tag "${donor_id}"
  label "long"
  publishDir "${params.out_dir}/${donor_id}", mode:"copy"
  input:
    tuple val(donor_id), path(tsv)
  output:
    tuple val(donor_id), path("*.report.tsv")
  script:
    """
    python3 ${params.scomatic}/GetCallableSites/GetAllCallableSites.py \
        --infile ${tsv} \
        --outfile ${donor_id}
    """
}

// Get callable sites on a cell level
// Mirror argument values based on earlier bamToTsv step for consistency
process callableSitesCell {
  tag "${donor_id}"
  label "week16core10gb"
  publishDir "${params.out_dir}/${donor_id}/cell_callable_sites", mode:"copy"
  input:
    tuple val(donor_id), val(celltype), path(bam), path(bai), path(fasta), path(fai), path(tsv)
  output:
    tuple val(donor_id), val(celltype), path("*.SitesPerCell.tsv")
  script:
    """
    mkdir -p temp
    python3 ${params.scomatic}/SitesPerCell/SitesPerCell.py --bam ${bam} \
        --infile ${tsv} \
        --ref ${fasta} \
        --min_bq ${params.min_bq} \
        --min_mq ${params.min_MQ} \
        --tmp_dir temp \
        --bin 1000000 \
        --nprocs ${task.cpus}
    rm -rf temp
    """
}

// Post-processing - convert the donor-celltype BAMs to per cell genotypes
// Up the bin parameter to decrease the number of temporary files made
// Also the LSF config for this uses scratch to not pollute the drive with those files
// (Even if they are just temporary and get removed at the end of the process)
process bamToGenotype {
  tag "${donor_id}"
  label "long16core10gb"
  publishDir "${params.out_dir}/${donor_id}-genotypes", mode:"copy"
  input:
    tuple val(donor_id), val(celltype), path(bam), path(bai), path(fasta), path(fai), path(allcelltypes), path(mutations)
  output:
    path("${donor_id}.${celltype}.single_cell_genotype.tsv", optional: true)
  script:
    """
    mkdir -p temp
    python3 ${params.scomatic}/SingleCellGenotype/SingleCellGenotype.py --bam ${bam} \
        --ref ${fasta} \
        --infile ${mutations}/${donor_id}/${donor_id}.calling.step2.intersect.tsv \
        --meta ${allcelltypes} \
        --outfile ${donor_id}.${celltype}.single_cell_genotype.tsv \
        --min_bq ${params.min_bq} \
        --min_mq ${params.min_MQ} \
        --tmp_dir temp \
        --bin 1000000 \
        --nprocs ${task.cpus}
    rm -rf temp
    """
}

// Perform the various data downloading and demultiplexing early in scomatic
// Yields demultiplexed cell type BAMs on a donor level
// Along with the master cell type file and the reference genome + index
// And the upcoming downstream process --min_mq (255 GEX, 30 ATAC)
// As both downstream uses want that to call the process appropriately
// The cell type BAMs are a good input for both normal scomatic
// And the optional secondary pipeline of single cell genotype calling
workflow STEP1 {
  main:

    // Make sure all the input files actually exist
    if (params.samplesheet == null) {
      error "Please provide a samplesheet CSV file via --samplesheet"
    }
    if (params.celltypes == null) {
      error "Please provide a cell types TSV file via --celltypes"
    }

    // Load the mapping as a CSV, skipping the one line of header for ease of use
    // Split the BAM file into its parent directory and the actual file name
    // As the parent directory is used both by irods and local for finding the BAI
    // It's fine to use file() here for the file name operations
    // Nextflow doesn't natively check the paths exist, allowing for irods splitting
    samplesheet = Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header: true)
      .map({row -> [row.sample_id, 
                    file(row.bam).getParent(), 
                    file(row.bam).getName(), 
                    row.donor_id]})
    allCelltypes = Channel.fromPath(params.celltypes, checkIfExists: true)
    fasta = Channel.fromPath(params.genome, checkIfExists: true)
    fai = Channel.fromPath(params.genome+".fai", checkIfExists: true)
    
    // Check the modality
    if (!(params.modality in ["GEX","ATAC"])) {
      error "Unknown modality, must be GEX or ATAC"
    }
    
    // Sanitise the samplesheet to just samples that show up in the cell types file
    // To do this, we need a list of sample IDs from the cell type file
    // Need to use splitCsv as regular splitText doesn't have a skip option
    // Skip the header line, then keep only the first column for the file
    // And split that string on _ and remove the last element from it
    // Then glue it back together with _ as the delimiter
    // The equivalent of "_".join(VARIABLE.split("_")[:-1]) in python
    celltypeSamples = allCelltypes
      .splitCsv(sep:"\t", skip:1)
      .map({it -> it[0].split("_")[0..-2].join("_")})
      .unique()
    
    // Doing a join here keeps the intersection of the sample IDs
    // That show up in both the samplesheet info and our parsed celltype sample list
    samplesheet = samplesheet.join(celltypeSamples)
    if (params.location == "irods") {
      // Download all the samplesheet from irods - yields a sample/donor/BAM/BAI tuple
      sampleBams = irods(samplesheet)
    }
    else if (params.location == "local") {
      // The samplesheet are available locally, make the same tuple
      // Symlink up the files so they're named like the downloads would have been
      sampleBams = local(samplesheet)
    }
    else {
      error "Unknown location, must be irods or local"
    }
    
    // Extract the samples as the first element of the tuple
    // Then combine it with the cell types file
    samples = sampleBams
      .map({it -> it[0]})
      .combine(allCelltypes)
    
    // Can now search the cell types file for each of the sample IDs
    sampleCelltypes = getSampleCelltypes(samples)
    
    // Join the cell type file with the BAM/BAI list from earlier
    sampleBamsCelltypes = sampleBams.join(sampleCelltypes)
    
    // optional pre-processing step: subset bams to regions of interest
    if (params.subset_bed == null) {
      bams_to_split_1 = sampleBamsCelltypes
    } else {
      bams_to_split_1 = subset_bams(sampleBamsCelltypes)
    }

    // optional pre-processing step: set base quality to maximum for all reads
    // to handle PacBio data where base quality scores are set to `*`
    if (params.ignore_base_quality) {
      bams_to_split_2 = ignore_base_quality(bams_to_split_1)
    } else {
      bams_to_split_2 = bams_to_split_1
    }
    
    // Step one of scomatic - split the BAMs to cell types on a per sample basis
    sampleFiles = splitBams(bams_to_split_2)
    // This outputs a "list of lists" of BAMs/BAIs, constructed for each sample
    // Flatten the files into a single list, and then index them as 
    // sample.donor.celltype
    // (the file names are sample.donor.celltype.bam[.bai])
    // Note that the dot needs to be escaped for the parsing to work, here and later
    bam = sampleFiles.bam
      .flatten()
      .map({file -> [file.getName().split('\\.bam')[0], file]})
    bai = sampleFiles.bai
      .flatten()
      .map({file -> [file.getName().split('\\.bam')[0], file]})
    
    // We can now safely join the two file lists on the index to combine correctly
    // Replace the index with two entries, one for donor, one for cell type
    // And then group up all the files from a given donor and cell type for merging
    sampleCelltypeBams = bam
      .join(bai)
      .map({name, file, ind -> [name.split('\\.')[1], name.split('\\.')[2], file, ind]})
      .groupTuple(by: [0,1])
    unindexedCelltypeBams = mergeCelltypeBams(sampleCelltypeBams)
    
    // Index the cell type BAMs, and then add the genome to the resulting tuple
    // All downstream uses of the BAMs want the genome present, so store that too
    indexedCelltypeBams = indexCelltypeBams(unindexedCelltypeBams)
      .combine(fasta)
      .combine(fai)
  
  emit:
    indexedCelltypeBams = indexedCelltypeBams
    allCelltypes = allCelltypes
    fasta = fasta
    fai = fai
}

// Turn an existing set of scomatic mutations, plus freshly re-generated cell type BAMs
// Into single cell level mutation calls
workflow genotypes {
    // A folder with prior scomatic mutations output
    // Needs DONOR-MODALITY subfolders present
    if (params.mutations == null) {
        error "Please provide a path to prior scomatic pipeline mutation output via --mutations"
    }
    mutations = Channel.fromPath(params.mutations, checkIfExists: true)
    // Run the first step of scomatic, getting the cell type BAMs
    STEP1()
    // Pull out the various outputs that we use
    indexedCelltypeBams = STEP1.out.indexedCelltypeBams
    allCelltypes = STEP1.out.allCelltypes
    // We have the donor cell type genotypes we're after here
    // Just need to include the master cell types file and the master mutations folder
    indexedCelltypeBams = indexedCelltypeBams
        .combine(allCelltypes)
        .combine(mutations)
    // Get the single cell genotypes
    cellGenotypes = bamToGenotype(indexedCelltypeBams)
}

// The main scomatic workflow
workflow {

  // Validate the parameters
  validateParameters()

  // Print a summary of supplied parameters
  log.info paramsSummaryLog(workflow)

  // Some inputs that do not matter to step one, and as such are absent there
  // But we care about them here, so load them up and check them
  pons = Channel.fromPath(params.pons, checkIfExists: true)
  // This is a GEX-only file
  if (params.modality == "GEX") {
    editing = Channel.fromPath(params.editing, checkIfExists: true)
  }
  bed = Channel.fromPath(params.bed, checkIfExists: true)
  // Run the first step of scomatic, getting the cell type BAMs
  STEP1()
  // Pull out the various outputs that we use
  indexedCelltypeBams = STEP1.out.indexedCelltypeBams
  fasta = STEP1.out.fasta
  fai = STEP1.out.fai
  // Step two of scomatic - convert the BAMs to scomatic's TSV format for each cell type
  celltypeTsvs = bamToTsv(indexedCelltypeBams)
  // This only stores output where the thing actually ran, don't have to filter on that
  // Have to rebuild donor info from the file name as tuples and optional don't work
  // Then group them up on donor for collapsing to a per-donor TSV
  donorCelltypeTsvs = celltypeTsvs
    .map({file -> [file.getName().split('\\.')[0], file]})
    .groupTuple()
  // Step three of scomatic - merge the cell type TSVs into a single one per donor
  // We need to also pass the genome to the next scomatic step, so combine it in
  masterTsvs = mergeTsvs(donorCelltypeTsvs)
    .combine(fasta)
    .combine(fai)
  // Step 4.1 of scomatic - perform the initial mutation calling
  unfilteredMutationsClean = callMutations(masterTsvs)
  // Combine in the PoNs for filtering, but keep a clean version for later
  unfilteredMutations = unfilteredMutationsClean.combine(pons)
  // Step 4.2 of scomatic - filter the called mutations based on editing/PoNs
  // This requires two separate processes for the modalities as ATAC has no editing
  // Regardless of process, combine in the BED for later
  if (params.modality == "GEX") {
    filteredMutations = filterMutationsGex(unfilteredMutations.combine(editing))
        .combine(bed)
  }
  else if (params.modality == "ATAC") {
    filteredMutations = filterMutationsAtac(unfilteredMutations)
        .combine(bed)
  }
  // Final filtering steps - intersect with a BED of interest
  intersectedMutations = intersectBed(filteredMutations)
  // Stringently keep just the PASS mutations
  passedMutations = passMutations(intersectedMutations)
  // Compute the number of callable sites
  // Start off with per cell type, can just feed it the step 4.1 output
  cellTypeSites = callableSitesCellType(unfilteredMutationsClean)
  // The cell level callable sites require the cell type BAMs from earlier
  // To which we can "non-unique join" the step 4.1 output
  // The "non-unique join" is accomplished via .combine(by: 0)
  cellSitesInput = indexedCelltypeBams.combine(unfilteredMutationsClean, by: 0)
  // And now we can get cell level callable sites
  cellSites = callableSitesCell(cellSitesInput)
}
