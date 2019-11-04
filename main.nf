#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnassembly
========================================================================================
 nf-core/rnassembly Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnassembly
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/rnassembly --bams '*bam' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:
      --genome                      Name of iGenomes reference

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --genome                      Name of iGenomes reference
      --hisat2_index                Path to HiSAT2 index
      --fasta                       Path to genome fasta file
      --gtf                         Path to GTF file
      --tx                          Transcriptome FASTA file
      --ucsc_gtf                    GTF files from UCSC alignment (custom made)
      --proteins                    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
      --saveReference               Save the generated reference files to the results directory
      --skipAlignment               Skip alignment altogether
      --saveAlignedIntermediates
      --skipTrimming
      --skipSpades
      --tx_assembled

    Strandedness:
      --forwardStranded             The library is forward stranded
      --reverseStranded             The library is reverse stranded
      --unStranded                  The default behaviour

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --maxMultiqcEmailFileSize     Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

//Check incompatible parameters
if (params.skipSpades && !params.tx_assembled){
    exit 1, "You need to provide a FASTA files with assembled transcriptome if you want to skip Spades."
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false

Channel.fromPath(params.fasta, checkIfExists: true)
    .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
    .into { ch_fasta_for_hisat_index; fasta }


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/*
 * Create a channel for input read files
 */

if(params.readPaths){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_trimming }
} else {
    Channel
        .fromFilePairs( params.reads, size: 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { ch_trimming }
}

Channel
    .fromPath(params.gtf, checkIfExists: true)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .into { gtf_stringtieFPKM; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_unify }

Channel
    .fromPath(params.tx, checkIfExists: true)
    .ifEmpty { exit 1, "Tx fasta file not found: ${params.tx}" }
    .into { tx_fasta }

Channel
    .fromPath(params.proteins, checkIfExists: true)
    .ifEmpty { exit 1, "Proteins fasta file not found: ${params.proteins}" }
    .into { prot_fasta; prot_fasta_ann; prot_fasta_qc }

Channel
    .fromPath(params.ucsc_gtf, checkIfExists: true)
    .ifEmpty { exit 1, "UCSC GTF file not found: ${params.ucsc}" }
    .into {ucsc2unify}

Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
       .into{ch_where_hisat2; ch_where_hisat2_sort}

if ( params.hisat2_index && !params.skipAlignment ){

    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*", checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }

}

forwardStranded = params.forwardStranded
reverseStranded = params.reverseStranded
unStranded = params.unStranded

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Fasta Ref']        = params.fasta
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-rnassembly-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rnassembly Workflow Summary'
    section_href: 'https://github.com/nf-core/rnassembly'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


if (!params.skipTrimming){

    /*
    * STEP N - Atropos
    */
    process trimming {
        label 'low_memory'
        

        input: 
        set val(name), file(reads) from ch_trimming

        output:
        set val(name), file("*fq.gz") into ch_fastq1, ch_fastq2

        script:
        if ( task.cpus > 1 ){
            threads = "--threads ${task.cpus}"
        } else 
        {
            threads = ""
        }
        """
        atropos $threads -a AGATCGGAAGAGC  -a AAAAAAAAA \\
                -A AGATCGGAAGAGC -A TTTTTTTTTTTT \\
                --minimum-length 50 \\
                -q 20 \\
                --trim-n \\
                -o $name.1.fq.gz -p $name.2.fq.gz \\
                -pe1 ${reads[0]} -pe2 ${reads[1]}

        """
    }

}else{
   ch_trimming
       .into {ch_fastq1; ch_fastq2}
}


ch_fastq1
    .map {sample -> sample[1][0]}
    .collect()
    .set { ch_fastq_r1}

ch_fastq2
    .map {sample -> sample[1][1]}
    .collect()
    .set { ch_fastq_r2}

process merge {
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input: 
    file(reads1) from ch_fastq_r1.collect()
    file(reads2) from ch_fastq_r2.collect()

    output:
    file "mreads*fq.gz" into ch_fq_merged, ch_fq_merged2bam
 
    script:
    """
    merge(){ 
        I=\$1 
        
        shift 
        if [ \$# -eq 1 ] 
        then  
            ln -s \$@ \$I  
        else  
            zcat \$@ | gzip > \$I  
        fi  
    }  
    merge mreads.1.fq.gz $reads1
    merge mreads.2.fq.gz $reads2
    """
}

/*
 * STEP N - SPAdes
 */
if (!params.skipSpades){
    process spades {
        label 'high_memory'
        publishDir "${params.outdir}", mode: 'copy'

        input:
        file(reads) from ch_fq_merged

        output:
        file "spades/*" into spades_output
        file "spades/transcripts.fasta" into spades2blastn_tx, spades2blastn2parser, spades2qc

        script:
        """
        spades.py --pe1-1 ${reads[0]} \\
                --pe1-2 ${reads[1]} \\
                --rna \\
                -t ${task.cpus} \\
                -o spades
        """
    }
} else {
    Channel
        .fromPath(params.tx_assembled, checkIfExists: true)
        .ifEmpty { exit 1, "tx assembled file not found: ${params.gtf}" }
        .into { spades2blastn_tx; spades2blastn2parser; spades2qc }

}


// # TODO: create process to get transripts from GTF  + genome to get all possible known RNA

/*
 * STEP N - Blastp-known
 */
process blastnt_known {
    label 'mid_memory'
    publishDir "${params.outdir}/blastn_known", mode: 'copy'
    
    input:
    file(cds) from spades2blastn_tx
    file(tx) from tx_fasta

    output:
    file "blastn.known.outfmt6" into blastn_known, blastn2qc

    script:
    """
    makeblastdb -in $tx -input_type fasta -dbtype nucl -parse_seqids -out known_tx
    blastn -query $cds  \\
    -db known_tx  -max_target_seqs 1 \\
    -outfmt 6 -evalue 1e-5 -num_threads ${task.cpus} > blastn.known.outfmt6
    """
}

/*
 * STEP N - Blastp-known
 */
process blastnt_parse {
    label 'mid_memory'
    publishDir "${params.outdir}/blastn_known", mode: 'copy'
    
    input:
    file(cds) from spades2blastn2parser
    file(outfmt) from blastn_known

    output:
    file "transcript_unknown.fa" into tx_unknown1, tx_unknown2

    script:
    """
    clean_blastnt.py $outfmt > transcript_unknown
    cat $cds | seqkit grep -v -f transcript_unknown > transcript_unknown.fa
    """
}

/*
 * STEP N - TransDecoder.longORF
 */
process transdec_longorf {
    label 'low_memory'
    publishDir "${params.outdir}/transdecoder", mode: 'copy'

    input:
    file(tx) from tx_unknown1

    output:
    file("transcripts.fasta.transdecoder_dir/*") into longorf_dir
    file "transcripts.fasta.transdecoder_dir/longest_orfs.pep" into longestorf_tx, tdlong2qc

    script:
    """
    ln -s $tx transcripts.fasta
    TransDecoder.LongOrfs -t transcripts.fasta
    """
}

/*
 * STEP N - Blastp
 */
process blastp {
    label 'mid_memory'
    publishDir "${params.outdir}/blastp", mode: 'copy'

    input:
    file(prot) from prot_fasta
    file(tx) from longestorf_tx

    output:
    file "blastp.outfmt6" into blastp, blastp2qc
    script:
    """
    makeblastdb -in $prot -input_type fasta -dbtype prot -parse_seqids -out known_prot
    blastp -query $tx  \
    -db known_prot  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads ${task.cpus} > blastp.outfmt6
    """
}

/*
 * STEP N - hmmscan
 */
// process hmmscan {
//     publishDir "${params.outdir}/hmmer", mode: 'copy'
    
//     input:
//     file(tx) from longestorf_tx


//     script:
//     """
//     # cpus 8
//     hmmscan --cpu ${task.cpus} --domtblout pfam.domtblout /path/to/Pfam-A.hmm $tx
//     """
// }

/*
* STEP N -  transdecode Predict
*/
process transdec_predict {
    label 'low_memory'
    publishDir "${params.outdir}/transdecoder", mode: 'copy'

    input:
    file(tx) from tx_unknown2
    file("predict/*") from longorf_dir.collect()
    file(fmt) from blastp
    
    output:
    file("*.transdecoder.*")
    file "transcript_unknown.fa.transdecoder.bed" into transdecoder_bed, tdpredict2qc
    file("transcript_with_function.fa") into transdecoder_prot
    script:
    """
    # --retain_pfam_hits pfam
    TransDecoder.Predict -t $tx  --retain_blastp_hits $fmt -O predict --no_refine_starts
    select_blastp.py transcript_unknown.fa.transdecoder.bed > matched.txt
    cat $tx | seqkit grep -f matched.txt > transcript_with_function.fa
    """
}

/*
 * STEP N - PASA
 */
process pasa {
    publishDir "${params.outdir}/PASA", mode: 'copy'
    label 'mid_memory'

    input:
    file(tx) from transdecoder_prot
    file(genome) from fasta

    output:
    file("*.gff3") into pasa_gff3

    script:
    prefix = tx.toString() - ~/(\.fa)?(\.fasta)?(\.gz)?$/
    """
    export PASAHOME=\$(dirname \$(which python))/../opt/pasa-2.3.3
    \$PASAHOME/scripts/run_spliced_aligners.pl --aligners blat \\
        --genome $genome \\
        --transcripts $tx -I 5000000 -N 1 --CPU ${task.cpus}
    mv blat.spliced_alignments.gff3 ${prefix}.gff3
    """
}

/*
 * STEP N - ANN PASA
 */
process annotation {
    publishDir "${params.outdir}/PASA", mode: 'copy'
    label 'low_memory'

    input:
    file(bed) from transdecoder_bed
    file(pasa) from pasa_gff3
    file(protein) from prot_fasta_ann

    output:
    file "*.gtf" into pasa2unify, pasa2qc

    script:
    prefix = pasa.toString() - ~/(\.fa)?(\.fasta)?(\.gz)?$/
    """
    gfftools.py $pasa $bed $protein > ${prefix}.gff3
    gffread --keep-genes -E -T --keep-exon-attrs -F ${prefix}.gff3 > ${prefix}.gtf
    """
}

/*
* PREPROCESSING - Build HISAT2 splice sites file
*/
if(!params.skipAlignment){
    process makeHisatSplicesites {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                    saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeHisatSplicesites

        output:
        file "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites, alignment_splicesites

        script:
        """
        hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
        """
    }
}

/*
* PREPROCESSING - Build HISAT2 index
*/
if(!params.skipAlignment && !params.hisat2_index && params.fasta){
    process makeHISATindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                    saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_hisat_index
        file indexing_splicesites from indexing_splicesites
        file gtf from gtf_makeHISATindex

        output:
        file "${fasta.baseName}.*.ht2*" into hs2_indices

        script:
        if( !task.memory ){
            log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
            avail_mem = 0
        } else {
            log.info "[HISAT2 index build] Available memory: ${task.memory}"
            avail_mem = task.memory.toGiga()
        }
        if( avail_mem > params.hisat_build_memory ){
            log.info "[HISAT2 index build] Over ${params.hisat_build_memory} GB available, so using splice sites and exons in HISAT2 index"
            extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
            ss = "--ss $indexing_splicesites"
            exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
        } else {
            log.info "[HISAT2 index build] Less than ${params.hisat_build_memory} GB available, so NOT using splice sites and exons in HISAT2 index."
            log.info "[HISAT2 index build] Use --hisat_build_memory [small number] to skip this check."
            extract_exons = ''
            ss = ''
            exon = ''
        }
        """
        $extract_exons
        hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
        """
    }
}

/*
* STEP N - align with HISAT2
*/
if(!params.skipAlignment){
    star_log = Channel.from(false)
    process hisat2Align {
        label 'high_memory'
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        file(reads) from ch_fq_merged2bam
        file hs2_indices from hs2_indices.collect()
        file alignment_splicesites from alignment_splicesites.collect()
        file wherearemyfiles from ch_where_hisat2.collect()

        output:
        file "${prefix}.bam" into hisat2_bam
        file "${prefix}.hisat2_summary.txt" into alignment_logs
        file "where_are_my_files.txt"
        file "unmapped.hisat2*" optional true

        script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2l?/
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        def rnastrandness = ''
        if (forwardStranded && !unStranded){
            rnastrandness = '--rna-strandness FR'
        } else if (reverseStranded && !unStranded){
            rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
        }
        

        unaligned = "--un-conc-gz unmapped.hisat2.gz" 
        """
        hisat2 -x $index_base \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                $rnastrandness \\
                --known-splicesite-infile $alignment_splicesites \\
                --no-mixed \\
                --no-discordant \\
                -p ${task.cpus} $unaligned\\
                --met-stderr \\
                --new-summary \\
                --dta \\
                --summary-file ${prefix}.hisat2_summary.txt \\
                | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
        """
    }

    process hisat2_sortOutput {
        label 'mid_memory'
        tag "${hisat2_bam.baseName}"
        publishDir "${params.outdir}/HISAT2", mode: 'copy'

        input:
        file hisat2_bam

        output:
        file "${hisat2_bam.baseName}.sorted.bam" into bam_stringtieFPKM
        file "${hisat2_bam.baseName}.sorted.bam.bai" into bam_index

        script:
        def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
        def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
        """
        samtools sort \\
            $hisat2_bam \\
            -@ ${task.cpus} ${avail_mem} \\
            -o ${hisat2_bam.baseName}.sorted.bam
        samtools index ${hisat2_bam.baseName}.sorted.bam
        """
    }

    /*
    * STEP N - stringtie FPKM
    */
    process stringtieFPKM {
        label 'mid_memory'
        publishDir "${params.outdir}/stringtieFPKM", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
                else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
                else "$filename"
            }


        input:
        file(bam) from bam_stringtieFPKM
        file gtf from gtf_stringtieFPKM.collect()

        output:
        file "*_transcripts.gtf"
        file "*_merged_transcripts.gtf" into stringtieGTF 
        file "*.gene_abund.txt"
        file "*.cov_refs.gtf"

        script:
        def st_direction = ''
        if (forwardStranded && !unStranded){
            st_direction = "--fr"
        } else if (reverseStranded && !unStranded){
            st_direction = "--rf"
        }
        name = bam.toString() - ~/(_R1)?(_trimmed)?(\.sorted\.bam)?(\.fq)?(\.fastq)?(\.gz)?$/

        """
        stringtie $bam \\
            $st_direction \\
            -o ${name}_transcripts.gtf \\
            -v \\
            -G $gtf \\
            -A ${name}.gene_abund.txt \\
            -C ${name}.cov_refs.gtf
        stringtie ${name}_transcripts.gtf --merge -G $gtf -o ${name}_merged_transcripts.gtf
        """
    }

}

/*
 * STEP N - Unify
 */
 process unify {
    label 'low_memory'
    publishDir "${params.outdir}/unify", mode: 'copy'

    input:
    file tx from gtf_unify
    // stringtie
    file stie from stringtieGTF
    // blat
    file spades from pasa2unify
    // ucsc
    file ucsc from ucsc2unify

    output:
    file "novel.gtf"
    file "complete_sorted.gtf"

    script:
    """
    cat $stie $spades $ucsc > full.gtf
    gffread -M --cluster-only -T  -F --keep-genes --keep-exon-attrs full.gtf > clustered.gtf
    unify.py --gtf clustered.gtf > novel.gtf
    grep "^#" novel.gtf > complete.gtf
    cat $tx novel.gtf >> complete.gtf
    bedtools sort -header -i complete.gtf > complete_sorted.gtf
    """
 }


/*
 * STEP N - QC
 */
process qc {
    publishDir "${params.outdir}/QC", mode: 'copy'
    label 'mid_memory'

    when:
    false

    input:
    file(tx) from spades2qc
    file(blastn) from blastn2qc
    file(txlong) from tdlong2qc
    file(blastp) from blastp2qc
    file(bed) from tdpredict2qc
    file(pasa) from pasa2qc
    file(protein) from prot_fasta_qc

    output:
    file "assembly_qc.db" into qc

    script:
    """
    make_db.py $prot_fasta_qc $spades2qc $blastn2qc $tdlong2qc $blastp2qc $tdpredict2qc $pasa2qc
    """
}


/*
 * STEP N - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    when:
    false

    input:
    file multiqc_config from ch_multiqc_config
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    //file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from create_workflow_summary(summary)
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config -m hisat2 .
    """
}


/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/rnassembly] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/rnassembly] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/rnassembly] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/rnassembly] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/rnassembly] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/rnassembly] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/rnassembly]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/rnassembly]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnassembly v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
