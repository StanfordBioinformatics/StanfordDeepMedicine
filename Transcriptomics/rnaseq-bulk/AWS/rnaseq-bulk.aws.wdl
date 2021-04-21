
######### RNA Seq Bulk Pipeline  ###############################################################
# Based on NextFlow Pipeline: https://nf-co.re/rnaseq
# See README.md for more information.
################################################################################################

# WDL Version 1.0
version 1.0

struct RNASample {
  File r1_fastq
  File r2_fastq
  String? strandedness
  String sample_name
}

workflow rnaseq {
  input {
    ############# Inputs ################
    Array[RNASample] input_samples

    # sortmerna input/options
    String? sortmerna_reference_path          # Will use reference_path if not provided
    Array[String]? sortmerna_reference_fasta  # Relative to reference_path

    # star-rsem inputs/options
    String star_rsem_index_base
    Array[String] star_rsem_index_files

    Boolean run_fastqc = true
    Boolean trim_reads = true
    Boolean remove_rrna = true
    #####################################
  }

  scatter( rna_sample in input_samples ) {
    if(run_fastqc) {
      call fastqc {
        input:
          r1_fastq = rna_sample.r1_fastq,
          r2_fastq = rna_sample.r2_fastq
      }
    }

    if(trim_reads) {
      call trim_galore {
        input:
          r1_fastq = rna_sample.r1_fastq,
          r2_fastq = rna_sample.r2_fastq
      }
    }

    ##### sortmerna #####
    if(remove_rrna) {
      # Use the base reference path if none was provided for sotmerna
      String sortmerna_rp = select_first([sortmerna_reference_path, ""])

      # This will prepend the reference path to the sortmerna fasta files
      Array[String] fasta_rel_paths = select_first([sortmerna_reference_fasta,[]])
      scatter(i in fasta_rel_paths) {
        File full_paths = sortmerna_rp + "/" + i
      }
      Array[File] sortmerna_reference_fasta_full = full_paths

      File sortmerna_r1_fastq = select_first([trim_galore.r1_trimmed_fastq, rna_sample.r1_fastq])
      File? sortmerna_r2_fastq = select_first([trim_galore.r2_trimmed_fastq, rna_sample.r2_fastq]) # Will this fail when no r2 is used?

      call sortmerna {
        input:
          r1_fastq = sortmerna_r1_fastq,
          r2_fastq = sortmerna_r2_fastq,
          sample_name = rna_sample.sample_name,
          reference_fasta = sortmerna_reference_fasta_full
      }
    }
    ########################

    ##### Star-RSEM #####
    scatter(i in star_rsem_index_files) {
      File full_star_rsem_index_path = star_rsem_index_base + "/" + i
    }
    Array[File] star_rsem_index_path = full_star_rsem_index_path # Gathers the scatter into an Array

    File r1_fastq_sr = select_first([sortmerna.r1_rrna_filtered_fastq, trim_galore.r1_trimmed_fastq, rna_sample.r1_fastq])
    File? r2_fastq_sr = select_first([sortmerna.r2_rrna_filtered_fastq, trim_galore.r2_trimmed_fastq, rna_sample.r2_fastq])

    String strandedness = select_first([rna_sample.strandedness, "none"])

    call star_rsem {
      input:
        r1_fastq = r1_fastq_sr,
        r2_fastq = r2_fastq_sr,
        sample_name = rna_sample.sample_name,
        strandedness = strandedness,
        star_rsem_index_files = star_rsem_index_path
    }
    ########################
  }

  call merge_counts {
    input:
      gene_counts = star_rsem.gene_counts,
      isoform_counts = star_rsem.transcript_counts
  }

  call multiqc {
    input:
      r1_fastqc_zip = fastqc.r1_fastqc_zip,
      r2_fastqc_zip = fastqc.r2_fastqc_zip,
      stats_cnt = star_rsem.stats_cnt,
      stats_model = star_rsem.stats_model,
      stats_theta = star_rsem.stats_theta,
      trim_html = trim_galore.htmls,
      trim_zip = trim_galore.zips
  }


  output {
    Array[File?] r1_fastqc_zip = fastqc.r1_fastqc_zip
    Array[File?] r1_fastqc_html = fastqc.r1_fastqc_html
    Array[File?] r2_fastqc_zip = fastqc.r2_fastqc_zip
    Array[File?] r2_fastqc_html = fastqc.r2_fastqc_html

    Array[Array[File]?] trim_html = trim_galore.htmls
    Array[Array[File]?] trim_zip = trim_galore.zips

    # RSEM STAR output
    Array[File] gene_counts = star_rsem.gene_counts
    Array[File] transcript_counts = star_rsem.transcript_counts
    Array[File] log = star_rsem.log
    Array[File] genome_bam = star_rsem.genome_bam
    Array[File] transcript_bam = star_rsem.transcript_bam

    Array[File] stats_cnt = star_rsem.stats_cnt
    Array[File] stats_model = star_rsem.stats_model
    Array[File] stats_theta = star_rsem.stats_theta

    # Merged files
    File merged_gene_counts = merge_counts.merged_gene_counts
    File merged_gene_tpm = merge_counts.merged_gene_tpm
    File merged_transcript_counts = merge_counts.merged_transcript_counts
    File merged_transcript_tpm = merge_counts.merged_transcript_tpm

    # MultiQC Report
    File multiqc_report = multiqc.multiqc_report
  }
}


###################################### Tasks #########################################
task multiqc {
  # Inputs are specific to this pipeline. More data can be added
  # to extend the MultiQC Report
  input {
    Array[File?] r1_fastqc_zip
    Array[File?] r2_fastqc_zip

    Array[File] stats_cnt
    Array[File] stats_model
    Array[File] stats_theta

    Array[Array[File]?] trim_html
    Array[Array[File]?] trim_zip

    String docker = "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    Int cpus = 1
    String memory_mb = "8000"
  }

  command {
    multiqc .
  }

  runtime {
    docker: docker
    cpu: cpus
    memory: memory_mb + " MB"
  }

  output  {
    File multiqc_report = "multiqc_report.html"
  }

}

task merge_counts {
  input {
    Array[File] gene_counts
    Array[File] isoform_counts

    String docker = "biocontainers/biocontainers:v1.2.0_cv1"
    Int cpus = 1
    String memory_mb = "2000"
  }

  command {
    set -e
    set -x

    cut -f 1,2 ~{gene_counts[0]} > gene_ids.txt
    mkdir -p tmp/genes
    for i in ~{sep=' ' gene_counts}; do
      filename=$(basename $i)
      sample_name=$(echo "$filename" | sed 's/.genes.results//')
      echo "$sample_name" > tmp/genes/$sample_name.counts.txt
      cut -f 5 $i | tail -n+2 >> tmp/genes/$sample_name.counts.txt

      echo "$sample_name" > tmp/genes/$sample_name.tpm.txt
      cut -f 6 $i | tail -n+2 >> tmp/genes/$sample_name.tpm.txt
    done

    cut -f 1,2 ~{isoform_counts[0]} > transcript_ids.txt
    mkdir -p tmp/isoforms
    for i in ~{sep=' ' isoform_counts}; do
      filename=$(basename $i)
      sample_name=$(echo "$filename" | sed 's/.isoforms.results//')
      echo "$sample_name" > tmp/isoforms/$sample_name.count.txt
      cut -f 5 $i | tail -n+2 >> tmp/isoforms/$sample_name.counts.txt

      echo "$sample_name" > tmp/isoforms/$sample_name.tpm.txt
      cut -f 6 $i | tail -n+2 >> tmp/isoforms/$sample_name.tpm.txt
    done

    paste gene_ids.txt tmp/genes/*.counts.txt > rsem.merged.gene_counts.tsv
    paste gene_ids.txt tmp/genes/*.tpm.txt > rsem.merged.gene_tpm.tsv
    paste transcript_ids.txt tmp/isoforms/*.counts.txt > rsem.merged.transcript_counts.tsv
    paste transcript_ids.txt tmp/isoforms/*.tpm.txt > rsem.merged.transcript_tpm.tsv
  }

  runtime {
    docker: docker
    cpu: cpus
    memory: memory_mb + " MB"
  }

  output {
    File merged_gene_counts = "rsem.merged.gene_counts.tsv"
    File merged_gene_tpm = "rsem.merged.gene_tpm.tsv"
    File merged_transcript_counts = "rsem.merged.transcript_counts.tsv"
    File merged_transcript_tpm = "rsem.merged.transcript_tpm.tsv"
  }
}


######################## STAR - RSEM ##############################################################
task star_rsem {
  # Will run rsem-calculate-expression with star aligner. Will take ~ 6 minutes to localize
  # reference data.
  input {
    File r1_fastq
    File? r2_fastq
    Array[File] star_rsem_index_files
    String sample_name
    String strandedness

    # If the index was built with provided workflow, this should hold.
    # Otherwise, will need to overwrite with whatever the reference name was used 
    # in the rsem-prepare-reference command
    String reference_name = "genome"
    String docker = "quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    Int cpus = 8
  }

  String other_opts = if defined(r2_fastq) then "--paired-end" else ""

  command {
    set -e

    star_rsem_index_basedir=$(dirname ~{star_rsem_index_files[0]})
    reference_path="$star_rsem_index_basedir/~{reference_name}"

    rsem-calculate-expression \
        --num-threads ~{cpus} \
        --temporary-folder ./tmp/ \
        ~{other_opts} \
        --strandedness ~{strandedness} \
        --star \
        --star-output-genome-bam \
        --star-gzipped-read-file \
        --estimate-rspd \
        --seed 1 \
        ~{r1_fastq} ~{r2_fastq} \
        $reference_path \
        ~{sample_name}
  }

  runtime {
    docker: docker
    cpu: cpus
    memory: 32000 + " MB"
  }

  output {
    File gene_counts = "~{sample_name}.genes.results"
    File transcript_counts = "~{sample_name}.isoforms.results"
    File log = "~{sample_name}.log"
    File genome_bam = "~{sample_name}.STAR.genome.bam"
    File transcript_bam = "~{sample_name}.transcript.bam"

    File stats_cnt = "~{sample_name}.stat/~{sample_name}.cnt"
    File stats_model = "~{sample_name}.stat/~{sample_name}.model"
    File stats_theta = "~{sample_name}.stat/~{sample_name}.theta"
  }
}

###############################################################################
task sortmerna {
  input {
    File r1_fastq
    File? r2_fastq
    Array[File] reference_fasta
    String sample_name

    String docker = "quay.io/biocontainers/sortmerna:4.2.0--0"
    Int cpus = 4
  }



  String other_opts = if defined(r2_fastq) then "--paired_in" else ""

  command {
    set -e
    sortmerna \
      ~{sep=' ' prefix("--ref ", reference_fasta)} \
      --reads ~{r1_fastq} ~{"--reads " + r2_fastq} \
      --threads ~{cpus} \
      --workdir . \
      --aligned rRNA_reads \
      --other non_rRNA_reads \
      ~{other_opts} \
      --out2 \
      --num_alignments 1 \
      --fastx \
      -v

      gzip -f non_rRNA_reads_fwd.fastq > ~{sample_name}_1_rrna_filtered.fastq.gz

      if [[ -e "non_rRNA_reads_rev.fastq" ]]
      then
        gzip -f non_rRNA_reads_rev.fastq > ~{sample_name}_2_rrna_filtered.fastq.gz
      fi

  }

  runtime {
    docker: docker
    cpu: cpus
    memory: 2000 + " MB"
  }

  output {
    File r1_rrna_filtered_fastq = "~{sample_name}_1_rrna_filtered.fastq.gz"
    File? r2_rrna_filtered_fastq = "~{sample_name}_2_rrna_filtered.fastq.gz"
  }
}

###############################################################################
task trim_galore {
  input {
    File r1_fastq
    File? r2_fastq
    
    String docker = "quay.io/biocontainers/trim-galore:0.6.6--0"
    Int cpus = 4
  }

  String r1_prefix_tmp = basename(r1_fastq, ".fastq.gz")
  String r1_prefix = basename(r1_prefix_tmp, ".fq.gz")

  String r2_filename_or_nothing = select_first([r2_fastq, ""])
  String r2_prefix_tmp = if defined( r2_fastq ) then basename(r2_filename_or_nothing, ".fastq.gz") else ""
  String r2_prefix = if defined( r2_fastq ) then basename(r2_prefix_tmp, ".fq.gz") else ""
  
  command {
    set -e
    trim_galore --fastqc --cores ~{cpus} --paired --gzip ~{r1_fastq} ~{r2_fastq}

    # Rename some files
    mv ~{r1_prefix}_val_1.fq.gz ~{r1_prefix}_trimmed.fastq.gz
    mv ~{r1_prefix}_val_1_fastqc.html ~{r1_prefix}_trimmed_fastqc.html
    mv ~{r1_prefix}_val_1_fastqc.zip ~{r1_prefix}_trimmed_fastqc.zip

    if [[ -e "~{r2_prefix}_val_2.fq.gz" ]]
    then
      mv ~{r2_prefix}_val_2.fq.gz ~{r2_prefix}_trimmed.fastq.gz
      mv ~{r2_prefix}_val_2_fastqc.html ~{r2_prefix}_trimmed_fastqc.html
      mv ~{r2_prefix}_val_2_fastqc.zip ~{r2_prefix}_trimmed_fastqc.zip
    fi
  }


  runtime {
    docker: docker
    cpu: cpus
    memory: 8000 + " MB"
  }

  output {
    File r1_trimmed_fastq = "~{r1_prefix}_trimmed.fastq.gz"
    File? r2_trimmed_fastq = "~{r2_prefix}_trimmed.fastq.gz"
    Array[File] htmls = glob("*.html")
    Array[File] zips = glob("*.zip")
  }
}

###############################################################################
task fastqc {
  input {
    File r1_fastq
    File? r2_fastq
    

    String docker = "quay.io/biocontainers/fastqc:0.11.9--0"
    Int cpus = 4
  }

  # Will work with either fastq.gz or fq.gz
  String r1_prefix_tmp = basename(r1_fastq, ".fastq.gz")
  String r1_prefix = basename(r1_prefix_tmp, ".fq.gz")

  # Will work with either fastq.gz or fq.gz
  String r2_filename_or_nothing = select_first([r2_fastq, ""])
  String r2_prefix_tmp = if defined( r2_fastq ) then basename(r2_filename_or_nothing, ".fastq.gz") else ""
  String r2_prefix = if defined( r2_fastq ) then basename(r2_prefix_tmp, ".fq.gz") else ""

  command {
    set -e
    fastqc --threads ~{cpus} --outdir ./ ~{r1_fastq} ~{r2_fastq}
  }

  runtime {
    docker: docker
    memory: 2000 + " MB"
    cpu: cpus
  }

  output {
    File r1_fastqc_zip  = "~{r1_prefix}_fastqc.zip"
    File r1_fastqc_html = "~{r1_prefix}_fastqc.html"

    File? r2_fastqc_zip  = "~{r2_prefix}_fastqc.zip"
    File? r2_fastqc_html = "~{r2_prefix}_fastqc.html"
  }
}
###############################################################################

