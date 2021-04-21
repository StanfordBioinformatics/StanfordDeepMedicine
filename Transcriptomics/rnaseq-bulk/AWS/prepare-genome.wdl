# © 2021 Amazon Web Services, Inc. or its affiliates. All Rights Reserved.
# 
# This AWS Content is provided subject to the terms of the AWS Customer Agreement
# available at http://aws.amazon.com/agreement or other written agreement between
# Customer and either Amazon Web Services, Inc. or Amazon Web Services EMEA SARL or both.

version 1.0

workflow RNASeqBulkPrepareGenome {
    input {
        File genes_gtf
        File genome_fa
    }

    call rsem_prepare_genome {
        input:
            genes_gtf = genes_gtf,
            genome_fa = genome_fa
    }

    call star_prepare_genome {
        input:
            genes_gtf = genes_gtf,
            genome_fa = genome_fa
    }

    output {
        File exon_ge_tab = star_prepare_genome.exon_ge_tab
        File gene_info_tab = star_prepare_genome.gene_info_tab
        File transcript_info_tab = star_prepare_genome.transcript_info_tab
        File exon_info_tab = star_prepare_genome.exon_info_tab
        File gtf_out_tab = star_prepare_genome.gtf_out_tab
        File chr_name = star_prepare_genome.chr_name
        File chr_start = star_prepare_genome.chr_start
        File chr_length = star_prepare_genome.chr_length
        File chr_name_length = star_prepare_genome.chr_name_length
        File sjdb_info = star_prepare_genome.sjdb_info
        File sjdb_list_tab = star_prepare_genome.sjdb_list_tab
        File genome_parameters = star_prepare_genome.genome_parameters
        File genome = star_prepare_genome.genome
        File sa = star_prepare_genome.sa
        File sa_index = star_prepare_genome.sa_index

        File rsem_grp = rsem_prepare_genome.rsem_grp
        File rsem_ti = rsem_prepare_genome.rsem_ti
        File rsem_transcripts_fa = rsem_prepare_genome.rsem_transcripts_fa
        File rsem_seq = rsem_prepare_genome.rsem_seq
        File rsem_chrlist = rsem_prepare_genome.rsem_chrlist
        File rsem_idx_fa = rsem_prepare_genome.rsem_idx_fa
        File rsem_n2g_idx_fa = rsem_prepare_genome.rsem_n2g_idx_fa

        File star_version = star_prepare_genome.star_version
        File rsem_version = rsem_prepare_genome.rsem_version
    }
}

task star_prepare_genome {
    input {
        File genome_fa
        File genes_gtf

        String docker = "quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
        Int cpus = 8
        Int memory_mb = 64000
    }

    command {
        STAR --runMode genomeGenerate \
        --genomeDir ./ \
        --genomeFastaFiles ~{genome_fa} \
        --sjdbGTFfile ~{genes_gtf} --runThreadN ~{cpus}

        STAR --version > star.version.txt
    }

    runtime {
        docker: docker
        cpu: cpus
        memory: memory_mb + " MB"
    }

    output {
        File exon_ge_tab = "exonGeTrInfo.tab"
        File gene_info_tab = "geneInfo.tab"
        File transcript_info_tab = "transcriptInfo.tab"
        File exon_info_tab = "exonInfo.tab"
        File gtf_out_tab = "sjdbList.fromGTF.out.tab"
        File chr_name = "chrName.txt"
        File chr_start = "chrStart.txt"
        File chr_length = "chrLength.txt"
        File chr_name_length = "chrNameLength.txt"
        File sjdb_info = "sjdbInfo.txt"
        File sjdb_list_tab = "sjdbList.out.tab"
        File genome_parameters = "genomeParameters.txt"
        File genome = "Genome"
        File sa = "SA"
        File sa_index = "SAindex"

        File star_version = "star.version.txt"
    }
}

task rsem_prepare_genome {
    input {
        File genome_fa
        File genes_gtf

        String docker = "quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
        Int cpus = 4
        Int memory_mb = 8000
    }

    command {
        set -e

        rsem-prepare-reference \
            --gtf ~{genes_gtf} \
            --num-threads ~{cpus} \
            ~{genome_fa} \
            genome

        rsem-calculate-expression --version > rsem.version.txt
    }

    runtime {
        docker: docker
        cpu: cpus
        memory: memory_mb + " MB"
    }

    output {
        File rsem_grp = "genome.grp"
        File rsem_ti = "genome.ti"
        File rsem_transcripts_fa = "genome.transcripts.fa"
        File rsem_seq = "genome.seq"
        File rsem_chrlist = "genome.chrlist"
        File rsem_idx_fa = "genome.idx.fa"
        File rsem_n2g_idx_fa = "genome.n2g.idx.fa"

        File rsem_version = "rsem.version.txt"
    }
}