
######### Mutect2 Tumor/Normal bam -> filtered VCF ###############################################################
## Adapted from Gatk Best Practice WDL found here:
## https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2.wdl
#
## Inputs:
# Mutect2.intervals - Intervals for Variant Calling. For exome, can use the set of Exome Targets. For WGS,
#  can use a set of callable regions. This will be used to create scatter groups. See:
#  https://gatk.broadinstitute.org/hc/en-us/articles/360035889551-When-should-I-restrict-my-analysis-to-specific-intervals-
#  This intervals files should be from the same reference used as alignment and also as input, in the GATK4 intervals
#  format.
# Mutect2.scatter_count - Number of scatter groups for Mutect2 jobs.
# Mutect2.ref_fasta, ref_fai, and ref_dict - Reference file. Should be same reference file used in alignment
#  for the input bam files.
# Mutect2.gatk_docker - Used "broadinstitute/gatk:4.1.9.0" in development.
# Mutect2.tumor_bam, tumor_bai, normal_bam and normal_bai - Aligned reads for tumor and normal samples.

## Optional:
# Mutect2.pon and pon_idx - Can include a panel of normals. This may help in removing alignment/assay artifacts,
#   but is not required.
# Mutect2.gnomad, gnomad_idx - Can be included to help remove germline specific SNPs (VCF). Can be retrieved from
#  here: http://gnomad.broadinstitute.org/downloads. Useful for tumor only calling, but may also be helpful in
#  tumor/normal calling as well.
# Mutect2.m2_extra_args - Any extra arguments to pass to the Mutect2 command directly.
# Mutect2.make_bamout - Will create an output file from Mutect2.
# Mutect2.run_ob_filter - This will run the Orientation Bias filter which is highly recommended for FFPE samples.
#  This may not be needed for other samples types (e.g. fresh/frozen) but should not be detrimental in terms of
#  mutation calls.
# Mutect2.compress - Will create compressed (vcf.gz) VCF files (and indexed with tabix).
#####################################################################################################################

# WDL Version 1.0
version 1.0

struct Runtime {
	String gatk_docker
	Int cpu
	Int command_mem
	Int machine_mem
}

workflow Mutect2 {
	input {
    File intervals
		Int scatter_count
		File ref_fasta
		File ref_fai
		File ref_dict
    String gatk_docker
    File tumor_bam
    File tumor_bai
    File normal_bam
    File normal_bai

    # Optional
    File? pon
    File? pon_idx
    File? gnomad
    File? gnomad_idx
    String? m2_extra_args
    Boolean? make_bamout
    Boolean? run_ob_filter
    Boolean? compress
  }

  Boolean make_bamout_or_default = select_first([make_bamout, false])
  Boolean run_ob_filter_or_default = select_first([run_ob_filter, false])
  Boolean compress_or_default = select_first([compress, false])

	# logic about output file names -- these are the names *without* .vcf extensions
  String output_basename = basename(basename(tumor_bam, ".bam"),".cram")  #hacky way to strip either .bam or .cram
  String unfiltered_name = output_basename + "-unfiltered"
  String filtered_name = output_basename + "-filtered"

  Runtime standard_runtime = {
     "gatk_docker": gatk_docker,
     "cpu": 2,
     "command_mem": 7000,
     "machine_mem": 8000
  }

  call SplitIntervals {
    input:
      intervals = intervals,
      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      ref_dict = ref_dict,
      scatter_count = scatter_count,
      runtime_params = standard_runtime
  }

    scatter (subintervals in SplitIntervals.interval_files ) {
        call M2 {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                pon = pon,
                pon_idx = pon_idx,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                m2_extra_args = m2_extra_args,
                make_bamout = make_bamout_or_default,
                run_ob_filter = run_ob_filter_or_default,
                compress = compress_or_default,
                runtime_params = standard_runtime
        }
    }

		if(run_ob_filter_or_default) {
			call LearnReadOrientationModel {
				input:
					f1r2_tar_gz = M2.f1r2_counts,
					runtime_params = standard_runtime,
			}
		}

		call MergeVCFs {
        input:
            input_vcfs = M2.unfiltered_vcf,
            input_vcf_indices = M2.unfiltered_vcf_idx,
            output_name = unfiltered_name,
            compress = compress_or_default,
            runtime_params = standard_runtime
    }

    call MergeStats {
      input:
        stats = M2.stats,
        runtime_params = standard_runtime
    }

		call Filter {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            unfiltered_vcf = MergeVCFs.merged_vcf,
            unfiltered_vcf_idx = MergeVCFs.merged_vcf_idx,
            output_name = filtered_name,
            compress = compress_or_default,
            mutect_stats = MergeStats.merged_stats,
            #contamination_table = CalculateContamination.contamination_table,
            #maf_segments = CalculateContamination.maf_segments,
            artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
            #m2_extra_filtering_args = m2_extra_filtering_args,
            runtime_params = standard_runtime,
    }

  output {
    File filtered_vcf = Filter.filtered_vcf
		File filtered_vcf_idx = Filter.filtered_vcf_idx
		File filtering_stats = Filter.filtering_stats
  }

}

task SplitIntervals {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count
      String? split_intervals_extra_args

      # runtime
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=/root/gatk.jar

        mkdir interval-files
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" SplitIntervals \
            -R ~{ref_fasta} \
            ~{"-L " + intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    }

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        cpu: runtime_params.cpu
    }

    output {
      Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File tumor_bam
      File tumor_bai
      File? normal_bam
      File? normal_bai
      File? pon
      File? pon_idx
      File? gnomad
      File? gnomad_idx
      String? m2_extra_args
      Boolean? make_bamout
      Boolean? run_ob_filter
      Boolean compress

			Runtime runtime_params
    }

    String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    command <<<
        set -e

        export GATK_LOCAL_JAR=/root/gatk.jar

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        touch f1r2.tar.gz
        echo "" > normal_name.txt

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{tumor_bam} -O tumor_name.txt -encode
        tumor_command_line="-I ~{tumor_bam} -tumor `cat tumor_name.txt`"

        if [[ ! -z "~{normal_bam}" ]]; then
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{normal_bam} -O normal_name.txt -encode
            normal_command_line="-I ~{normal_bam} -normal `cat normal_name.txt`"
        fi

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" Mutect2 \
            -R ~{ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            -O "~{output_vcf}" \
            ~{true='--bam-output bamout.bam' false='' make_bamout} \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{m2_extra_args} \
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        cpu: runtime_params.cpu
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
    }
}

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter 
task LearnReadOrientationModel {
    input {
      Array[File] f1r2_tar_gz
      Runtime runtime_params
      Int? mem  #override memory
    }

    Int machine_mem = select_first([mem, runtime_params.machine_mem])
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx~{command_mem}m" LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    }

    runtime {
        docker: runtime_params.gatk_docker
        memory: machine_mem + " MB"
        cpu: runtime_params.cpu
    }

    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }

}

task MergeVCFs {
    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      String output_name
      Boolean compress
      Runtime runtime_params
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    }

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}
task MergeStats {
    input {
      Array[File]+ stats
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"


        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats
    }

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        cpu: runtime_params.cpu
    }

    output {
        File merged_stats = "merged.stats"
    }
}

task Filter {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File unfiltered_vcf
      File unfiltered_vcf_idx
      String output_name
      Boolean compress
      File? mutect_stats
      File? artifact_priors_tar_gz
      File? contamination_table
      File? maf_segments
      String? m2_extra_filtering_args

      Runtime runtime_params
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    command {
        set -e

        export GATK_LOCAL_JAR="/root/gatk.jar"

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{ref_fasta} \
            -O ~{output_vcf} \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            ~{"--ob-priors " + artifact_priors_tar_gz} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats filtering.stats \
            ~{m2_extra_filtering_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "filtering.stats"
    }
}
