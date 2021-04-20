# Bulk RNASeq Workflow

## Overview
Adapted from [nf-core/rnaseq](https://github.com/nf-core/rnaseq) workflow. Performs the following steps:
1. FastQC - Quality metrics for input fastq files.
2. TrimGalore - Adapter and low quality trimming (also runs FastQC post-trimming)
3. SortMeRNA - Filter rRNA Reads
4. STAR - Align Reads to reference
5. RSEM - Quantification of reads aligned to genes/transcripts.


There is one auxiliarry workflow provided in addition to the Bulk RNASeq workflow:
* prepare-genome - This pipeline will create RSEM indexes from publically available STAR index files/reference.

## Preparing Data
Input genome reference files can be retrieved from [iGenome](https://ewels.github.io/AWS-iGenomes/). The prepare-genome workflow will take the genome.fa and genes.gtf and create the STAR and RSEM index files which can be used as reference for the Bulk RNASeq workflow. *Note:* There is a STARIndex directory which contains STARIndex files for the genomes, but this may not work with newer versions of STAR/RSEM, so rebuilding these before running with the same docker containers used to process data is required. After this is run, the STAR and RSEM index files should be stored together in a bucket (under the same prefix/directory).

There are also SortMeRna reference files which can be used if filtering out rRNA reads. The example set used in development were the same files used in the [nf-core/rnaseq](https://github.com/nf-core/rnaseq). The manifest file for fasta used can be found [here](https://github.com/nf-core/rnaseq/blob/master/assets/rrna-db-defaults.txt). The files were retrieved from github, uploaded to an S3 bucket (which the AWS Genomics Workflow core EC2 instances would have access to) and configured in the rnaseq-bulk.aws.json. This is optional and only required if the rnaseq.remove_rrna option is set to True.

## RNASeq Inputs
An example set of inputs is provided in rnaseq-bulk.aws.json. 

| Parameter | Type | Description |
| --- | --- | --- |
| rnaseq.input_samples | Array of maps | Map keys: r1_fastq, r2_fastq, sample_name, strandedness (optional). FastQ should be gzipped fastq with extensions of fq.gz or fastq.gz |
| rnaseq.run_fastqc | Boolean | True to run FastQ on input fastq |
| rnaseq.trim_reads | Boolean | True to trim reads with trim_galore (and also fastqc on trimmed reads) |
| rnaseq.remove_rrna | Boolean | True to remove rRNA reads using SortMeRNA (requires database fasta files described above) |
| rnaseq.star_rsem_index_base | String (s3 bucket and prefix) | Bucket prefix where STAR and RSEM index files are stored (generated from prepare-genome workflow). These should all be in the same bucket, under the same prefix. |
| rnaseq.star_rsem_index_files | Array of Strings | List of STAR/RSEM reference index file names relative to star_rsem_index_base |
| rnaseq.sortmerna_reference_path | String (s3 bucket and prefix) | Bucket path where SortMeRNA fasta files are located. Only required if running sortmerna (remove_rrna=true) |
| rnaseq.sortmerna_reference_fasta | Array of Strings | List of SortMeRNA fasta file names relative to sortmerna_reference_path |

## RNASeq Outputs
| Task | Type | Output Description |
| --- | --- | --- |
| FastQC | HTML and zip | Files describing sequencing quality of input reads |
| TrimGalore | FastQ | Trimmed reads and FastQC output post-trimming |
| SortMeRNA | FastQ | FastQ with rRNA filtered out |
| STAR-RSEM | BAM | transcript.bam and genome.bam for each input sample|
| | Stat | Stats output directory from the RSEM run for each sample |
| | Gene Counts | Raw counts for aligned fragments per gene for each samples |
| | Transcript Counts | Raw counts for aligned fragments per transcript/isoform for each sample |
| | merged transcript tpm | TPM (transcripts per million) measurement for each transcript/isoform, includes all samples |
| | merged gene tpm | TPM mesaure for each gene, includes all samples |
| | merged transcript counts | Raw counts of aligned fragment per transcript/isoform, includes all samples |
| | merged gene counts | Raw counts of aligned fragment per gene, includes all samples |

## Running Cromwell on AWS
1. Download cromwell-55.jar from https://github.com/broadinstitute/cromwell/releases/tag/55
2. Install the aws-cli by following the instructions listed in https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html
3. Execute ```aws configure```. Follow this guide to understand the parameters to be provided https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html. In the region enter ```us-east-2```. This is because the Cromwell account has been set up in us-east-2 and the aws.conf file also points to us-east-2
4. Download the wdl, json and conf file(s) from this repo to your local machine
5. Do not change the aws.conf file
6. First you need to execute the prepare-genome workflow. Edit the prepare-genome.options.json file to add your bucket where you want the output files to be stored. These output files will need to be referenced while executing the next workflow, which is the rnaseq-bulk workflow
7. Then execute ```java -Dconfig.file=aws.conf -jar /home/utsab/Cromwell/cromwell-55.jar run prepare-genome.wdl -i prepare-genome.json -o prepare-genome.options.json```
   If you get a bucket not found error during execution, first check whether you have access to the buckets listed in prepare-genome.json by executing ```aws s3 ls <bucket-name>```
   If you do have access, copy the files listed in prepare-genome.json to your S3 bucket and provide the path in the prepare-genome.json file
   Also check the region you've configured the aws cli in, and it should match with the region listed in aws.conf(us-east-2)
8. 
