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

## RNASeq Inputs
An example set of inputs is provided in rnaseq-bulk.json. 

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

## Preparing Data
All the data needs to be transfered from S3 buckets to Google Cloud Storage buckets. 
1. Create a VM on Google Cloud
2. Login to the VM and install aws-cli
3. Execute ```aws configure```
4. In the AWS directory of [rnaseq-bulk](https://github.com/StanfordBioinformatics/StanfordDeepMedicine/tree/main/Transcriptomics/rnaseq-bulk/AWS) look at the file paths of the input files in prepare-genome.json

   Execute ```aws s3 cp s3://<bucket_path> .```

   This will copy all the input files from S3 to your Google Cloud VM.

   Now from the VM execute ```gsutil cp <input_files> gs://<gcp_storage_bucket_path>```

   Now you have all the input files that are needed for executing the workflow

## Running Cromwell on Google Cloud
1. For setup instructions refer to https://cloud.google.com/life-sciences/docs/tutorials/gatk or the instructions given in HW2. Both are the same however the instructions in HW2 are a bit more verbose
2. Edit the prepare-genome.json file to correctly reflect the file paths of the input files. Keep in mind that you'll have to copy over the files from S3 buckets to Google Cloud Storage buckets, and the instructions to do so are given in the previous section
3. Execute the prepare-genome workflow using the lifesciences API. Below is a sample command that I used for execution:

   ```gcloud beta lifesciences pipelines run --pipeline-file wdl_pipeline.yaml --location us-central1 --regions us-central1 --inputs-from-file WDL=prepare-genome.wdl,WORKFLOW_INPUTS=prepare-genome.json,WORKFLOW_OPTIONS=prepare-genome.options.json --env-vars WORKSPACE=gs://utsab-test-bucket/cromwell-out-prepare-genome/work,OUTPUTS=gs://utsab-test-bucket/cromwell-out-prepare-genome/out --logging gs://utsab-test-bucket/cromwell-out-prepare-genome/logging```
4. Once the prepare-genome workflow is done executing, take a note of the output directory which you'll need to add in to the ```rnaseq.star_rsem_index_base``` parameter in rnaseq-bulk.json.    
5. Execute the rnaseq-bulk workflow using the lifesciences API. Below is a sample command that I used for execution:

   ```gcloud beta lifesciences pipelines run --pipeline-file wdl_pipeline.yaml --location us-central1 --regions us-central1 --inputs-from-file WDL=rnaseq-bulk.wdl,WORKFLOW_INPUTS=rnaseq-bulk.json,WORKFLOW_OPTIONS=rnaseq-bulk.options.json --env-vars WORKSPACE=gs://utsab-test-bucket/cromwell-out-rnaseq-bulk/work,OUTPUTS=gs://utsab-test-bucket/cromwell-out-rnaseq-bulk/out --logging gs://utsab-test-bucket/cromwell-out-rnaseq-bulk/logging```




