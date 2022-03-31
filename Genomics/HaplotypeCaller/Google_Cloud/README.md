# HaplotypeCaller on Google Cloud

## Getting Started on Google Cloud
Have [Google Cloud SDK](https://cloud.google.com/sdk/docs/quickstarts) installed and run:
```
gcloud init
```
This will set up your default project and grant credentials to the Google Cloud SDK. Also, provide credentials so that dsub can call Google APIs:
```
gcloud auth application-default login
```
Install [dsub](https://github.com/DataBiosphere/dsub) using the instructions made in their github page. Try to familiarize yourself with the various parameters required to execute a job on Google Cloud using dsub, as it is the job scheduler we will be using.

## Gathering Input Files
The input FASTQ files can be found at gs://genomics-public-data/platinum-genomes/fastq. Remember to use paired end FASTQ files. For example if you are using ERR194146_1.fastq.gz, you must also use ERR194146_2.fastq.gz. The reference files can be found at ```gs://bwa-reference-files```.

To create a storage bucket to store the FASTQ files and the reference files, you can follow the instructions at https://cloud.google.com/storage/docs/creating-buckets#storage-create-bucket-console.

## Downsampling
Since FASTQ files are large in size, it becomes very time consuming and expensive to work with them. As a workaround, it is possible to execute pipelines with downsampled FASTQ files. Downsampling is the concept of reducing the size of a FASTQ file to ensure faster execution of pipelines.

For downsampling we will be using the [seqtk](https://github.com/lh3/seqtk) tool. 

Create a VM by following the instructions listed at https://cloud.google.com/compute/docs/instances/create-start-instance. To connect to the instance you just created, you can look at https://cloud.google.com/compute/docs/instances/connecting-to-instance.

For the VM creation it would be easier if you choose an Ubuntu image. Once you've connected to the VM you created, execute ```sudo apt-get install build-essential``` and then follow the instructions given on the seqtk github page to install seqtk.

Once seqtk is installed, download the FASTQ files from the bucket to the VM. Execute
```
seqtk sample -s100 file1.fastq.gz 10000 > downsampled_file1.fastq
seqtk sample -s100 file2.fastq.gz 10000 > downsampled_file2.fastq
```
100 is the random seed using which seqtk downsamples the fastq file, and 10000 is the number of reads that will be present in the downsampled file. You can change the number of reads to anything less than the number of reads in the whole input(the original FASTQ file).

Once the files have been downsampled, upload the downsampled files to your bucket.

Now you are ready to execute the HaplotypeCaller pipeline.

## HaplotypeCaller
The HaplotypeCaller pipeline consists of multiple stages. We will walk through each of these stages and how to execute them using dsub, which is the job scheduler we will be using for Google Cloud. The exact dsub commands we used to run the pipeline will be provided for each stage. You will need to edit out the parts that would be different in your case(project id, zone, input bucket path, output bucket path, reference files bucket path, minimum RAM, minimum cores, etc).

Chaining of outputs is an important thing to keep in mind while executing pipelines like these. The output of one stage is more often than not used as the input to the next stage. For this reason names of output files need to be noted very carefully to ensure that the correct file name is given as input to subsequent stages.

1. [BWA](http://bio-bwa.sourceforge.net/): BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome
   A prototype dsub command for BWA would look like
   
   ```dsub --provider google-v2 --project <insert project id> --zones "<insert zone>" --logging <insert logging bucket path> --input-recursive FASTQ_INPUT=<insert FASTQ input files path> --input-recursive REFERENCE=<insert reference files path> --output OUTPUT_FILE=<insert output bucket path, with /* at the end> --min-ram <RAM of VM you want BWA to be executed on> --min-cores <vCPUs of VM you want BWA to be executed on> --image pegi3s/bwa --command 'bwa mem -t 4 -M -R "@RG\\tID:0\\tLB:Library\\tPL:Illumina\\tSM:" "${REFERENCE}"/GRCh37-lite.fa "${FASTQ_INPUT}"/<insert filename of first fastq file> "${FASTQ_INPUT}"/<insert filename of second fastq file> > "$(dirname ${OUTPUT_FILE})"/bwa-sam.sam'```
   
   The command we used for execution looked like 
   
   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input-recursive FASTQ_INPUT=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/FASTQ_INPUT --input-recursive REFERENCE=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/REFERENCE_GRCH37 --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image pegi3s/bwa --command 'bwa mem -t 4 -M -R "@RG\\tID:0\\tLB:Library\\tPL:Illumina\\tSM:" "${REFERENCE}"/GRCh37-lite.fa "${FASTQ_INPUT}"/ERR194159_1_seqtk_707646.fastq "${FASTQ_INPUT}"/ERR194159_2_seqtk_707646.fastq > "$(dirname ${OUTPUT_FILE})"/bwa-sam.sam'```
   
   After generating the SAM file you will need to convert the SAM file to a BAM file. The dsub command to do so would look like
   
   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input SAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/bwa-sam.sam --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'samtools view -bS "${SAM_INPUT}" > "$(dirname ${OUTPUT_FILE})"/output_bwa.bam'```
   
   Once the SAM file has been converted to a BAM file, the BAM file needs to be reheadered. This can be done using a command like
   
   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/output_bwa.bam --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'samtools view -H "${BAM_INPUT}" | sed -e 's/SM/SM:ERR194159/' | samtools reheader - "${BAM_INPUT}" > "$(dirname ${OUTPUT_FILE})"/output_bwa_reheadered.bam'```
   
   You would need to change the ERR194159 part of the above command with the prefix of your FASTQ file. For example if your FASTQ files are named ERR194146_1.fastq.gz and ERR194146_2.fastq.gz, you would put ERR194146 instead of ERR194159
2. [SortSam](https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard): The dsub command for this stage would look like
   
   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/output_bwa_reheadered.bam --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'gatk SortSam INPUT="${BAM_INPUT}" OUTPUT="$(dirname ${OUTPUT_FILE})"/output_picard_sorted.bam SORT_ORDER="queryname" CREATE_MD5_FILE=true CREATE_INDEX=true TMP_DIR=`pwd`/tmp'```

3. [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard): The dsub command for this stage would look like

   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/output_picard_sorted.bam --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'gatk MarkDuplicates -I "${BAM_INPUT}" -M "$(dirname ${OUTPUT_FILE})"/metrics_md -O "$(dirname ${OUTPUT_FILE})"/output_dedup_bam.bam --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --TMP_DIR `pwd`/tmp'```
   
   The dedup bam that we obtain needs to be sorted again using a command like
   
   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/output_dedup_bam.bam --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'gatk SortSam -INPUT "${BAM_INPUT}" -OUTPUT "$(dirname ${OUTPUT_FILE})"/output_dedup_sorted.bam -SORT_ORDER coordinate -CREATE_MD5_FILE true -TMP_DIR `pwd`/tmp -CREATE_INDEX true'```
   
4. [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator): The dsub command for this stage would look like

   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/output_dedup_sorted.bam --input-recursive REFERENCE_GR=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/REFERENCE_GRCH37 --input-recursive REFERENCE_DB=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/REFERENCE_DBSNP --input-recursive REFERENCE_MI=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/REFERENCE_MILLS --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'gatk BaseRecalibrator -I "${BAM_INPUT}" -R "${REFERENCE_GR}"/GRCh37-lite.fa --known-sites "${REFERENCE_DB}"/dbSNP.b150.GRCh37p13.All_20170710.vcf.gz --known-sites "${REFERENCE_MI}"/Mills_and_1000G_gold_standard.indels.b37.vcf -O "$(dirname ${OUTPUT_FILE})"/recal_data.table'```
   
5. [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR): The dsub command for this stage would look like

   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/output_dedup_sorted.bam --input RECAL=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/recal_data.table --input-recursive REFERENCE_GR=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/REFERENCE_GRCH37 --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --machine-type n1-standard-4 --image broadinstitute/gatk --boot-disk-size 20 --command 'gatk ApplyBQSR -I "${BAM_INPUT}" -R "${REFERENCE_GR}"/GRCh37-lite.fa --bqsr-recal-file "${RECAL}" -O "$(dirname ${OUTPUT_FILE})"/output_recal_bam.bam'```
   
6. [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller): The dsub command for this stage would look like

   ```dsub --provider google-v2 --project gbsc-gcp-project-cba --zones "us-east1-*" --logging gs://gbsc-gcp-project-cba_user-uray/Logging --input-recursive BAM_INPUT=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA --input-recursive REFERENCE_GR=gs://gbsc-gcp-project-cba_user-uray/Input_Bucket/REFERENCE_GRCH37 --output OUTPUT_FILE=gs://gbsc-gcp-project-cba_user-uray/Output/OUTPUT_BWA/* --min-ram 15 --min-cores 4 --image broadinstitute/gatk --boot-disk-size 20 --command 'gatk HaplotypeCaller -R "${REFERENCE_GR}"/GRCh37-lite.fa -I "${BAM_INPUT}"/output_recal_bam.bam -O "$(dirname ${OUTPUT_FILE})"/output_vcf.vcf'```
   
The output_vcf.vcf file is the final desired output of the HaplotypeCaller.
