# AWS Environment Setup

## Getting Started on AWS
Have [aws cli](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) installed and run:
```
aws configure
```
Also, create a key pair using the instructions on this [website](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html)

## Cromwell
We will be using [Cromwell](https://cromwell.readthedocs.io/en/stable/) as the job scheduler. Follow this [tutorial](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) to get an idea of how to configure Cromwell on your local machine.

Once you have Cromwell up and running on your local machine we will have to configure AWS Batch such that you can execute workflows on it using Cromwell.

This [website](https://docs.opendata.aws/genomics-workflows/) gives an overview of how to configure the AWS environment properly. Go to the Cromwell [subsection](https://docs.opendata.aws/genomics-workflows/orchestration/cromwell/cromwell-overview/) of the website and launch the three CloudFormation templates provided. These CloudFormation templates will help in launch and configuring all the necessary resources required to execute workflows using Cromwell on AWS. Keep in mind the stackname you provide during launching the genomics workflow core has to be provided as a reference to the CloudFormation template used for launching Cromwell.

Now you have to create a configuration file which will be used by Cromwell during execution. Under subsection "Configuring Cromwell to use AWS Batch" on this [website](https://docs.opendata.aws/genomics-workflows/orchestration/cromwell/cromwell-overview/) you will find a sample configuration file. You can remove the database section of the configuration file as we will not be using that capability of Cromwell for this pipeline. Another sample configuration file can be found [here](https://github.com/StanfordBioinformatics/GENE222/blob/main/Genomics/aws.conf).

Fill in the region, the root bucket, queue ARN and script bucket. The queue ARN will be of the batch queue which was created previously using the CloudFormation template. Go to the Batch section of the AWS console to find the previously created queue. Keep the script bucket to be one level higher than the root bucket. For example if your root bucket is s3://course/cromwell-execution, keep the script bucket as s3://course. Do not use the s3 prefix while entering the script bucket in the configuration file.

Now we are ready to execute pipelines on AWS using Cromwell.

# Google Cloud Environment Setup

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

For executing WDL files on Google Cloud follow the tutorial at: https://cloud.google.com/life-sciences/docs/tutorials/gatk
