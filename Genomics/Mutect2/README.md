# [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) on Google Cloud

Executing wdl files on Google Cloud requires usage of the Google Cloud Life Sciences API. For setup plpease refer to [this link](https://cloud.google.com/life-sciences/docs/tutorials/gatk). The input files that we will be using for this pipeline are all publicly available at ```gs://gatk-tutorials/workshop_2002/3-somatic```. The exact path of the input files that we require can be found in the Mutect2.gcp.json file in this directory.

The following command illustrates how to run Mutect2 on Google Cloud. The current directory in which the command below was executed is wdl-runner/wdl-runner. The Mutect2 workflow files were stored under that directory in a folder named Mutect2.

```gcloud beta lifesciences pipelines run --pipeline-file wdl_pipeline.yaml --location us-central1 --regions us-central1 --inputs-from-file WDL=./Mutect2/Mutect2.gcp.wdl,WORKFLOW_INPUTS=./Mutect2/Mutect2.gcp.json,WORKFLOW_OPTIONS=./Mutect2/Mutect2.gcp.options.json --env-vars WORKSPACE=gs://utsab-test-bucket/cromwell-out1/work,OUTPUTS=gs://utsab-test-bucket/cromwell-out1 --logging gs://utsab-test-bucket/logging```

The final output files will be available at ```gs://utsab-test-bucket/cromwell-out1```
