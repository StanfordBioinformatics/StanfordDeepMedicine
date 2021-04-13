## Code repository
1. Fork this GitHub repository
2. Clone the repository to your local computer 
  
    ```git clone {your repo URL}```
## Configure Cloud Build and Functions settings: 
2. Enable the Cloud Build API (https://cloud.google.com/build/docs/automating-builds/create-manage-triggers)
3. Enable the Cloud Functions API (https://cloud.google.com/build/docs/deploying-builds/deploy-functions#yaml_1)
4. Set the Cloud Build the Cloud Functions Developer role to Enabled
5. Connect your forked GitHub repository
## Create Cloud Build triggers
6. Create triggers for importVCF and annotateAPC
  1. trigger names:
      1. gcf-importVCF
      2. gcf-annotateAPC
  2. Set the included files filters to only build the functions when changes to relevant files are made
     ```
     Serverless/importVCF/*
     Serverless/annnotateAPC/*
     ```
  3. Specify the path to your cloudbuild.yaml file
## Create Cloud Storage bucket for you data object
7. Create (2) storage buckets (https://console.cloud.google.com/storage/browser)
  1. Create `gene222-serverless-variants-{YOUR NAME}` bucket
  2. Create `gene222-serverless-annotations-{YOUR NAME}` bucket
## Create a BigQuery dataset
8. Create a "1000g" BigQuery dataset (https://console.cloud.google.com/bigquery)
## Validate that the importVCF function works
9. Copy the data object into your bucket
  1. Open Cloud Shell using the console (>_) button in the top right of your screen, next to the question mark
  2. In the Cloud Shell, run the following command:
    
      ```
      gsutil cp gs://gene222_datasets/1000g_APC.csv gs://gene222-serverless-variants-{YOUR NAME}
      ```
10. Navigate to the BigQuery console and check that a new "1000g_APC" table has been created in the "1000g" dataset (https://console.cloud.google.com/bigquery)
## Create a logging sink for BigQuery events to trigger annotateAPC
11. Create a Logging sink for BigQuery data (https://console.cloud.google.com/logs/router)
  1. Click on the "CREATE SINK" button at the top of the page, in blue
  2. Add sink details
      1. Name: bigquery-insert-events
      2. Description: Send metadata describing BigQuery insertion events to a Pub/Sub topic
  3. Specify sink destination
      1. Sink service: Cloud Pub/Sub topic
      2. Select Cloud Pub/Sub topic: new topic (name: bigquery-insert-events)
  4. Create the query to get logs specific to BigQuery insertion events
      1. Click on "PREVIEW LOGS"
      2. Select "BigQuery Dataset" under "Resource Type" to just see BigQuery logs
      3. Select the bottom-most entry and select the following labels (click them and them select "Show matching entries")
          1. bigquery.googleapis.com
          2. projects/gbsc-gcp-class-gene222-spr21/datasets/1000g/tables/1000g_APC
          3. google.cloud.bigquery.v2.JobService.InsertJob
      4. Copy the query and paste it back into the "Create logs routing sink" page
      5. Click the "PREVIEW LOGS" button again to validate that your query is getting the right logs
      6. Go back to the "Create logs routing sink" page and smash that "CREATE SINK" button (like & subscribe)
## Validate that the importVCF function works
12. Copy the data object into your bucket
  1. Open Cloud Shell using the console (>_) button in the top right of your screen, next to the question mark
  2. In the Cloud Shell, run the following command:
    
      ```
      gsutil cp gs://gene222_datasets/1000g_APC.csv gs://gene222-serverless-variants-{YOUR NAME}
      ```
13. Navigate to the `gene222-serverless-annotations-{YOUR NAME}` and check for the annotation object
