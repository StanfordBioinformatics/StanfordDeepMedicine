## Code repository
1. Fork this GitHub repository
## Configure Cloud Build and Functions settings: https://cloud.google.com/build/docs/deploying-builds/deploy-functions#yaml_1
2. Enable the Cloud Build API 
3. Enable the Cloud Functions API
4. Set the Cloud Build the Cloud Functions Developer role to Enabled
5. Connect your forked GitHub repository
## Create Cloud Build triggers
6. Create triggers for importVCF and annotateAPC
  1. Set the included files filters to only build the functions when changes to relevant files are made
  2. Specify the path to your cloudbuild.yaml file
## Create Cloud Storage bucket for you data object
7. Create the `gene222-serverless-demo-{YOUR NAME}` bucket (https://console.cloud.google.com/storage/browser)
## Create a BigQuery dataset
8. Create a "1000g" BigQuery dataset (https://console.cloud.google.com/bigquery)
## Validate that the importVCF function works
9. Copy the data object into your bucket
  1. Open Cloud Shell using the console (>_) button in the top right of your screen, next to the question mark
  2. In the Cloud Shell, run the following command:
    ```gsutil cp gs://gene222_datasets/1000g_APC.csv gs://{YOUR BUCKET}```
10. Navigate to the BigQuery console and check that a new "1000g_APC" table has been created in the "1000g" dataset (https://console.cloud.google.com/bigquery)
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
    2.
