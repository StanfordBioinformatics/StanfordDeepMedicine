import os
from google.cloud import pubsub_v1
from google.cloud import storage
from google.cloud import bigquery

# Instantiates a Pub/Sub client
publisher = pubsub_v1.PublisherClient()

# TODO(developer): Set dataset_id to the ID of the dataset.
dataset_id = '1000g'

# TODO(developer): Set project_id to the ID of the GCP project.
project_id = 'gbsc-gcp-class-gene222-spr21'

# TODO(developer): Set topic_id to the ID of the PubSub topic.
topic_id = 'projects/gbsc-gcp-class-gene222-spr21/topics/test'


def main(event, context):
    """Triggered by a change to a Cloud Storage bucket.
    Args:
         event (dict): Event payload.
         context (google.cloud.functions.Context): Metadata for the event.
    """

    print('Event ID: {}'.format(context.event_id))
    print('Event type: {}'.format(context.event_type))
    print('Bucket: {}'.format(event['bucket']))
    print('File: {}'.format(event['name']))

    try:
        client = bigquery.Client()

        # TODO(developer): Set table_id to the ID of the table to create.
        #table_id = "gbsc-gcp-class-gene222-spr21.1000g.sampleVCF"
        
#        table_id = '{}.{}.{}'.format(project_id,dataset_id,'sampleVCF')
        table_id = '{}.{}.{}'.format(project_id,dataset_id,os.path.splitext(event['name'])[0])

        #uri = "gs://gene222_datasets/sample.csv"
        uri = 'gs://{}/{}'.format(event['bucket'],event['name'])

        job_config = bigquery.LoadJobConfig(
            schema=[
                bigquery.SchemaField("chrm", "STRING"),
                bigquery.SchemaField("start_position", "INTEGER"),
                bigquery.SchemaField("end_position", "INTEGER"),
                bigquery.SchemaField("reference_bases", "STRING"),
                bigquery.SchemaField("alternate_bases", "STRING"),
                bigquery.SchemaField("rsID", "STRING"),
                bigquery.SchemaField("qual", "STRING"),
                bigquery.SchemaField("filter", "STRING"), 
                bigquery.SchemaField("info", "STRING"),                    
            ],
            skip_leading_rows=1
        )


        load_job = client.load_table_from_uri(
            uri, table_id, job_config=job_config
        )  # Make an API request.

        # Check whether table exists and create if not
        try: 
            table = client.get_table(table_id)
        except:
            table = bigquery.Table(table_id, schema=schema)

        load_job.result()  # Wait for the job to complete.

        table = client.get_table(table_id)
        print("Loaded {} rows to table {}".format(table.num_rows, table_id))

        print(f"Push Message")
        tableAttribute = table_id 
        print(topic_id)
        publisher.publish(topic_id, b'', tableAttribute=tableAttribute)
        
        return f"OK"
    except Exception as e:
        print(e)
        return (e, 500)    



#    client = storage.Client()
#    vcf_bucket = client.get_bucket(event['bucket'])
#    print(vcf_bucket)
#    vcf_blob = vcf_bucket.get_blob(event['name'])
#    print(vcf_blob)

#    data_string = vcf_blob.download_as_string()
# Construct a BigQuery client object.
