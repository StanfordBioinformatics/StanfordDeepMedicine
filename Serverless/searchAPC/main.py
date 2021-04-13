import base64
from google.cloud import bigquery
from google.cloud import storage

def main(event, context):
    """Triggered from a message on a Cloud Pub/Sub topic.
    Args:
         event (dict): Event payload.
         context (google.cloud.functions.Context): Metadata for the event.
    """
    
    print(f"This Function was triggered by messageId {context.event_id} published at {context.timestamp}.")

    #if 'data' in event:
    data = base64.b64decode(event['data']).decode('utf-8')
    print(f"> Pubsub data: {data}.")
    print(f"> Context: {context}.")

    data = json.loads(data)

    # Example: "projects/gbsc-gcp-class-gene222-spr21/datasets/1000g/tables/1000g_APC_PBR"
    resource_name = data["protoPayload"]["resourceName"]
    elements = resouce_name.split("/")
    project_id = elements[1]
    dataset_id = elements[3]
    table_name = elements[5]

    table_id = f"{project_id}.{dataset_id}.{table_name}" 

    client = bigquery.Client()
    # Perform a query.
    QUERY = (
      'SELECT rsID ' +
      f"FROM `{table_id}` " +
      'WHERE start_position > 112073558 AND start_position < 112181934 ' +
      'LIMIT 100')
    query_job = client.query(QUERY)  # API request
    rows = query_job.result()  # Waits for query to finish
    
    rows_string = ""
    for row in rows:
      rows_string += row + "\n"

    storage_client = storage.Client()
    bucket = storage_client.get_bucket(STORAGE_BUCKET)
    composed_blob = bucket.blob("apc-gene-rsids.txt")
    composed_blob.upload_from_string(rows_string)