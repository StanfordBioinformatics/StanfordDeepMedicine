import base64
from google.cloud import bigquery


def main(event, context):
    """Triggered from a message on a Cloud Pub/Sub topic.
    Args:
         event (dict): Event payload.
         context (google.cloud.functions.Context): Metadata for the event.
    """
    
    print(f"This Function was triggered by messageId {context.event_id} published at {context.timestamp}.")

    if 'data' in event:
        data = base64.b64decode(event['data']).decode('utf-8')
    print(f"> Pubsub data: {data}.")
    print(f"> Context: {context}.")  

    #print(event['data'])
    #client = bigquery.Client()
    # Perform a query.
    QUERY = (
      'SELECT name FROM `bigquery-public-data.usa_names.usa_1910_2013` '
      'WHERE state = "TX" '
      'LIMIT 100')
#    query_job = client.query(QUERY)  # API request
#    rows = query_job.result()  # Waits for query to finish
    
#    for row in rows:
#      print(row.name)