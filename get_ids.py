#!/usr/bin/env python3

import requests

# Base URL for the modern UniProt Knowledgebase search endpoint
url = "https://rest.uniprot.org/uniprotkb/stream"

# Package advanced search parameters as a clean dictionary
query_parameters = {
    "query": 'protein_name:"cytochrome p450" AND reviewed:true',
    "format": "list"
}

print("Connecting to UniProt REST API...")

try:
    # Execute the GET request with the dynamic parameters
    response = requests.get(url, params=query_parameters)

    # Check if the connection succeeded
    response.raise_for_status()

    # Clean up the response text and split it into a list of individual lines
    uniprot_ids = [line.strip() for line in response.text.strip().split("\n") if line.strip()]
    
    # Fetching 1000 CYP450 IDs
    uniprot_ids = uniprot_ids[:1001]
    print(f"Successfully fetched {len(uniprot_ids)} unique protein IDs.")

    # Save the IDs to a text file for your downstream AlphaFold download loop
    with open("cytochrome_ids.txt", "w") as out_file:
        for up_id in uniprot_ids:
            out_file.write(f"{up_id}\n")

    print("-> Target IDs have been saved to 'cytochrome_ids.txt'")
    print(f"   First 5 samples: {uniprot_ids[:5]}")

except requests.exceptions.RequestException as e:
    print(f"An error occurred while querying the UniProt API: {e}")
