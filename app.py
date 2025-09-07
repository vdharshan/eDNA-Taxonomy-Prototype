import streamlit as st
import requests
import time
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
import io

# --- Configuration ---
st.set_page_config(
    page_title="eDNA Taxonomy Prototype",
    page_icon="🧬",
    layout="wide"
)

# --- Helper Functions ---

def parse_fasta(uploaded_file):
    """
    Parses an uploaded FASTA file and returns the first sequence.
    Returns None if the file is invalid or empty.
    """
    if uploaded_file is None:
        return None
    try:
        # To read the uploaded file, we need to decode it to a string
        stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        # Use Biopython to parse the FASTA file
        records = list(SeqIO.parse(stringio, "fasta"))
        if not records:
            st.error("The uploaded FASTA file is empty or invalid.")
            return None
        return records[0]
    except Exception as e:
        st.error(f"Error parsing FASTA file: {e}")
        return None

def run_blast(sequence):
    """
    Submits a sequence to NCBI BLAST, polls for results, and returns the top 5 hits.
    """
    base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # 1. Submit the BLAST job (CMD=Put)
    put_payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'DATABASE': 'nt',
        'QUERY': sequence.seq
    }
    try:
        with st.spinner("Submitting your sequence to NCBI BLAST..."):
            req = requests.post(base_url, data=put_payload)
            req.raise_for_status() # Raise an exception for bad status codes
        
        # Extract the Request ID (RID)
        rid_search = [line for line in req.text.split('\n') if 'RID =' in line]
        if not rid_search:
            st.error("Could not find RID in BLAST submission response.")
            return None
        rid = rid_search[0].split('=')[1].strip()
        st.info(f"BLAST job submitted successfully. Request ID (RID): **{rid}**")

        # 2. Poll for job status until it's ready (CMD=Get)
        with st.spinner(f"Waiting for BLAST results (RID: {rid}). This may take a moment..."):
            while True:
                time.sleep(5) # Wait 5 seconds between checks
                get_payload = {
                    'CMD': 'Get',
                    'RID': rid,
                    'FORMAT_OBJECT': 'SearchInfo'
                }
                status_req = requests.post(base_url, data=get_payload)
                status_req.raise_for_status()
                if 'Status=READY' in status_req.text:
                    break
        st.success("BLAST search complete! Fetching results...")
        
        # 3. Fetch the results in XML format
        with st.spinner("Parsing results..."):
            fetch_payload = {
                'CMD': 'Get',
                'RID': rid,
                'FORMAT_TYPE': 'XML'
            }
            results_req = requests.post(base_url, data=fetch_payload)
            results_req.raise_for_status()
            
            # 4. Parse the XML results
            blast_record = NCBIXML.read(io.StringIO(results_req.text))
            
            hits = []
            for i, alignment in enumerate(blast_record.alignments):
                if i >= 5: # Limit to top 5 hits
                    break
                hsp = alignment.hsps[0] # High-scoring Pair
                percent_identity = (hsp.identities / hsp.align_length) * 100
                hit_data = {
                    'Accession': alignment.accession,
                    'Title': alignment.title,
                    'E-value': f"{hsp.expect:.2e}",
                    'Percent Identity': f"{percent_identity:.2f}%"
                }
                hits.append(hit_data)
        
        if not hits:
            st.warning("No significant hits were found for your sequence.")
            return None
            
        return pd.DataFrame(hits)

    except requests.exceptions.RequestException as e:
        st.error(f"A network error occurred: {e}")
        return None
    except Exception as e:
        st.error(f"An unexpected error occurred during the BLAST process: {e}")
        return None

# --- Streamlit UI ---

st.title("🧬 eDNA Taxonomy Prototype")
st.write(
    "Upload a FASTA file containing a single DNA sequence to identify its taxonomy "
    "using the NCBI BLAST API."
)

# File uploader
uploaded_file = st.file_uploader(
    "Choose a FASTA file (.fasta)",
    type=['fasta', 'fa', 'fna']
)

if uploaded_file is not None:
    # Parse the sequence and store it in session state
    sequence_record = parse_fasta(uploaded_file)
    
    if sequence_record:
        st.session_state.sequence_record = sequence_record
        
        st.subheader("Sequence Information")
        st.text_area(
            "First Sequence in File",
            f">{sequence_record.id} {sequence_record.description}\n{sequence_record.seq}",
            height=150
        )

# Run BLAST button - only active if a sequence is loaded
if 'sequence_record' in st.session_state:
    if st.button("Run BLAST", type="primary"):
        results_df = run_blast(st.session_state.sequence_record)
        
        if results_df is not None:
            st.session_state.results_df = results_df

# Display results if they exist in session state
if 'results_df' in st.session_state:
    results_df = st.session_state.results_df
    st.subheader("Top 5 BLAST Hits")
    st.table(results_df)

    st.subheader("Percent Identity of Top Hits")
    
    # Prepare data for charting (convert Percent Identity back to float)
    chart_df = results_df[['Accession', 'Percent Identity']].copy()
    chart_df['Percent Identity'] = chart_df['Percent Identity'].str.replace('%', '').astype(float)
    chart_df = chart_df.set_index('Accession')
    
    st.bar_chart(chart_df)
