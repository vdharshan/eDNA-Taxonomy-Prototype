import streamlit as st  # This creates web apps
import requests  # This helps us talk to websites
import time  # This helps us wait between actions
import pandas as pd  # This helps us work with data tables
from Bio import SeqIO  # This reads DNA files
import io  # This helps us work with file data

# HELPER FUNCTIONS - These do the main backend work

def send_dna_to_ncbi_blast(dna_sequence_text):
    """
    This function sends a DNA sequence to NCBI's computer to find matches.
    
    Think of it like sending a letter to a DNA detective service.
    You give them a DNA sequence, and they give you back a tracking number.
    
    What it does:
    - Takes your DNA sequence as input
    - Sends it to NCBI's BLAST service on the internet
    - Gets back a tracking number (called RID) to check results later
    
    Returns:
    - The tracking number if successful
    - None if something went wrong
    """
    # This is NCBI's web address for BLAST searches
    ncbi_blast_website = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # These are the instructions we send to NCBI
    search_instructions = {
        'CMD': 'Put',  # This means "start a new search"
        'PROGRAM': 'blastn',  # This means "compare DNA sequences"
        'DATABASE': 'nt',  # This means "search in the main DNA database"
        'QUERY': dna_sequence_text  # This is our DNA sequence to search for
    }
    
    try:
        # Send our DNA sequence to NCBI (like mailing a letter)
        ncbi_response = requests.post(ncbi_blast_website, data=search_instructions, timeout=30)
        
        # Check if NCBI received our request properly
        ncbi_response.raise_for_status()
        
        # Look through NCBI's response to find our tracking number
        response_lines = ncbi_response.text.split('\n')
        for single_line in response_lines:
            if 'RID =' in single_line:
                # Found the tracking number! Extract it from the line
                tracking_number = single_line.strip().split('RID =')[1].strip()
                return tracking_number
        
        # If we get here, we couldn't find the tracking number
        st.error("Could not find tracking number in NCBI's response.")
        return None
        
    except requests.exceptions.RequestException as connection_error:
        # Something went wrong with the internet connection
        st.error(f"Error sending DNA to NCBI: {connection_error}")
        return None


def check_if_blast_search_is_done(tracking_number):
    """
    This function checks if NCBI has finished analyzing our DNA sequence.
    
    Think of it like calling a restaurant to ask "Is my food ready yet?"
    We keep calling until they say "Yes, it's ready!" or "Sorry, we messed up."
    
    What it does:
    - Uses the tracking number to check on our search
    - Keeps checking every 10 seconds until it's done
    - Stops after 5 minutes if it takes too long
    
    Returns:
    - True if the search finished successfully
    - False if it failed or took too long
    """
    # NCBI's web address for checking search status
    ncbi_status_website = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # Instructions to check our search status
    status_check_instructions = {
        'CMD': 'Get',  # This means "get information"
        'FORMAT_OBJECT': 'SearchInfo',  # This means "tell me the status"
        'RID': tracking_number  # This is our tracking number
    }
    
    # We'll try checking for 5 minutes maximum (30 attempts × 10 seconds each)
    maximum_check_attempts = 30
    current_attempt = 0
    
    while current_attempt < maximum_check_attempts:
        try:
            # Ask NCBI about our search status
            status_response = requests.get(ncbi_status_website, params=status_check_instructions, timeout=30)
            status_response.raise_for_status()
            
            # Check what NCBI told us
            if "Status=READY" in status_response.text:
                # Great! Our search is done
                return True
            elif "Status=FAILED" in status_response.text:
                # Oh no! Something went wrong with our search
                st.error("BLAST search failed on NCBI's server. Please try a different DNA sequence.")
                return False
            
            # Search is still running, so we wait 10 seconds and try again
            time.sleep(10)
            current_attempt += 1
            
        except requests.exceptions.RequestException as connection_error:
            # Something went wrong with checking the status
            st.error(f"Error checking search status: {connection_error}")
            return False
    
    # If we get here, we waited too long (5 minutes)
    st.error("Search took too long. Please try again.")
    return False


def download_blast_results_from_ncbi(tracking_number):
    """
    This function downloads the final results from NCBI after the search is done.
    
    Think of it like picking up your photos from a photo lab.
    You give them your receipt (tracking number) and they give you the photos (DNA matches).
    
    What it does:
    - Uses the tracking number to get our results
    - Downloads the results in a special computer format (XML)
    - Returns the results as text
    
    Returns:
    - The search results as text if successful
    - None if something went wrong
    """
    # NCBI's web address for getting results
    ncbi_results_website = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # Instructions to download our results
    results_download_instructions = {
        'CMD': 'Get',  # This means "get information" 
        'FORMAT_TYPE': 'XML',  # This means "give me results in XML format"
        'RID': tracking_number  # This is our tracking number
    }
    
    try:
        # Download our results from NCBI
        results_response = requests.get(ncbi_results_website, params=results_download_instructions, timeout=60)
        results_response.raise_for_status()
        
        # Return the results as text
        return results_response.text
        
    except requests.exceptions.RequestException as connection_error:
        # Something went wrong with downloading
        st.error(f"Error downloading results: {connection_error}")
        return None


def convert_xml_results_to_readable_table(xml_results_text):
    """
    This function converts NCBI's computer results into a human-readable table.
    
    Think of it like translating a technical manual into plain English.
    NCBI gives us results in XML (computer language), and we turn it into a nice table.
    
    What it does:
    - Reads the XML results from NCBI
    - Extracts information about the top 5 DNA matches
    - Creates a table with species names, match quality, and links
    
    Returns:
    - A pandas table with the top 5 species matches
    - An empty table if no matches were found
    """
    # Import the special tool for reading BLAST results
    from Bio.Blast import NCBIXML
    
    # This will hold our table data
    species_matches_data = []
    
    try:
        # Convert the text results into something we can read
        xml_file_like_object = io.StringIO(xml_results_text)
        blast_results_reader = NCBIXML.parse(xml_file_like_object)
        first_blast_record = next(blast_results_reader)
        
        # Look at each species match (we only want the top 5)
        for match_number, species_alignment in enumerate(first_blast_record.alignments):
            if match_number >= 5:  # Stop after 5 matches
                break
                
            # Each match can have multiple comparisons, we want the best one
            if species_alignment.hsps:
                best_comparison = species_alignment.hsps[0]
                
                # Calculate how similar the DNA sequences are (as a percentage)
                similarity_percentage = (best_comparison.identities / best_comparison.align_length) * 100
                
                # Get the species database ID (accession number)
                species_database_id = species_alignment.accession if hasattr(species_alignment, 'accession') else 'Unknown'
                
                # Create a link to see more info about this species
                species_info_link = f"https://www.ncbi.nlm.nih.gov/nuccore/{species_database_id}"
                
                # Clean up the species description
                species_description = species_alignment.title
                if ' ' in species_description:
                    # Remove the first part (usually just the ID) to make it more readable
                    species_description = species_description.split(' ', 1)[1]
                
                # Add this match to our data
                species_matches_data.append({
                    "Database ID": species_database_id,
                    "Species Description": species_description,
                    "E-value": best_comparison.expect,  # Lower numbers = better matches
                    "DNA Similarity %": round(similarity_percentage, 2),
                    "More Info Link": species_info_link
                })
                
    except StopIteration:
        # This means NCBI didn't find any matches
        st.warning("No DNA matches found in the results.")
        
    except Exception as parsing_error:
        # Something went wrong while reading the results
        st.error(f"Error reading the results: {parsing_error}")
        
    # Convert our data into a pandas table and return it
    return pd.DataFrame(species_matches_data)


# WEB APP SETUP - This creates the user interface

# Set up the web page
st.set_page_config(
    page_title="eDNA Species Finder",  # Name that appears in browser tab
    layout="wide",  # Use the full width of the screen
    initial_sidebar_state="collapsed"  # Start with sidebar hidden
)

# Create storage for information that persists as the user interacts with the app
# (Think of this like the app's memory)
if 'user_dna_sequence' not in st.session_state:
    st.session_state.user_dna_sequence = ""  # The DNA sequence from uploaded file
    
if 'species_results_table' not in st.session_state:
    st.session_state.species_results_table = pd.DataFrame()  # Table of species matches
    
if 'should_run_blast_search' not in st.session_state:
    st.session_state.should_run_blast_search = False  # Whether to start a new search


# WEB PAGE LAYOUT - This creates what the user sees

# Main title at the top of the page
st.title("🧬 DeepBioDiv Prototype")
st.markdown("eDNA Taxonomy Identifier: Upload a DNA sequence file and find out what species it might be from!")

# Sidebar with instructions (appears on the left)
st.sidebar.header("How This App Works")
st.sidebar.info(
    """
    **Step-by-step guide:**
    
    1. **Upload your DNA file** - Choose a FASTA file containing one DNA sequence
    
    2. **Review your sequence** - The app shows you what DNA it found in your file
    
    3. **Run the search** - Click the button to send your DNA to NCBI's database
    
    4. **Wait for results** - NCBI compares your DNA to millions of known species
    
    5. **See the matches** - View the top 5 species that match your DNA
    
    **What is eDNA?**
    Environmental DNA (eDNA) is genetic material found in water, soil, or air that can identify what species live in an area.

    **What is Taxonomy?**
    Taxonomy is the system scientists use to name, group, and organize living things.

    **What is a FASTA File?**
    A FASTA file is a simple text file that stores DNA, RNA, or protein sequences.
    
    **What is BLAST?**
    BLAST (Basic Local Alignment Search Tool) is a computer program that compares a DNA or protein sequence to a database and finds what it matches.
    """
)

# SECTION 1: FILE UPLOAD

# Create a box for the file upload section
file_upload_section = st.container(border=True)

with file_upload_section:
    st.header("Upload Your DNA Sequence File")
    
    # Create the file uploader widget
    uploaded_dna_file = st.file_uploader(
        "Choose a FASTA file (.fasta, .fa, .fna)",
        type=['fasta', 'fa', 'fna'],
        help="FASTA files contain DNA sequences in a specific text format",
        # When a new file is uploaded, reset the app's memory
        on_change=lambda: st.session_state.update(
            should_run_blast_search=False, 
            species_results_table=pd.DataFrame()
        )
    )

    # If user uploaded a file, try to read the DNA sequence from it
    if uploaded_dna_file:
        try:
            # Convert the uploaded file into text we can read
            file_text_content = io.StringIO(uploaded_dna_file.getvalue().decode("utf-8"))
            
            # Read all DNA sequences from the file
            dna_sequences_in_file = list(SeqIO.parse(file_text_content, "fasta"))
            
            if dna_sequences_in_file:
                # We found at least one DNA sequence! Use the first one
                first_dna_sequence = dna_sequences_in_file[0]
                st.session_state.user_dna_sequence = str(first_dna_sequence.seq)
                
                # Show success message
                st.success(f"✅ Successfully loaded DNA sequence: **{first_dna_sequence.id}**")
                st.info(f"📏 Your DNA sequence is {len(st.session_state.user_dna_sequence)} letters long")
            else:
                # The file was empty or not formatted correctly
                st.warning("⚠️ The FASTA file appears to be empty or incorrectly formatted.")
                st.session_state.user_dna_sequence = ""
                
        except Exception as file_reading_error:
            # Something went wrong reading the file
            st.error(f"❌ Error reading your FASTA file: {file_reading_error}")
            st.session_state.user_dna_sequence = ""


# SECTION 2: SEQUENCE REVIEW AND SEARCH

# Only show this section if we have a DNA sequence
if st.session_state.user_dna_sequence:
    
    sequence_review_section = st.container(border=True)
    
    with sequence_review_section:
        st.header("Review Your DNA Sequence")
        
        # Show some basic information about the sequence
        dna_length = len(st.session_state.user_dna_sequence)
        st.write(f"**Sequence length:** {dna_length:,} base pairs (DNA letters)")
        
        # Count each type of DNA letter
        a_count = st.session_state.user_dna_sequence.count('A')
        t_count = st.session_state.user_dna_sequence.count('T') 
        g_count = st.session_state.user_dna_sequence.count('G')
        c_count = st.session_state.user_dna_sequence.count('C')
        
        # Show the counts in columns
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("A's (Adenine)", a_count)
        with col2:
            st.metric("T's (Thymine)", t_count)
        with col3:
            st.metric("G's (Guanine)", g_count)
        with col4:
            st.metric("C's (Cytosine)", c_count)
        
        # Let user view the actual DNA sequence (hidden by default to save space)
        with st.expander("🔍 Click here to view your DNA sequence"):
            st.code(st.session_state.user_dna_sequence, language="text")
            st.caption("This is your DNA sequence made up of A, T, G, and C letters")
        
        # Big button to start the search
        if st.button("🚀 Find Species Matches", type="primary", use_container_width=True):
            st.session_state.should_run_blast_search = True
            st.session_state.species_results_table = pd.DataFrame()  # Clear old results

# SECTION 3: RUNNING THE BLAST SEARCH

# Only run the search if the user clicked the button
if st.session_state.should_run_blast_search:
    
    st.header("Searching for Species Matches")
    
    # STEP 1: Send DNA to NCBI
    with st.spinner("📤 Sending your DNA sequence to NCBI's database..."):
        tracking_number = send_dna_to_ncbi_blast(st.session_state.user_dna_sequence)
    
    # If we got a tracking number, continue with the search
    if tracking_number:
        st.success(f"✅ DNA sequence sent successfully! Tracking number: `{tracking_number}`")
        
        # STEP 2: Wait for NCBI to finish the search
        with st.spinner("⏳ NCBI is comparing your DNA to millions of species... This usually takes 1-2 minutes."):
            search_completed_successfully = check_if_blast_search_is_done(tracking_number)
        
        # STEP 3: If search completed, download the results
        if search_completed_successfully:
            st.success("🎉 Search completed! Downloading your results...")
            
            with st.spinner("📥 Getting your species matches..."):
                xml_results = download_blast_results_from_ncbi(tracking_number)
                
                if xml_results:
                    # Convert results to a readable table
                    st.session_state.species_results_table = convert_xml_results_to_readable_table(xml_results)
                    st.success("✅ Results downloaded and processed!")
                    
                    # Stop the search process
                    st.session_state.should_run_blast_search = False

# SECTION 4: DISPLAY RESULTS

# Only show results if we have some
if not st.session_state.species_results_table.empty:
    
    results_section = st.container(border=True)
    
    with results_section:
        st.header("📊 Your Species Match Results")
        
        # Show how many matches we found
        num_matches = len(st.session_state.species_results_table)
        st.subheader(f"🔍 Found {num_matches} potential species matches")
        
        # Explain what the results mean
        with st.expander("ℹ️ How to read these results"):
            st.markdown("""
            **Understanding your results:**
            
            - **Database ID**: NCBI's identification number for this species
            - **Species Description**: The name and description of the matching species
            - **E-value**: How likely this match happened by chance (smaller numbers = better matches)
            - **DNA Similarity %**: How similar your DNA is to this species (higher % = better match)
            - **More Info Link**: Click to see detailed information about this species on NCBI's website
            
            **What makes a good match?**
            - DNA Similarity should be above 95% for a confident species identification
            - E-value should be very small (like 0.0 or with many zeros)
            - Look at the species description to see if it makes sense for your sample
            """)
        
        # Show the results table
        st.dataframe(
            st.session_state.species_results_table,
            column_config={
                "More Info Link": st.column_config.LinkColumn(
                    "More Info", 
                    display_text="🔗 View Details"
                ),
                "DNA Similarity %": st.column_config.ProgressColumn(
                    "DNA Similarity %",
                    format="%.2f%%",
                    min_value=0,
                    max_value=100,
                ),
            },
            hide_index=True,  # Don't show row numbers
            use_container_width=True  # Use full width
        )
        
        # Show a chart if we have data
        if len(st.session_state.species_results_table) > 0:
            st.subheader("📈 DNA Similarity Comparison")
            
            # Create chart data
            chart_data = st.session_state.species_results_table.set_index("Database ID")[["DNA Similarity %"]]
            st.bar_chart(chart_data, use_container_width=True)
            
            # Show the best match prominently
            best_match = st.session_state.species_results_table.iloc[0]  # First row is best match
            
            st.success(f"""
            **🏆 Best Match:** {best_match['Species Description']}
            
            - **DNA Similarity:** {best_match['DNA Similarity %']}%
            - **Reliability (E-value):** {best_match['E-value']}
            """)
        
        # Option to download results
        if st.button("💾 Download Results as CSV", use_container_width=True):
            csv_data = st.session_state.species_results_table.to_csv(index=False)
            st.download_button(
                label="Click here to download",
                data=csv_data,
                file_name="dna_species_matches.csv",
                mime="text/csv"
            )