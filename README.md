# 🧬 DeepBioDiv Prototype

A user-friendly eDNA Species Finder web app that identifies species from environmental DNA sequences using NCBI BLAST search.

![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/streamlit-v1.28+-red.svg)
![BioPython](https://img.shields.io/badge/biopython-v1.80+-green.svg)

## 🎯 What Does This App Do?

Upload a DNA sequence file, and the app will:
1. **Send your DNA** to NCBI's massive species database
2. **Find matches** by comparing your DNA to millions of known species
3. **Show results** with the top 5 most likely species matches
4. **Display similarity** percentages and confidence scores

Perfect for researchers, students, and anyone curious about identifying species from DNA samples!

## 🚀 Quick Start

### Installation
```bash
# Clone the repository
git clone https://github.com/yourusername/edna-species-finder.git
cd edna-species-finder

# Install required packages
pip install streamlit biopython pandas requests

# Run the app
streamlit run app.py
```

### Usage
1. **Upload** a FASTA file containing your DNA sequence
2. **Review** your sequence (length, composition, etc.)
3. **Click** "Find Species Matches" to start the search
4. **Wait** 1-2 minutes for NCBI to process your request
5. **View** your results with species names and similarity scores

## 📁 File Requirements

Your DNA file should be in **FASTA format** (.fasta, .fa, .fna):
```
>Sample_DNA_Sequence_001
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

## 🧬 Scientific Terms Explained

### **Environmental DNA (eDNA)**
Think of it like genetic fingerprints left behind by living things. When animals, plants, or microbes live in water or soil, they leave tiny bits of their DNA behind. Scientists can collect these "genetic breadcrumbs" to find out what species live in an area without actually seeing them!

### **FASTA Format**
A simple text format for storing DNA sequences:
- First line starts with `>` and has a description
- Following lines contain the actual DNA letters (A, T, G, C)
- Like a labeled recipe card for genetic information

### **NCBI BLAST**
**NCBI** = National Center for Biotechnology Information (a huge government database)
**BLAST** = Basic Local Alignment Search Tool (a DNA comparison program)

Think of BLAST like a super-powered Google search, but instead of searching websites, it searches through millions of DNA sequences to find matches to your sample.

### **DNA Base Pairs**
The building blocks of DNA:
- **A** (Adenine) - Always pairs with T
- **T** (Thymine) - Always pairs with A  
- **G** (Guanine) - Always pairs with C
- **C** (Cytosine) - Always pairs with G

Like letters in a 4-letter alphabet that spells out the instructions for life!

### **E-value**
A measure of how likely your DNA match happened by pure chance:
- **Small numbers** (like 0.0001) = Great match, very unlikely to be random
- **Big numbers** (like 0.5) = Poor match, could easily be random
- Think of it like a "confidence score" - smaller is better!

### **Percent Identity**
How similar your DNA is to the database sequence:
- **95-100%** = Excellent match, likely the same species
- **85-95%** = Good match, possibly related species  
- **Below 85%** = Distant relationship or different species entirely

## 🛠️ Technical Features

- **Real-time progress tracking** during NCBI searches
- **Automatic result parsing** from XML to readable tables
- **Interactive data visualization** with charts and progress bars
- **Error handling** for network issues and malformed files
- **Responsive design** that works on desktop and mobile

## 🔬 Use Cases

- **Biodiversity surveys** - Identify species in water or soil samples
- **Food authenticity** - Verify what's actually in processed foods
- **Conservation research** - Monitor endangered species populations
- **Educational projects** - Learn about DNA analysis and bioinformatics
- **Quality control** - Confirm lab sample identities

## 📊 Example Results

The app shows you:
| Species                               | DNA Similarity | Confidence      | Database Link  |
|---------------------------------------|----------------|-----------------|----------------|
| *Salmo trutta* (Brown trout)          | 98.5%          | E-value: 0.0    | [View Details] |
| *Oncorhynchus mykiss* (Rainbow trout) | 92.3%          | E-value: 2e-150 | [View Details] |

## ⚠️ Important Notes

- **Internet required** - App communicates with NCBI's online database
- **Processing time** - Searches typically take 1-3 minutes
- **File size limits** - Keep DNA sequences under 10,000 base pairs for best performance
- **Rate limits** - NCBI may limit frequent searches from the same IP address

## 🙋‍♀️ Support

Having trouble? Check out:
- [NCBI BLAST Documentation](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs)
- [BioPython Tutorial](https://biopython.org/DIST/docs/tutorial/Tutorial.html)
- [Streamlit Documentation](https://docs.streamlit.io/)
