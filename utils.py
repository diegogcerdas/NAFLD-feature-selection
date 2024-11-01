from Bio import Entrez

Entrez.email = "diego.gcerdas@gmail.com"


# Function to convert GenBank IDs to gene symbols
def genbank_to_gene_symbol(genbank_id):
    try:
        # Use Entrez efetch to get gene information
        handle = Entrez.efetch(
            db="nucleotide", id=genbank_id, rettype="gb", retmode="text"
        )
        record = handle.read()
        handle.close()
        # Extract the gene symbol from the GenBank record
        for line in record.splitlines():
            if "/gene=" in line:
                gene_symbol = line.split("=")[1].strip().replace('"', "")
                return gene_symbol
        else:
            return None  # If no gene symbol found
    except Exception as e:
        return None
