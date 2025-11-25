from Bio import Entrez

# --- CONFIGURATION ---
Entrez.email = "gararavi172@gmail.com"  # Must be valid!
# This is the exact string your agent is currently sending:
bad_query = "Patient: Adults with Type 2 Diabetes, Intervention: Metformin, Comparison: Lifestyle changes only, Outcome: Cardiovascular mortality"

def test_search(query):
    print(f"\nTesting Query: '{query}'")
    try:
        # Clean the query manually to test the fix
        clean_query = query.replace("Patient:", "").replace("Intervention:", "")\
                           .replace("Comparison:", "").replace("Outcome:", "")\
                           .replace(",", " ")
        
        print(f"Sending to PubMed as: '{clean_query.strip()}'")
        
        handle = Entrez.esearch(db="pubmed", term=clean_query, retmax=5)
        record = Entrez.read(handle)
        count = int(record["Count"])
        print(f"✅ Found {count} papers.")
        
        if count == 0:
            print("❌ ZERO RESULTS. The query is too specific or contains stopwords.")
            
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    test_search(bad_query)