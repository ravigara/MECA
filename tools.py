import re
from Bio import Entrez
from langchain_community.tools import DuckDuckGoSearchResults

# 1. SETUP EMAIL (CRITICAL)
Entrez.email = "gararavi172@gmail.com" 

def search_pubmed(query, max_results=5):
    """
    Smart PubMed Query Engine.
    1. Tries full PICO search.
    2. If 0 results, falls back to "Intervention + Outcome" broad search.
    """
    
    # --- HELPER: Parse PICO from the string ---
    # We extract the specific inputs to build a better fallback query
    i_match = re.search(r"Intervention: (.*?)(,|$| Comparison:)", query, re.IGNORECASE)
    o_match = re.search(r"Outcome: (.*?)(,|$)", query, re.IGNORECASE)
    
    intervention = i_match.group(1).strip() if i_match else ""
    outcome = o_match.group(1).strip() if o_match else ""

    # --- STRATEGY 1: Full Query Cleaned ---
    # Remove labels and commas
    clean_query = query.replace("Patient:", "").replace("Intervention:", "")\
                       .replace("Comparison:", "").replace("Outcome:", "")\
                       .replace(",", " ")
    
    # Run Search 1
    results = run_entrez_search(clean_query, max_results)
    
    # --- STRATEGY 2: Fallback (Broad Search) ---
    if "No peer-reviewed papers" in results and intervention and outcome:
        print(f"⚠️ Strict search failed. Retrying with broad search: {intervention} AND {outcome}...")
        broad_query = f"{intervention} AND {outcome}"
        results = run_entrez_search(broad_query, max_results)
        
    return results

def run_entrez_search(search_term, max_results):
    """
    Executes the actual API call to NCBI.
    """
    try:
        # Search for IDs
        # We append (hasabstract[text]) to ensure quality
        term_with_filter = f"{search_term} AND (hasabstract[text])"
        
        handle = Entrez.esearch(db="pubmed", term=term_with_filter, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list:
            return "No peer-reviewed papers found on PubMed for this query."

        # Fetch details
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        papers_raw = handle.read()
        handle.close()
        
        # Parse output
        summary_list = []
        papers = papers_raw.split("\n\n")
        
        for i, paper_text in enumerate(papers):
            if len(paper_text) < 50: continue
            
            # Robust Regex Extraction
            title_search = re.search(r'TI  - (.*)', paper_text)
            abstract_search = re.search(r'AB  - (.*)', paper_text, re.DOTALL)
            source_search = re.search(r'SO  - (.*)', paper_text)
            pmid_search = re.search(r'PMID- (.*)', paper_text)
            
            title = title_search.group(1) if title_search else "No Title"
            # Clean newlines in abstract
            abstract = abstract_search.group(1).replace('\n      ', ' ') if abstract_search else "No Abstract"
            source = source_search.group(1) if source_search else "Unknown Journal"
            pmid = pmid_search.group(1).strip() if pmid_search else ""
            
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "#"
            
            entry = (
                f"PAPER {i+1}:\n"
                f"TITLE: {title}\n"
                f"JOURNAL: {source}\n"
                f"LINK: {link}\n"
                f"ABSTRACT: {abstract[:2500]}\n" # Increased limit slightly
                f"--------------------------------------------------\n"
            )
            summary_list.append(entry)
            
        return "\n".join(summary_list) if summary_list else "No papers found."

    except Exception as e:
        return f"PubMed API Error: {str(e)}"

def check_for_bias_keywords(text_content):
    """
    Deterministic Clinical Rules Engine.
    """
    flags = []
    text_lower = text_content.lower()
    
    if "meta-analysis" in text_lower or "systematic review" in text_lower:
        pass 
    elif "observational" in text_lower or "cohort" in text_lower:
        flags.append("Notice: Observational Study (Correlation ≠ Causation)")
    elif "case report" in text_lower:
        flags.append("Warning: Low Level of Evidence (Case Report)")
        
    n_matches = re.findall(r'\bn\s*=\s*(\d+)', text_lower)
    if n_matches:
        small_n = [int(x) for x in n_matches if int(x) < 30]
        if small_n:
            flags.append(f"Critical: Small sample size detected (n={small_n[0]})")

    if "mice" in text_lower or "rats" in text_lower:
        flags.append("Critical: Animal Study (Pre-clinical data)")

    if "funded by" in text_lower:
        flags.append("Note: Funding source disclosed - check for industry bias.")

    if not flags:
        return "No specific keyword flags detected."
    
    return "; ".join(flags)