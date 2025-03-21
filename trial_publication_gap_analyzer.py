import psycopg2
from Bio import Entrez
from xml.etree import ElementTree as ET
import litellm
import os
import json
from typing import List, Dict

# Configuration
Entrez.email = "abhay.saini@grazitti.com"  # Required for PubMed API
AACT_CREDS = {
    "host": "aact-db.ctti-clinicaltrials.org",
    "database": "aact",
    "user": "pmadan13",
    "password": "Apple227kgi"
}

def get_search_terms(user_query: str) -> Dict:
    """Convert layperson's query to MeSH terms and keywords using LLM"""
    prompt = f"""Convert this medical query into:
    1. Top 3 MeSH terms (format: ["term1", "term2"...])
    2. 5-7 keywords/synonyms (format: ["kw1", "kw2"...])

    Query: {user_query}
    Example response for "early cancer detection via blood":
    {{"mesh": ["Early Detection of Cancer", "Liquid Biopsy", "Neoplasms"], 
     "keywords": ["liquid biopsy", "ctDNA", "circulating tumor cells", "blood marker", "early diagnosis"]}}"""
    
    response = litellm.completion(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3
    )
    return json.loads(response.choices[0].message.content)

def fetch_pubmed_papers(terms: Dict, start_year: int, end_year: int, max_results: int = 100) -> List[Dict]:
    """Search PubMed with combined MeSH and keyword strategy"""
    search_query = (
        "(" + " OR ".join([term + "[MeSH]" for term in terms['mesh']]) + ") " +
        "AND (" + " OR ".join(["\"" + kw + "\"[Title/Abstract]" for kw in terms['keywords']]) + ") " +
        "AND (" + str(start_year) + ":" + str(end_year) + "[PDAT])"
    )

    try:
        handle = Entrez.esearch(db="pubmed", term=search_query, retmax=max_results)
        record = Entrez.read(handle)
        pmids = record["IdList"]
        
        papers = []
        if pmids:
            handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml")
            tree = ET.parse(handle)
            for article in tree.findall(".//PubmedArticle"):
                title = article.find(".//ArticleTitle").text
                abstract = " ".join([t.text for t in article.findall(".//AbstractText") if t.text])
                papers.append({
                    "title": title,
                    "abstract": abstract,
                    "keywords": [m.find("DescriptorName").text 
                                for m in article.findall(".//MeshHeading")]
                })
        return papers
    
    except Exception as e:
        print(f"PubMed Error: {e}")
        return []

def fetch_aact_trials(terms: Dict, start_year: int, end_year: int, max_trials: int = 100) -> List[Dict]:
    """Fetch clinical trials matching terms in conditions or outcomes"""
    try:
        conn = psycopg2.connect(**AACT_CREDS)
        cursor = conn.cursor()
        query = f"""
        SELECT 
            s.nct_id, s.brief_title, s.phase, s.enrollment,
            s.start_date, 
            COALESCE(o.title, s.official_title) AS outcome_or_title,
            string_agg(DISTINCT sp.name, '; ') AS sponsors
        FROM studies s
        LEFT JOIN browse_conditions bc ON s.nct_id = bc.nct_id
        LEFT JOIN outcomes o ON s.nct_id = o.nct_id
        LEFT JOIN sponsors sp ON s.nct_id = sp.nct_id
        WHERE (
            bc.mesh_term ILIKE ANY(%s) 
            OR o.title ILIKE ANY(%s)
            OR s.brief_title ILIKE ANY(%s)
            OR s.official_title ILIKE ANY(%s)
        )

        AND s.overall_status NOT IN ('Withdrawn', 'Terminated')
        AND s.start_date BETWEEN '{start_year}-01-01' AND '{end_year}-12-31'
        GROUP BY s.nct_id, s.brief_title, s.phase, s.enrollment, 
                 s.start_date, o.title, s.official_title
        LIMIT {max_trials};
        """

        # Generate patterns for all search fields
        patterns = [f"%{term}%" for term in terms['mesh'] + terms['keywords']]
        cursor.execute(query, (patterns, patterns, patterns, patterns))
        columns = [desc[0] for desc in cursor.description]
        return [dict(zip(columns, row)) for row in cursor.fetchall()]
    
    except Exception as e:
        print(f"AACT Error: {e}")
        return []
    finally:
        if conn: conn.close()

def calculate_research_gaps(papers, trials, keywords):
    gaps = []
    for kw in keywords:
        paper_count = sum(
            1 for p in papers
            if kw.lower() in (p.get("abstract", "") or "").lower() 
            or (p.get("title") and kw.lower() in p["title"].lower())
        )
        trial_count = sum(
            1 for t in trials if kw.lower() in (
                (t.get("outcome_or_title") or "").lower() +
                (t.get("brief_title") or "").lower() +
                (t.get("sponsors") or "").lower()
            )
        )
        if paper_count > trial_count:
            gaps.append({
                "topic": kw,
                "papers": paper_count,
                "trials": trial_count,
                "gap_score": paper_count - trial_count
            })
    return sorted(gaps, key=lambda x: x["gap_score"], reverse=True)

if __name__ == "__main__":
    # Hardcoded user input
    user_query = "early cancer detection via blood tests"
    start_year = 2002
    end_year = 2023
    OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
    os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY

    # Generate search terms
    terms = get_search_terms(user_query)
    
    # Fetch data
    papers = fetch_pubmed_papers(terms, start_year, end_year)

    trials = fetch_aact_trials(terms, start_year, end_year)
    # Analyze gaps
    gaps = calculate_research_gaps(papers, trials, terms["keywords"])
    
    # Generate output
    output = {
        "query": user_query,
        "timeframe": f"{start_year}-{end_year}",
        "total_papers": len(papers),
        "total_trials": len(trials),
        "top_gaps": gaps[:5],
        "interpretation": (
            f"Research papers outnumber clinical trials {len(papers):,}:{len(trials):,} "
            f"({len(papers)/max(len(trials),1):.1f}x ratio)\n"
            "Largest opportunities in: " 
            + ", ".join([g["topic"] for g in gaps[:3]])
        )
    }
    
    print(json.dumps(output, indent=2))
