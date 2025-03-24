import psycopg2
from Bio import Entrez
from xml.etree import ElementTree as ET
import litellm
import os
import json
from typing import List, Dict
litellm.suppress_debug_info = True
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

def fetch_pubmed_papers(terms: Dict, start_year: int, end_year: int, max_results: int = 5) -> List[Dict]:
    """Search PubMed with combined MeSH and keyword strategy"""
    search_query = (
        "(" + " OR ".join([term + "[MeSH]" for term in terms['mesh']]) + ") " +
        "AND (" + " OR ".join(["\"" + kw + "\"[Title/Abstract]" for kw in terms['keywords']]) + ") " +
        "AND (" + str(start_year) + ":" + str(end_year) + "[PDAT])"
    )
    print("search query is")
    print(search_query)
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

def fetch_aact_trials(terms: Dict, start_year: int, end_year: int, max_trials: int = 5) -> List[Dict]:
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

def analyze_entity(text: str, entity_type: str) -> list:
    """Use LLM to extract specific entities from text"""
    prompt = f"""Extract {entity_type} from this text. Return as JSON dictionary.
    Entities should be concise (2-5 words max). If no entities are found for a category, return an empty list for that key.
    Types:
    1. Demographics (age groups, populations)
    2. Technologies (methods, tools)
    3. Problems (medical challenges addressed)

    Text: {text[:3000]}  # Truncate for token limits
    Example output for 'early detection in elderly': {{"{entity_type}": ["elderly patients"]}}"""

    response = litellm.completion(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.0
    )
    return json.loads(response.choices[0].message.content).get(entity_type, [])

def generate_contextual_insights(metrics, papers, trials) -> str:
    """Generate rich insights using metrics context and structured prompt"""
    # Build metrics summary
    metrics_context = "## Research Focus Areas\n"
    for category in ['demographics', 'technologies', 'problems']:
        metrics_context += (
            f"### {category.capitalize()}\n"
            f"- **Paper Focus**: {', '.join(set(metrics['papers'][category][:5]))}\n"
            f"- **Trial Focus**: {', '.join(set(metrics['trials'][category][:5]))}\n\n"
        )
    
    # Add content examples
    content_samples = (
        "## Sample Context\n"
        "**Papers**\n" + "\n".join(
            f"- {p['title'][:100]}..." for p in papers[:3]
        ) + "\n\n**Trials**\n" + "\n".join(
            f"- {t['brief_title'][:100]}: {t['outcome_or_title'][:200]}..." 
            for t in trials[:3]
        )
    )

    prompt = f"""Analyze these medical research components and provide structured insights:
    {metrics_context}
    {content_samples}
    Generate a report with these sections:

    1. **Emerging Patterns** - Novel paper concepts lacking trials
    2. **Translation Opportunities** - Ready-for-trial technologies
    3. **Demographic Alignment** - Population focus comparisons  
    4. **Innovation Pipeline** - Research-to-trial maturity spectrum
    5. **Commercial Landscape** - Industry participation trends

    For each section:
    - 3 concise bullet points
    - Specific examples from metrics/content
    - Strategic recommendations

    Format with markdown headers and bullet points."""
    response = litellm.completion(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.7,
        max_tokens=1500
    )
    return response.choices[0].message.content

