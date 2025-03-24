from typing import Dict
import os, json
from workflows.logs.log import setup_logger
from trial_publication_gap_analyzer import get_search_terms, fetch_pubmed_papers, fetch_aact_trials, analyze_entity, generate_contextual_insights
from dotenv import load_dotenv
load_dotenv()
class QueryProcessor:
    def __init__(self, api_key: str):
        """Initialize the QueryProcessor with OpenAI API key."""
        os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
        self.logger = setup_logger()
        self.logger.info("Initialized QueryProcessor")

    def process_query(self, query_text: str) -> Dict:
        # Hardcoded user input
        user_query = "early cancer detection via blood tests"
        start_year = 2002
        end_year = 2023
        OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
        os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY
        user_query = "early cancer detection via blood tests"
        # Generate search terms
        terms = get_search_terms(user_query)
        # Fetch data
        papers = fetch_pubmed_papers(terms, start_year, end_year, 5)
        trials = fetch_aact_trials(terms, start_year, end_year, 5)
        # Analyze gaps
        user_query = "early cancer detection via blood tests"
        metrics = {
            'papers': {'demographics': [], 'technologies': [], 'problems': []},
            'trials': {'demographics': [], 'technologies': [], 'problems': []}
        }
        # Process papers
        for p in papers:
            for category in metrics['papers']:
                metrics['papers'][category].extend(analyze_entity(p['abstract'], category))
        
        # Process trials
        for t in trials:
            text = f"{t['brief_title']} {t['outcome_or_title']}"
            for category in metrics['trials']:
                metrics['trials'][category].extend(analyze_entity(text, category))
        
        # Generate insights
        insights = generate_contextual_insights(metrics, papers, trials)
        # Final output
        output = {
            "query": user_query,
            "timeframe": f"{start_year}-{end_year}",
            "total_papers": len(papers),
            "total_trials": len(trials),
            "metrics_summary": {
                "papers": {k: list(set(v))[:5] for k,v in metrics['papers'].items()},
                "trials": {k: list(set(v))[:5] for k,v in metrics['trials'].items()}
            },
            "analysis": insights
        }

        return output

 
