from typing import Dict
import os
from workflows.logs.log import setup_logger
from trial_publication_gap_analyzer import get_search_terms, fetch_pubmed_papers, fetch_aact_trials, calculate_research_gaps, generate_insights
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
        # Generate search terms
        terms = get_search_terms(user_query)
        # Fetch data
        papers = fetch_pubmed_papers(terms, start_year, end_year)
        trials = fetch_aact_trials(terms, start_year, end_year)
        # Analyze gaps
        gaps = calculate_research_gaps(papers, trials, terms["keywords"])
        insights = generate_insights(papers, trials)
        # Generate output
        output = {
            "query": user_query,
            "timeframe": f"{start_year}-{end_year}",
            "total_papers": len(papers),
            "total_trials": len(trials),
            "insights" : insights
        }

 
