from typing import Dict
import os
from workflows.logs.log import setup_logger
from trial_publication_gap_analyzer import get_search_terms, fetch_pubmed_papers, fetch_aact_trials, calculate_research_gaps

class QueryProcessor:
    def __init__(self, api_key: str):
        """Initialize the QueryProcessor with OpenAI API key."""
        os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
        self.logger = setup_logger()
        self.logger.info("Initialized QueryProcessor")

    def process_query(self, query_text: str) -> Dict:
        """Process the research gap analysis query and return response."""
        self.logger.info(f"Processing query: {query_text}")
        
        # Set time range (you can make these configurable if needed)
        start_year = 2002
        end_year = 2023
        # Generate search terms
        terms = get_search_terms(query_text)
        # Fetch data
        papers = fetch_pubmed_papers(terms, start_year, end_year)
        trials = fetch_aact_trials(terms, start_year, end_year)
        # Analyze gaps
        gaps = calculate_research_gaps(papers, trials, terms["keywords"])
        gaps_string = ""
        #for gap in gaps[:5]:
        #    gaps_string += "Topic " + gap['topic'] + " "
        #    gaps_string += "Research Papers " + str(gap['papers']) + " "
        #    gaps_string += "Clinical Trials " + str(gap['trials']) + " "
        #    gaps_string += "Research Gap Score " + str(gap['gap_score']) + " "

        # Create output without any indentation or template literals
        #output = "Research Analysis Summary "
        #output += "Query " + query_text + " "
        #output += "Time Period " + str(start_year) + " to " + str(end_year) + " "
        #output += "Publication Statistics "
        #output += "Total Research Papers " + str(len(papers)) + " "
        #output += "Total Clinical Trials " + str(len(trials)) + " "
        #output += "Paper to Trial Ratio " + str(round(len(papers) / max(len(trials), 1), 1)) + " to 1 "
        #output += "Research Gap Analysis "
        #output += gaps_string
        #output += "Key Findings "
        #output += "Research papers outnumber clinical trials " + str(len(papers)) + " to " + str(len(trials)) + " with a ratio of " + str(round(len(papers) / max(len(trials), 1), 1)) + " times "
        #output += "The largest research opportunities are in " + " ".join([g['topic'] for g in gaps[:3]]) + " "
        #output += "This analysis suggests areas where additional clinical trials might be beneficial to bridge the gap between research and clinical implementation "
        print('433')
        papers_to_trials = str(len(papers)) + "to" + str(len(trials))
        ratio = str(len(papers) / len(trials))
        print('fdf')
        #print(output)
        self.logger.info("Successfully processed query")
        return {
            "query": query_text,
            "time_period": f"{start_year} to {end_year}",
            "total_papers": len(papers),
            "total_trials": len(trials),
            "paper_to_trial_ratio": round(len(papers) / max(len(trials), 1), 1),
            "top_research_gaps": [gap['topic'] for gap in gaps[:5]],
            "analysis": f"Research papers outnumber clinical trials with a ratio of {ratio}"
        }


 