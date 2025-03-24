# Clinical Trial Trend Detector

An AI-powered tool that analyzes the relationship between medical research papers and clinical trials to identify trends, gaps, and opportunities in medical research.

## Features

- **Smart Analysis**: Uses GPT-4 to analyze research papers and clinical trials
- **Gap Detection**: Identifies areas where research hasn't translated to trials
- **Trend Insights**: Provides comprehensive insights about emerging trends
- **Interactive Interface**: User-friendly web interface for easy analysis

## Tech Stack

- Backend: Python (FastAPI)
- Frontend: HTML, JavaScript, CSS
- AI: OpenAI GPT-4
- Databases: AACT Clinical Trials Database
- APIs: PubMed API

## Prerequisites

- Python 3.8+
- OpenAI API key
- AACT Database credentials
- PubMed API access

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/clinical-trial-trend-detector.git
cd clinical-trial-trend-detector
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Create a `.env` file in the root directory:
```env
OPENAI_API_KEY=your_openai_api_key
AACT_USER=your_aact_username
AACT_PASSWORD=your_aact_password
```

## Local Development

1. Start the FastAPI server:
```bash
uvicorn app:app --reload
```

2. Open `http://localhost:8000` in your browser

## Deployment on Render

1. Create a new Web Service on Render
2. Connect your GitHub repository
3. Configure the following:
   - Build Command: `pip install -r requirements.txt`
   - Start Command: `uvicorn app:app --host 0.0.0.0 --port $PORT`
4. Add environment variables:
   - `OPENAI_API_KEY`
   - `AACT_USER`
   - `AACT_PASSWORD`
5. Deploy

## API Documentation

The API has the following endpoints:

- `GET /`: Home page
- `GET /chat`: Analysis interface
- `POST /process-query`: Process analysis request

## Usage Example

1. Visit the home page
2. Click "Start Analyzing Trends"
3. Enter your research topic (e.g., "early cancer detection via blood tests")
4. Review the comprehensive analysis report

## Contributing

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- AACT Database for clinical trials data
- PubMed for research paper access
- OpenAI for GPT-4 API

## Support

For support, please open an issue in the GitHub repository or contact [your-email]. 