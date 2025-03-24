from fastapi import FastAPI, HTTPException, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import os
from workflows.jobs.query_processor import QueryProcessor

app = FastAPI()

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")

# Initialize query processor
query_processor = QueryProcessor(api_key=os.getenv('OPENAI_API_KEY'))

class QueryRequest(BaseModel):
    query: str

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/analyze", response_class=HTMLResponse)
async def analyze(request: Request):
    return templates.TemplateResponse("chat.html", {"request": request})

@app.post("/process-query")
async def process_query(query_request: QueryRequest):
    if not query_request.query:
        raise HTTPException(status_code=400, detail="No query provided")
    try:
        response = query_processor.process_query(query_request.query)
        print("Response:", response)  # Debug print
        return response
    except Exception as e:
        print("Error:", str(e))  # Debug print
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000) 
