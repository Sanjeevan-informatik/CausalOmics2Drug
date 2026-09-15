from pathlib import Path
import json
from fastapi import FastAPI, HTTPException, Request
from fastapi.staticfiles import StaticFiles
from pydantic import ValidationError
from .analysis import analyze
from .models import AnalysisRequest, Dataset, StrictModel
from .importers import import_table

ROOT = Path(__file__).resolve().parents[2]
app = FastAPI(title='CausalOmics2Drug', version='0.1.0', description='Explainable research target prioritization')

@app.middleware('http')
async def limit_body(request: Request, call_next):
    if request.method == 'POST':
        size = 0
        chunks = []
        async for chunk in request.stream():
            size += len(chunk)
            if size > 8_000_000:
                from fastapi.responses import JSONResponse
                return JSONResponse({'detail':'Request exceeds 8 MB limit'}, status_code=413)
            chunks.append(chunk)
        request._body = b''.join(chunks)
    return await call_next(request)

@app.get('/api/health')
def health():
    return {'status':'ok', 'version':'0.1.0'}

@app.get('/api/demo', response_model=Dataset)
def demo():
    return Dataset.model_validate_json((ROOT / 'data/demo.json').read_text())

@app.post('/api/analyze')
def run_analysis(request: AnalysisRequest):
    return analyze(request)

@app.post('/api/validate', response_model=Dataset)
def validate_dataset(dataset: Dataset):
    return dataset

class ImportRequest(StrictModel):
    text: str
    disease: str
    tissue: str
    name: str

@app.post('/api/import', response_model=Dataset)
def upload(request: ImportRequest):
    try:
        return import_table(request.text, request.disease, request.tissue, request.name)
    except (ValueError, TypeError, ValidationError) as exc:
        raise HTTPException(422, detail=str(exc)) from exc

dist = ROOT / 'frontend/dist'
if dist.exists():
    app.mount('/', StaticFiles(directory=dist, html=True), name='frontend')
