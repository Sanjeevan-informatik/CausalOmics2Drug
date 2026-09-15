from typing import Literal
from pydantic import BaseModel, ConfigDict, Field, model_validator

LAYERS = ('genomics', 'epigenomics', 'bulk_rna', 'single_cell', 'spatial', 'proteomics', 'metabolomics')
Layer = Literal['genomics', 'epigenomics', 'bulk_rna', 'single_cell', 'spatial', 'proteomics', 'metabolomics']

class StrictModel(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False, str_strip_whitespace=True)

class Evidence(StrictModel):
    id: str = Field(min_length=1, max_length=120)
    gene: str = Field(pattern=r'^[A-Z][A-Z0-9.-]{0,39}$')
    layer: Layer
    disease: str = Field(min_length=1, max_length=120)
    tissue: str = Field(min_length=1, max_length=120)
    study: str = Field(min_length=1, max_length=120)
    source: str = Field(min_length=1, max_length=500)
    observation: str = Field(min_length=1, max_length=500)
    strength: float = Field(ge=0, le=1)
    effect: float | None = None
    unit: str = Field(default='context', max_length=60)
    q_value: float | None = Field(default=None, ge=0, le=1)
    cell_type: str = Field(default='all', min_length=1, max_length=120)
    kind: Literal['association', 'fine_mapping_coloc', 'perturbation', 'mendelian_randomization'] = 'association'
    action: Literal['inhibit', 'activate', 'unknown'] = 'unknown'
    synthetic: bool = False
    x: float | None = None
    y: float | None = None

    @model_validator(mode='after')
    def coordinates(self):
        if (self.x is None) != (self.y is None):
            raise ValueError('Spatial coordinates require both x and y')
        return self

class Pathway(StrictModel):
    id: str = Field(min_length=1, max_length=120)
    name: str = Field(min_length=1, max_length=200)
    genes: list[str] = Field(min_length=1, max_length=5000)
    source: str = Field(min_length=1, max_length=500)

class Dataset(StrictModel):
    name: str = Field(min_length=1, max_length=200)
    evidence: list[Evidence] = Field(min_length=1, max_length=20000)
    pathways: list[Pathway] = Field(default_factory=list, max_length=1000)
    universe: list[str] = Field(default_factory=list, max_length=60000)
    notes: list[str] = Field(default_factory=list, max_length=100)

    @model_validator(mode='after')
    def identifiers(self):
        ids = [r.id for r in self.evidence]
        if len(ids) != len(set(ids)):
            raise ValueError('Evidence ids must be unique')
        pids = [p.id for p in self.pathways]
        if len(pids) != len(set(pids)):
            raise ValueError('Pathway ids must be unique')
        if self.universe and not {r.gene for r in self.evidence} <= set(self.universe):
            raise ValueError('Measured universe must contain every evidence gene')
        return self

class AnalysisRequest(StrictModel):
    dataset: Dataset
    disease: str = Field(min_length=1)
    tissue: str = Field(min_length=1)
    cell_type: str = 'all'
    min_strength: float = Field(default=0.5, ge=0, le=1)
    max_q: float = Field(default=0.05, ge=0, le=1)
    weights: dict[Layer, float] = Field(default_factory=lambda: {k: 1.0 for k in LAYERS})

    @model_validator(mode='after')
    def valid_weights(self):
        import math
        if not self.weights or any(not math.isfinite(v) or v < 0 or v > 10 for v in self.weights.values()) or sum(self.weights.values()) <= 0:
            raise ValueError('Weights must be finite, between 0 and 10, with positive total')
        return self
