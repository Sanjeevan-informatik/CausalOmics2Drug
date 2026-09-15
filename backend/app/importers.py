"""Import normalized evidence or the inspected CellOmics expression format."""
import csv
import io
import math
from .models import Dataset, Evidence

def import_table(text: str, disease: str, tissue: str, name: str) -> Dataset:
    delimiter = '\t' if '\t' in text.partition('\n')[0] else ','
    reader = csv.DictReader(io.StringIO(text.lstrip('\ufeff')), delimiter=delimiter)
    rows = list(reader)
    if not rows or len(rows) > 20000:
        raise ValueError('Provide 1–20,000 evidence rows')
    if {'gene', 'tumor_tpm', 'control_tpm'} <= set(reader.fieldnames or []):
        evidence = []
        for i, row in enumerate(rows):
            tumor, control = float(row['tumor_tpm']), float(row['control_tpm'])
            if not all(math.isfinite(v) and v >= 0 for v in (tumor, control)):
                raise ValueError('TPM values must be finite and non-negative')
            effect = math.log2((tumor + 1) / (control + 1))
            evidence.append(Evidence(id=f'import-{i}', gene=row['gene'], layer='bulk_rna', disease=disease,
                tissue=tissue, study=name, source=name, observation=f'Tumor {tumor:g} TPM; control {control:g} TPM',
                effect=effect, unit='log2((TPM+1)/(control+1))', strength=min(abs(effect)/4, 1)))
        return Dataset(name=name, evidence=evidence, notes=['Imported descriptive TPM contrast; no replicates or inferential statistics. Strength = min(abs(log2 ratio with +1 pseudocount)/4, 1).'])
    evidence = []
    for row in rows:
        clean = {k:v for k,v in row.items() if v not in ('', None)}
        if None in clean:
            raise ValueError('A row contains more fields than the header')
        evidence.append(Evidence.model_validate(clean))
    return Dataset(name=name, evidence=evidence, notes=['Normalized user evidence; strengths and statistical results supplied by uploader.'])
