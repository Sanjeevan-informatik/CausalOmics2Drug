"""Reproducible teaching dataset. No invented patient observations or p-values."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'https://github.com/Sanjeevan-informatik/cell-omics-explorer/blob/main/data/'
rows = []

def add(gene, layer, observation, strength, effect=None, unit='context', source=None, **extra):
    rows.append(dict(id=f'demo-{len(rows)+1}', gene=gene, layer=layer, disease='Illustrative EGFR-high tumor',
        tissue='Unspecified tumor (teaching)', study='CellOmics teaching example', source=source or 'Generated teaching extension; not measured',
        observation=observation, strength=strength, effect=effect, unit=unit, synthetic=True, **extra))

for r in csv.DictReader((ROOT/'data/cell_omics/example_expression.csv').open()):
    effect = math.log2((float(r['tumor_tpm'])+1)/(float(r['control_tpm'])+1))
    add(r['gene'], 'bulk_rna', f"Tumor {r['tumor_tpm']} TPM; control {r['control_tpm']} TPM", min(abs(effect)/4, 1),
        effect, 'log2 ratio (+1 TPM)', SOURCE+'example_expression.csv')
add('EGFR','genomics','Amplification CN approximately 6; no variant-to-gene causal analysis',0.75,6,'copy number',SOURCE+'example_multiomics.tsv')
add('EGFR','epigenomics','High promoter accessibility; qualitative context',0.5,source=SOURCE+'example_multiomics.tsv')
add('EGFR','proteomics','Protein abundance 2.8x control; pEGFR Y1068 3.4x control',0.65,2.8,'fold change',SOURCE+'example_multiomics.tsv')
for r in csv.DictReader((ROOT/'data/cell_omics/example_methylation.tsv').open(),delimiter='\t'):
    add(r['locus'].split('_')[0], 'epigenomics',r['locus']+' methylation; no matched control or regulatory inference',0.25,
        float(r['beta_value']),'beta value',SOURCE+'example_methylation.tsv')
# These extensions demonstrate alignment, not discoveries from the source repository.
add('EGFR','single_cell','Simulated tumor epithelial pseudobulk contrast',0.7,1.8,'log2 fold change',cell_type='Tumor epithelial')
for x,y,v in [(0,0,0.2),(1,0,0.4),(2,0,0.65),(0,1,0.35),(1,1,0.9),(2,1,0.8),(0,2,0.2),(1,2,0.5),(2,2,0.6)]:
    add('EGFR','spatial','Simulated spatial EGFR signal; coordinates in grid units',v,v,'relative signal',cell_type='Tumor epithelial',x=x,y=y)
add('LDHA','metabolomics','Lactate 2.2x control; LDHA mapping is a teaching extension, not evidence of LDHA activity',0.3,2.2,'fold change',
    source=SOURCE+'example_multiomics.tsv (lactate value only; gene mapping added)')
for gene in ['KRAS','BRAF','MAPK1','PIK3CA','AKT1','MTOR']:
    add(gene,'proteomics','Generated low-strength pathway context; no measurement in source',0.15)

data = dict(name='EGFR: from multi-omics observations to a testable hypothesis',evidence=rows,
    pathways=[dict(id='teaching-egfr-mapk',name='EGFR → MAPK / growth (teaching subset)',genes=['EGFR','KRAS','BRAF','MAPK1','MYC'],source='https://reactome.org/content/detail/R-HSA-177929'),
              dict(id='teaching-pi3k',name='PI3K–AKT / survival (teaching subset)',genes=['EGFR','PIK3CA','AKT1','MTOR'],source='https://reactome.org/content/detail/R-HSA-2219528')],
    universe=[],notes=['All observations are synthetic or teaching examples. Source data has no patient identifiers, disease subtype, replication or q-values.',
    'Strengths are explicitly hand-set illustration weights except expression: min(abs(log2((tumor TPM+1)/(control TPM+1)))/4,1).',
    'No causal evidence was supplied. The example intentionally ranks associations without claiming causality.',
    'The original six-cell immune-marker table is preserved under data/cell_omics; it has no EGFR or donor identifiers and is not joined to the tumor example.',
    'Pathways are small explanatory subsets, not complete Reactome gene sets. No enrichment p-values are calculated for the demo.'])
(ROOT/'data/demo.json').write_text(json.dumps(data,indent=2)+'\n')
# All columns for normalized uploads; the downloaded demo also documents optional fields.
fields=['id','gene','layer','disease','tissue','study','source','observation','strength','effect','unit','q_value','cell_type','kind','action','synthetic','x','y']
with (ROOT/'data/evidence_template.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    for r in rows[:5]:writer.writerow(r)
print(f'Wrote {len(rows)} teaching evidence records')
