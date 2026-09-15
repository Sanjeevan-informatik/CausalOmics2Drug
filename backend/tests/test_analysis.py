import math
import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError
from backend.app.main import app, demo
from backend.app.analysis import analyze, bh_adjust
from backend.app.models import Evidence, Dataset, AnalysisRequest
from backend.app.importers import import_table

client = TestClient(app)

def row(**kw):
    data=dict(id='1',gene='EGFR',layer='bulk_rna',disease='D',tissue='T',study='study',source='fixture',observation='test',strength=0.8)
    return Evidence(**(data|kw))

def run(rows, **kw):
    return analyze(AnalysisRequest(dataset=Dataset(name='test',evidence=rows),disease='D',tissue='T',**kw))

def test_demo_has_seven_layers_but_no_causal_or_statistical_claims():
    d=demo(); r=analyze(AnalysisRequest(dataset=d,disease=d.evidence[0].disease,tissue=d.evidence[0].tissue))
    assert len({e.layer for e in d.evidence})==7
    assert r['targets'][0]['gene']=='EGFR'
    assert all(t['causal_score']==0 for t in r['targets'])
    assert all(p['q_value'] is None for p in r['pathways'])

def test_duplicate_observations_do_not_inflate_score():
    a=run([row()]); b=run([row(),row(id='2')])
    assert a['targets'][0]['score']==b['targets'][0]['score']
    assert a['targets'][0]['score']==round(80/7,2)

def test_missing_evidence_is_explicit():
    t=run([row()])['targets'][0]
    assert t['coverage']==1 and len(t['missing'])==6 and t['action']=='unknown'

def test_context_isolation():
    r=run([row(),row(id='2',gene='TP53',disease='other'),row(id='3',gene='MYC',tissue='other')])
    assert [t['gene'] for t in r['targets']]==['EGFR']

def test_cell_filter_retains_general_evidence():
    r=run([row(),row(id='2',gene='TP53',cell_type='B'),row(id='3',gene='MYC',cell_type='T')],cell_type='T')
    assert {t['gene'] for t in r['targets']}=={'EGFR','MYC'}

def test_conflicting_causal_directions():
    r=run([row(kind='perturbation',action='inhibit'),row(id='2',kind='fine_mapping_coloc',action='activate')])
    assert r['targets'][0]['action']=='conflicting'
    assert r['targets'][0]['causal_score']==round(160/3,2)

def test_failed_fdr_does_not_score():
    t=run([row(q_value=0.5,kind='perturbation',action='inhibit')])['targets'][0]
    assert t['score']==0 and t['causal_score']==0 and t['action']=='unknown'

@pytest.mark.parametrize('value',[float('nan'),float('inf'),-0.1,1.1])
def test_invalid_strength(value):
    with pytest.raises(ValidationError): row(strength=value)

def test_duplicate_ids_rejected():
    with pytest.raises(ValidationError): Dataset(name='x',evidence=[row(),row()])

def test_universe_must_cover_observed_genes():
    with pytest.raises(ValidationError): Dataset(name='x',evidence=[row()],universe=['TP53'])

def test_bh_known_vector():
    assert bh_adjust([0.01,0.04,0.03])==pytest.approx([0.03,0.04,0.04])

def test_pathway_enrichment_known_hypergeometric():
    d=Dataset(name='x',evidence=[row(q_value=0.01)],universe=['EGFR','TP53','MYC','BRCA1'],
      pathways=[dict(id='p',name='p',genes=['EGFR','TP53'],source='test')])
    p=analyze(AnalysisRequest(dataset=d,disease='D',tissue='T'))['pathways'][0]
    assert p['p_value']==pytest.approx(0.5) and p['q_value']==pytest.approx(0.5)

def test_import_source_format():
    d=import_table('gene,tumor_tpm,control_tpm\nEGFR,320,30\n','D','T','x')
    assert d.evidence[0].effect==pytest.approx(math.log2(321/31))
    assert d.evidence[0].q_value is None

@pytest.mark.parametrize('value',['nan','inf','-1'])
def test_bad_expression_rejected(value):
    with pytest.raises(ValueError): import_table(f'gene,tumor_tpm,control_tpm\nEGFR,{value},30\n','D','T','x')

def test_api_roundtrip_and_fingerprint():
    d=client.get('/api/demo').json()
    payload=dict(dataset=d,disease=d['evidence'][0]['disease'],tissue=d['evidence'][0]['tissue'])
    r=client.post('/api/analyze',json=payload)
    assert r.status_code==200
    assert r.json()['fingerprint']==client.post('/api/analyze',json=payload).json()['fingerprint']
    assert client.get('/api/health').json()['status']=='ok'

def test_api_rejects_bad_upload():
    assert client.post('/api/import',json=dict(text='a,b\n1,2',disease='D',tissue='T',name='x')).status_code==422

def test_zero_weights_rejected():
    with pytest.raises(ValidationError): run([row()],weights={'bulk_rna':0})

def test_json_validation_expands_optional_defaults():
    r=client.post('/api/validate',json={'name':'minimal','evidence':[row().model_dump(exclude_defaults=True)]})
    assert r.status_code==200
    assert r.json()['notes']==[] and r.json()['pathways']==[]
    assert r.json()['evidence'][0]['x'] is None
