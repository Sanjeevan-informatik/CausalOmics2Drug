"""Transparent evidence prioritization, not a fitted causal discovery model."""
from collections import defaultdict
from hashlib import sha256
import json
from scipy.stats import hypergeom
from .models import AnalysisRequest, LAYERS

VERSION = '0.1.0'

def bh_adjust(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [1.0] * len(values)
    running = 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        running = min(running, values[index] * len(values) / rank)
        result[index] = running
    return result

def analyze(request: AnalysisRequest):
    rows = [r for r in request.dataset.evidence if r.disease == request.disease and r.tissue == request.tissue
            and (request.cell_type == 'all' or r.cell_type in ('all', request.cell_type))]
    grouped = defaultdict(list)
    for row in rows:
        grouped[row.gene].append(row)
    targets = []
    denominator = sum(request.weights.values())
    for gene, evidence in grouped.items():
        # Max per layer prevents repeated evidence and correlated studies inflating scores.
        qualified = [r for r in evidence if r.q_value is None or r.q_value <= request.max_q]
        layers = {layer: max((r.strength for r in qualified if r.layer == layer), default=0) for layer in LAYERS}
        observed = sorted({r.layer for r in evidence})
        score = sum(layers[k] * request.weights.get(k, 0) for k in LAYERS) / denominator
        causal = {kind: max((r.strength for r in qualified if r.kind == kind), default=0)
                  for kind in ('fine_mapping_coloc', 'perturbation', 'mendelian_randomization')}
        causal_score = sum(causal.values()) / 3
        actions = {r.action for r in qualified if r.kind != 'association' and r.strength >= request.min_strength and r.action != 'unknown'}
        action = next(iter(actions)) if len(actions) == 1 else 'conflicting' if actions else 'unknown'
        ncausal = sum(v >= request.min_strength and v > 0 for v in causal.values())
        tier = 'Convergent causal support' if ncausal >= 2 else 'Single causal evidence type' if ncausal else 'Association only'
        if all(r.synthetic for r in evidence):
            tier = 'Illustrative · ' + tier
        targets.append(dict(gene=gene, score=round(100 * score, 2), causal_score=round(100 * causal_score, 2),
            causal_components=causal, tier=tier, action=action, layers=layers, coverage=len(observed),
            missing=[k for k in LAYERS if k not in observed], evidence=[r.model_dump() for r in evidence],
            untested=sum(r.q_value is None for r in evidence), synthetic=any(r.synthetic for r in evidence)))
    targets.sort(key=lambda t: (-t['causal_score'], -t['score'], t['gene']))
    for i, target in enumerate(targets, 1):
        target['rank'] = i
    # Gene set ORA uses only externally supplied FDR-qualified evidence, never teaching strengths alone.
    selected = {r.gene for r in rows if r.q_value is not None and r.q_value <= request.max_q and r.strength >= request.min_strength}
    universe = set(request.dataset.universe)
    pathways = []
    for pathway in request.dataset.pathways:
        genes = set(pathway.genes)
        observed = genes & set(grouped)
        hits = genes & selected
        background = genes & universe
        p = float(hypergeom.sf(len(hits) - 1, len(universe), len(background), len(selected))) if universe and selected and background else None
        pathways.append(dict(id=pathway.id, name=pathway.name, source=pathway.source, measured_genes=sorted(observed),
                             hits=sorted(hits), p_value=p, q_value=None, universe_size=len(universe),
                             status='Exploratory ORA' if p is not None else 'Descriptive overlap; no statistical test'))
    tested = [p for p in pathways if p['p_value'] is not None]
    for p, q in zip(tested, bh_adjust([p['p_value'] for p in tested])):
        p['q_value'] = q
    warnings = ['Scores are heuristic evidence summaries, not probabilities of causality or treatment benefit.',
                'Missing layers reduce coverage and are not evidence against a target. Correlated evidence is not independent replication.',
                'Disease and tissue labels select evidence; they do not establish specificity against other diseases.']
    if any(r.synthetic for r in rows): warnings.append('Synthetic teaching observations are included. No patient-level conclusions can be drawn.')
    if any(r.q_value is None for r in rows): warnings.append('Some observations have no supplied q-value. They remain descriptive and cannot create enrichment hits.')
    if not rows: warnings.append('No evidence matches the selected context.')
    payload = request.model_dump(mode='json')
    return dict(version=VERSION, dataset=request.dataset.name, context=dict(disease=request.disease, tissue=request.tissue, cell_type=request.cell_type),
                fingerprint=sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(),
                parameters={k:v for k,v in payload.items() if k != 'dataset'}, targets=targets, pathways=pathways,
                selected_genes=sorted(selected), warnings=warnings, evidence_count=len(rows))
