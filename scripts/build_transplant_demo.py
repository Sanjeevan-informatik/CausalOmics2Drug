"""Generate fully synthetic transplant research fixtures; no patient data or predictor execution."""
import json,csv
from pathlib import Path
R=Path(__file__).resolve().parents[1]
out=R/'data/transplant';out.mkdir(exist_ok=True)
source='Synthetic teaching fixture; not a measurement or clinical result'
pairs=[dict(id='TX-001',diagnosis='AML · synthetic longitudinal case',donor_id='D-001',recipient_id='R-001',build='GRCh38',donor_hla=['HLA-A*02:01','HLA-B*07:02','HLA-C*07:02'],recipient_hla=['HLA-A*02:01','HLA-B*07:02','HLA-C*07:02'],germline_source='nonhematopoietic',source=source),dict(id='TX-002',diagnosis='MDS/AML · synthetic comparison',donor_id='D-002',recipient_id='R-002',build='GRCh38',donor_hla=['HLA-A*02:01','HLA-B*07:02'],recipient_hla=['HLA-A*02:01','HLA-B*07:02'],germline_source='unknown',source=source)]
samples=[]
for pair in pairs:
 for compartment,assay,days,total in [('blood','bulk_tcr',[0,30,90,180,365],10000),('marrow','single_cell',[30,90,180,365],500)]:
  for day in days:
   samples.append(dict(id=f'{pair["id"]}-{compartment}-{day}',pair_id=pair['id'],day=day,compartment=compartment,assay=assay,unit='cells' if assay=='single_cell' else 'templates',total=total,origin='donor' if day<365 else 'mixed',origin_method='Illustrative upstream genotype assignment',source=source))
candidates=[]
def cand(i,gene,**changes):
 d=dict(id=f'MIHA-{i:02}',pair_id='TX-001',gene=gene,chrom='19',pos=1000000+i*100,ref='G',alt='A',transcript=f'TEACHING_TRANSCRIPT_{i}',donor_gt='0/0',recipient_gt='0/1',donor_dp=55,recipient_dp=72,donor_gq=99,recipient_gq=99,target_allele='ALT',peptide_ref='ALRDFLLEA',peptide_alt='ALHDFLLEA',variant_index=3,hla='HLA-A*02:01',rank_el=0.35,affinity_nm=95,predictor='Synthetic rank fixture; no predictor executed',leukemia_tpm=22,healthy_hematopoietic_tpm=6,nonhematopoietic_max_tpm=0.2,immunopeptidomics='not_tested',source=source)
 d.update(changes);candidates.append(d)
cand(1,'ARHGAP45',immunopeptidomics='detected')
cand(2,'MYO1G',rank_el=0.9,leukemia_tpm=11,peptide_ref='YIGEVLVSA',peptide_alt='YIGEVLVSM',variant_index=9)
cand(3,'PTK2B',donor_gt='./.',donor_dp=4,donor_gq=None)
cand(4,'BCL2',donor_gt='1/1',recipient_gt='0/0')
cand(5,'UGT2B17',nonhematopoietic_max_tpm=28,rank_el=0.12)
cand(6,'DDX3X',donor_gt='0/1',recipient_gt='0/1')
cand(7,'RPS4X',donor_dp=8,donor_gq=12)
cand(8,'HIVEP1',rank_el=None,affinity_nm=None,leukemia_tpm=None,nonhematopoietic_max_tpm=None)
cand(9,'NISCH',rank_el=4.2,nonhematopoietic_max_tpm=2.5)
cand(10,'ARHGAP45',pair_id='TX-002')
clones=[]
seqs=[('CAVRDSNYQLIW','CASSLGQETQYF'),('CAVSDRGSTLGRLYF','CASSQAGTDTQYF'),('CAVNAGGTSYGKLTF','CASSFSGANTGELFF'),('CAVKDTDKLIF','CASSPGQGYEQYF'),(None,'CASSLGQETQYF'),('CAVRDSNYQLIW','CASSQDRDTQYF')]
for i,(a,b) in enumerate(seqs,1):
 clones.append(dict(id=f'TCR-{i:02}',pair_id='TX-001',tra_v='TRAV12-2' if a else '',tra_j='TRAJ33' if a else '',tra_cdr3=a,trb_v='TRBV7-9',trb_j='TRBJ2-7',trb_cdr3=b,pairing='beta_only' if a is None else 'ambiguous' if i==6 else 'paired',productive=True,source=source))
clones.append(dict(clones[0],id='TCR-07',pair_id='TX-002'))
observations=[]
counts=[[20,180,950,650,95],[10,25,100,310,440],[500,400,270,160,100],[0,15,70,35,10],[8,60,180,90,25],[2,15,25,20,5]]
for s in samples:
 cs=clones[:6] if s['pair_id']=='TX-001' else clones[6:]
 for i,c in enumerate(cs):
  dayidx=[0,30,90,180,365].index(s['day']);v=counts[i][dayidx]
  if s['compartment']=='marrow':v=round(v/20)
  if c['id']=='TCR-02' and s['day']==180 and s['compartment']=='marrow':continue
  observations.append(dict(sample_id=s['id'],clone_id=c['id'],count=v,phenotype='memory' if s['day']<=30 else 'effector' if s['day']<=90 else 'dysfunctional',source=source))
assays=[]
def assay(cand,clone,kind,result='positive',replicates=3,controls='adequate',**kw):
 a=dict(id=f'EXP-{len(assays)+1:03}',pair_id='TX-001',candidate_id=cand,clone_id=clone,sample_id='TX-001-blood-90',kind=kind,result=result,replicates=replicates,controls=controls,positive_events=None,total_events=None,control_positive=None,control_total=None,source=source,notes='Example only; no experiment was performed.')
 a.update(kw);assays.append(a)
assay('MIHA-01','TCR-01','tetramer',positive_events=820,total_events=10000,control_positive=12,control_total=10000)
for kind in ['reexpression','killing','hla_blocking']:assay('MIHA-01','TCR-01',kind)
assay('MIHA-01','TCR-01','normal_tissue',result='negative')
assay('MIHA-01','TCR-01','cross_reactivity',result='not_done',replicates=0,controls='not_recorded')
assay('MIHA-02','TCR-02','tetramer',positive_events=95,total_events=10000,control_positive=8,control_total=10000)
assay('MIHA-02','TCR-02','reexpression',replicates=1,controls='incomplete')
assay('MIHA-05','TCR-03','tetramer',positive_events=370,total_events=10000,control_positive=22,control_total=10000)
assay('MIHA-05','TCR-03','normal_tissue',result='positive')
assay('MIHA-01','TCR-05','tetramer',replicates=1)
escape=[]
for s in samples:
 if s['compartment']!='marrow':continue
 i=[30,90,180,365].index(s['day'])
 escape.append(dict(id='ESC-'+s['id'],sample_id=s['id'],blast_percent=[1.2,0.4,4.5,27][i],donor_chimerism=[99.2,99.8,94,73][i],hla_status=['retained','retained','reduced_expression','loss_supported'][i],hla='HLA-A*02:01',antigen_tpm=[20,18,9,2][i],somatic_clone='LEUK-A',vaf=[0.012,0.004,0.05,0.22][i],source=source,notes='HLA status is a supplied example annotation; VAF is not purity-corrected and cannot establish clone phylogeny.'))
p=dict(schema_version='1.0',name='AlloTrace · synthetic transplant study',synthetic=True,pairs=pairs,samples=samples,candidates=candidates,clones=clones,observations=observations,assays=assays,escape=escape,notes=['All genotypes, coordinates, peptide sequences, expression, binding ranks, receptors and experimental results are synthetic teaching fixtures, not real miHA candidates or patient data.','No HLA binding prediction or germline annotation tool has been run. Values exercise software behavior only.','The missing TCR-02 marrow observation at day +180 is deliberately missing, not a zero.'])
(out/'demo.json').write_text(json.dumps(p,indent=2)+'\n')
for name in ['pairs','samples','candidates','clones','observations','assays','escape']:
 with (out/f'{name}.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=p[name][0],lineterminator='\n');writer.writeheader()
  writer.writerows({k:';'.join(v) if isinstance(v,list) else str(v).lower() if isinstance(v,bool) else v for k,v in row.items()} for row in p[name])
vcf='##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD-001\tR-001\n19\t1000100\tSYNTHETIC\tG\tA\t99\tPASS\tGENE=ARHGAP45;TRANSCRIPT=TEACHING_TRANSCRIPT_1;PEP_REF=ALRDFLLEA;PEP_ALT=ALHDFLLEA;PEP_POS=3;HLA=HLA-A*02:01;RANK_EL=0.35;LEUK_TPM=22;NORMAL_TPM=0.2;SOURCE=Synthetic_fixture\tGT:DP:GQ\t0/0:55:99\t0/1:72:99\n'
(out/'annotated-example.vcf').write_text(vcf)
with (out/'airr-example.tsv').open('w',newline='') as f:
 fields=['sequence_id','cell_id','sample_id','locus','v_call','j_call','junction_aa','productive'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for i in range(4):
  for locus in ['TRA','TRB']:
   c=clones[i%2];stem=locus.lower();w.writerow(dict(sequence_id=f'CELL{i}-{locus}',cell_id=f'CELL{i}',sample_id='TX-001-marrow-90',locus=locus,v_call=c[stem+'_v'],j_call=c[stem+'_j'],junction_aa=c[stem+'_cdr3'],productive='T'))
print('Generated transplant teaching dataset and 9 input templates')
