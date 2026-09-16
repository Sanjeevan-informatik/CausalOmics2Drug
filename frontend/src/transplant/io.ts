import Papa from 'papaparse';
import {Project,TableKey,projectSchema,candidateSchema,Candidate,Clone} from './schema';

export function parseRows(text:string){const r=Papa.parse<Record<string,string>>(text.replace(/^\uFEFF/,''),{header:true,skipEmptyLines:'greedy',delimiter:text.split('\n')[0].includes('\t')?'\t':','});if(r.errors.length)throw new Error(r.errors.map(e=>`Row ${e.row??'?'}: ${e.message}`).join('; '));if(!r.data.length)throw new Error('The file contains no records');return r.data;}
const numeric=new Set(['day','total','pos','donor_dp','recipient_dp','donor_gq','recipient_gq','variant_index','rank_el','affinity_nm','leukemia_tpm','healthy_hematopoietic_tpm','nonhematopoietic_max_tpm','count','replicates','positive_events','total_events','control_positive','control_total','blast_percent','donor_chimerism','antigen_tpm','vaf']);
export function replaceTable(p:Project,key:TableKey,text:string){
 const rows=parseRows(text).map(r=>Object.fromEntries(Object.entries(r).map(([k,v])=>{
  if(['donor_hla','recipient_hla'].includes(k))return [k,v.split(';').map(s=>s.trim()).filter(Boolean)];
  if(['productive'].includes(k)){if(!['true','false','T','F','1','0'].includes(v))throw new Error(`${k} requires true/false`);return [k,['true','T','1'].includes(v)];}
  if(numeric.has(k))return [k,v.trim()===''?null:Number(v)];
  if(['tra_cdr3','hla'].includes(k)&&!v)return [k,null];
  return [k,v];
 })));
 return projectSchema.parse({...p,[key]:rows,synthetic:false});
}
export function csvFor(rows:unknown[]){return Papa.unparse(rows.map(r=>Object.fromEntries(Object.entries(r as object).map(([k,v])=>[k,Array.isArray(v)?v.join(';'):v]))),{escapeFormulae:true});}
export function saveFile(name:string,content:string,type='application/json'){const url=URL.createObjectURL(new Blob([content],{type}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
export function saveJson(name:string,value:unknown){saveFile(name,JSON.stringify(value,null,2));}

export function importAnnotatedVcf(text:string,pairId:string,donor:string,recipient:string){
 const header=text.split(/\r?\n/).find(l=>l.startsWith('#CHROM\t'));if(!header)throw new Error('VCF needs a #CHROM header');
 const columns=header.split('\t'),di=columns.indexOf(donor),ri=columns.indexOf(recipient);
 if(di<9||ri<9||di===ri)throw new Error('Select two different, existing donor and recipient sample columns');
 if(!/##reference=.*(?:GRCh38|hg38)/i.test(text))throw new Error('Declare ##reference=GRCh38; mixed or unknown genome builds are not supported');
 const rows:Candidate[]=[],skipped:{line:number;reason:string}[]=[];
 const number=(value:string|undefined)=>value&&value!=='.'?Number(value):null;
 text.split(/\r?\n/).forEach((line,index)=>{
  if(!line||line.startsWith('#'))return;
  const c=line.split('\t');
  if(c.length!==columns.length)throw new Error(`VCF line ${index+1}: incorrect column count`);
  if(c[6]!=='PASS'){skipped.push({line:index+1,reason:'Variant filter is not PASS'});return;}
  if(!/^[ACGT]$/.test(c[3])||!/^[ACGT]$/.test(c[4])){skipped.push({line:index+1,reason:'Only normalized biallelic SNVs are supported'});return;}
  const info=Object.fromEntries(c[7].split(';').map(x=>{const i=x.indexOf('=');return i<0?[x,'']:[x.slice(0,i),decodeURIComponent(x.slice(i+1))];}));
  const required=['GENE','TRANSCRIPT','PEP_REF','PEP_ALT','PEP_POS','HLA'];
  if(required.some(k=>!info[k])){skipped.push({line:index+1,reason:'Missing missense peptide / HLA annotation; annotate upstream'});return;}
  const format=c[8].split(':'),call=(col:number)=>Object.fromEntries(format.map((k,j)=>[k,c[col].split(':')[j]]));
  const d=call(di),r=call(ri);let target='ALT';
  if(info.TARGET)target=info.TARGET;else if(r.GT&&d.GT&&r.GT.split(/[|/]/).includes('0')&&!d.GT.split(/[|/]/).includes('0')&&!d.GT.includes('.'))target='REF';
  rows.push(candidateSchema.parse({id:`VCF-${pairId}-${index+1}`,pair_id:pairId,gene:info.GENE,chrom:c[0].replace(/^chr/,''),pos:Number(c[1]),ref:c[3],alt:c[4],transcript:info.TRANSCRIPT,
   donor_gt:d.GT||'./.',recipient_gt:r.GT||'./.',donor_dp:number(d.DP),recipient_dp:number(r.DP),donor_gq:number(d.GQ),recipient_gq:number(r.GQ),target_allele:target,
   peptide_ref:info.PEP_REF,peptide_alt:info.PEP_ALT,variant_index:Number(info.PEP_POS),hla:info.HLA,rank_el:number(info.RANK_EL),affinity_nm:number(info.AFF_NM),predictor:info.PREDICTOR||'not supplied',
   leukemia_tpm:number(info.LEUK_TPM),healthy_hematopoietic_tpm:number(info.HEM_TPM),nonhematopoietic_max_tpm:number(info.NORMAL_TPM),immunopeptidomics:'not_tested',source:info.SOURCE||`Annotated VCF line ${index+1}`}));
 });
 if(!rows.length)throw new Error(`No supported records. ${skipped.length} skipped; ${skipped[0]?.reason||'empty input'}`);
 return {rows,skipped};
}

export function importAirr(p:Project,text:string){
 const rows=parseRows(text),required=['sequence_id','cell_id','sample_id','locus','v_call','j_call','junction_aa','productive'];
 if(required.some(k=>!(k in rows[0])))throw new Error('Single-cell AIRR subset requires '+required.join(', '));
 if(rows.length>100000)throw new Error('Maximum 100,000 AIRR rows');
 const cells=new Map<string,typeof rows>(),seen=new Set<string>();let excluded=0;
 for(const r of rows){
  if(!r.cell_id)throw new Error('Bulk AIRR has no supported pairing here; import normalized clone/count tables instead');
  const sample=p.samples.find(s=>s.id===r.sample_id);if(!sample||sample.assay!=='single_cell')throw new Error(`Unknown or non-single-cell sample ${r.sample_id}`);
  const seqKey=r.sample_id+'|'+r.sequence_id;if(seen.has(seqKey))throw new Error('Duplicate sequence_id within sample');seen.add(seqKey);
  if(!['T','F','true','false','1','0'].includes(r.productive))throw new Error('Invalid AIRR productive flag');
  if(!['T','true','1'].includes(r.productive)||!['TRA','TRB'].includes(r.locus)){excluded++;continue;}
  if(!/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(r.junction_aa))throw new Error('Invalid AIRR junction_aa');
  const k=r.sample_id+'|'+r.cell_id;cells.set(k,[...(cells.get(k)||[]),r]);
 }
 const clones:Clone[]=[],observations:Project['observations']=[],keys=new Map<string,string>();let orphan=0,ambiguous=0;
 for(const chains of cells.values()){
  const alpha=chains.filter(c=>c.locus==='TRA').sort((a,b)=>a.sequence_id.localeCompare(b.sequence_id)),beta=chains.filter(c=>c.locus==='TRB').sort((a,b)=>a.sequence_id.localeCompare(b.sequence_id));
  if(!beta.length){orphan++;continue;}
  const a=alpha[0],b=beta[0],sample=p.samples.find(s=>s.id===b.sample_id)!;
  const pairing=alpha.length>1||beta.length>1?'ambiguous':a?'paired':'beta_only';if(pairing==='ambiguous')ambiguous++;
  // Store a representative for review; ambiguous cells are excluded from the paired shortlist.
  const key=JSON.stringify([sample.pair_id,a?.v_call||'',a?.j_call||'',a?.junction_aa||null,b.v_call,b.j_call,b.junction_aa,pairing]);
  let cloneId=keys.get(key);if(!cloneId){cloneId=`AIRR-${keys.size+1}`;keys.set(key,cloneId);clones.push({id:cloneId,pair_id:sample.pair_id,tra_v:a?.v_call||'',tra_j:a?.j_call||'',tra_cdr3:a?.junction_aa||null,trb_v:b.v_call,trb_j:b.j_call,trb_cdr3:b.junction_aa,pairing,productive:true,source:'Imported AIRR subset; exact paired amino-acid + V/J definition'});}
  const obs=observations.find(o=>o.sample_id===sample.id&&o.clone_id===cloneId);if(obs)obs.count++;else observations.push({sample_id:sample.id,clone_id:cloneId,count:1,phenotype:'unknown',source:'Unique AIRR cell_id within sample'});
 }
 if(!clones.length)throw new Error('No productive TRB-containing cells');
 // Explicit replacement prevents old assays being silently assigned to newly numbered clones.
 const project=projectSchema.parse({...p,clones,observations,assays:[],synthetic:false});
 return {project,summary:{cells:cells.size,clones:clones.length,excluded,orphan,ambiguous,clearedAssays:p.assays.length}};
}
