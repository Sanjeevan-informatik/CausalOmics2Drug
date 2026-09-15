/** Browser deployment of the reference Python algorithm. Tested against SciPy. */
import {z} from 'zod';
import Papa from 'papaparse';
import demo from '../../data/demo.json';
import {layers} from './types';

const text=(max:number)=>z.string().trim().min(1).max(max);
const finite=z.number().finite();
const strength=finite.min(0).max(1);
const evidence=z.object({
 id:text(120),gene:z.string().trim().regex(/^[A-Z][A-Z0-9.-]{0,39}$/),layer:z.enum(layers),
 disease:text(120),tissue:text(120),study:text(120),source:text(500),observation:text(500),strength,
 effect:finite.nullable().default(null),unit:z.string().trim().max(60).default('context'),
 q_value:strength.nullable().default(null),cell_type:text(120).default('all'),
 kind:z.enum(['association','fine_mapping_coloc','perturbation','mendelian_randomization']).default('association'),
 action:z.enum(['inhibit','activate','unknown']).default('unknown'),synthetic:z.boolean().default(false),
 x:finite.nullable().default(null),y:finite.nullable().default(null),
}).strict().refine(r=>(r.x===null)===(r.y===null),'Spatial coordinates require both x and y');
const dataset=z.object({name:text(200),evidence:z.array(evidence).min(1).max(20000),
 pathways:z.array(z.object({id:text(120),name:text(200),genes:z.array(z.string()).min(1).max(5000),source:text(500)}).strict()).max(1000).default([]),
 universe:z.array(z.string()).max(60000).default([]),notes:z.array(z.string()).max(100).default([]),
}).strict().superRefine((d,ctx)=>{
 if(new Set(d.evidence.map(r=>r.id)).size!==d.evidence.length)ctx.addIssue({code:'custom',message:'Evidence ids must be unique'});
 if(new Set(d.pathways.map(p=>p.id)).size!==d.pathways.length)ctx.addIssue({code:'custom',message:'Pathway ids must be unique'});
 const universe=new Set(d.universe);
 if(universe.size&&d.evidence.some(r=>!universe.has(r.gene)))ctx.addIssue({code:'custom',message:'Measured universe must contain every evidence gene'});
});
const requestSchema=z.object({dataset,disease:text(120),tissue:text(120),cell_type:z.string().default('all'),
 min_strength:strength.default(0.5),max_q:strength.default(0.05),
 weights:z.record(z.string(),finite.min(0).max(10)).default(Object.fromEntries(layers.map(l=>[l,1])))
}).strict().refine(r=>Object.keys(r.weights).every(k=>layers.includes(k as typeof layers[number]))&&Object.values(r.weights).reduce((a,b)=>a+b,0)>0,'Weights require known layers and a positive total');

export function validateDataset(input:unknown){return dataset.parse(input);}
export function bhAdjust(values:number[]){
 const order=values.map((_,i)=>i).sort((a,b)=>values[a]-values[b]);const result=values.map(()=>1);let running=1;
 for(let rank=order.length;rank>0;rank--){const i=order[rank-1];running=Math.min(running,values[i]*order.length/rank);result[i]=running;}
 return result;
}
// Log-space summation avoids factorial overflow and tiny-probability cancellation.
export function hypergeomTail(hits:number,N:number,K:number,n:number){
 const lower=Math.max(0,n-(N-K)),upper=Math.min(K,n);
 if(hits<=lower)return 1;if(hits>upper)return 0;
 const logFact=new Float64Array(N+1);for(let i=2;i<=N;i++)logFact[i]=logFact[i-1]+Math.log(i);
 const choose=(a:number,b:number)=>b<0||b>a?-Infinity:logFact[a]-logFact[b]-logFact[a-b];
 const base=choose(N,n);let max=-Infinity;const logs=[];
 for(let k=hits;k<=upper;k++){const v=choose(K,k)+choose(N-K,n-k)-base;logs.push(v);max=Math.max(max,v);}
 return Math.min(1,Math.exp(max)*logs.reduce((a,v)=>a+Math.exp(v-max),0));
}
const sorted=(items:Iterable<string>)=>[...new Set(items)].sort();
const round=(v:number)=>Math.round((v+Number.EPSILON)*100)/100;
function canonical(value:unknown):string {
 if(Array.isArray(value))return '['+value.map(canonical).join(',')+']';
 if(value!==null&&typeof value==='object')return '{'+Object.entries(value).sort(([a],[b])=>a<b?-1:a>b?1:0).map(([k,v])=>JSON.stringify(k)+':'+canonical(v)).join(',')+'}';
 return JSON.stringify(value);
}
export async function analyze(input:unknown){
 const r=requestSchema.parse(input);
 const rows=r.dataset.evidence.filter(e=>e.disease===r.disease&&e.tissue===r.tissue&&(r.cell_type==='all'||e.cell_type==='all'||e.cell_type===r.cell_type));
 const grouped=new Map<string,typeof rows>();for(const e of rows)grouped.set(e.gene,[...(grouped.get(e.gene)||[]),e]);
 const targets=[...grouped].map(([gene,records])=>{
  const qualified=records.filter(e=>e.q_value===null||e.q_value<=r.max_q);
  const layerValues=Object.fromEntries(layers.map(l=>[l,Math.max(0,...qualified.filter(e=>e.layer===l).map(e=>e.strength))]));
  const observed=new Set(records.map(e=>e.layer));
  const causal=Object.fromEntries(['fine_mapping_coloc','perturbation','mendelian_randomization'].map(k=>[k,Math.max(0,...qualified.filter(e=>e.kind===k).map(e=>e.strength))]));
  const causalScore=Object.values(causal).reduce((a,b)=>a+b,0)/3;
  const actions=new Set(qualified.filter(e=>e.kind!=='association'&&e.strength>=r.min_strength&&e.action!=='unknown').map(e=>e.action));
  const n=Object.values(causal).filter(v=>v>=r.min_strength&&v>0).length;
  let tier=n>=2?'Convergent causal support':n?'Single causal evidence type':'Association only';
  if(records.every(e=>e.synthetic))tier='Illustrative · '+tier;
  return {gene,rank:0,score:round(100*layers.reduce((a,l)=>a+layerValues[l]*(r.weights[l]||0),0)/Object.values(r.weights).reduce((a,b)=>a+b,0)),
   causal_score:round(causalScore*100),causal_components:causal,tier,action:actions.size===1?[...actions][0]:actions.size?'conflicting':'unknown',
   layers:layerValues,coverage:observed.size,missing:layers.filter(l=>!observed.has(l)),evidence:records,untested:records.filter(e=>e.q_value===null).length,synthetic:records.some(e=>e.synthetic)};
 });
 targets.sort((a,b)=>b.causal_score-a.causal_score||b.score-a.score||(a.gene<b.gene?-1:a.gene>b.gene?1:0));targets.forEach((t,i)=>t.rank=i+1);
 const selected=new Set(rows.filter(e=>e.q_value!==null&&e.q_value<=r.max_q&&e.strength>=r.min_strength).map(e=>e.gene));
 const universe=new Set(r.dataset.universe);
 const pathways=r.dataset.pathways.map(p=>{
  const genes=sorted(p.genes),hits=genes.filter(g=>selected.has(g)),background=genes.filter(g=>universe.has(g));
  const pValue=universe.size&&selected.size&&background.length?hypergeomTail(hits.length,universe.size,background.length,selected.size):null;
  return {id:p.id,name:p.name,source:p.source,measured_genes:genes.filter(g=>grouped.has(g)),hits,p_value:pValue,q_value:null as number|null,universe_size:universe.size,status:pValue===null?'Descriptive overlap; no statistical test':'Exploratory ORA'};
 });
 const tested=pathways.filter(p=>p.p_value!==null);bhAdjust(tested.map(p=>p.p_value!)).forEach((q,i)=>tested[i].q_value=q);
 const warnings=['Scores are heuristic evidence summaries, not probabilities of causality or treatment benefit.',
  'Missing layers reduce coverage and are not evidence against a target. Correlated evidence is not independent replication.',
  'Disease and tissue labels select evidence; they do not establish specificity against other diseases.'];
 if(rows.some(e=>e.synthetic))warnings.push('Synthetic teaching observations are included. No patient-level conclusions can be drawn.');
 if(rows.some(e=>e.q_value===null))warnings.push('Some observations have no supplied q-value. They remain descriptive and cannot create enrichment hits.');
 if(!rows.length)warnings.push('No evidence matches the selected context.');
 const digest=await crypto.subtle.digest('SHA-256',new TextEncoder().encode(canonical(r)));
 const {dataset:_,...parameters}=r;
 return {version:'0.1.0',runtime:'browser-1',fingerprint_format:'sorted-json-v1',dataset:r.dataset.name,context:{disease:r.disease,tissue:r.tissue,cell_type:r.cell_type},
  fingerprint:Array.from(new Uint8Array(digest),b=>b.toString(16).padStart(2,'0')).join(''),parameters,targets,pathways,selected_genes:sorted(selected),warnings,evidence_count:rows.length};
}
export function importTable(input:unknown){
 const {text:contents,disease,tissue,name}=z.object({text:z.string(),disease:z.string(),tissue:z.string(),name:z.string()}).strict().parse(input);
 const parsed=Papa.parse<Record<string,string>>(contents.replace(/^\uFEFF/,''),{header:true,skipEmptyLines:'greedy',delimiter:contents.split('\n')[0].includes('\t')?'\t':','});
 if(parsed.errors.length)throw new Error(parsed.errors.map(e=>e.message).join('; '));
 if(!parsed.data.length||parsed.data.length>20000)throw new Error('Provide 1–20,000 evidence rows');
 const expression=['gene','tumor_tpm','control_tpm'].every(k=>parsed.meta.fields?.includes(k));
 const records=parsed.data.map((raw,i)=>{
  if(expression){
   if(!raw.tumor_tpm?.trim()||!raw.control_tpm?.trim())throw new Error('TPM values are required');
   const tumor=Number(raw.tumor_tpm),control=Number(raw.control_tpm);
   if(!Number.isFinite(tumor)||!Number.isFinite(control)||tumor<0||control<0)throw new Error('TPM values must be finite and non-negative');
   const effect=Math.log2((tumor+1)/(control+1));
   return {id:`import-${i}`,gene:raw.gene,layer:'bulk_rna',disease,tissue,study:name,source:name,observation:`Tumor ${tumor} TPM; control ${control} TPM`,effect,unit:'log2((TPM+1)/(control+1))',strength:Math.min(Math.abs(effect)/4,1)};
  }
  const clean:Record<string,unknown>=Object.fromEntries(Object.entries(raw).filter(([,v])=>v!==''&&v!==null));
  for(const key of ['strength','effect','q_value','x','y'])if(key in clean){if(!String(clean[key]).trim())throw new Error(`${key} must be numeric`);clean[key]=Number(clean[key]);}
  if('synthetic' in clean){const v=String(clean.synthetic).toLowerCase();if(!['true','false','1','0','yes','no'].includes(v))throw new Error('synthetic must be boolean');clean.synthetic=['true','1','yes'].includes(v);}
  return clean;
 });
 return validateDataset({name,evidence:records,notes:[expression?'Imported descriptive TPM contrast; no replicates or inferential statistics. Strength = min(abs(log2 ratio with +1 pseudocount)/4, 1).':'Normalized user evidence; strengths and statistical results supplied by uploader.']});
}
export async function browserApi(path:string,body?:unknown){
 if(path==='demo')return validateDataset(demo);
 if(path==='validate')return validateDataset(body);
 if(path==='analyze')return analyze(body);
 if(path==='import')return importTable(body);
 if(path==='health')return {status:'ok',runtime:'browser-1'};
 throw new Error('Unknown analysis operation');
}
