import {Project,Candidate,Pair,Assay,Clone} from './schema';
export type Settings={minDepth:number;minGq:number;maxRank:number;minExpression:number;offTissueThreshold:number};
export const defaults:Settings={minDepth:20,minGq:30,maxRank:2,minExpression:1,offTissueThreshold:1};
export function directionalMismatch(c:Candidate){
 if(c.donor_gt==='./.'||c.recipient_gt==='./.')return 'unknown';
 const target=c.target_allele==='ALT'?'1':'0';
 const d=c.donor_gt.split(/[|/]/),r=c.recipient_gt.split(/[|/]/);
 return r.includes(target)&&!d.includes(target)?'recipient_only':d.includes(target)&&!r.includes(target)?'donor_only':'shared_or_absent';
}
const supported=(a:Assay)=>a.result==='positive'&&a.controls==='adequate'&&a.replicates>=2;
export function candidateResult(c:Candidate,pair:Pair,assays:Assay[],settings=defaults,clones:Clone[]=[]){
 const mismatch=directionalMismatch(c);
 const missingQc=[c.donor_dp,c.recipient_dp,c.donor_gq,c.recipient_gq].some(v=>v===null);
 const quality=!missingQc&&c.donor_dp!>=settings.minDepth&&c.recipient_dp!>=settings.minDepth&&c.donor_gq!>=settings.minGq&&c.recipient_gq!>=settings.minGq;
 const sharedHla=pair.donor_hla.includes(c.hla)&&pair.recipient_hla.includes(c.hla);
 const germline=pair.germline_source!=='unknown';
 const eligible=quality&&mismatch==='recipient_only'&&sharedHla&&germline;
 const matches=assays.filter(a=>a.candidate_id===c.id);
 const binding=c.rank_el!==null&&c.rank_el<=settings.maxRank;
 const expression=c.leukemia_tpm!==null&&c.leukemia_tpm>=settings.minExpression;
 const tetramer=matches.some(a=>a.kind==='tetramer'&&supported(a));
 // Functional evidence must converge on the same receptor, not unrelated positive assays.
 const functional=clones.filter(cl=>cl.pairing==='paired'&&cl.productive).some(cl=>['reexpression','killing','hla_blocking'].every(kind=>{
  const adequate=matches.filter(a=>a.clone_id===cl.id&&a.kind===kind&&a.controls==='adequate'&&a.replicates>=2);
  return adequate.some(supported)&&!adequate.some(a=>a.result==='negative');
 }));
 const offTissue=c.nonhematopoietic_max_tpm===null?'unknown':c.nonhematopoietic_max_tpm>settings.offTissueThreshold?'flag':'low_observed';
 const safetySignal=matches.some(a=>['normal_tissue','cross_reactivity'].includes(a.kind)&&a.result==='positive');
 const reasons:string[]=[];
 const conflicting=[...new Set(matches.map(a=>a.clone_id+'|'+a.kind))].some(key=>{
  const records=matches.filter(a=>a.clone_id+'|'+a.kind===key&&a.controls==='adequate'&&a.replicates>=2);
  return records.some(a=>a.result==='positive')&&records.some(a=>a.result==='negative');
 });
 if(conflicting)reasons.push('Conflicting controlled assay outcomes require review');
 if(mismatch!=='recipient_only')reasons.push(mismatch==='unknown'?'No-call: donor absence is unproven':mismatch==='donor_only'?'Reverse mismatch (host-versus-graft direction)':'Target allele shared or absent in recipient');
 if(!quality)reasons.push(missingQc?'Coverage or genotype quality missing':'Coverage or genotype quality below threshold');
 if(!sharedHla)reasons.push('Restricting HLA is not shared by donor and recipient');
 if(!germline)reasons.push('Recipient germline provenance is unresolved');
 if(!binding)reasons.push(c.rank_el===null?'No imported HLA presentation prediction':'EL percentile above selected threshold');
 if(!expression)reasons.push(c.leukemia_tpm===null?'Leukemia expression missing':'Leukemia expression below threshold');
 if(offTissue!=='low_observed')reasons.push(offTissue==='flag'?'Nonhematopoietic expression concern':'Normal tissue expression not assessed');
 if(safetySignal)reasons.push('Positive normal-tissue / cross-reactivity assay requires review');
 const presentation=c.immunopeptidomics==='detected';
 const stage= !eligible?'review':functional?'functional':tetramer?'binding':presentation?'presented':binding&&expression?'predicted':'incomplete';
 const evidenceCount=[quality,mismatch==='recipient_only',sharedHla,binding,expression,presentation,tetramer,functional].filter(Boolean).length;
 return {candidate:c,mismatch,quality,sharedHla,germline,eligible,binding,expression,presentation,tetramer,functional,offTissue,safetySignal,reasons,stage,evidenceCount};
}
export function rankedCandidates(p:Project,pairId:string,s=defaults){
 const pair=p.pairs.find(x=>x.id===pairId)!;
 const order:Record<string,number>={functional:5,binding:4,presented:3,predicted:2,incomplete:1,review:0};
 return p.candidates.filter(c=>c.pair_id===pairId).map(c=>candidateResult(c,pair,p.assays,s,p.clones)).sort((a,b)=>Number(b.eligible)-Number(a.eligible)||Number(a.safetySignal)-Number(b.safetySignal)||order[b.stage]-order[a.stage]||b.evidenceCount-a.evidenceCount||(a.candidate.rank_el??101)-(b.candidate.rank_el??101)||a.candidate.id.localeCompare(b.candidate.id));
}
export function wilson(k:number,n:number):[number,number]|null{if(n<=0||k<0||k>n)return null;const z=1.95996398454,p=k/n,d=1+z*z/n,center=(p+z*z/(2*n))/d,half=z*Math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return [Math.max(0,center-half),Math.min(1,center+half)];}
export function cloneSeries(p:Project,cloneId:string,compartment:string,assay:string){
 const clone=p.clones.find(c=>c.id===cloneId);if(!clone)return [];
 return p.samples.filter(s=>s.pair_id===clone.pair_id&&s.compartment===compartment&&s.assay===assay).sort((a,b)=>a.day-b.day||a.id.localeCompare(b.id)).map(s=>{
  const o=p.observations.find(o=>o.sample_id===s.id&&o.clone_id===cloneId);
  // Absent rows are missing, never zero. Explicit count=0 is a measured nondetection.
  return {sample:s,count:o?.count??null,frequency:o&&s.total>0?o.count/s.total:null,ci:o?wilson(o.count,s.total):null,phenotype:o?.phenotype??'unknown'};
 });
}
export function cloneChange(series:ReturnType<typeof cloneSeries>){const measured=series.filter(x=>x.frequency!==null);if(measured.length<2)return {fold:null,delta:null,label:'insufficient'};const a=measured[0],b=measured[measured.length-1];if(a.sample.day===b.sample.day)return {fold:null,delta:null,label:'same_day'};return {fold:a.frequency!>0?b.frequency!/a.frequency!:null,delta:(b.frequency!-a.frequency!)*100,label:a.frequency===0?'baseline_zero':'descriptive'};}
export function receptorEvidence(p:Project,clone:Clone,candidateId:string){
 const assays=p.assays.filter(a=>a.clone_id===clone.id&&a.candidate_id===candidateId);
 const steps=['tetramer','reexpression','killing','hla_blocking','normal_tissue','cross_reactivity'] as const;
 const status=Object.fromEntries(steps.map(kind=>{
  const records=assays.filter(a=>a.kind===kind);const adequate=records.filter(a=>a.controls==='adequate'&&a.replicates>=2&&['positive','negative'].includes(a.result));
  const pos=adequate.some(a=>a.result==='positive'),neg=adequate.some(a=>a.result==='negative');
  return [kind,pos&&neg?'conflicting':pos?'positive':neg?'negative':records.some(a=>a.result!=='not_done')?'preliminary':'missing'];
 }));
 const positiveSafety=assays.some(a=>['normal_tissue','cross_reactivity'].includes(a.kind)&&a.result==='positive');
 const linkedSamples=assays.map(a=>p.samples.find(s=>s.id===a.sample_id)!);
 const donorSupported=linkedSamples.length>0&&linkedSamples.every(s=>s.origin==='donor'&&s.origin_method!=='unknown');
 const functional=['reexpression','killing','hla_blocking'].every(k=>status[k]==='positive');
 return {status,assays,positiveSafety,donorSupported,functional,paired:clone.pairing==='paired'&&clone.productive,
  next:positiveSafety?'Resolve off-tissue or cross-reactivity signal':clone.pairing!=='paired'?'Resolve alpha/beta pairing':!donorSupported?'Confirm donor origin':status.tetramer!=='positive'?'Establish pMHC binding with controls':!functional?'Complete functional and HLA-dependence evidence':status.normal_tissue!=='negative'||status.cross_reactivity!=='negative'?'Broaden specificity and normal-tissue testing':'Independent replication and translational review'};
}
export function escapeSeries(p:Project,pairId:string,compartment:string){return p.escape.map(e=>({...e,sample:p.samples.find(s=>s.id===e.sample_id)!})).filter(e=>e.sample.pair_id===pairId&&e.sample.compartment===compartment).sort((a,b)=>a.sample.day-b.sample.day);}
export async function fingerprint(value:unknown){if(!crypto.subtle)throw new Error('Report fingerprinting requires HTTPS or localhost. Study JSON can still be exported from Data & provenance.');function sorted(v:unknown):unknown{if(Array.isArray(v))return v.map(sorted);if(v&&typeof v==='object')return Object.fromEntries(Object.entries(v).sort(([a],[b])=>a.localeCompare(b)).map(([k,v])=>[k,sorted(v)]));return v;}const bytes=await crypto.subtle.digest('SHA-256',new TextEncoder().encode(JSON.stringify(sorted(value))));return [...new Uint8Array(bytes)].map(b=>b.toString(16).padStart(2,'0')).join('');}
