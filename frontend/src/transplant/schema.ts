import {z} from 'zod';

const id=z.string().trim().min(1).max(100).regex(/^[A-Za-z0-9_.:-]+$/,'Use a coded identifier without names or spaces');
const label=z.string().trim().min(1).max(500);
const num=z.number().finite();
const nonnegative=num.min(0);
const pct=nonnegative.max(100);
const nullable=nonnegative.nullable().default(null);
const aa=z.string().regex(/^[ACDEFGHIKLMNPQRSTVWY]+$/,'Use canonical amino-acid letters');
const hla=z.string().regex(/^HLA-[ABC]\*\d{2,3}:\d{2,3}$/,'Use two-field class-I HLA, e.g. HLA-A*02:01');
const gt=z.enum(['0/0','0/1','1/0','1/1','0|0','0|1','1|0','1|1','./.']);
export const pairSchema=z.object({id,diagnosis:label,donor_id:id,recipient_id:id,build:z.literal('GRCh38'),donor_hla:z.array(hla).min(1),recipient_hla:z.array(hla).min(1),germline_source:z.enum(['nonhematopoietic','pretransplant_remission','unknown']),source:label}).strict().refine(p=>p.donor_id!==p.recipient_id,'Donor and recipient identifiers must differ');
export const sampleSchema=z.object({id,pair_id:id,day:num.int().min(-3650).max(3650),compartment:z.enum(['blood','marrow','graft','skin']),assay:z.enum(['bulk_tcr','single_cell']),unit:z.enum(['templates','cells']),total:nonnegative.int().max(1e9),origin:z.enum(['donor','recipient','mixed','unknown']),origin_method:label,source:label}).strict().refine(s=>s.assay!=='single_cell'||s.unit==='cells','Single-cell counts must use cells');
export const candidateSchema=z.object({id,pair_id:id,gene:label,chrom:z.string().regex(/^(?:[1-9]|1[0-9]|2[0-2]|X|Y)$/),pos:nonnegative.int().min(1),ref:z.string().regex(/^[ACGT]$/),alt:z.string().regex(/^[ACGT]$/),transcript:label,
 donor_gt:gt,recipient_gt:gt,donor_dp:nonnegative.int().nullable().default(null),recipient_dp:nonnegative.int().nullable().default(null),donor_gq:nullable,recipient_gq:nullable,target_allele:z.enum(['REF','ALT']),
 peptide_ref:aa.min(8).max(11),peptide_alt:aa.min(8).max(11),variant_index:nonnegative.int().min(1).max(11),hla,
 rank_el:nullable,affinity_nm:nullable,predictor:label,leukemia_tpm:nullable,healthy_hematopoietic_tpm:nullable,nonhematopoietic_max_tpm:nullable,
 immunopeptidomics:z.enum(['detected','not_detected','not_tested']),source:label,
}).strict().superRefine((c,ctx)=>{
 if(c.ref===c.alt)ctx.addIssue({code:'custom',message:'Reference and alternate bases must differ'});
 if(c.peptide_ref.length!==c.peptide_alt.length||c.variant_index>c.peptide_ref.length)ctx.addIssue({code:'custom',message:'Peptide lengths and variant index must match'});
 const changes=[...c.peptide_ref].flatMap((x,i)=>x!==c.peptide_alt[i]?[i+1]:[]);
 if(changes.length!==1||changes[0]!==c.variant_index)ctx.addIssue({code:'custom',message:'This prototype requires one annotated missense change at variant_index'});
 if(c.rank_el!==null&&c.rank_el>100)ctx.addIssue({code:'custom',message:'EL percentile rank must be 0–100'});
});
export const cloneSchema=z.object({id,pair_id:id,tra_v:z.string().max(60),tra_j:z.string().max(60),tra_cdr3:aa.nullable(),trb_v:label,trb_j:label,trb_cdr3:aa,
 pairing:z.enum(['paired','beta_only','ambiguous']),productive:z.boolean(),source:label}).strict().superRefine((c,ctx)=>{
 if(c.pairing==='beta_only'&&(c.tra_cdr3||c.tra_v||c.tra_j))ctx.addIssue({code:'custom',message:'Beta-only receptors cannot contain an alpha chain'});
 if(c.pairing==='paired'&&(!c.tra_cdr3||!c.tra_v||!c.tra_j))ctx.addIssue({code:'custom',message:'Paired receptors require TRA and TRB sequences and V/J calls'});
});
export const observationSchema=z.object({sample_id:id,clone_id:id,count:nonnegative.int(),phenotype:z.enum(['naive','memory','effector','dysfunctional','unknown']),source:label}).strict();
export const assaySchema=z.object({id,pair_id:id,candidate_id:id,clone_id:id,sample_id:id,kind:z.enum(['tetramer','reexpression','killing','normal_tissue','hla_blocking','cross_reactivity']),
 result:z.enum(['positive','negative','inconclusive','not_done']),replicates:nonnegative.int().max(1000),controls:z.enum(['adequate','incomplete','not_recorded']),
 positive_events:nonnegative.int().nullable().default(null),total_events:nonnegative.int().nullable().default(null),control_positive:nonnegative.int().nullable().default(null),control_total:nonnegative.int().nullable().default(null),
 source:label,notes:z.string().max(1500).default('')}).strict().superRefine((a,ctx)=>{
 for(const [p,t] of [[a.positive_events,a.total_events],[a.control_positive,a.control_total]])if((p===null)!==(t===null)||(p!==null&&t!==null&&(t===0||p>t)))ctx.addIssue({code:'custom',message:'Event numerator and nonzero denominator must be paired, with positive ≤ total'});
});
export const escapeSchema=z.object({id,sample_id:id,blast_percent:pct.nullable(),donor_chimerism:pct.nullable(),hla_status:z.enum(['retained','loss_supported','reduced_expression','unknown']),hla:hla.nullable(),antigen_tpm:nullable,
 somatic_clone:id,vaf:nonnegative.max(1).nullable(),source:label,notes:z.string().max(1000).default('')}).strict();
export const projectSchema=z.object({schema_version:z.literal('1.0'),name:label,synthetic:z.boolean(),pairs:z.array(pairSchema).min(1).max(500),samples:z.array(sampleSchema).max(5000),candidates:z.array(candidateSchema).max(20000),
 clones:z.array(cloneSchema).max(20000),observations:z.array(observationSchema).max(100000),assays:z.array(assaySchema).max(20000),escape:z.array(escapeSchema).max(10000),notes:z.array(z.string().max(2000)).max(50).default([])}).strict().superRefine((p,ctx)=>{
 const fail=(message:string)=>ctx.addIssue({code:'custom',message});
 for(const key of ['pairs','samples','candidates','clones','assays','escape'] as const)if(new Set(p[key].map(x=>x.id)).size!==p[key].length)fail(`Duplicate identifiers in ${key}`);
 const pairs=new Map(p.pairs.map(x=>[x.id,x])),samples=new Map(p.samples.map(x=>[x.id,x])),clones=new Map(p.clones.map(x=>[x.id,x])),candidates=new Map(p.candidates.map(x=>[x.id,x]));
 for(const x of [...p.samples,...p.clones,...p.candidates,...p.assays])if(!pairs.has(x.pair_id))fail(`Unknown pair ${x.pair_id}`);
 const seen=new Set<string>(),sums=new Map<string,number>();
 for(const o of p.observations){const s=samples.get(o.sample_id),c=clones.get(o.clone_id);if(!s||!c||s.pair_id!==c.pair_id)fail(`Observation crosses or misses pair/sample/clone: ${o.sample_id}/${o.clone_id}`);const k=o.sample_id+'|'+o.clone_id;if(seen.has(k))fail(`Duplicate observation ${k}`);seen.add(k);sums.set(o.sample_id,(sums.get(o.sample_id)||0)+o.count);}
 for(const [s,n] of sums)if(n>(samples.get(s)?.total??0))fail(`Clone counts exceed total for ${s}`);
 for(const a of p.assays){if(candidates.get(a.candidate_id)?.pair_id!==a.pair_id||clones.get(a.clone_id)?.pair_id!==a.pair_id||samples.get(a.sample_id)?.pair_id!==a.pair_id)fail(`Assay ${a.id} has incompatible references`);}
 for(const e of p.escape)if(!samples.has(e.sample_id))fail(`Unknown escape sample ${e.sample_id}`);
 for(const c of p.candidates)if(!pairs.get(c.pair_id)?.recipient_hla.includes(c.hla))fail(`Candidate ${c.id} HLA is absent from recipient typing`);
 const receptorKeys=new Set<string>();
 for(const c of p.clones){const key=JSON.stringify([c.pair_id,c.tra_v,c.tra_j,c.tra_cdr3,c.trb_v,c.trb_j,c.trb_cdr3,c.pairing]);if(receptorKeys.has(key))fail('Identical receptor definitions need one clone ID per pair');receptorKeys.add(key);}
});
export type Project=z.infer<typeof projectSchema>;
export type Candidate=z.infer<typeof candidateSchema>;
export type Clone=z.infer<typeof cloneSchema>;
export type Assay=z.infer<typeof assaySchema>;
export type Sample=z.infer<typeof sampleSchema>;
export type Pair=z.infer<typeof pairSchema>;
export type TableKey='pairs'|'samples'|'candidates'|'clones'|'observations'|'assays'|'escape';
export const tableSchemas={pairs:pairSchema,samples:sampleSchema,candidates:candidateSchema,clones:cloneSchema,observations:observationSchema,assays:assaySchema,escape:escapeSchema};
