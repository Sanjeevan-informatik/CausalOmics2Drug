import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
import {analyze,hypergeomTail,importTable,validateDataset} from '../src/browser-engine';
const demo=JSON.parse(readFileSync('../data/demo.json','utf8'));
const requests=[{dataset:demo,disease:demo.evidence[0].disease,tissue:demo.evidence[0].tissue}];
for(let i=0;i<8;i++){
 const genes=Array.from({length:40},(_,j)=>'G'+j);
 const evidence=genes.slice(0,12).map((gene,j)=>({id:String(j),gene,layer:['genomics','bulk_rna','proteomics'][j%3],disease:'D',tissue:'T',study:'S',source:'test',observation:'test',strength:(j+1)/12,q_value:j%2?0.01:0.3,kind:j%3?'association':'perturbation',action:j%3?'unknown':'inhibit'}));
 requests.push({dataset:{name:'parity',evidence,universe:genes,pathways:[{id:'p',name:'p',source:'test',genes:genes.slice(i,i+12)}]},disease:'D',tissue:'T'} as typeof requests[number]);
}
const reference=JSON.parse(execFileSync('python',['-c',`import json,sys\nfrom backend.app.models import AnalysisRequest\nfrom backend.app.analysis import analyze\nprint(json.dumps([analyze(AnalysisRequest.model_validate(r)) for r in json.load(sys.stdin)]))`],{cwd:'..',input:JSON.stringify(requests),encoding:'utf8'}));
for(let i=0;i<requests.length;i++){
 const actual=await analyze(requests[i]);const expected=reference[i];
 assert.deepEqual(actual.targets,expected.targets);
 assert.deepEqual(actual.selected_genes,expected.selected_genes);
 actual.pathways.forEach((p,j)=>{const q=expected.pathways[j];assert.deepEqual(p.hits,q.hits);if(p.p_value!==null)assert.ok(Math.abs(p.p_value-q.p_value)<1e-10);else assert.equal(q.p_value,null);});
}
assert.equal(hypergeomTail(0,100,20,10),1);
assert.throws(()=>validateDataset({name:'bad',evidence:[]}));
assert.throws(()=>importTable({name:'bad',disease:'D',tissue:'T',text:'gene,tumor_tpm,control_tpm\nEGFR,nan,20'}));
assert.equal(importTable({name:'csv',disease:'D',tissue:'T',text:'gene,tumor_tpm,control_tpm\nEGFR,320,30'}).evidence[0].gene,'EGFR');
console.log('Browser/Python parity: 9 datasets; rankings, evidence and pathway probabilities agree. Import validation checks passed.');
