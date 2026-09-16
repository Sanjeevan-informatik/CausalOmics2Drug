import {spawnSync} from 'node:child_process';
import {cpSync,rmSync,mkdirSync} from 'node:fs';
const built=spawnSync('npm',['run','build','--prefix','frontend'],{stdio:'inherit',env:{...process.env,VITE_ANALYSIS_RUNTIME:'browser'}});
if(built.status!==0)process.exit(built.status||1);
rmSync('dist',{recursive:true,force:true});
cpSync('frontend/dist','dist',{recursive:true});
mkdirSync('dist/documentation',{recursive:true});
for(const file of ['INPUTS.md','SCIENCE.md','PROVENANCE.md'])cpSync('docs/'+file,'dist/documentation/'+file);

cpSync('docs/transplant','dist/documentation/transplant',{recursive:true});
