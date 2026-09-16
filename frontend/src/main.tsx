import React from 'react';
import {createRoot} from 'react-dom/client';
import TransplantApp from './transplant/TransplantApp';
const Legacy=React.lazy(()=>import('./App'));
function Root(){const [legacy,setLegacy]=React.useState(location.hash==='#omics');React.useEffect(()=>{const change=()=>setLegacy(location.hash==='#omics');window.addEventListener('hashchange',change);return()=>window.removeEventListener('hashchange',change);},[]);return legacy?<React.Suspense fallback={<p>Loading omics workspace…</p>}><a href="#transplant" style={{display:'block',padding:12}}>← AlloTrace transplant workspace</a><div className="legacy"><Legacy/></div></React.Suspense>:<TransplantApp/>;}
createRoot(document.getElementById('root')!).render(<React.StrictMode><Root/></React.StrictMode>);
