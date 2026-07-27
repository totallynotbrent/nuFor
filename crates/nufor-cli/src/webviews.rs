//! the web ui view: a plain app shell with one panel per feature.
//!
//! intentionally unstyled and functional; an external design tool restyles it
//! against the same data api, so the markup is the contract.

pub fn app_html() -> String {
    r##"<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><title>nuFor</title>
<style>
body{font-family:system-ui,sans-serif;margin:0;color:#111}
header{display:flex;align-items:baseline;gap:1rem;padding:.6rem 1rem;border-bottom:1px solid #ddd}
header h1{font-size:1.1rem;margin:0}
nav{display:flex;gap:.4rem;padding:.6rem 1rem;border-bottom:1px solid #eee}
nav button{background:#f5f5f5;border:1px solid #ddd;border-radius:6px;padding:.3rem .7rem;cursor:pointer}
nav button.active{background:#0b5;color:#fff;border-color:#0b5}
main{padding:1rem}
section{display:none}
section.active{display:block}
.row{display:flex;gap:1rem;flex-wrap:wrap;align-items:center;margin:.4rem 0}
label{font-size:.85rem}
input{width:5rem;padding:.2rem .4rem}
svg{background:#111;border-radius:6px;display:block}
table{border-collapse:collapse;font-size:.85rem}
td,th{border:1px solid #ddd;padding:.2rem .6rem;text-align:right}
.stat{font-size:.85rem;color:#444}
#msg{font-size:.85rem;min-height:1.2em}
</style></head>
<body>
<header><h1>nuFor</h1><span class="stat">1D Euler — local dashboard</span></header>
<nav id="nav">
  <button data-v="solve" class="active">Solve</button>
  <button data-v="case">Case</button>
  <button data-v="verify">Verify</button>
  <button data-v="output">Output</button>
  <button data-v="bench">Benchmark</button>
  <button data-v="history">History</button>
  <button data-v="image">2D field</button>
  <button data-v="probe">Probe</button>
  <button data-v="compare">Compare</button>
</nav>
<main>

<section id="v-solve" class="active">
  <div class="row">
    <button id="runBtn">Run</button>
    <span class="stat">field:</span>
    <select id="field">
      <option value="rho">density</option>
      <option value="u">velocity</option>
      <option value="p">pressure</option>
      <option value="mach">mach</option>
      <option value="m">momentum</option>
      <option value="e">energy</option>
    </select>
    <label><input type="checkbox" id="exactChk"> exact overlay</label>
    <span id="msg"></span>
  </div>
  <svg id="plot" width="760" height="260"></svg>
  <div class="row stat" id="statline"></div>
</section>

<section id="v-case">
  <h2>Case</h2>
  <div class="row"><label>kind <select id="cfg-kind"><option>sod</option><option>lax</option></select></label>
    <label>cells n <input id="cfg-n" type="number" value="200" min="2"></label>
    <label>t_end <input id="cfg-t" type="number" value="0.2" step="0.01"></label>
    <label>gamma <input id="cfg-gamma" type="number" value="1.4" step="0.1"></label>
    <label>cfl <input id="cfg-cfl" type="number" value="0.5" step="0.05"></label>
    <label>boundary <select id="cfg-bc"><option>transmissive</option><option>reflective</option></select></label>
  </div>
  <div class="row stat">the run button in the Solve panel uses these values.</div>
</section>

<section id="v-verify">
  <h2>Verification</h2>
  <p class="stat">Runs against the exact 1D Riemann solution. Toggle the overlay in the
  Solve panel to see it; the L1 error is shown below.</p>
  <div class="row"><button id="verifyBtn">Run + verify</button><span id="verifyline" class="stat"></span></div>
</section>

<section id="v-output">
  <h2>Output & export</h2>
  <p class="stat">Writes the current snapshot and lets you download it.</p>
  <div class="row">
    <a id="dl-csv" class="stat" href="/api/export?format=csv" download="snapshot.csv">download CSV</a>
    <a id="dl-vtk" class="stat" href="/api/export?format=vtk" download="snapshot.vtk">download VTK</a>
    <a id="dl-h5" class="stat" href="/api/export?format=h5" download="snapshot.h5">download HDF5</a>
  </div>
</section>

<section id="v-bench">
  <h2>Benchmark</h2>
  <div class="row"><button id="benchBtn">Run sweep</button><span id="benchline" class="stat"></span></div>
  <table id="benchtable"></table>
</section>

<section id="v-history">
  <h2>Run history</h2>
  <table id="historytable"></table>
</section>

<section id="v-image">
  <h2>2D field</h2>
  <p class="stat">A 2D blast wave solved and rendered to an image; pick the scalar
  field to view. Reloads on change.</p>
  <div class="row">
    <label>field <select id="imgField">
      <option value="rho">density</option>
      <option value="mach">mach</option>
      <option value="p">pressure</option>
    </select></label>
  </div>
  <img id="img2d" src="/api/image?n=128&field=rho" alt="2D blast field" height="420">
  <script>document.getElementById('imgField').addEventListener('change', function(){
    document.getElementById('img2d').src='/api/image?n=128&field='+this.value;
  });</script>
  <div class="row stat">density colormap; blue is ambient, red is compressed.</div>
</section>

<section id="v-probe">
  <h2>Line probe</h2>
  <p class="stat">Sample a 2D field along a straight line and chart it. Pick a line
  preset or give explicit endpoints, then Sample. The curve is value vs distance
  along the line.</p>
  <div class="row">
    <label>field <select id="prField">
      <option value="rho">density</option>
      <option value="mach">mach</option>
      <option value="p">pressure</option>
    </select></label>
    <label>line <select id="prLine">
      <option value="H">horizontal (mid y)</option>
      <option value="V">vertical (mid x)</option>
      <option value="D">diagonal</option>
    </select></label>
    <label>samples <input id="prN" type="number" value="40" min="2" max="200"></label>
    <button id="prBtn">Sample</button><span id="prline" class="stat"></span>
  </div>
  <svg id="prPlot" width="760" height="220"></svg>
</section>

<section id="v-compare">
  <h2>Comparison dashboard</h2>
  <p class="stat">Overlay the same probe line for several scalars and mesh
  resolutions on one chart, to compare profiles across a cut. Pick fields and
  resolutions, then Draw.</p>
  <div class="row">
    <label>fields
      <label><input type="checkbox" class="cpField" value="rho" checked> rho</label>
      <label><input type="checkbox" class="cpField" value="mach" checked> mach</label>
      <label><input type="checkbox" class="cpField" value="p"> p</label>
    </label>
    <label>resolutions <input id="cpN" type="text" value="32,64,128"></label>
    <label>line <select id="cpLine"><option value="H">horizontal</option><option value="V">vertical</option><option value="D">diagonal</option></select></label>
    <label>samples <input id="cpSamples" type="number" value="48" min="2" max="200"></label>
    <button id="cpBtn">Draw</button><span id="cpline" class="stat"></span>
  </div>
  <div class="row stat" id="cpLegend"></div>
  <svg id="cpPlot" width="760" height="240"></svg>
</section>

</main>
<script>
async function j(url){const r=await fetch(url);return r.json();}
async function runSolve(){
  const kind=document.getElementById('cfg-kind').value;
  const n=+document.getElementById('cfg-n').value;
  const t=+document.getElementById('cfg-t').value;
  const gamma=+document.getElementById('cfg-gamma').value;
  const cfl=+document.getElementById('cfg-cfl').value;
  const bc=document.getElementById('cfg-bc').value;
  const url=`/api/run?kind=${kind}&n=${n}&t=${t}&gamma=${gamma}&cfl=${cfl}&bc=${bc}`;
  try{const d=await j(url);afterRun(d);return d;}
  catch(e){document.getElementById('msg').textContent='run failed: '+e;}
}
function plot(field,exact){
  const s=window.snap; if(!s)return;
  const x=s.centers, y=s[field];
  const w=760,h=260,mx=Math.max(...y),mn=Math.min(...y);
  const px=v=>((v-x[0])/(x[x.length-1]-x[0]))*w;
  const py=v=>(1-(v-mn)/(mx-mn))*h;
  let html='<polyline points="'+y.map((v,i)=>px(x[i])+','+py(v)).join(' ')+'" fill="none" stroke="#4f8" stroke-width="1.5"/>';
  if(exact && window.exact && window.exact[field]){const ey=window.exact[field];
    html+='<polyline points="'+ey.map((v,i)=>px(x[i])+','+py(v)).join(' ')+'" fill="none" stroke="#f84" stroke-width="1.5" stroke-dasharray="4 3"/>';}
  document.getElementById('plot').innerHTML=html;
}
function afterRun(d){
  window.snap=d.snapshot; window.exact=d.exact;
  document.getElementById('statline').textContent=
    `n=${d.snapshot.n} t=${d.snapshot.time.toFixed(4)} steps=${d.steps} residual=${d.residual.toExponential(2)} ${d.reason}`;
  if(d.exact && d.exact.l1!==undefined){const ev=document.getElementById('verifyline');
    if(ev)ev.textContent='exact overlay ready; L1 rho error: '+d.exact.l1.toExponential(2);}
  plot(document.getElementById('field').value,document.getElementById('exactChk').checked);
  loadHistory();
}
document.getElementById('runBtn').onclick=runSolve;
document.getElementById('verifyBtn').onclick=async()=>{const d=await runSolve();if(d)document.getElementById('exactChk').checked=true;};
document.getElementById('field').onchange=()=>plot(document.getElementById('field').value,document.getElementById('exactChk').checked);
document.getElementById('exactChk').onchange=()=>plot(document.getElementById('field').value,document.getElementById('exactChk').checked);
async function loadHistory(){try{
  const h=await j('/api/history');
  const tb=document.getElementById('historytable');
  tb.innerHTML='<tr><th>#</th><th>kind</th><th>n</th><th>steps</th><th>time</th><th>residual</th><th>reason</th></tr>'+
   h.map(r=>`<tr><td>${r.index}</td><td>${r.kind}</td><td>${r.n}</td><td>${r.steps}</td><td>${r.time.toFixed(4)}</td><td>${r.residual.toExponential(2)}</td><td>${r.reason}</td></tr>`).join('');
}catch(e){}}
document.getElementById('benchBtn').onclick=async()=>{
  const d=await j('/api/benchmark?steps=1500');
  const simd=await j('/api/simd');
  document.getElementById('benchline').textContent=`us/step/cell (single core); SIMD: ${simd}; lower is better.`;
  const tb=document.getElementById('benchtable');
  tb.innerHTML='<tr><th>cells</th><th>us/step/cell</th><th>cell-steps/s</th></tr>'+
   d.map(r=>`<tr><td>${r.cells}</td><td>${r.us_per_step_per_cell.toFixed(3)}</td><td>${Number(r.cell_steps_per_second).toLocaleString()}</td></tr>`).join('');
};
document.getElementById('nav').onclick=e=>{
  const b=e.target.closest('button'); if(!b)return;
  document.querySelectorAll('#nav button').forEach(x=>x.classList.remove('active'));
  document.querySelectorAll('section').forEach(s=>s.classList.remove('active'));
  b.classList.add('active');
  document.getElementById('v-'+b.dataset.v).classList.add('active');
};
async function runProbe(){
  const field=document.getElementById('prField').value;
  const line=document.getElementById('prLine').value;
  const samples=+document.getElementById('prN').value;
  const L={H:[0,0.5,1,0.5], V:[0.5,0,0.5,1], D:[0.1,0.1,0.9,0.9]}[line];
  const url=`/api/probe?n=128&field=${field}&x0=${L[0]}&y0=${L[1]}&x1=${L[2]}&y1=${L[3]}&samples=${samples}`;
  try{
    const d=await j(url);
    const pts=d.samples, w=760, h=220;
    const mx=Math.max(...pts.map(p=>p.v)), mn=Math.min(...pts.map(p=>p.v));
    const xs=Math.max(...pts.map(p=>p.s))||1;
    const px=v=>v/xs*w, py=v=>h-(v-mn)/((mx-mn)||1)*h;
    document.getElementById('prPlot').innerHTML=
      '<polyline points="'+pts.map(p=>px(p.s)+','+py(p.v)).join(' ')+'" fill="none" stroke="#4f8" stroke-width="1.5"/>';
    document.getElementById('prline').textContent=
      `field=${d.field} samples=${pts.length} min=${mn.toFixed(3)} max=${mx.toFixed(3)}`;
  }catch(e){document.getElementById('prline').textContent='probe failed: '+e;}
}
document.getElementById('prBtn').onclick=runProbe;
document.getElementById('prField').onchange=runProbe;
document.getElementById('prLine').onchange=runProbe;
async function runCompare(){
  const fields=[...document.querySelectorAll('.cpField')].filter(c=>c.checked).map(c=>c.value);
  const res=document.getElementById('cpN').value.trim()||'32,64';
  const line=document.getElementById('cpLine').value;
  const samples=+document.getElementById('cpSamples').value;
  const url=`/api/compare?fields=${fields.join(',')}&n=${res}&line=${line}&samples=${samples}`;
  const colors=['#4f8','#f84','#4cf','#fc4','#f4c','#8f4'];
  try{
    const d=await j(url);
    const series=d.series, w=760, h=240;
    const all=d.series.flatMap(s=>s.points.map(p=>p.v));
    const mx=Math.max(...all), mn=Math.min(...all);
    const xs=Math.max(...d.series.flatMap(s=>s.points.map(p=>p.s)))||1;
    const px=v=>v/xs*w, py=v=>h-(v-mn)/((mx-mn)||1)*h;
    let polys='';let legend='';
    series.forEach((s,i)=>{const c=colors[i%colors.length];
      polys+='<polyline points="'+s.points.map(p=>px(p.s)+','+py(p.v)).join(' ')+'" fill="none" stroke="'+c+'" stroke-width="1.4"/>';
      legend+='<span style="color:'+c+';margin-right:.8rem">'+s.label+'</span>';});
    document.getElementById('cpPlot').innerHTML=polys;
    document.getElementById('cpLegend').innerHTML=legend;
    document.getElementById('cpline').textContent=`${series.length} series · ${d.line}-line · ${d.series[0].points.length} pts each`;
  }catch(e){document.getElementById('cpline').textContent='compare failed: '+e;}
}
document.getElementById('cpBtn').onclick=runCompare;
document.getElementById('cpLine').onchange=runCompare;
(async()=>{try{const d=await j('/api/result?auto');afterRun(d);}catch(e){}})().then(runProbe);
</script>
</body></html>"##
        .to_string()
}
