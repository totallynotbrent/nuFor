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
  <p class="stat">A 2D blast wave (high-pressure disc) solved and rendered to an
  image. Reload the panel to recompute it.</p>
  <img src="/api/image?n=128" alt="2D blast density" height="420">
  <div class="row stat">density colormap; blue is ambient, red is compressed.</div>
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
  document.getElementById('benchline').textContent='us/step/cell (single core); lower is better.';
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
(async()=>{try{const d=await j('/api/result?auto');afterRun(d);}catch(e){}})();
</script>
</body></html>"##
        .to_string()
}
