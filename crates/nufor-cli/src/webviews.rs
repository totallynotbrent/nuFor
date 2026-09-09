//! the web ui: a cfd-software style single workspace.
//!
//! one central flow viewport with a data/display sidebar on the left, a
//! playback/color control panel, and a diagnostics dock below. the markup mirrors
//! how parapview and tecplot lay out: the plot is always visible and everything
//! else docks around it.

pub fn app_html() -> String {
    r##"<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>nuFor</title>
<style>
:root{--bg:#0f1114;--panel:#16181d;--panel2:#1b1e24;--line:#262a32;--text:#cdd3dc;--dim:#8b93a1;--accent:#7b97aa;--accent2:#84a59d;--mono:ui-monospace,SFMono-Regular,Menlo,monospace}
*{box-sizing:border-box}
html,body{height:100%}
body{margin:0;font-family:system-ui,sans-serif;background:var(--bg);color:var(--text);display:flex;flex-direction:column;height:100vh;overflow:hidden}
header{display:flex;align-items:center;gap:.8rem;padding:.45rem .9rem;border-bottom:1px solid var(--line);background:var(--panel);font-size:.8rem}
header .brand{font-weight:700;letter-spacing:.02em}
header .ver{color:var(--dim)}
header .spacer{flex:1}
#app{flex:1;display:flex;min-height:0}
aside{width:248px;flex:0 0 248px;background:var(--panel);border-right:1px solid var(--line);overflow-y:auto;padding:.4rem .6rem .8rem;font-size:.78rem}
aside h4{margin:.55rem 0 .3rem;font-size:.68rem;text-transform:uppercase;letter-spacing:.08em;color:var(--dim);border-bottom:1px solid var(--line);padding-bottom:.2rem}
aside label{display:flex;align-items:center;gap:.4rem;margin:.22rem 0;color:var(--dim)}
aside label span{flex:1}
aside input[type=number],aside select,aside button{background:var(--panel2);border:1px solid var(--line);color:var(--text);border-radius:4px;padding:.22rem .4rem;font-size:.76rem;font-family:var(--mono)}
aside input[type=number]{width:4.6rem}
aside input[type=checkbox]{accent-color:var(--accent)}
aside .runs{width:100%;margin-top:.35rem;background:var(--accent);color:#0b0d10;border:none;font-weight:600;padding:.34rem;border-radius:4px;cursor:pointer}
aside .stat{color:var(--dim);margin-top:.3rem;font-family:var(--mono);white-space:pre-line}
#viewport{flex:1;position:relative;background:#0c0d10;display:flex;align-items:center;justify-content:center;min-width:0;overflow:hidden;padding:8px}
#stage{position:relative;box-shadow:0 0 0 1px var(--line)}
#fimage{display:block;max-width:100%;max-height:100%;background:#0c0d10}
#meshOv{position:absolute;top:0;left:0;pointer-events:none}
#colorbar{position:absolute;right:2px;top:2px;bottom:2px;width:20px;pointer-events:none}
#viewhead{position:absolute;top:6px;left:8px;font-size:.72rem;color:var(--dim);font-family:var(--mono);pointer-events:none}
#viewt{position:absolute;bottom:34px;left:8px;font-size:.7rem;color:var(--dim);font-family:var(--mono);pointer-events:none}
#playbar{position:absolute;left:0;right:0;bottom:0;height:30px;background:rgba(20,22,27,.9);border-top:1px solid var(--line);display:flex;align-items:center;gap:.5rem;padding:0 .6rem;font-size:.74rem}
#playbar button{background:var(--panel2);border:1px solid var(--line);color:var(--text);border-radius:4px;padding:.14rem .5rem;cursor:pointer}
#tSlide{flex:1;accent-color:var(--accent)}
#tRead{font-family:var(--mono);color:var(--dim);min-width:3.4em;text-align:right}
#dock{border-top:1px solid var(--line);background:var(--panel);height:172px;display:flex;flex-direction:column}
#docktabs{display:flex;gap:.15rem;padding:.3rem .6rem 0;border-bottom:1px solid var(--line);font-size:.76rem}
#docktabs button{background:none;border:none;color:var(--dim);padding:.28rem .7rem;cursor:pointer;border-bottom:2px solid transparent}
#docktabs button.active{color:var(--text);border-bottom-color:var(--accent)}
#dockbodies{flex:1;min-height:0;padding:.4rem .6rem;font-size:.76rem}
#dockbodies>div{display:none;height:100%}
#dockbodies>div.active{display:block}
svg{display:block}
table{border-collapse:collapse;font-size:.74rem;margin-top:.3rem}
td,th{border:1px solid var(--line);padding:.16rem .5rem;text-align:right;font-family:var(--mono)}
.moncol{display:flex;gap:.8rem;flex-wrap:wrap;align-items:flex-start}
.monstat{color:var(--dim);font-family:var(--mono);white-space:pre-line;min-width:9rem}
#monStat2{background:var(--panel2);border:1px solid var(--line);border-radius:4px;padding:.45rem .7rem}
button{cursor:pointer}
</style></head>
<body>
<header>
  <span class="brand">nuFor</span><span class="ver">__VER__</span>
  <span class="stat">finite-volume CFD</span>
  <span class="spacer"></span>
  <span class="stat" id="simdStat"></span>
</header>
<div id="app">
  <aside id="data">
    <h4>Simulation</h4>
    <label><span>case</span><select id="cfg-kind"><option value="sod">sod</option><option value="lax">lax</option></select></label>
    <label><span>cells</span><input id="cfg-n" type="number" value="200" min="2"></label>
    <label><span>t&nbsp;end</span><input id="cfg-t" type="number" value="0.2" step="0.01"></label>
    <label><span>gamma</span><input id="cfg-gamma" type="number" value="1.4" step="0.1"></label>
    <label><span>cfl</span><input id="cfg-cfl" type="number" value="0.5" step="0.05"></label>
    <button class="runs" id="runBtn">Run</button>
    <div class="stat" id="monStat"></div>

    <h4>Display</h4>
    <label><span>dim</span><select id="dim">
      <option value="2d">2D field</option>
      <option value="3d">3D slice</option>
    </select></label>
    <label><span>field</span><select id="imgField">
      <option value="rho">density</option>
      <option value="mach">mach</option>
      <option value="p">pressure</option>
    </select></label>
    <label style="display:none" id="dimAX"><span>axis</span><select id="ax3d">
      <option value="z">xy (z const)</option>
      <option value="y">xz (y const)</option>
      <option value="x">yz (x const)</option>
    </select></label>
    <label style="display:none" id="dimU"><span>slice</span><input id="u3d" type="range" min="0" max="1" step="0.05" value="0.5"></label>
    <label><span>layers</span><span></span></label>
    <label style="padding-left:.5rem"><input type="checkbox" id="fieldOn" checked><span>field color</span></label>
    <label style="padding-left:.5rem"><input type="checkbox" id="meshOn"><span>mesh</span></label>
    <label><span>&nbsp;mesh cell/edge</span><select id="meshN">
      <option value="8">8</option>
      <option value="16" selected>16</option>
      <option value="32">32</option>
      <option value="64">64</option>
    </select></label>

    <h4>Playback</h4>
    <div class="stat" id="playCtrl">
      <button id="playBtn">&#9654; play</button> &nbsp;
      <button id="toolBtn">&#9198;</button>
    </div>
    <label><span>t</span><input type="range" id="tSlide" min="0.02" max="0.35" step="0.005" value="0.10" style="flex:1"></label>
    <div class="stat" id="tRead">t = 0.100</div>
    <div class="stat" id="valRead"></div>
  </aside>

  <div id="viewport">
    <div id="stage">
      <img id="fimage" src="/api/image?n=128&field=rho&t=0.10" alt="flow field">
      <svg id="meshOv" width="520" height="520" style="display:none"></svg>
      <svg id="colorbar"></svg>
      <div id="viewhead">blast · density</div>
      <div id="viewt"></div>
    </div>
    <div id="playbar">
      <button id="playBtn2">&#9654;</button>
      <input id="tSlide2" type="range" min="0.02" max="0.35" step="0.005" value="0.10">
      <span id="tRead2">0.100</span>
    </div>
  </div>
</div>

<div id="dock">
  <div id="docktabs">
    <button data-d="mon" class="active">Monitor</button>
    <button data-d="probe">Probe</button>
    <button data-d="compare">Compare</button>
    <button data-d="history">History</button>
  </div>
  <div id="dockbodies">
    <div id="d-mon" class="active">
      <div class="moncol">
        <svg id="monPlot" width="560" height="120"></svg>
        <div class="monstat" id="monStat2"></div>
      </div>
    </div>
    <div id="d-probe">
      <div class="row" style="display:flex;gap:.5rem;align-items:center;margin-bottom:.3rem">
        <label style="color:var(--dim)">field
          <select id="prField"><option>rho</option><option>mach</option><option>p</option></select></label>
        <label style="color:var(--dim)">line
          <select id="prLine"><option value="H">horizontal</option><option value="V">vertical</option><option value="D">diagonal</option></select></label>
        <label style="color:var(--dim)">samples <input id="prN" type="number" value="48" min="2" max="200" style="width:4rem"></label>
        <button id="prBtn">Sample</button>
      </div>
      <svg id="prPlot" width="860" height="96"></svg>
    </div>
    <div id="d-compare">
      <div class="row" style="display:flex;gap:.5rem;align-items:center;margin-bottom:.3rem">
        <label style="color:var(--dim)">fields
          <label style="color:var(--dim)"><input type="checkbox" class="cpField" value="rho" checked> rho</label>
          <label style="color:var(--dim)"><input type="checkbox" class="cpField" value="mach" checked> mach</label>
          <label style="color:var(--dim)"><input type="checkbox" class="cpField" value="p"> p</label>
        </label>
        <label style="color:var(--dim)">n <input id="cpN" type="text" value="32,64,128" style="width:6rem"></label>
        <button id="cpBtn">Draw</button>
        <span class="stat" id="cpLegend" style="color:var(--dim)"></span>
      </div>
      <svg id="cpPlot" width="860" height="96"></svg>
    </div>
    <div id="d-history"><table id="historyTable"></table></div>
  </div>
</div>

<script>
var simd='';
fetch('/api/simd').then(function(r){return r.text();}).then(function(x){document.getElementById('simdStat').textContent='SIMD avx2'+x;}).catch(function(){});
var imgTimer=null, imgT=0.10, viz=520;
function fitViz(){
  var vp=document.getElementById('viewport');
  var s=Math.floor(Math.min(vp.clientWidth-24, vp.clientHeight-42));
  if(s<96)s=96; viz=s;
  document.getElementById('fimage').style.width=s+'px';
  document.getElementById('fimage').style.height=s+'px';
  var ov=document.getElementById('meshOv'); ov.setAttribute('width',s); ov.setAttribute('height',s);
  var cb=document.getElementById('colorbar'); cb.setAttribute('height',s);
  redrawMesh(); colorbar(fld());
}
function fld(){return document.getElementById('imgField').value;}
function rangeFor(f){return f==='p'?[1,5]:(f==='mach'?[0,3]:[1,3]);}
function colorbar(f){
  var el=document.getElementById('colorbar'),h=viz;
  var grad='<defs><linearGradient id="cbg" x1="0" y1="1" x2="0" y2="0">'+
    '<stop offset="0" stop-color="#302278"/><stop offset="0.2" stop-color="#3092d9"/>'+
    '<stop offset="0.4" stop-color="#21d6b5"/><stop offset="0.6" stop-color="#92eb40"/>'+
    '<stop offset="0.8" stop-color="#fda823"/><stop offset="1" stop-color="#d61c1c"/></linearGradient></defs>';
  var lo=rangeFor(f)[0],hi=rangeFor(f)[1];
  el.innerHTML=grad+'<rect x="0" y="0" width="16" height="'+h+'" fill="url(#cbg)"/>'+
    '<text x="0" y="12" fill="#8b93a1" font-size="14">'+hi+'</text>'+
    '<text x="0" y="'+(h-6)+'" fill="#8b93a1" font-size="14">'+lo+'</text>';
}
function loadFrame(t){
  document.getElementById('imgTRef'); 
  var tgt=document.getElementById('tSlide'),tgt2=document.getElementById('tSlide2');
  tgt.value=t; tgt2.value=t;
  document.getElementById('tRead').textContent='t = '+t.toFixed(3);
  document.getElementById('tRead2').textContent=t.toFixed(3);
  document.getElementById('fimage').src='/api/image?n=128&field='+fld()+'&t='+t.toFixed(3);
  document.getElementById('viewt').textContent='t = '+t.toFixed(3);
}
function redrawMesh(){
  var el=document.getElementById('meshOv'),show=document.getElementById('meshOn').checked&&document.getElementById('fieldOn').checked;
  el.style.display=show?'block':'none';
  if(!show){return;}
  var n=+document.getElementById('meshN').value,W=viz,H=viz,step=W/n,s='',i;
  for(i=0;i<=n;i++){var p=i*step;s+='<line x1="'+p+'" y1="0" x2="'+p+'" y2="'+H+'" stroke="rgba(255,255,255,0.32)"/>';
    s+='<line x1="0" y1="'+p+'" x2="'+W+'" y2="'+p+'" stroke="rgba(255,255,255,0.32)"/>';}
  el.innerHTML='<rect width="'+W+'" height="'+H+'" fill="none"/>'+s;
}
function dim3d(){return document.getElementById('dim').value==='3d';}
function refreshViz(){
  if(dim3d()){
    if(imgTimer){clearInterval(imgTimer);imgTimer=null;playLabel(false);}
    document.getElementById('dimAX').style.display='flex';
    document.getElementById('dimU').style.display='flex';
    var ax=document.getElementById('ax3d').value,u=document.getElementById('u3d').value;
    document.getElementById('fimage').src='/api/image3d?n=40&field='+fld()+'&axis='+ax+'&u='+u;
    document.getElementById('viewhead').textContent='blast3d · '+fld()+' · '+ax+' slice';
    document.getElementById('viewt').textContent='u = '+u;
    return;
  }
  document.getElementById('dimAX').style.display='none';
  document.getElementById('dimU').style.display='none';
  document.getElementById('viewhead').textContent='blast · '+fld();
  loadFrame(imgT);
}
function playLabel(playing){document.getElementById('playBtn').textContent=playing?'\u23f8\u23f8':'\u25b6 play';document.getElementById('playBtn2').textContent=playing?'\u23f8\u23f8':'\u25b6';}
function togglePlay(){
  if(dim3d()){stepSlice();return;}
  var btn=document.getElementById('playBtn'),btn2=document.getElementById('playBtn2');
  if(imgTimer){clearInterval(imgTimer);imgTimer=null;playLabel(false);}
  else{playLabel(true);
    imgTimer=setInterval(function(){imgT+=0.005;if(imgT>0.35){imgT=0.02;}loadFrame(imgT);},80);}
}
function stepT(){if(dim3d()){stepSlice();return;}imgT+=0.01;if(imgT>0.35){imgT=0.02;}loadFrame(imgT);}
function stepSlice(){
  var s=document.getElementById('u3d'),v=parseFloat(s.value)+0.05;
  if(v>1){v=0;}s.value=v;refreshViz();
}
function load3d(){refreshViz();}
function syncVol(){var v=+document.getElementById('tSlide').value;imgT=v;loadFrame(v);}
document.getElementById('playBtn').onclick=togglePlay;
document.getElementById('playBtn2').onclick=togglePlay;
document.getElementById('toolBtn').onclick=stepT;
document.getElementById('tSlide').oninput=syncVol;
document.getElementById('tSlide2').oninput=syncVol;
document.getElementById('dim').onchange=refreshViz;
document.getElementById('ax3d').onchange=refreshViz;
document.getElementById('u3d').oninput=refreshViz;
document.getElementById('meshOn').onchange=redrawMesh;
document.getElementById('fieldOn').onchange=function(){var f=document.getElementById('fieldOn').checked;
  document.getElementById('fimage').style.visibility=f?'visible':'hidden';redrawMesh();};
document.getElementById('meshN').onchange=redrawMesh;
document.getElementById('imgField').onchange=function(){colorbar(fld());refreshViz();};
colorbar('rho');fitViz();window.addEventListener('resize',fitViz);refreshViz();

// monitor (1d solve envelope)
async function j(url){var r=await fetch(url);return r.json();}
async function runCase(){
  var dim=document.getElementById('dim').value;
  var cfl=(+document.getElementById('cfg-cfl').value)||0.4;
  var t=(+document.getElementById('cfg-t').value)||0.15;
  // 1d monitor case (sod/lax) keeps the Monitor + history populated
  try{
    var d1=await j('/api/run?kind='+document.getElementById('cfg-kind').value+'&n='+document.getElementById('cfg-n').value+
      '&t='+t+'&gamma='+document.getElementById('cfg-gamma').value+'&cfl='+cfl);
    monitor(d1);loadHistory();
  }catch(e){}
  // the viewport case: solve the configured dim to t, save a results vtk, report convergence
  try{
    var n=Math.max(20,+document.getElementById('cfg-n').value||96);
    var nn=(dim==='3d')?Math.min(40,n):n;
    var d=await j('/api/run-case?dim='+dim+'&n='+nn+'&t='+t+'&cfl='+cfl+'&name=case-'+dim);
    document.getElementById('monStat').textContent=dim.toUpperCase()+' blast '+nn+(dim==='3d'?'^3':'\u00d7'+nn)+
      '\nsteps='+d.steps+'\ntime='+d.time.toFixed(4)+'\nsolve '+d.seconds.toFixed(2)+'s  \u2192 results/case-'+dim+'.vtk';
    refreshViz();
  }catch(e){}
}
function monitor(d){
  var s=d.snapshot,w=560,h=120,f='rho',y=s[f],x=s.centers;
  var mx=Math.max.apply(0,y),mn=Math.min.apply(0,y),W=w,H=h;
  var px=function(v){return(v-x[0])/(x[x.length-1]-x[0])*W;},py=function(v){return H-(v-mn)/((mx-mn)||1)*H;};
  document.getElementById('monPlot').innerHTML='<polyline points="'+y.map(function(v,i){return px(x[i])+','+py(v);}).join(' ')+'" fill="none" stroke="#4f8" stroke-width="1.4"/>';
  document.getElementById('monStat2').textContent='steps='+d.steps+'\ntime='+s.time.toFixed(4)+'\nresidual='+d.residual.toExponential(2)+'\nreason='+d.reason;
}
document.getElementById('runBtn').onclick=runCase;

// probe
async function runProbe(){
  var f=document.getElementById('prField').value,l=document.getElementById('prLine').value,
      n=+document.getElementById('prN').value,L={H:[0,0.5,1,0.5],V:[0.5,0,0.5,1],D:[0.1,0.1,0.9,0.9]}[l];
  var d=await j('/api/probe?n=128&field='+f+'&x0='+L[0]+'&y0='+L[1]+'&x1='+L[2]+'&y1='+L[3]+'&samples='+n);
  var pts=d.samples,w=860,h=96,mx=Math.max.apply(0,pts.map(function(p){return p.v;})),
      mn=Math.min.apply(0,pts.map(function(p){return p.v;})),xs=Math.max.apply(0,pts.map(function(p){return p.s;}))||1;
  document.getElementById('prPlot').innerHTML='<polyline points="'+pts.map(function(p){return p.s/xs*w+','+(h-(p.v-mn)/((mx-mn)||1)*h);}).join(' ')+'" fill="none" stroke="#4f8" stroke-width="1.4"/>';
}
document.getElementById('prBtn').onclick=runProbe;
document.getElementById('prField').onchange=runProbe;
document.getElementById('prLine').onchange=runProbe;

// compare
async function runCompare(){
  var fields=[].slice.call(document.querySelectorAll('.cpField')).filter(function(c){return c.checked;}).map(function(c){return c.value;});
  var d=await j('/api/compare?fields='+fields.join(',')+'&n='+(document.getElementById('cpN').value||'64')+'&line=H&samples=48');
  var w=860,h=96,all=[],i;
  d.series.forEach(function(s){all=all.concat(s.points.map(function(p){return p.v;}));});
  var mx=Math.max.apply(0,all),mn=Math.min.apply(0,all),xs=Math.max.apply(0,d.series[0].points.map(function(p){return p.s;}))||1;
  var colors=['#4f8','#f84','#4cf','#fc4','#f4c','#8f4'],html='';
  d.series.forEach(function(s,k){html+='<polyline points="'+s.points.map(function(p){return p.s/xs*w+','+(h-(p.v-mn)/((mx-mn)||1)*h);}).join(' ')+'" fill="none" stroke="'+colors[k%6]+'" stroke-width="1.3"/>';});
  document.getElementById('cpPlot').innerHTML=html;
  document.getElementById('cpLegend').textContent=d.series.map(function(s,k){return colors[k%6]+' '+s.label;}).join('  ');
}
document.getElementById('cpBtn').onclick=runCompare;

// history
async function loadHistory(){try{
  var h=await j('/api/history');
  document.getElementById('historyTable').innerHTML='<tr><th>#</th><th>kind</th><th>n</th><th>steps</th><th>time</th><th>residual</th></tr>'+
    h.map(function(r){return '<tr><td>'+r.index+'</td><td>'+r.kind+'</td><td>'+r.n+'</td><td>'+r.steps+'</td><td>'+r.time.toFixed(4)+'</td><td>'+r.residual.toExponential(2)+'</td></tr>';}).join('');
}catch(e){}
}
document.getElementById('docktabs').onclick=function(e){
  var b=e.target.closest('button');if(!b)return;
  [].slice.call(document.querySelectorAll('#docktabs button')).forEach(function(x){x.classList.remove('active');});
  [].slice.call(document.querySelectorAll('#dockbodies>div')).forEach(function(x){x.classList.remove('active');});
  b.classList.add('active');document.getElementById('d-'+b.dataset.d).classList.add('active');
};
runCase();loadHistory();runProbe();runCompare();
</script>
</body></html>"##
        .replace("__VER__", env!("CARGO_PKG_VERSION"))
}
