//! the web ui: a cfd-software style single workspace.
//!
//! one central flow viewport with a data/display sidebar on the left, a
//! playback/color control panel, and a diagnostics dock below. the markup mirrors
//! how parapview and tecplot lay out: the plot is always visible and everything
//! else docks around it. the visual language follows the house style: warm
//! charcoal, a serif display brand, square corners, flat borders, and the
//! amber accent reserved for the primary action.

pub fn app_html() -> String {
    r##"<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>nuFor</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,500&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
:root{--bg:#141414;--panel:#191919;--panel2:#1c1c22;--line:#2b2b2b;--line2:#3a3a42;--text:#d8d8d8;--dim:#8b8b8b;--dimmer:#5a5a5a;--accent:#e0a458;--ok:#7fbf7f;--mono:'JetBrains Mono',ui-monospace,Menlo,monospace}
*{box-sizing:border-box;border-radius:0 !important}
html,body{height:100%}
body{margin:0;font-family:'Inter',system-ui,sans-serif;background:var(--bg);color:var(--text);display:flex;flex-direction:column;height:100vh;overflow:hidden;font-size:14px;line-height:1.5}
::selection{background:#3a3020}
::-webkit-scrollbar{width:8px}::-webkit-scrollbar-thumb{background:#3a3a3a}::-webkit-scrollbar-track{background:transparent}
header{display:flex;align-items:center;gap:.8rem;padding:.45rem .9rem;border-bottom:1px solid var(--line);background:var(--panel);font-size:.8rem}
header .brand{font-family:'Newsreader',Georgia,serif;font-weight:500;letter-spacing:-.02em;font-size:1.05rem}
header .ver{color:var(--dimmer);font-family:var(--mono);font-size:.7rem}
header .spacer{flex:1}
#app{flex:1;display:flex;min-height:0}
aside{width:248px;flex:0 0 248px;background:var(--panel);border-right:1px solid var(--line);overflow-y:auto;padding:.4rem .6rem .8rem;font-size:.78rem}
aside h4{margin:.55rem 0 .3rem;font-size:.68rem;text-transform:lowercase;letter-spacing:.5px;color:var(--dim);border-bottom:1px solid var(--line);padding-bottom:.2rem;font-weight:600}
aside label{display:flex;align-items:center;gap:.4rem;margin:.22rem 0;color:var(--dim)}
aside label span{flex:1}
aside input[type=number],aside select,aside input[type=text],aside button{background:var(--panel2);border:1px solid var(--line2);color:#c8c8d0;padding:.22rem .4rem;font-size:.76rem;font-family:var(--mono)}
aside input[type=number]{width:4.6rem}
aside input[type=checkbox]{accent-color:var(--accent)}
aside .runs{width:100%;margin-top:.35rem;background:var(--accent);color:#141414;border:1px solid var(--accent);font-weight:600;padding:.34rem;cursor:pointer;font-family:'Inter',sans-serif}
aside .runs:hover{background:#eab76d;border-color:#eab76d}
aside .stat{color:var(--dim);margin-top:.3rem;font-family:var(--mono);white-space:pre-line;font-size:.72rem}
#viewport{flex:1;position:relative;background:#101010;display:flex;flex-direction:column;align-items:stretch;justify-content:center;min-width:0;overflow:hidden;padding:8px}
#viewwrap{display:flex;gap:.6rem;align-items:flex-start;justify-content:center;flex:1;min-height:0}
#meshpane{width:240px;background:var(--panel);border:1px solid var(--line);display:flex;flex-direction:column;align-self:stretch}
#meshpane svg{width:100%;flex:1}
#meshhead{color:var(--dim);font-size:.7rem;padding:.3rem .4rem 0}
#meshStat{color:var(--dim);font-size:.68rem;font-family:var(--mono);white-space:pre-line;padding:.2rem .4rem .4rem}
#meshview{flex:1;display:flex;gap:.6rem;align-items:stretch;justify-content:center;padding:.6rem;min-height:0}
.mvpane{background:var(--panel);border:1px solid var(--line);display:flex;flex-direction:column;flex:1;min-width:0}
.mvpane svg,.mvpane img,.mvpane canvas{width:100%;flex:1;object-fit:contain;min-height:0;background:#101010}
.mvhead{color:var(--dim);font-size:.7rem;padding:.3rem .4rem 0}
.mvstat{color:var(--dim);font-size:.68rem;font-family:var(--mono);white-space:pre-line;padding:.2rem .4rem .4rem}
#stage{position:relative;box-shadow:0 0 0 1px var(--line);align-self:center}
#fimage{display:block;max-width:100%;max-height:100%;background:#101010}
#meshOv{position:absolute;top:0;left:0;pointer-events:none}
#colorbar{position:absolute;right:2px;top:2px;bottom:2px;width:20px;pointer-events:none}
#viewhead{position:absolute;top:6px;left:8px;font-size:.72rem;color:var(--dim);font-family:var(--mono);pointer-events:none}
#viewt{position:absolute;bottom:34px;left:8px;font-size:.7rem;color:var(--dim);font-family:var(--mono);pointer-events:none}
#playbar{height:32px;background:#141414;border-top:1px solid var(--line);display:flex;align-items:center;gap:.5rem;padding:0 .6rem;font-size:.74rem}
#playbar button{background:var(--panel2);border:1px solid var(--line2);color:#c8c8d0;padding:.14rem .5rem;cursor:pointer}
#playbar button:hover{border-color:var(--accent);color:var(--accent)}
#tSlide{flex:1;accent-color:var(--accent)}
#tRead{font-family:var(--mono);color:var(--dim);min-width:3.4em;text-align:right}
#dock{border-top:1px solid var(--line);background:var(--panel);height:172px;display:flex;flex-direction:column}
#docktabs{display:flex;gap:.15rem;padding:.3rem .6rem 0;border-bottom:1px solid var(--line);font-size:.76rem}
#docktabs button{background:none;border:none;color:var(--dim);padding:.28rem .7rem;cursor:pointer;border-bottom:2px solid transparent;font-family:'Inter',sans-serif}
#docktabs button:hover{color:var(--text)}
#docktabs button.active{color:var(--accent);border-bottom-color:var(--accent)}
#dockbodies{flex:1;min-height:0;padding:.4rem .6rem;font-size:.76rem}
#dockbodies>div{display:none;height:100%}
#dockbodies>div.active{display:block}
svg{display:block}
table{border-collapse:collapse;font-size:.74rem;margin-top:.3rem}
td,th{border:1px solid var(--line);padding:.16rem .5rem;text-align:right;font-family:var(--mono)}
th{color:var(--dim);font-weight:500}
.moncol{display:flex;gap:.8rem;flex-wrap:wrap;align-items:flex-start}
.monstat{color:var(--dim);font-family:var(--mono);white-space:pre-line;min-width:9rem}
#monStat2{background:#101010;border:1px solid var(--line);padding:.45rem .7rem}
button{cursor:pointer}
button:hover{border-color:var(--accent);color:var(--accent)}
aside .runs:hover{color:#141414}
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
    <h4>Mesh view</h4>
    <label style="display:block"><span>upload mesh (vtk / coords / msh)</span>
      <input type="file" id="meshFile" style="width:100%;font-size:.7rem"></label>
    <div class="stat" id="meshUpStat">2d + 3d rectilinear, gmsh 2d</div>
    <label><span>&nbsp;field</span><select id="mvField">
      <option value="rho">density</option>
      <option value="mach">mach</option>
      <option value="p">pressure</option>
    </select></label>
    <label id="mvAxisRow" style="display:none"><span>&nbsp;axis</span><select id="mvAxis">
      <option value="z">xy</option>
      <option value="y">xz</option>
      <option value="x">zy</option>
    </select></label>
    <label id="mvSliceRow" style="display:none"><span>&nbsp;slice</span><input type="range" id="mvSlice" min="0" max="1" step="0.05" value="0.5" style="flex:1"></label>
    <label style="padding-left:.5rem"><span>t end</span><input id="mvT" type="number" value="0.08" step="0.01" style="width:3.4rem"></label>
    <button id="mvSolveBtn" class="runs">solve on mesh</button>
    <div id="mvbar" style="height:4px;background:var(--panel2);margin:.35rem 0;display:none"><div style="height:100%;width:0;background:var(--accent)"></div></div>
    <div class="stat" id="mvStat"></div>
    <h4>Gallery</h4>
    <label><span>case</span><select id="gallery">
      <option value="blast">2d blast</option>
      <option value="laminar">laminar plate</option>
      <option value="turbulent">turbulent plate</option>
    </select></label>
    <div class="stat" id="galStat"></div>
    <h4>Simulation</h4>
    <label><span>case</span><select id="cfg-kind"><option value="sod">sod</option><option value="lax">lax</option></select></label>
    <label><span>cells</span><input id="cfg-n" type="number" value="200" min="2"></label>
    <label><span>t&nbsp;end</span><input id="cfg-t" type="number" value="0.2" step="0.01"></label>
    <label><span>gamma</span><input id="cfg-gamma" type="number" value="1.4" step="0.1"></label>
    <label><span>cfl</span><input id="cfg-cfl" type="number" value="0.5" step="0.05"></label>
    <button class="runs" id="runBtn">Run</button>
    <div id="runbar"><div></div></div>
    <div id="runmsg"></div>
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
    <label><span>mesh file</span><input id="meshPath" type="text" placeholder="/path/grid.vtk" style="flex:1"></label>
    <label style="padding-left:.5rem"><button id="meshLoad" class="runs">load mesh</button></label>

    <div class="stat" id="valRead"></div>
  </aside>

  <div id="viewport">
    <div id="meshview" style="display:none">
      <div class="mvpane" id="mvLeft">
        <div class="mvhead">mesh</div>
        <svg id="mvMeshSvg" viewBox="0 0 460 360" preserveAspectRatio="xMidYMid meet"></svg>
        <div class="mvstat" id="mvMeshStat"></div>
      </div>
      <div class="mvpane" id="mvRight">
        <div class="mvhead">result · <span id="mvFieldLbl">density</span></div>
        <img id="mvImg" alt="mesh result" style="display:none">
        <canvas id="mvSpin" width="460" height="360" style="display:none;cursor:grab"></canvas>
        <div class="mvstat" id="mvResStat"></div>
      </div>
    </div>
    <div id="viewwrap">
      <div id="stage">
        <img id="fimage" src="/api/image?n=128&field=rho&t=0.10" alt="flow field">
        <svg id="meshOv" width="520" height="520" style="display:none"></svg>
        <svg id="colorbar"></svg>
        <div id="viewhead">blast · density</div>
        <div id="viewt"></div>
      </div>
      <div id="meshpane">
        <div id="meshhead">mesh</div>
        <svg id="meshSvg" viewBox="0 0 240 200" preserveAspectRatio="none"></svg>
        <div id="meshStat"></div>
      </div>
    </div>
    <div id="playbar">
      <button id="playBtn2">&#9654; play</button>
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
    <button data-d="valid">Validation</button>
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
    <div id="d-valid">
      <div class="row" style="display:flex;gap:.5rem;align-items:center;margin-bottom:.3rem">
        <label style="color:var(--dim)">plate
          <select id="vdPlate"><option value="laminar">laminar blasius</option><option value="turbulent">turbulent sa</option></select></label>
        <span class="stat" id="vdLegend" style="color:var(--dim)"></span>
        <span class="stat" id="vdProgress" style="color:var(--accent)"></span>
      </div>
      <svg id="vdPlot" width="860" height="96"></svg>
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
    '<text x="0" y="12" fill="#8b8b8b" font-size="14">'+hi+'</text>'+
    '<text x="0" y="'+(h-6)+'" fill="#8b8b8b" font-size="14">'+lo+'</text>';
}
function loadFrame(t){
  document.getElementById('imgTRef'); 
  var tgt2=document.getElementById('tSlide2');
  tgt2.value=t;
  document.getElementById('tRead2').textContent=t.toFixed(3);
  document.getElementById('fimage').src='/api/frames?t='+t.toFixed(3);
  document.getElementById('viewt').textContent='t = '+t.toFixed(3);
}
var loadedMesh=null;
function redrawMesh(){
  var el=document.getElementById('meshOv'),show=document.getElementById('meshOn').checked&&document.getElementById('fieldOn').checked;
  el.style.display=show?'block':'none';
  if(!show){return;}
  var W=viz,H=viz,s='',i;
  if(loadedMesh&&loadedMesh.x&&loadedMesh.y){
    var xmin=loadedMesh.x[0],xmax=loadedMesh.x[loadedMesh.x.length-1];
    var ymin=loadedMesh.y[0],ymax=loadedMesh.y[loadedMesh.y.length-1];
    var sx=W/(xmax-xmin),sy=H/(ymax-ymin);
    var px=function(v){return(v-xmin)*sx;},py=function(v){return H-(v-ymin)*sy;};
    for(i=0;i<loadedMesh.x.length;i++){var p=px(loadedMesh.x[i]);s+='<line x1="'+p+'" y1="0" x2="'+p+'" y2="'+H+'" stroke="rgba(120,200,255,0.6)"/>';}
    for(i=0;i<loadedMesh.y.length;i++){var q=py(loadedMesh.y[i]);s+='<line x1="0" y1="'+q+'" x2="'+W+'" y2="'+q+'" stroke="rgba(120,200,255,0.6)"/>';}
  }else{
    var n=+document.getElementById('meshN').value,step=W/n;
    for(i=0;i<=n;i++){var p=i*step;s+='<line x1="'+p+'" y1="0" x2="'+p+'" y2="'+H+'" stroke="rgba(255,255,255,0.32)"/>';
      s+='<line x1="0" y1="'+p+'" x2="'+W+'" y2="'+p+'" stroke="rgba(255,255,255,0.32)"/>';}
  }
  el.innerHTML='<rect width="'+W+'" height="'+H+'" fill="none"/>'+s;
}
async function drawMeshPane(){
  var kind=galKind||'blast';
  try{
    var m=await j('/api/mesh-faces?kind='+kind);
    var el=document.getElementById('meshSvg');
    var xmin=m.xmin,xmax=m.xmax,ymin=m.ymin,ymax=m.ymax,W=240,H=200;
    var sx=W/(xmax-xmin),sy=H/(ymax-ymin);
    var px=function(v){return (v-xmin)*sx;},py=function(v){return H-(v-ymin)*sy;};
    var s='',i;
    for(i=0;i<m.x.length;i++){var p=px(m.x[i]);s+='<line x1="'+p+'" y1="0" x2="'+p+'" y2="'+H+'" stroke="rgba(255,255,255,0.18)"/>';}
    for(i=0;i<m.y.length;i++){var q=py(m.y[i]);s+='<line x1="0" y1="'+q+'" x2="'+W+'" y2="'+q+'" stroke="rgba(255,255,255,0.18)"/>';}
    el.innerHTML=s;
    document.getElementById('meshStat').textContent=(kind==='blast'?'uniform ':'clustered ')+m.nx+'\u00d7'+m.ny+
      '\nx ['+xmin.toFixed(2)+', '+xmax.toFixed(2)+']\ny ['+ymin.toFixed(2)+', '+ymax.toFixed(2)+']';
  }catch(e){document.getElementById('meshStat').textContent='mesh unavailable';}
}
async function loadImportedMesh(){
  var p=document.getElementById('meshPath').value.trim();
  if(!p){return;}
  try{
    loadedMesh=await j('/api/mesh?path='+encodeURIComponent(p));
    document.getElementById('meshOn').checked=true;
    document.getElementById('fieldOn').checked=false;
    document.getElementById('fimage').style.visibility='hidden';
    redrawMesh();
    document.getElementById('monStat').textContent='mesh '+loadedMesh.x.length+'x'+loadedMesh.y.length+' faces loaded';
  }catch(e){document.getElementById('monStat').textContent='mesh load failed: '+e;loadedMesh=null;}
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
function playLabel(playing){document.getElementById('playBtn2').textContent=playing?'\u23f8 pause':'\u25b6 play';}
function togglePlay(){
  if(dim3d()){stepSlice();return;}
  if(imgTimer){clearInterval(imgTimer);imgTimer=null;playLabel(false);return;}
  playLabel(true);
  imgTimer=setInterval(function(){imgT+=0.01;if(imgT>0.35){imgT=0.02;}loadFrame(imgT);},120);
}
function stepSlice(){
  var s=document.getElementById('u3d'),v=parseFloat(s.value)+0.05;
  if(v>1){v=0;}s.value=v;refreshViz();
}
function load3d(){refreshViz();}
function syncVol(){var v=+document.getElementById('tSlide2').value;imgT=v;loadFrame(v);}
document.getElementById('playBtn2').onclick=togglePlay;
document.getElementById('tSlide2').oninput=syncVol;
drawMeshPane();
document.getElementById('dim').onchange=refreshViz;
document.getElementById('ax3d').onchange=refreshViz;
document.getElementById('u3d').oninput=refreshViz;
document.getElementById('meshOn').onchange=redrawMesh;
document.getElementById('fieldOn').onchange=function(){var f=document.getElementById('fieldOn').checked;
  document.getElementById('fimage').style.visibility=f?'visible':'hidden';redrawMesh();};
document.getElementById('meshN').onchange=redrawMesh;
document.getElementById('meshLoad').onclick=loadImportedMesh;
document.getElementById('imgField').onchange=function(){colorbar(fld());refreshViz();};
colorbar('rho');fitViz();window.addEventListener('resize',fitViz);refreshViz();

// monitor (1d solve envelope)
async function j(url){var r=await fetch(url);return r.json();}
var runTimer=null;
function runShow(pct,msg,err){
  var bar=document.getElementById('runbar'),fill=bar.firstElementChild,el=document.getElementById('runmsg');
  bar.style.display=err?'none':'block';el.style.display='block';
  fill.style.width=(pct||0)+'%';el.textContent=msg||'';
}
async function runCase(){
  var dim=document.getElementById('dim').value;
  var cfl=(+document.getElementById('cfg-cfl').value)||0.4;
  var t=(+document.getElementById('cfg-t').value)||0.15;
  var n=Math.max(20,+document.getElementById('cfg-n').value||96);
  var nn=(dim==='3d')?Math.min(40,n):n;
  // 1d monitor case (sod/lax) keeps the Monitor + history populated
  try{
    var d1=await j('/api/run?kind='+document.getElementById('cfg-kind').value+'&n='+n+
      '&t='+t+'&gamma='+document.getElementById('cfg-gamma').value+'&cfl='+cfl);
    monitor(d1);loadHistory();
  }catch(e){}
  // the viewport case: background march, live progress
  try{
    runShow(0,dim+' blast '+nn+' starting\u2026');
    await j('/api/run-case?dim='+dim+'&n='+nn+'&t='+t+'&cfl='+cfl+'&name=case-'+dim);
    if(runTimer){clearInterval(runTimer);}
    runTimer=setInterval(async function(){
      try{
        var s=await j('/api/run-status');
        if(s.failed){runShow(0,'run failed: '+s.error,true);clearInterval(runTimer);runTimer=null;return;}
        if(s.done){runShow(100,dim+' blast '+nn+' done \u00b7 '+s.steps+' steps \u00b7 t='+s.t.toFixed(3)+' \u2192 results/case-'+dim+'.vtk');
          clearInterval(runTimer);runTimer=null;refreshViz();return;}
        runShow(s.pct,dim+' '+s.pct+'% \u00b7 step '+s.steps+' \u00b7 t='+s.t.toFixed(3)+' / '+s.t_end.toFixed(2));
      }catch(e){}
    },700);
  }catch(e){runShow(0,'run failed to start',true);}
}
function monitor(d){
  var s=d.snapshot,w=560,h=120,f='rho',y=s[f],x=s.centers;
  var mx=Math.max.apply(0,y),mn=Math.min.apply(0,y),W=w,H=h;
  var px=function(v){return(v-x[0])/(x[x.length-1]-x[0])*W;},py=function(v){return H-(v-mn)/((mx-mn)||1)*H;};
  document.getElementById('monPlot').innerHTML='<polyline points="'+y.map(function(v,i){return px(x[i])+','+py(v);}).join(' ')+'" fill="none" stroke="#e0a458" stroke-width="1.4"/>';
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
  document.getElementById('prPlot').innerHTML='<polyline points="'+pts.map(function(p){return p.s/xs*w+','+(h-(p.v-mn)/((mx-mn)||1)*h);}).join(' ')+'" fill="none" stroke="#e0a458" stroke-width="1.4"/>';
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
  var colors=['#e0a458','#9dc183','#6fc9b4','#6fa8c9','#e06c5a','#b8a0d8'],html='';
  d.series.forEach(function(s,k){html+='<polyline points="'+s.points.map(function(p){return p.s/xs*w+','+(h-(p.v-mn)/((mx-mn)||1)*h);}).join(' ')+'" fill="none" stroke="'+colors[k%6]+'" stroke-width="1.3"/>';});
  document.getElementById('cpPlot').innerHTML=html;
  document.getElementById('cpLegend').textContent=d.series.map(function(s,k){return colors[k%6]+' '+s.label;}).join('  ');
}
document.getElementById('cpBtn').onclick=runCompare;

// gallery + validation
var galKind='blast';
function isPlate(){return galKind==='laminar'||galKind==='turbulent';}
function galLoad(kind){
  galKind=kind;
  var head=document.getElementById('viewhead');
  var img=document.getElementById('fimage');
  document.getElementById('playbar').style.display=isPlate()?'none':'flex';
  document.getElementById('dim').style.display=isPlate()?'none':'block';
  document.getElementById('dimAX').style.display='none';
  document.getElementById('dimU').style.display='none';
  if(kind==='laminar'){
    head.textContent='laminar blasius plate \u00b7 '+fld();
    img.src='/api/plates/laminar/image?field='+fld();
    document.getElementById('galStat').textContent='re_x 40..2000\nlayer held vs exact\ncf vs 0.664/sqrt(re_x)';
  }else if(kind==='turbulent'){
    head.textContent='turbulent sa plate \u00b7 '+fld();
    img.src='/api/plates/turbulent/image?field='+fld();
    document.getElementById('galStat').textContent='re_x 5e4..5e5\nsa layer vs power law\nwarming in background\u2026';
    pollTurb();
  }else{
    head.textContent='blast \u00b7 '+fld();
    document.getElementById('galStat').textContent='';
    refreshViz();
  }
}
function pollTurb(){
  j('/api/plates/progress').then(function(p){
    var el=document.getElementById('vdProgress');
    if(p.turbulent_pct<100&&galKind==='turbulent'){
      el.textContent='turbulent solve '+p.turbulent_pct+'%';
      setTimeout(pollTurb,4000);
    }else{el.textContent='';}
  }).catch(function(){});
}
document.getElementById('gallery').onchange=function(e){galLoad(e.target.value);drawMeshPane();};
document.getElementById('imgField').onchange=function(){colorbar(fld());isPlate()?galLoad(galKind):refreshViz();};

async function drawValidation(){
  var kind=document.getElementById('vdPlate').value;
  try{
    var d=await j('/api/plates/'+kind+'/validation');
    var w=860,h=96,html='',i;
    var xs=d.stations.map(function(s){return s.re_x;}),ys=d.stations.map(function(s){return s.cf;});
    var mx=Math.max.apply(0,ys),mn=0,xmax=Math.max.apply(0,xs);
    var px=function(v){return v/xmax*w;},py=function(v){return h-v/mx*h;};
    html+='<polyline points="'+d.stations.map(function(s){return px(s.re_x)+','+py(s.cf);}).join(' ')+'" fill="none" stroke="#e0a458" stroke-width="1.4"/>';
    var refKey=d.turbulent?'corr':'blasius';
    var refs=d.stations.filter(function(s){return s[refKey]!==''&&s.re_x>0;});
    if(refs.length){
      html+='<polyline points="'+refs.map(function(s){return px(s.re_x)+','+py(parseFloat(s[refKey]));}).join(' ')+'" fill="none" stroke="#6fa8c9" stroke-width="1.2" stroke-dasharray="4 3"/>';
    }
    document.getElementById('vdPlot').innerHTML=html;
    document.getElementById('vdLegend').textContent=(d.turbulent?'cf vs schlichting power law':'cf vs 0.664/sqrt(re_x)')+' \u00b7 stations '+d.stations.length;
  }catch(e){document.getElementById('vdLegend').textContent='validation: plate not ready';}
}
document.getElementById('vdPlate').onchange=drawValidation;
document.getElementById('docktabs').addEventListener('click',function(e){
  var b=e.target.closest('button');if(!b)return;
  if(b.dataset.d==='valid'){drawValidation();}
});

// mesh view: upload, wireframe, solve, slices
var mvData=null, mvTimer=null, mvDrag=null;
function mvShow(on){
  document.getElementById('meshview').style.display=on?'flex':'none';
  document.getElementById('viewwrap').style.display=on?'none':'flex';
  document.getElementById('playbar').style.display=on?'none':'flex';
}
document.getElementById('meshFile').onchange=async function(){
  var f=this.files[0]; if(!f){return;}
  var stat=document.getElementById('meshUpStat');
  stat.textContent='parsing '+f.name+'\u2026';
  try{
    var buf=await f.arrayBuffer();
    var r=await fetch('/api/upload?name='+encodeURIComponent(f.name),{method:'POST',body:buf});
    var d=await r.json();
    if(!r.ok){stat.textContent='upload failed: '+(d&&d.error?d:'bad mesh');return;}
    mvData=d;
    stat.textContent=d.kind==='gmsh2d'?(d.cells+' nodes'):(d.kind+' '+d.nx+'\u00d7'+d.ny+(d.nz>1?'\u00d7'+d.nz:'')+'\n'+d.cells+' cells');
    document.getElementById('mvStat').textContent='mesh loaded \u00b7 '+(d.kind==='rect3d'?'3d spin view':d.kind==='gmsh2d'?'wireframe':'2d field');
    document.getElementById('mvAxisRow').style.display=d.kind==='rect3d'?'flex':'none';
    document.getElementById('mvSliceRow').style.display=d.kind==='rect3d'?'flex':'none';
    mvShow(true); drawMvWire();
  }catch(e){stat.textContent='upload error: '+e;}
};
async function drawMvWire(){
  if(!mvData){return;}
  try{
    var w=await j('/api/meshview/wire');
    if(w.kind==='rect2d'){
      drawMvRect2d(w);
    }else if(w.kind==='rect3d'){
      drawMvSpin(w);
    }else{
      drawMvGmsh(w);
    }
  }catch(e){document.getElementById('mvMeshStat').textContent='wireframe error';}
}
function drawMvRect2d(w){
  var svg=document.getElementById('mvMeshSvg'),cv=document.getElementById('mvSpin');
  svg.style.display='block';cv.style.display='none';
  var W=460,H=360,x0=w.x[0],x1=w.x[w.x.length-1],y0=w.y[0],y1=w.y[w.y.length-1];
  var s=Math.min(W/(x1-x0),H/(y1-y0));
  var px=function(v){return (v-x0)*s+(W-(x1-x0)*s)/2;},py=function(v){return H-(v-y0)*s-(H-(y1-y0)*s)/2;};
  var out='',i;
  for(i=0;i<w.x.length;i++){var p=px(w.x[i]);out+='<line x1="'+p+'" y1="'+py(y0)+'" x2="'+p+'" y2="'+py(y1)+'" stroke="rgba(255,255,255,0.25)"/>';}
  for(i=0;i<w.y.length;i++){var q=py(w.y[i]);out+='<line x1="'+px(x0)+'" y1="'+q+'" x2="'+px(x1)+'" y2="'+q+'" stroke="rgba(255,255,255,0.25)"/>';}
  svg.innerHTML=out;
  document.getElementById('mvMeshStat').textContent=w.x.length-1+'\u00d7'+(w.y.length-1)+' cells';
}
function drawMvGmsh(w){
  var svg=document.getElementById('mvMeshSvg'),cv=document.getElementById('mvSpin');
  svg.style.display='block';cv.style.display='none';
  var W=460,H=360,i;
  var xs=w.nodes.map(function(n){return n[0];}),ys=w.nodes.map(function(n){return n[1];});
  var x0=Math.min.apply(0,xs),x1=Math.max.apply(0,xs),y0=Math.min.apply(0,ys),y1=Math.max.apply(0,ys);
  var s=Math.min(W/(x1-x0),H/(y1-y0));
  var px=function(v){return (v-x0)*s+(W-(x1-x0)*s)/2;},py=function(v){return H-(v-y0)*s-(H-(y1-y0)*s)/2;};
  var out='';
  for(i=0;i<w.edges.length;i++){
    var a=w.nodes[w.edges[i][0]],b=w.nodes[w.edges[i][1]];
    out+='<line x1="'+px(a[0])+'" y1="'+py(a[1])+'" x2="'+px(b[0])+'" y2="'+py(b[1])+'" stroke="rgba(224,164,88,0.5)" stroke-width="0.6"/>';
  }
  svg.innerHTML=out;
  document.getElementById('mvMeshStat').textContent=w.nodes.length+' nodes\n'+w.edges.length+' edges';
}
// the 3d spin view: draw the segs json with a rotation matrix, drag to orbit.
var spin={ry:-0.6,rx:0.45,px:0,py:0};
function drawMvSpin(w){
  var svg=document.getElementById('mvMeshSvg'),cv=document.getElementById('mvSpin');
  svg.style.display='none';cv.style.display='block';
  var ctx=cv.getContext('2d'),W=cv.width,H=cv.height,i;
  var bbox=w.bbox,cx=(bbox[0]+bbox[3])/2,cy=(bbox[1]+bbox[4])/2,cz=(bbox[2]+bbox[5])/2;
  var span=Math.max(bbox[3]-bbox[0],bbox[4]-bbox[1],bbox[5]-bbox[2]);
  var scale=0.62*Math.min(W,H)/span*2;
  var cy_=(Math.cos(spin.ry)),sy=(Math.sin(spin.ry)),cx_=(Math.cos(spin.rx)),sx=(Math.sin(spin.rx));
  ctx.fillStyle='#101010';ctx.fillRect(0,0,W,H);
  // painter's sort by depth
  var pts=[];
  for(i=0;i<w.segs.length;i++){
    var g=w.segs[i];
    var p1=proj(g[0],g[1],g[2]),p2=proj(g[3],g[4],g[5]);
    pts.push({d:(p1.z+p2.z),a:p1,b:p2,w:g[6]});
  }
  function proj(x,y,z){
    x-=cx;y-=cy;z-=cz;
    var x1=x*cy_+z*sy, z1=-x*sy+z*cy_;
    var y2=y*cx_-z1*sx, z2=y*sx+z1*cx_;
    var d=1/(1+z2/span*0.9);
    return {x:W/2+x1*scale*d,y:H/2-y2*scale*d,z:z2};
  }
  pts.sort(function(a,b){return b.d-a.d;});
  for(i=0;i<pts.length;i++){
    var s=pts[i];
    var al=s.w>=1?0.9:0.14+s.w*0.2;
    ctx.strokeStyle=s.w>=1?'rgba(224,164,88,'+al+')':'rgba(216,216,216,'+al+')';
    ctx.lineWidth=s.w>=1?1.4:0.7;
    ctx.beginPath();ctx.moveTo(s.a.x,s.a.y);ctx.lineTo(s.b.x,s.b.y);ctx.stroke();
  }
  document.getElementById('mvMeshStat').textContent='drag to spin\n'+w.segs.length+' wire segments';
}
(function(){
  var cv=document.getElementById('mvSpin');
  cv.addEventListener('mousedown',function(e){mvDrag={x:e.clientX,y:e.clientY};cv.style.cursor='grabbing';});
  window.addEventListener('mouseup',function(){mvDrag=null;cv.style.cursor='grab';});
  window.addEventListener('mousemove',function(e){
    if(!mvDrag){return;}
    spin.ry+=(e.clientX-mvDrag.x)*0.01;
    spin.rx=Math.max(-1.4,Math.min(1.4,spin.rx+(e.clientY-mvDrag.y)*0.01));
    mvDrag={x:e.clientX,y:e.clientY};
    if(mvData&&mvData.kind==='rect3d'){drawMvWire();}
  });
})();
document.getElementById('mvSolveBtn').onclick=async function(){
  if(!mvData){document.getElementById('mvStat').textContent='upload a mesh first';return;}
  try{
    var t=(+document.getElementById('mvT').value)||0.08;
    await j('/api/meshview/solve?t='+t);
    if(mvTimer){clearInterval(mvTimer);}
    var bar=document.getElementById('mvbar');
    bar.style.display='block';
    mvTimer=setInterval(async function(){
      try{
        var s=await j('/api/meshview/status');
        bar.firstElementChild.style.width=(s.pct||0)+'%';
        if(s.failed){clearInterval(mvTimer);mvTimer=null;document.getElementById('mvStat').textContent='solve failed: '+s.error;return;}
        if(s.done){clearInterval(mvTimer);mvTimer=null;document.getElementById('mvStat').textContent='solve done \u00b7 '+s.steps+' steps \u00b7 t='+s.t.toFixed(3);mvShowResult();return;}
        document.getElementById('mvStat').textContent='solving '+s.pct+'% \u00b7 step '+s.steps;
      }catch(e){}
    },600);
    document.getElementById('mvStat').textContent='marching on the mesh\u2026';
  }catch(e){document.getElementById('mvStat').textContent='solve error: '+e;}
};
function mvShowResult(){
  document.getElementById('mvFieldLbl').textContent=document.getElementById('mvField').selectedOptions[0].text;
  var img=document.getElementById('mvImg'),cv=document.getElementById('mvSpin'),svg=document.getElementById('mvMeshSvg');
  var fld=document.getElementById('mvField').value;
  var url;
  if(mvData.kind==='rect3d'){
    var ax=document.getElementById('mvAxis').value,u=document.getElementById('mvSlice').value;
    url='/api/meshview/slice?field='+fld+'&axis='+ax+'&u='+u;
    cv.style.display='none';svg.style.display='none';img.style.display='block';
  }else if(mvData.kind==='rect2d'){
    url='/api/meshview/slice?field='+fld;
    cv.style.display='none';svg.style.display='none';img.style.display='block';
  }else{
    img.style.display='none';
    return; // gmsh: solve rejected, keep the wireframe
  }
  img.src=url;
  document.getElementById('mvResStat').textContent='field: '+fld+(mvData.kind==='rect3d'?' \u00b7 axis '+document.getElementById('mvAxis').selectedOptions[0].text+' \u00b7 slice '+document.getElementById('mvSlice').value:'');
}
document.getElementById('mvField').onchange=mvShowResult;
document.getElementById('mvAxis').onchange=mvShowResult;
document.getElementById('mvSlice').oninput=mvShowResult;
document.getElementById('gallery').onchange=function(e){galLoad(e.target.value);drawMeshPane();mvShow(false);};

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
