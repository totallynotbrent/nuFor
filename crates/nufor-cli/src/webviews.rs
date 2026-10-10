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
aside{width:248px;flex:0 0 248px;background:var(--panel);border-right:1px solid var(--line);overflow-y:auto;overflow-x:hidden;padding:.4rem .6rem .8rem;font-size:.78rem;box-sizing:border-box}
aside h4{margin:.55rem 0 .3rem;font-size:.68rem;text-transform:lowercase;letter-spacing:.5px;color:var(--dim);border-bottom:1px solid var(--line);padding-bottom:.2rem;font-weight:600}
aside label{display:flex;align-items:center;gap:.4rem;margin:.22rem 0;color:var(--dim);min-width:0;max-width:100%}
aside label span{flex:1;min-width:0}
aside input[type=number],aside select,aside input[type=text],aside button{background:var(--panel2);border:1px solid var(--line2);color:#c8c8d0;padding:.22rem .4rem;font-size:.76rem;font-family:var(--mono);max-width:100%;min-width:0;width:100%;box-sizing:border-box}
aside input[type=number]{width:4.6rem;min-width:0;flex:0 0 auto}
aside input[type=checkbox]{accent-color:var(--accent)}
aside .runs{width:100%;margin-top:.35rem;background:var(--accent);color:#141414;border:1px solid var(--accent);font-weight:600;padding:.34rem;cursor:pointer;font-family:'Inter',sans-serif}
aside .runs:hover{background:#eab76d;border-color:#eab76d}
aside .stat{color:var(--dim);margin-top:.3rem;font-family:var(--mono);white-space:pre-line;font-size:.72rem}
#setupTree details{margin:.2rem 0 .2rem .1rem}
#setupTree summary{cursor:pointer;color:var(--dim);font-size:.72rem;letter-spacing:.4px;padding:.15rem 0;user-select:none}
#setupTree summary:hover{color:var(--text)}
#setupTree details label{margin:.12rem 0;min-width:0}
#setupTree{max-height:44vh;overflow-y:auto;overflow-x:hidden}
#viewport{flex:1;position:relative;background:#101010;display:flex;flex-direction:column;align-items:stretch;justify-content:center;min-width:0;overflow:hidden;padding:8px}
#viewwrap{display:flex;flex-direction:column;flex:1;min-height:0}
#modebar{display:flex;gap:.15rem;padding:.35rem .6rem 0;border-bottom:1px solid var(--line);font-size:.76rem}
#modebar .mbtn{background:none;border:none;color:var(--dim);padding:.3rem .9rem;cursor:pointer;border-bottom:2px solid transparent;font-family:'Inter',sans-serif}
#modebar .mbtn:hover{color:var(--text)}
#modebar .mbtn.active{color:var(--accent);border-bottom-color:var(--accent)}
#stage{flex:1;position:relative;background:#101010;display:flex;align-items:center;justify-content:center;overflow:hidden;min-height:0}
#stage>img,#stage>canvas{max-width:100%;max-height:100%}
#meshCanvas{position:relative;cursor:grab}
#meshpane{width:240px;background:var(--panel);border:1px solid var(--line);display:flex;flex-direction:column;align-self:stretch}
#meshpane svg{width:100%;flex:1}
#meshhead{color:var(--dim);font-size:.7rem;padding:.3rem .4rem 0}
#meshStat{color:var(--dim);font-size:.68rem;font-family:var(--mono);white-space:pre-line;padding:.2rem .4rem .4rem}
.mvpane{background:var(--panel);border:1px solid var(--line);display:flex;flex-direction:column;flex:1;min-width:0}
.mvpane svg,.mvpane img,.mvpane canvas{width:100%;flex:1;object-fit:contain;min-height:0;background:#101010}
.mvhead{color:var(--dim);font-size:.7rem;padding:.3rem .4rem 0}
.mvstat{color:var(--dim);font-size:.68rem;font-family:var(--mono);white-space:pre-line;padding:.2rem .4rem .4rem}
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
.hbtn{background:var(--panel2);border:1px solid var(--line2);color:#c8c8d0;padding:.3rem .8rem;font-size:.78rem;font-family:'Inter',sans-serif}
.hbtn:hover{border-color:var(--accent);color:var(--accent)}
.hbtn.run{background:var(--accent);color:#141414;border-color:var(--accent);font-weight:600}
.hbtn.run:hover{background:#eab76d;color:#141414}
aside .runs:hover{color:#141414}
</style></head>
<body>
<header>
  <span class="brand">nuFor</span><span class="ver">__VER__</span>
  <span class="stat">finite-volume CFD</span>
  <button id="loadCaseBtn" class="hbtn">Load case</button>
  <button id="saveCaseBtn" class="hbtn">Save</button>
  <button id="newCaseBtn" class="hbtn">New</button>
  <span class="spacer"></span>
  <span class="stat" id="caseName" style="color:var(--accent)"></span>
  <button id="runBtn2" class="hbtn run">Run</button>
  <div id="runBarWrap" style="display:none;flex:0 0 140px;height:6px;border:1px solid #2b2b2b;background:#141414;margin:0 10px;align-self:center">
    <div id="runBarFill" style="width:0%;height:100%;background:#e0a458"></div>
  </div>
  <span class="stat" id="caseStat"></span>
  <span class="stat" id="simdStat"></span>
</header>
<div id="app">
  <aside id="data">
    <h4>Setup</h4>
    <div id="setupTree"><div class="stat">no case loaded — use Load case</div></div>
    <label style="margin-top:.3rem"><button id="rawToggle" class="runs">Edit raw TOML</button></label>
    <label id="rawRow" style="display:none;color:var(--dim);font-size:.72rem">case.toml
      <textarea id="caseText" spellcheck="false" style="width:calc(100% - 1rem);height:150px;font-family:'JetBrains Mono',monospace;font-size:.66rem;background:var(--panel2);color:var(--text);border:1px solid var(--line);padding:.3rem;margin-top:.2rem"></textarea>
    </label>

    <h4>Display</h4>
    <label><span>field</span><select id="imgField">
      <option value="rho">density</option>
      <option value="mach">mach</option>
      <option value="p">pressure</option>
    </select></label>
    <label style="padding-left:.5rem"><input type="checkbox" id="fieldOn" checked><span>field color</span></label>
    <label style="padding-left:.5rem"><input type="checkbox" id="meshOn"><span>mesh overlay</span></label>
    <label><span>mesh lines</span><select id="meshN">
      <option value="8">8</option>
      <option value="16" selected>16</option>
      <option value="32">32</option>
      <option value="64">64</option>
    </select></label>

    <h4>Results</h4>
    <div class="stat" id="resStat">run a case to see results</div>
    <svg id="csPlot" width="236" height="96"></svg>
    <div class="stat" id="csCA" style="color:var(--accent)"></div>

    <details>
    <summary style="color:var(--dim);font-size:.72rem;letter-spacing:.4px;cursor:pointer">Import mesh</summary>
    <label style="display:block"><span>upload mesh (vtk / coords / msh)</span>
      <input type="file" id="meshFile" style="width:100%;font-size:.7rem"></label>
    <div class="stat" id="meshUpStat">2d + 3d rectilinear, gmsh 2d</div>
    <label><span>mesh path (server)</span><input id="meshPath" type="text" placeholder="/path/grid.vtk" style="flex:1"></label>
    <button id="meshLoad" class="runs">load mesh by path</button>
    <div id="stlSliceRow" style="display:none">
      <label><span>slice z</span><input id="stlZ" type="number" step="0.01" value="0.05" style="width:4.5rem"></label>
      <label><span>aoa</span><input id="stlAoa" type="number" step="0.5" value="0" style="width:4rem"></label>
      <label><span>scale</span><input id="stlScale" type="number" step="0.1" value="1" style="width:4rem"></label>
      <label><span>place x</span><input id="stlDx" type="number" step="0.1" placeholder="nose x" style="width:4.5rem"></label>
      <label><span>y</span><input id="stlDy" type="number" step="0.1" placeholder="nose y" style="width:4rem"></label>
      <button id="stlCutBtn" class="runs">slice to body</button>
      <button id="stlApplyBtn" class="runs" style="display:none">apply to case</button>
      <div class="stat" id="stlStat"></div>
    </div>
    <div class="stat" id="mvStat"></div>
    <button id="mvSolveBtn" class="runs" style="display:none">solve on mesh</button>
    <div id="mvbar" style="height:4px;background:var(--panel2);margin:.35rem 0;display:none"><div style="height:100%;width:0;background:var(--accent)"></div></div>
    <label id="mvFieldRow" style="display:none"><span>field</span><select id="mvField">
      <option value="rho">density</option>
      <option value="mach">mach</option>
      <option value="p">pressure</option>
    </select></label>
    <label id="mvAxisRow" style="display:none"><span>axis</span><select id="mvAxis">
      <option value="z">xy</option>
      <option value="y">xz</option>
      <option value="x">zy</option>
    </select></label>
    <label id="mvSliceRow" style="display:none"><span>slice</span><input type="range" id="mvSlice" min="0" max="1" step="0.05" value="0.5" style="flex:1"></label>
    <label id="mvTRow" style="display:none"><span>t end</span><input id="mvT" type="number" value="0.08" step="0.01" style="width:3.4rem"></label>
    </details>

    <div class="stat" id="valRead"></div>
  </aside>

  <div id="viewport">
    <div id="viewwrap">
      <div id="modebar">
        <button id="modeMesh" class="mbtn active">Mesh</button>
        <button id="modeResults" class="mbtn">Results</button>
      </div>
      <div id="stage">
        <canvas id="meshCanvas" width="1600" height="800" style="display:block;width:100%;height:100%"></canvas>
        <img id="fimage" src="" alt="flow field" style="display:none;position:absolute;top:0;left:0">
        <div id="noFrameMsg" style="display:none;position:absolute;top:0;left:0;width:100%;height:100%;align-items:center;justify-content:center;color:#8a8a8b;font-size:.85rem;font-family:'Inter',sans-serif"></div>
        <svg id="meshOv" width="520" height="520" style="display:none"></svg>
        <svg id="axesOv" width="520" height="520" style="position:absolute;top:0;left:0;pointer-events:none"></svg>
        <svg id="colorbar"></svg>
        <div id="viewhead">load a case</div>
        <div id="viewt"></div>
      </div>
    </div>
    <div id="playbar">
      <button id="playBtn2">&#9654; play</button>
      <span style="color:var(--dim);font-size:.72rem">time</span>
      <input id="tSlide2" type="range" min="0" max="1" step="0.01" value="0" disabled>
      <span id="tRead2" style="color:var(--dim);font-family:var(--mono);font-size:.74rem">t</span>
    </div>
  </div>
</div>

<script>
var simd='';
fetch('/api/simd').then(function(r){return r.text();}).then(function(x){document.getElementById('simdStat').textContent='SIMD avx2'+x;}).catch(function(){});
var imgTimer=null, imgT=0.10, viz=520, vizW=520, vizH=260;
function fitViz(){
  // size the stage to the viewport, preserving the loaded image's
  // aspect ratio instead of forcing a square.
  var vp=document.getElementById('viewport');
  var img=document.getElementById('fimage');
  var ar=(img.naturalWidth&&img.naturalHeight)?(img.naturalWidth/img.naturalHeight):2.0;
  var aw=vp.clientWidth-48, ah=vp.clientHeight-96;
  var w,h;
  if(aw/ah>ar){h=ah;w=h*ar;}else{w=aw;h=w/ar;}
  if(w<96){w=96;h=w/ar;}
  if(h<96){h=96;w=h*ar;}
  vizW=Math.round(w);vizH=Math.round(h);
  img.style.width=vizW+'px';
  img.style.height=vizH+'px';
  var st=document.getElementById('stage').getBoundingClientRect();
  img.style.left=Math.round((st.width-vizW)/2)+'px';
  img.style.top=Math.round((st.height-vizH)/2)+'px';
  var ov=document.getElementById('meshOv'); ov.setAttribute('width',vizW); ov.setAttribute('height',vizH);
  var cb=document.getElementById('colorbar'); cb.setAttribute('height',vizH);
  drawAxes(); redrawMesh(); colorbar(fld());
}
// coordinate axes around the field: extents from the case's grid when
// known, else 0..1.
function drawAxes(){
  var el=document.getElementById('axesOv');
  var m=caseMeshData;
  var x0=m?m.xmin:0, x1=m?m.xmax:1, y0=m?m.ymin:0, y1=m?m.ymax:1;
  var W=vizW,H=vizH,f='{font:10px JetBrains Mono,monospace;fill:#8b8b8b}';
  el.setAttribute('width',W); el.setAttribute('height',H+16);
  el.innerHTML='<line x1="0" y1="'+(H-0.5)+'" x2="'+W+'" y2="'+(H-0.5)+'" stroke="#3a3a42"/>'+
    '<line x1="0.5" y1="0" x2="0.5" y2="'+H+'" stroke="#3a3a42"/>'+
    '<text x="2" y="'+(H-4)+'" '+f+'>'+y1.toFixed(2)+'</text>'+
    '<text x="2" y="11" '+f+'>'+y0.toFixed(2)+'</text>'+
    '<text x="'+(W-30)+'" y="'+(H+12)+'" '+f+'>'+x1.toFixed(2)+'</text>'+
    '<text x="2" y="'+(H+12)+'" '+f+'>'+x0.toFixed(2)+'</text>';
}
document.getElementById('fimage').addEventListener('load',function(){
  document.getElementById('noFrameMsg').style.display='none';
  document.getElementById('fimage').style.display='block';
  fitViz();
});
document.getElementById('fimage').addEventListener('error',function(){
  // no frame at that time (a case that never ran under this build):
  // show an honest message instead of the broken-image icon.
  document.getElementById('fimage').style.display='none';
  document.getElementById('fimage').removeAttribute('src');
  var m=document.getElementById('noFrameMsg');
  m.style.display='flex';
  m.textContent='no frames for this case yet \u00b7 press Run';
});
function fld(){return document.getElementById('imgField').value;}
function rangeFor(f){return f==='p'?[0,6]:(f==='mach'?[0,2]:[0.85,3.2]);}
function colorbar(f){
  var el=document.getElementById('colorbar'),h=vizH;
  var grad='<defs><linearGradient id="cbg" x1="0" y1="1" x2="0" y2="0">'+
    '<stop offset="0" stop-color="#302278"/><stop offset="0.2" stop-color="#3092d9"/>'+
    '<stop offset="0.4" stop-color="#21d6b5"/><stop offset="0.6" stop-color="#92eb40"/>'+
    '<stop offset="0.8" stop-color="#fda823"/><stop offset="1" stop-color="#d61c1c"/></linearGradient></defs>';
  var lo=rangeFor(f)[0],hi=rangeFor(f)[1];
  el.innerHTML=grad+'<rect x="0" y="0" width="16" height="'+h+'" fill="url(#cbg)"/>'+
    '<text x="0" y="12" fill="#8b8b8b" font-size="14">'+hi+'</text>'+
    '<text x="0" y="'+(h-6)+'" fill="#8b8b8b" font-size="14">'+lo+'</text>';
}
// what owns the viewport right now: 'case' (a run's baked frames),
// 'gallery' (plates/blast demo), or '3d' (demo slice).
function loadFrame(t){
  var tgt2=document.getElementById('tSlide2');
  tgt2.value=t;
  document.getElementById('tRead2').textContent=t.toFixed(3);
  document.getElementById('fimage').src='/api/frames?t='+t.toFixed(3)+'&field='+fld();
  document.getElementById('viewt').textContent='t = '+t.toFixed(3);
  var cn=currentCase||'case';
  document.getElementById('viewhead').textContent=cn+' \u00b7 '+fld();
  var pb=document.getElementById('playbar');
  if(pb){pb.style.display='flex';}
}

// keep the field image and its placeholder in sync: an img with no
// src must never be visible (the browser paints a broken-icon box).
function syncFieldImg(){
  var img=document.getElementById('fimage'),m=document.getElementById('noFrameMsg');
  var has=(img.getAttribute('src')||'').length>0;
  img.style.display=has?'block':'none';
  if(m){m.style.display=has?'none':'flex';}
}

var loadedMesh=null;
var caseMeshData=null;
function redrawMesh(){
  var el=document.getElementById('meshOv'),show=document.getElementById('meshOn').checked&&document.getElementById('fieldOn').checked;
  el.style.display=show?'block':'none';
  if(!show){return;}
  var W=vizW,H=vizH,s='',i;
  if(viewMode==='case'&&caseMeshData&&caseMeshData.x){
    var m=caseMeshData;
    var xmin=m.xmin,xmax=m.xmax,ymin=m.ymin,ymax=m.ymax;
    var sx=W/(xmax-xmin),sy=H/(ymax-ymin);
    var px=function(v){return (v-xmin)*sx;},py=function(v){return H-(v-ymin)*sy;};
    var xstep=Math.ceil(m.x.length/40),ystep=Math.ceil(m.y.length/40);
    for(i=0;i<m.x.length;i+=xstep){var p=px(m.x[i]);s+='<line x1="'+p+'" y1="0" x2="'+p+'" y2="'+H+'" stroke="rgba(120,200,255,0.45)"/>';}
    for(i=0;i<m.y.length;i+=ystep){var q=py(m.y[i]);s+='<line x1="0" y1="'+q+'" x2="'+W+'" y2="'+q+'" stroke="rgba(120,200,255,0.45)"/>';}
    if(m.body){
      var pts=m.body.map(function(r){return px(r[0]).toFixed(1)+','+py(r[1]).toFixed(1);}).join(' ');
      s+='<polyline points="'+pts+'" fill="none" stroke="#e0a458" stroke-width="1.6"/>';
    }
    el.innerHTML='<rect width="'+W+'" height="'+H+'" fill="none"/>'+s;
    return;
  }
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
async function loadImportedMesh(){
  var p=document.getElementById('meshPath').value.trim();
  if(!p){return;}
  try{
    loadedMesh=await j('/api/mesh?path='+encodeURIComponent(p));
    document.getElementById('meshOn').checked=true;
    document.getElementById('fieldOn').checked=false;
    document.getElementById('fimage').style.visibility='hidden';
    redrawMesh();
    document.getElementById('meshUpStat').textContent='mesh '+loadedMesh.x.length+'x'+loadedMesh.y.length+' faces loaded';
  }catch(e){document.getElementById('meshUpStat').textContent='mesh load failed: '+e;loadedMesh=null;}
}
// playback over the current run's baked frames
function playLabel(playing){document.getElementById('playBtn2').textContent=playing?'\u23f8 pause':'\u25b6 play';}
function togglePlay(){
  if(imgTimer){clearInterval(imgTimer);imgTimer=null;playLabel(false);return;}
  playLabel(true);
  var sl=document.getElementById('tSlide2');
  imgTimer=setInterval(function(){
    imgT+=0.01;
    if(imgT>+sl.max){imgT=+sl.min;}
    loadFrame(imgT);
  },120);
}
function syncVol(){var v=+document.getElementById('tSlide2').value;imgT=v;loadFrame(v);}
document.getElementById('playBtn2').onclick=togglePlay;
document.getElementById('tSlide2').oninput=syncVol;
document.getElementById('meshOn').onchange=redrawMesh;
document.getElementById('fieldOn').onchange=function(){var f=document.getElementById('fieldOn').checked;
  document.getElementById('fimage').style.visibility=f?'visible':'hidden';redrawMesh();};
document.getElementById('meshN').onchange=redrawMesh;
document.getElementById('meshLoad').onclick=loadImportedMesh;
document.getElementById('imgField').onchange=function(){
  colorbar(fld());
  loadFrame(imgT);
};
colorbar('rho');fitViz();window.addEventListener('resize',fitViz);

// the shared json fetch helper
async function j(url){var r=await fetch(url);return r.json();}
// the Load-case picker: a dropdown under the header button listing the
// case directory. nothing hardcoded — the menu is the /api/cases list.
async function loadCasePicker(){
  var d;
  try{d=await j('/api/cases');}catch(e){return;}
  var btn=document.getElementById('loadCaseBtn');
  var old=document.getElementById('caseMenu');
  if(old){old.remove();return;}
  var menu=document.createElement('div');
  menu.id='caseMenu';
  menu.style.cssText='position:fixed;z-index:50;background:#191919;border:1px solid #2b2b2b;padding:.3rem;min-width:200px;box-shadow:0 6px 18px rgba(0,0,0,.5)';
  var r=btn.getBoundingClientRect();
  menu.style.left=r.left+'px';menu.style.top=(r.bottom+4)+'px';
  d.cases.forEach(function(n){
    var it=document.createElement('div');
    it.textContent=n;
    it.style.cssText='padding:.32rem .5rem;cursor:pointer;color:#d8d8d8;font-size:.8rem';
    it.onmouseenter=function(){it.style.background='#232323';};
    it.onmouseleave=function(){it.style.background='none';};
    it.onclick=function(){menu.remove();openCase(n);};
    menu.appendChild(it);
  });
  document.body.appendChild(menu);
  setTimeout(function(){
    document.addEventListener('click',function close(e){
      if(!menu.contains(e.target)&&e.target!==btn){menu.remove();document.removeEventListener('click',close);}
    });
  },10);
}
async function openCase(n){
  try{
    var r=await fetch('/api/cases/load?name='+encodeURIComponent(n));
    caseTOML=await r.text();
    document.getElementById('caseText').value=caseTOML;
    document.getElementById('caseName').textContent=n;
    document.getElementById('caseStat').textContent='';
    currentCase=n;
    buildTree();
    drawMeshViewer();
    restoreCaseResults(n);
  }catch(e){document.getElementById('caseStat').textContent='load failed';}
}
// a case loaded from disk may already have a finished vtk (the cli
// runner writes one per run). restore it so the results panel works
// without re-running: frames, surface cp, and c_a populate from the
// file. no vtk: quietly leave the panel empty.
async function restoreCaseResults(n){
  try{
    var r=await fetch('/api/cases/restore?name='+encodeURIComponent(n));
    if(!r.ok){return;}
    var times=await j('/api/frames/times');
    var tlist=(times&&times.times)?times.times:times;
    if(!tlist||!tlist.length){return;}
    var sl=document.getElementById('tSlide2');
    sl.disabled=false;
    sl.min=tlist[0];sl.max=tlist[tlist.length-1];
    sl.step=(+sl.max-+sl.min)/Math.max(tlist.length-1,1);
    sl.value=+sl.max;imgT=+sl.max;
    document.getElementById('viewhead').textContent=n+' \u00b7 '+fld()+' \u00b7 restored';
    document.getElementById('caseStat').textContent=n+' restored from disk';
    loadFrame(imgT);
    drawCaseSurface();
    setViewMode('results');
  }catch(e){}
}

var currentCase=null;

// the mesh viewer: the case's grid in the main graphics window, like a
// commercial mesher — cell edges, boundary zones colored by type, the
// body geometry filled, wheel zoom + drag pan, and a stats line.
var mvZoom=1.0, mvPanX=0, mvPanY=0, mvGrid=null;
function drawMeshViewer(){
  var n=currentCase;
  if(!n){return;}
  j('/api/cases/mesh?name='+encodeURIComponent(n)).then(function(m){
    if(m.error){return;}
    mvGrid=m;
    renderMeshViewer();
  }).catch(function(){});
}
function renderMeshViewer(){
  var cv=document.getElementById('meshCanvas');
  if(!mvGrid){
    var ctx=cv.getContext('2d');
    ctx.fillStyle='#101010';ctx.fillRect(0,0,cv.width,cv.height);
    ctx.fillStyle='#6f6f78';ctx.font='15px Inter,sans-serif';
    ctx.fillText('load a case to view its mesh',cv.width/2-150,cv.height/2);
    return;
  }
  var rect=cv.getBoundingClientRect();
  if(rect.width>0){
    var dpr=window.devicePixelRatio||1;
    cv.width=Math.round(rect.width*dpr);
    cv.height=Math.round(rect.height*dpr);
  }
  var ctx=cv.getContext('2d');
  var W=cv.width,H=cv.height;
  ctx.fillStyle='#101010';ctx.fillRect(0,0,W,H);
  var m=mvGrid;
  var xmin=m.xmin,xmax=m.xmax,ymin=m.ymin,ymax=m.ymax;
  // fit the domain in PHYSICAL units, then scale to pixels
  var ar=(xmax-xmin)/(ymax-ymin);
  var vw,vh;
  if(W/H>ar){vh=(ymax-ymin)*1.05;vw=vh*ar;}else{vw=(xmax-xmin)*1.05;vh=vw/ar;}
  vw/=mvZoom;vh/=mvZoom;
  var cx=(xmin+xmax)/2+mvPanX*(xmax-xmin), cy=(ymin+ymax)/2-mvPanY*(ymax-ymin);
  var x0=cx-vw/2,x1=cx+vw/2,y0=cy-vh/2,y1=cy+vh/2;
  var sx=W/(x1-x0),sy=H/(y1-y0);
  var px=function(v){return (v-x0)*sx;},py=function(v){return H-(v-y0)*sy;};
  var kinds={left:'#4f8fd9',right:'#d97f4f',bottom:'#7fd97f',top:'#b98fd9'};
  var tree=parseTOML(caseTOML);
  var b=tree.boundaries||{};
  var zones=[['left',b.left],['right',b.right],['bottom',b.bottom],['top',b.top]];
  var lw=Math.max(1.0,1.4/mvZoom);
  zones.forEach(function(z){
    if(!z[1]){return;}
    ctx.strokeStyle=kinds[z[0]]||'#555';
    ctx.lineWidth=lw*2.4;
    ctx.beginPath();
    if(z[0]==='left'){ctx.moveTo(px(xmin),py(y0));ctx.lineTo(px(xmin),py(y1));}
    if(z[0]==='right'){ctx.moveTo(px(xmax),py(y0));ctx.lineTo(px(xmax),py(y1));}
    if(z[0]==='bottom'){ctx.moveTo(px(x0),py(ymin));ctx.lineTo(px(x1),py(ymin));}
    if(z[0]==='top'){ctx.moveTo(px(x0),py(ymax));ctx.lineTo(px(x1),py(ymax));}
    ctx.stroke();
  });
  var xstep=Math.max(1,Math.round(m.x.length/240)),ystep=Math.max(1,Math.round(m.y.length/240));
  ctx.strokeStyle='rgba(140,140,150,0.5)';
  ctx.lineWidth=lw;
  ctx.beginPath();
  for(var i=0;i<m.x.length;i+=xstep){
    var p=px(m.x[i]);ctx.moveTo(p,0);ctx.lineTo(p,H);
  }
  for(var i2=0;i2<m.y.length;i2+=ystep){
    var q=py(m.y[i2]);ctx.moveTo(0,q);ctx.lineTo(W,q);
  }
  ctx.stroke();
  if(m.body){
    ctx.beginPath();
    m.body.forEach(function(pt,k){
      var X=px(pt[0]),Y=py(pt[1]);
      if(k===0){ctx.moveTo(X,Y);}else{ctx.lineTo(X,Y);}
    });
    ctx.closePath();
    ctx.fillStyle='rgba(224,164,88,0.3)';
    ctx.fill();
    ctx.strokeStyle='#e0a458';
    ctx.lineWidth=lw*2.4;
    ctx.stroke();
  }
  ctx.font='13px JetBrains Mono,monospace';
  var lx=12;
  zones.forEach(function(z){
    if(!z[1]){return;}
    ctx.fillStyle=kinds[z[0]];
    ctx.fillRect(lx,H-28,12,12);
    ctx.fillStyle='#8b8b8b';
    ctx.fillText(z[0]+' = '+z[1],lx+16,H-18);
    lx+=16+ctx.measureText(z[0]+' = '+z[1]).width+18;
  });
  var stats=m.nx+'\u00d7'+m.ny+' = '+(m.nx*m.ny)+' cells';
  ctx.fillStyle='#8b8b8b';
  ctx.fillText(stats,W-ctx.measureText(stats).width-12,H-18);
}
function meshViewerInput(cv){
  cv.addEventListener('wheel',function(e){
    e.preventDefault();
    mvZoom=Math.max(0.3,Math.min(40,mvZoom*(e.deltaY<0?1.15:1/1.15)));
    renderMeshViewer();
  },{passive:false});
  var drag=null;
  cv.addEventListener('mousedown',function(e){drag={x:e.clientX,y:e.clientY};});
  window.addEventListener('mouseup',function(){drag=null;});
  window.addEventListener('mousemove',function(e){
    if(!drag){return;}
    mvPanX+=(e.clientX-drag.x)/cv.clientWidth;
    mvPanY-=(e.clientY-drag.y)/cv.clientHeight;
    drag={x:e.clientX,y:e.clientY};
    renderMeshViewer();
  });
}

// cases: list, edit, save, run, surface table
var caseTimer=null;var casePollN=0;
async function caseSave(asNew){
  var n=currentCase;
  var text=document.getElementById('rawRow').style.display!=='none'?document.getElementById('caseText').value:caseTOML;
  if(asNew){
    n=prompt('new case name (letters, digits, dash, underscore):');
    if(!n)return;
    var r=await fetch('/api/cases/create?name='+encodeURIComponent(n),{method:'POST',body:text});
    if(!r.ok){document.getElementById('caseStat').textContent='create failed: '+(await r.text());return;}
    currentCase=n;
    document.getElementById('caseName').textContent=n;
    document.getElementById('caseStat').textContent='created '+n;
    caseTOML=text;
    buildTree();
    drawMeshViewer();
  }else{
    if(!n){document.getElementById('caseStat').textContent='no case selected';return;}
    var r=await fetch('/api/cases/save?name='+encodeURIComponent(n),{method:'POST',body:text});
    if(!r.ok){document.getElementById('caseStat').textContent='save failed: '+(await r.text());return;}
    document.getElementById('caseStat').textContent='saved '+n;
    caseTOML=text;
    buildTree();
    drawMeshViewer();
  }
}
async function caseRun(){
  var n=currentCase;
  if(!n){document.getElementById('caseStat').textContent='no case loaded';return;}
  try{
    var d=await j('/api/cases/run?name='+encodeURIComponent(n));
    document.getElementById('caseStat').textContent='marching '+n+'\u2026';
    document.getElementById('runBarWrap').style.display='flex';
    document.getElementById('runBarFill').style.width='0%';
    if(caseTimer){clearInterval(caseTimer);}
    caseTimer=setInterval(async function(){
      try{
        var s=await j('/api/cases/run-status');
        if(s.failed){document.getElementById('caseStat').textContent='run failed: '+s.error;document.getElementById('runBarWrap').style.display='none';clearInterval(caseTimer);caseTimer=null;return;}
        document.getElementById('runBarFill').style.width=Math.min(s.pct,100)+'%';
        var doneYet=false;
        if(s.done){
          doneYet=true;
          clearInterval(caseTimer);caseTimer=null;
          document.getElementById('runBarWrap').style.display='none';
          document.getElementById('caseStat').textContent=n+' done \u00b7 '+s.steps+' steps \u00b7 t='+s.t.toFixed(3);
          drawCaseSurface();
          setViewMode('results');
          var times=await j('/api/frames/times');
          var tlist=(times&&times.times)?times.times:times;
          if(tlist&&tlist.length){
            var sl=document.getElementById('tSlide2');
            sl.disabled=false;
            sl.min=tlist[0];sl.max=tlist[tlist.length-1];sl.step=(+sl.max-+sl.min)/Math.max(tlist.length-1,1);
            sl.value=sl.max;
            document.getElementById('viewhead').textContent=n+' \u00b7 '+fld();
            imgT=+sl.max;
            loadFrame(imgT);
          }
          return;
        }
        document.getElementById('caseStat').textContent=n+' '+s.pct+'% \u00b7 step '+s.steps+' \u00b7 t='+s.t.toFixed(3)+' / '+s.t_end.toFixed(2);
        if(!doneYet&&(casePollN%5===4)){
          var times=await j('/api/frames/times');
          var tlist=(times&&times.times)?times.times:times;
          if(tlist&&tlist.length){
            var sl=document.getElementById('tSlide2');
            sl.disabled=false;
            sl.min=tlist[0];sl.max=tlist[tlist.length-1];sl.step=(+sl.max-+sl.min)/Math.max(tlist.length-1,1);
            document.getElementById('viewhead').textContent=n+' \u00b7 '+fld()+' \u00b7 live t='+s.t.toFixed(3);
            if(imgT>=+sl.max-1e-12){imgT=+sl.max;loadFrame(imgT);}
            if(viewMode!=='results'){setViewMode('results');}
          }
        }
        casePollN++;
      }catch(e){}
    },900);
  }catch(e){document.getElementById('caseStat').textContent='run failed to start';}
}
async function drawCaseSurface(){
  try{
    var d=await j('/api/cases/surface');
    var el=document.getElementById('resStat'),ca=document.getElementById('csCA');
    if(d.error){
      var nm=currentCase||'run';
      el.textContent=nm+' finished — line plot in the viewport';
      ca.textContent='';
      return;
    }
    el.textContent=(d.name||'case')+' \u00b7 surface cp, '+d.cp.length+' stations';
    ca.textContent='C_A = '+d.c_a.toFixed(4);
    var w=860,h=96,html='';
    if(d.cp.length>1){
      var xs=d.cp.map(function(r){return r[0];}),ys=d.cp.map(function(r){return r[1];});
      var x0=Math.min.apply(0,xs),x1=Math.max.apply(0,xs);
      var y0=Math.min(0,Math.min.apply(0,ys)),y1=Math.max.apply(0,ys);
      var px=function(v){return (v-x0)/Math.max(x1-x0,1e-9)*w;};
      var py=function(v){return h-(v-y0)/Math.max(y1-y0,1e-9)*h;};
      html+='<polyline points="'+d.cp.map(function(r){return px(r[0]).toFixed(1)+','+py(r[1]).toFixed(1);}).join(' ')+'" fill="none" stroke="#e0a458" stroke-width="1.4"/>';
      if(y1>0){html+='<line x1="0" y1="'+py(0)+'" x2="'+w+'" y2="'+py(0)+'" stroke="#2b2b2b" stroke-width="1"/>';}
    }
    document.getElementById('csPlot').innerHTML=html;
  }catch(e){document.getElementById('resStat').textContent='run a case to see results';}
}
// the setup tree over the loaded case.toml
var caseTOML='';
function parseTOML(text){
  // a pragmatic case.toml reader: sections, key = value lines, arrays,
  // and one level of dotted sub-tables. the editor only needs what the
  // case format uses.
  var root={},sec=root,lines=text.split('\n');
  for(var li=0;li<lines.length;li++){
    var line=lines[li].trim();
    if(!line||line[0]==='#'){continue;}
    var m=line.match(/^\[([^\]]+)\]$/);
    if(m){
      var path=m[1].split('.');
      sec=root;
      for(var pi=0;pi<path.length;pi++){
        if(!sec[path[pi]]||typeof sec[path[pi]]!=='object'){sec[path[pi]]={};}
        sec=sec[path[pi]];
      }
      continue;
    }
    m=line.match(/^([A-Za-z0-9_]+)\s*=\s*(.+)$/);
    if(m){
      var raw=m[2].trim();
      var v;
      if(raw[0]==='['){
        v=raw.replace(/"/g,'').replace(/[\[\]]/g,'').split(',').map(function(s){return s.trim();});
      }else if(raw[0]==='"'){
        v=raw.replace(/^"|"$/g,'');
      }else if(raw==='true'||raw==='false'){
        v=raw==='true';
      }else{
        v=parseFloat(raw);
        if(isNaN(v)){v=raw;}
      }
      sec[m[1]]=v;
    }
  }
  return root;
}
function tomlVal(v){
  if(Array.isArray(v)){return '['+v.map(function(x){return typeof x==='string'?'"'+x+'"':x;}).join(', ')+']';}
  if(typeof v==='string'){return '"'+v+'"';}
  return v;
}
// per-key option lists for the enum-valued fields
function enumOptions(key){
  var map={
    'physics.equations':['euler_1d','euler_2d','euler_3d','euler_axi','rans_2d_sa'],
    'physics.eos':['perfect','eqair'],
    'initial_condition.type':['uniform','two_state','blast'],
    'mesh.source':['uniform','file'],
    'body.type':['sphere_cone','sphere','polygon']
  };
  return map[key]||null;
}
function buildTree(){
  var tree=parseTOML(caseTOML);
  var el=document.getElementById('setupTree');
  var html='',order=['mesh','body','physics','initial_condition','boundaries','numerics','time','output','metadata'];
  var keys=Object.keys(tree).filter(function(k){return tree[k]&&typeof tree[k]==='object';});
  order.forEach(function(name){
    var i=keys.indexOf(name);
    if(i<0){return;}
    keys.splice(i,1);
    html+=treeNode(name,tree[name]);
  });
  keys.forEach(function(name){html+=treeNode(name,tree[name]);});
  el.innerHTML=html;
  el.querySelectorAll('input[data-path]').forEach(function(inp){
    inp.onchange=treeEdit;
  });
  el.querySelectorAll('select[data-path]').forEach(function(s){
    s.onchange=treeEdit;
  });
}
function treeNode(name,node){
  var rows=Object.keys(node).map(function(k){
    var v=node[k],path=name+'.'+k;
    if(v===null||v===undefined){return '';}
    if(v&&typeof v==='object'&&!Array.isArray(v)){
      return treeNode(path,v);
    }
    var ctl;
    if(typeof v==='number'){
      ctl='<input data-path="'+path+'" type="number" step="any" value="'+v+'" style="width:5.4rem">';
    }else if(typeof v==='boolean'){
      ctl='<input data-path="'+path+'" type="checkbox"'+(v?' checked':'')+'>';
    }else if(Array.isArray(v)){
      ctl='<input data-path="'+path+'" type="text" value="'+v.join(', ')+'" style="flex:1">';
    }else{
      var opts=enumOptions(path);
      if(opts){
        ctl='<select data-path="'+path+'">'+opts.map(function(o){return '<option value="'+o+'"'+(String(v)===o?' selected':'')+'>'+o+'</option>';}).join('')+'</select>';
      }else{
        ctl='<input data-path="'+path+'" type="text" value="'+String(v).replace(/"/g,'"')+'" style="flex:1">';
      }
    }
    return '<label><span>'+k+'</span>'+ctl+'</label>';
  }).join('');
  return '<details open><summary>'+name+'</summary>'+rows+'</details>';
}
function treeEdit(){
  var path=this.dataset.path.split('.');
  var tree=parseTOML(caseTOML);
  var sec=tree;
  for(var i=0;i<path.length-1;i++){sec=sec[path[i]];}
  var k=path[path.length-1];
  if(this.type==='checkbox'){sec[k]=this.checked;}
  else if(this.type==='number'){sec[k]=parseFloat(this.value)||0;}
  else if(this.tagName==='SELECT'){sec[k]=this.value;}
  else{
    var t=this.value.trim();
    if(t[0]==='['){sec[k]=t.replace(/"/g,'').replace(/[\[\]]/g,'').split(',').map(function(s){return s.trim();});}
    else{sec[k]=t;}
  }
  caseTOML=writeTOML(tree);
  document.getElementById('caseText').value=caseTOML;
  document.getElementById('caseStat').textContent='edited: save to apply';
}
function writeTOML(tree){
  // emit sections in the canonical order with scalar tables only
  var order=['metadata','physics','mesh','initial_condition','boundaries','body','numerics','time','output'];
  var out='schema_version = '+(tree.schema_version||1)+'\n\n';
  var keys=Object.keys(tree).filter(function(k){return k!=='schema_version'&&typeof tree[k]==='object';});
  order.forEach(function(name){
    var i=keys.indexOf(name);
    if(i<0){return;}
    keys.splice(i,1);
    out+=secTOML(name,tree[name]);
  });
  keys.forEach(function(name){out+=secTOML(name,tree[name]);});
  return out;
}
function secTOML(name,node,depth){
  var ind=depth?'.'.repeat(depth):'';
  var out='['+ind+name+']\n';
  Object.keys(node).forEach(function(k){
    var v=node[k];
    if(v&&typeof v==='object'&&!Array.isArray(v)){
      out+=secTOML(name+'.'+k,v,depth+1);
    }else{
      out+=k+' = '+tomlVal(v)+'\n';
    }
  });
  return out+'\n';
}
document.getElementById('loadCaseBtn').onclick=loadCasePicker;
document.getElementById('runBtn2').onclick=function(){caseRun();};
document.getElementById('saveCaseBtn').onclick=function(){caseSave(false);};
document.getElementById('newCaseBtn').onclick=function(){caseSave(true);};
document.getElementById('modeMesh').onclick=function(){setViewMode('mesh');};
document.getElementById('modeResults').onclick=function(){setViewMode('results');};
function setViewMode(m){
  var mc=document.getElementById('meshCanvas'),img=document.getElementById('fimage');
  var meshOn=m==='mesh';
  mc.style.display=meshOn?'block':'none';
  img.style.display=meshOn?'none':'block';
  syncFieldImg();
  document.getElementById('modeMesh').classList.toggle('active',meshOn);
  document.getElementById('modeResults').classList.toggle('active',!meshOn);
  document.getElementById('colorbar').style.display=(meshOn||!img.src)?'none':'block';
  var ax=document.getElementById('axesOv');
  if(ax){ax.style.display=meshOn?'none':'block';}
  var mo=document.getElementById('meshOv');
  if(mo){mo.style.display=(!meshOn&&document.getElementById('meshOn').checked&&img.src)?'block':'none';}
  var pb=document.getElementById('playbar');
  if(pb){pb.style.display=(!meshOn&&img.src&&img.src.indexOf('/api/frames')>=0)?'flex':'none';}
  if(meshOn){renderMeshViewer();}
  else{fitViz();}
}
meshViewerInput(document.getElementById('meshCanvas'));

document.getElementById('rawToggle').onclick=function(){
  var row=document.getElementById('rawRow');
  var on=row.style.display==='none';
  row.style.display=on?'block':'none';
  document.getElementById('rawToggle').textContent=on?'Hide raw TOML':'Edit raw TOML';
  if(on){document.getElementById('caseText').value=caseTOML;}
  else{
    caseTOML=document.getElementById('caseText').value;
    buildTree();
  }
};
setViewMode('mesh');

// mesh import: an uploaded grid renders in the mesh viewer canvas
var mvData=null;
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
    stat.textContent=d.kind==='gmsh2d'?(d.cells+' nodes'):(d.kind==='stl'?(d.cells+' tris \u00b7 wireframe preview'):(d.kind+' '+d.nx+'\u00d7'+d.ny+(d.nz>1?'\u00d7'+d.nz:'')+' \u00b7 '+d.cells+' cells'));
    document.getElementById('stlSliceRow').style.display=d.kind==='stl'?'block':'none';
    document.getElementById('mvStat').textContent='mesh loaded';
    showUploadedMesh();
  }catch(e){stat.textContent='upload error: '+e;}
};
document.getElementById('stlCutBtn').onclick=async function(){
  var q='z='+encodeURIComponent(document.getElementById('stlZ').value||'0.05');
  q+='&aoa='+encodeURIComponent(document.getElementById('stlAoa').value||'0');
  q+='&scale='+encodeURIComponent(document.getElementById('stlScale').value||'1');
  if(document.getElementById('stlDx').value!==''){q+='&dx='+encodeURIComponent(document.getElementById('stlDx').value);}
  if(document.getElementById('stlDy').value!==''){q+='&dy='+encodeURIComponent(document.getElementById('stlDy').value);}
  var stat=document.getElementById('stlStat');
  stat.textContent='slicing\u2026';
  try{
    var r=await fetch('/api/stl/verts?'+q);
    var d=await r.json();
    if(!r.ok){stat.textContent='slice failed: '+(d&&d.error?d:'bad plane');return;}
    window._stlVerts=d.verts;
    stat.textContent=d.points+' points \u2014 apply to case below';
    document.getElementById('stlApplyBtn').style.display='inline-block';
  }catch(e){stat.textContent='slice failed: '+e;}
};
document.getElementById('stlApplyBtn').onclick=function(){
  var v=window._stlVerts;
  if(!v){return;}
  if(caseTOML.indexOf('verts = [')<0){
    if(caseTOML.indexOf('[body]')<0){
      caseTOML+='\n[body]\ntype = "polygon"\n'+v+'\n';
    }else{
      caseTOML=caseTOML+'\n'+v+'\n';
    }
  }else{
    var head=caseTOML.slice(0,caseTOML.indexOf('verts = ['));
    var rest=caseTOML.slice(caseTOML.indexOf('verts = ['));
    var eol=rest.indexOf('\n]');
    var tail=eol<0?'':rest.slice(eol+2);
    caseTOML=head+v+'\n'+tail;
  }
  document.getElementById('caseText').value=caseTOML;
  document.getElementById('rawRow').style.display='';
  document.getElementById('stlStat').textContent='verts in the case \u2014 save + run';
};
document.getElementById('meshLoad').onclick=async function(){
  var path=document.getElementById('meshPath').value;
  if(!path){return;}
  var stat=document.getElementById('mvStat');
  stat.textContent='loading '+path+'\u2026';
  try{
    var r=await fetch('/api/meshview/load?path='+encodeURIComponent(path));
    var d=await r.json();
    if(!r.ok){stat.textContent='load failed: '+(d&&d.error?d:'bad mesh');return;}
    mvData=d;
    stat.textContent=d.kind==='gmsh2d'?(d.cells+' nodes'):(d.kind==='stl'?(d.cells+' tris \u00b7 wireframe preview'):(d.kind+' '+d.nx+'\u00d7'+d.ny+(d.nz>1?'\u00d7'+d.nz:'')+' \u00b7 '+d.cells+' cells'));
    document.getElementById('stlSliceRow').style.display=d.kind==='stl'?'block':'none';
    showUploadedMesh();
  }catch(e){stat.textContent='load error: '+e;}
};
async function showUploadedMesh(){
  // render the uploaded grid in the mesh viewer canvas
  setViewMode('mesh');
  try{
    var w=await j('/api/meshview/wire');
    if(!mvGrid||!mvGrid.uploaded){
      mvGrid={uploaded:true,x:[],y:[],xmin:w.xmin,xmax:w.xmax,ymin:w.ymin,ymax:w.ymax,nx:0,ny:0,body:null};
    }
    mvGrid.xmin=w.xmin;mvGrid.xmax=w.xmax;mvGrid.ymin=w.ymin;mvGrid.ymax=w.ymax;
    if(w.kind==='rect2d'&&w.x&&w.x.length){
      mvGrid.x=w.x;mvGrid.y=w.y;mvGrid.nx=w.x.length;mvGrid.ny=w.y.length;
    }
    mvZoom=1;mvPanX=0;mvPanY=0;
    renderMeshViewer();
  }catch(e){document.getElementById('mvStat').textContent='render error';}
}
document.getElementById('mvSolveBtn').onclick=async function(){
  if(!mvData){return;}
  var stat=document.getElementById('mvStat');
  stat.textContent='marching on the mesh\u2026';
  document.getElementById('mvSolveBtn').style.display='none';
  try{
    var r=await fetch('/api/meshview/run?field=rho&t='+encodeURIComponent(document.getElementById('mvT').value));
    var d=await r.json();
    var iv=setInterval(async function(){
      try{
        var s=await j('/api/run-status');
        if(s.failed){stat.textContent='run failed: '+s.error;clearInterval(iv);return;}
        if(s.done){
          clearInterval(iv);
          stat.textContent='done \u00b7 '+s.steps+' steps';
          var img=document.getElementById('fimage');
          img.src='/api/meshview/slice?field='+document.getElementById('mvField').value;
          setViewMode('results');
          fitViz();
        }
      }catch(e){}
    },600);
    document.getElementById('mvStat').textContent='marching on the mesh\u2026';
  }catch(e){stat.textContent='solve error: '+e;}
};
</script>
</body></html>"##
        .replace("__VER__", env!("CARGO_PKG_VERSION"))
}
