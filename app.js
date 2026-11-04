(function () {
  'use strict';
  const data = window.UNIVERSITY_LAB_DEMO_DATA;
  if (!data || !data.genes) return;
  const $ = id => document.getElementById(id);
  const fmt = x => Number(x).toLocaleString('en-US', {maximumFractionDigits:3});
  const names = {above:'stable above', sensitive:'reference-sensitive', below:'stable below'};
  const colors = {above:'#1f7a8c', sensitive:'#fff', below:'#b1b1b1'};
  let marker = 'SNX8', k = 3, selected = null;

  function stats(g, s) {
    const zs = g.leave_one_out.map(ref => (s.fc - ref.mean) / ref.sd);
    const n = zs.filter(z => z > k).length;
    return {z:(s.fc-g.mean)/g.sd, min:Math.min(...zs), max:Math.max(...zs), n,
      category:n===45?'above':n===0?'below':'sensitive'};
  }
  function element(tag, attrs, text) {
    const el = document.createElementNS('http://www.w3.org/2000/svg',tag);
    Object.entries(attrs || {}).forEach(([a,v])=>el.setAttribute(a,v));
    if (text!==undefined) el.textContent = text;
    return el;
  }
  function jitter(id) {
    let n=2166136261; for(const ch of id) n=Math.imul(n^ch.charCodeAt(0),16777619);
    return ((n>>>0)%1000)/1000-.5;
  }
  function showSample(g, id) {
    const s = [...g.cases,...g.controls].find(x=>x.id===id);
    if (!s) return;
    selected=id; $('sample-select').value=id;
    const isControl=g.controls.some(x=>x.id===id), st=stats(g,s);
    $('sample-detail').textContent = isControl
      ? `${id} · at-risk control · ${fmt(s.fc)} gene-level fold change. Its own held-out z score is ${fmt(s.heldout_z)} (calculated without this control). This is a calibration specimen, not an independent validation sample.`
      : `${id} · MIS-C · ${fmt(s.fc)} gene-level fold change · all-control z ${fmt(st.z)} · leave-one-out z ${fmt(st.min)}–${fmt(st.max)} · above ${st.n}/45 omitted-control cutoffs (${names[st.category]}). Treatment status and assay confirmation are unknown here.`;
    $('hero-viz').querySelectorAll('circle[data-sample]').forEach(p=>p.classList.toggle('selected',p.dataset.sample===id));
  }
  function draw(g, setting) {
    const svg=$('hero-viz'); svg.replaceChildren();
    const W=600, left=54, right=580, top=32, bottom=309;
    const max=Math.max(...g.cases.map(x=>x.fc),...g.controls.map(x=>x.fc),setting.loo_max);
    const yMax=Math.ceil(Math.log10(1+max)*2)/2;
    const y=v=>bottom-Math.log10(1+v)/yMax*(bottom-top);
    const title=element('title',{},`${g.gene}: ${setting.full_reference_above} of 199 MIS-C specimens exceed mean plus ${k} control SD. ${setting.sensitive} calls vary across omitted-control references.`);
    svg.append(title);
    for(let t=0;t<=yMax+1e-9;t+=.5){
      const yy=bottom-t/yMax*(bottom-top), val=10**t-1;
      svg.append(element('line',{x1:left,y1:yy,x2:right,y2:yy,stroke:'#e9e9e9'}));
      svg.append(element('text',{x:left-8,y:yy+4,'text-anchor':'end'},val<10?fmt(val):Math.round(val).toLocaleString('en-US')));
    }
    svg.append(element('text',{x:left,y:16},'Gene-level fold change / mock IP (log scale)'));
    svg.append(element('rect',{x:left,y:y(setting.loo_max),width:right-left,height:Math.max(1,y(setting.loo_min)-y(setting.loo_max)),fill:'#d7d7d7',opacity:.65}));
    svg.append(element('line',{x1:left,y1:y(setting.threshold),x2:right,y2:y(setting.threshold),stroke:'#111','stroke-dasharray':'5 4','stroke-width':1.2}));
    const groups=[{items:g.controls,center:165,label:'At-risk controls · 45',control:true},{items:g.cases,center:420,label:'MIS-C · 199',control:false}];
    for(const group of groups){
      svg.append(element('text',{x:group.center,y:339,'text-anchor':'middle'},group.label));
      for(const s of group.items){
        const st=stats(g,s);
        const circle=element('circle',{cx:group.center+jitter(s.id)*(group.control?100:220),cy:y(s.fc),r:group.control?3.2:3.5,
          fill:group.control?'#777':colors[st.category],stroke:group.control?'#fff':st.category==='sensitive'?'#111':'#fff','stroke-width':.8,
          opacity:group.control?.75:.85,'data-sample':s.id});
        circle.append(element('title',{},`${s.id}: ${fmt(s.fc)}; ${group.control?'calibration control':names[st.category]}`));
        circle.addEventListener('click',()=>{document.querySelector('.specimen-detail').open=true;showSample(g,s.id);});
        svg.append(circle);
      }
    }
    svg.setAttribute('aria-label',title.textContent);
  }
  function comparison() {
    const wrap=$('comparison');wrap.replaceChildren();
    for(const g of data.genes){
      const s=g.settings[String(k)]; const b=document.createElement('button'); b.type='button';
      b.setAttribute('aria-pressed',String(g.gene===marker));
      b.setAttribute('aria-label',`${g.gene}: ${s.stable_above} stable above, ${s.sensitive} reference-sensitive, ${s.stable_below} stable below. Select antigen.`);
      const name=document.createElement('span');name.textContent=g.gene;
      const track=document.createElement('span');track.className='comparison-track';
      for(const [cat,n] of [['above',s.stable_above],['sensitive',s.sensitive],['below',s.stable_below]]){
        const seg=document.createElement('span');seg.className=`comparison-segment ${cat}`;seg.style.width=`${n/199*100}%`;seg.title=`${n} ${names[cat]}`;track.append(seg);
      }
      const counts=document.createElement('span');counts.textContent=`${s.stable_above} / ${s.sensitive} / ${s.stable_below}`;
      b.append(name,track,counts);b.addEventListener('click',()=>{marker=g.gene;render();});wrap.append(b);
    }
  }
  function candidates(g) {
    const cat=$('category').value;
    const rows=g.cases.map(s=>({...s,...stats(g,s)})).filter(s=>s.category===cat).sort((a,b)=>Math.abs(a.z-k)-Math.abs(b.z-k)||a.id.localeCompare(b.id));
    $('candidate-count').textContent=`${g.gene} · ${rows.length} ${names[cat]} specimens at mean + ${k} SD. ${rows.length?'Showing '+Math.min(6,rows.length)+'.':'Choose another group to inspect examples.'}`;
    const tbody=$('candidates');tbody.replaceChildren();
    for(const s of rows.slice(0,6)){
      const tr=document.createElement('tr');
      [s.id,fmt(s.fc),fmt(s.z),`${fmt(s.min)}–${fmt(s.max)}`].forEach(v=>{const td=document.createElement('td');td.textContent=v;tr.append(td);});tbody.append(tr);
    }
  }
  function render() {
    const g=data.genes.find(x=>x.gene===marker), s=g.settings[String(k)];
    $('target-select').value=marker;draw(g,s);comparison();candidates(g);
    $('metric-one').textContent=`${s.full_reference_above}/199`;$('metric-two').textContent=s.sensitive;
    $('readout').textContent=`${marker}: the all-control cutoff is ${fmt(s.threshold)}. Omitting one control moves it from ${fmt(s.loo_min)} to ${fmt(s.loo_max)}. ${s.stable_above} MIS-C specimens stay above all 45 cutoffs; ${s.sensitive} change classification.${s.threshold<1?' The cutoff is below 1: above-reference does not mean confirmed enrichment or binding.':''}`;
    $('influence-note').textContent=`Largest reference shift: omit ${s.largest_cutoff_shift_id}, and the cutoff ${s.largest_cutoff_shift<0?'decreases':'increases'} by ${fmt(Math.abs(s.largest_cutoff_shift))} (${fmt(Math.abs(s.largest_cutoff_shift)/s.threshold*100)}%). This is a sensitivity check, not a reason to discard that control. Population SD (ddof = 0) matches the inspected upstream function.`;
    $('sample-select').replaceChildren();
    for(const [label,arr] of [['MIS-C',g.cases],['At-risk control',g.controls]]){
      const group=document.createElement('optgroup');group.label=label;
      for(const sample of arr){const opt=document.createElement('option');opt.value=sample.id;opt.textContent=sample.id;group.append(opt);}$('sample-select').append(group);
    }
    showSample(g,selected||g.cases.find(x=>stats(g,x).category==='sensitive')?.id||g.cases[0].id);
  }
  for(const gene of data.genes){const b=document.createElement('button');b.type='button';b.textContent=gene.gene;b.addEventListener('click',()=>{marker=gene.gene;render();});$('suggestions').append(b);}
  $('target-form').addEventListener('submit',e=>e.preventDefault());
  $('target-select').addEventListener('change',e=>{marker=e.target.value;render();});
  $('cutoff').addEventListener('change',e=>{k=Number(e.target.value);render();});
  $('category').addEventListener('change',()=>candidates(data.genes.find(x=>x.gene===marker)));
  $('sample-select').addEventListener('change',e=>showSample(data.genes.find(x=>x.gene===marker),e.target.value));
  render();
}());
