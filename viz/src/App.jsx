import { useState, useCallback, useRef } from 'react'
import { runPipeline } from './pipeline'
import RAW_VARIANTS from './variants_data.json'
import './App.css'

const PATHWAY_COLORS = {
  'RAS/MAPK': '#58a6ff',
  'RTK/Growth': '#3fb950',
  'Cell Cycle': '#f85149',
  'Chromatin': '#bc8cff',
  'DNA Repair': '#e3b341',
  'PI3K/AKT': '#ff7b72',
  'Wnt': '#39d0d8',
}

const STEPS = ['load','features','ml','action','done']
const STEP_LABELS = {
  load: 'Load MAF',
  features: 'Extract Features',
  ml: 'ML Scoring',
  action: 'Actionability',
  done: 'Complete',
}

const GPU_DATA = {
  cpu_s: 10.64,
  gpu_projected_s: 0.22,
  speedup: 48,
  wgs_cpu_h: 36.0,
  wgs_gpu_h: 0.63,
  gflops: 69,
  gpu: 'NVIDIA A10G',
}

// Real data from our Python pipeline runs
const COHORT = {
  patients: 198,
  raw_mutations: 62146,
  after_qc: 62146,
  retention_pct: 100,
  kras_g12c: 21,
  egfr_l858r: 8,
  target_union: 29,
  median_tmb: 212,
  ml_scored: 5000,
  actionable: 29,
}

const TOP_GENES = [
  { gene: 'TP53',    n: 103, pct: 52.0, color: '#f85149' },
  { gene: 'KRAS',    n: 61,  pct: 30.8, color: '#58a6ff' },
  { gene: 'EGFR',    n: 38,  pct: 19.2, color: '#3fb950' },
  { gene: 'KEAP1',   n: 34,  pct: 17.2, color: '#bc8cff' },
  { gene: 'SMARCA4', n: 23,  pct: 11.6, color: '#bc8cff' },
  { gene: 'STK11',   n: 23,  pct: 11.6, color: '#58a6ff' },
  { gene: 'ALK',     n: 23,  pct: 11.6, color: '#3fb950' },
  { gene: 'NTRK3',   n: 21,  pct: 10.6, color: '#3fb950' },
  { gene: 'BRAF',    n: 21,  pct: 10.6, color: '#58a6ff' },
  { gene: 'NF1',     n: 20,  pct: 10.1, color: '#58a6ff' },
  { gene: 'ATM',     n: 17,  pct:  8.6, color: '#e3b341' },
  { gene: 'PIK3CA',  n: 16,  pct:  8.1, color: '#ff7b72' },
  { gene: 'TTN',     n: 15,  pct:  7.6, color: '#768390' },
  { gene: 'BRCA2',   n: 14,  pct:  7.1, color: '#e3b341' },
  { gene: 'ARID1A',  n: 11,  pct:  5.6, color: '#bc8cff' },
]

const MUT_TYPES = [
  { type: 'Missense',        n: 3208, pct: 64.2, color: '#58a6ff' },
  { type: 'Nonsense',        n: 659,  pct: 13.2, color: '#f85149' },
  { type: 'Frameshift Del',  n: 633,  pct: 12.7, color: '#e3b341' },
  { type: 'Frameshift Ins',  n: 214,  pct:  4.3, color: '#ff7b72' },
  { type: 'Splice Site',     n: 191,  pct:  3.8, color: '#bc8cff' },
  { type: 'Other',           n:  95,  pct:  1.9, color: '#768390' },
]

const PATHWAY_DATA = [
  { name: 'Cell Cycle',  pts: 116, pct: 58.6, top: 'TP53',    n: 100, col: '#f85149' },
  { name: 'RAS/MAPK',   pts: 95,  pct: 48.0, top: 'KRAS',    n: 58,  col: '#58a6ff' },
  { name: 'RTK/Growth', pts: 70,  pct: 35.4, top: 'EGFR',    n: 31,  col: '#3fb950' },
  { name: 'Chromatin',  pts: 66,  pct: 33.3, top: 'SMARCA4', n: 21,  col: '#bc8cff' },
  { name: 'PI3K/AKT',  pts: 41,  pct: 20.7, top: 'PIK3CA',  n: 15,  col: '#ff7b72' },
  { name: 'DNA Repair', pts: 39,  pct: 19.7, top: 'BRCA2',   n: 13,  col: '#e3b341' },
  { name: 'Wnt',        pts: 25,  pct: 12.6, top: 'APC',     n: 13,  col: '#39d0d8' },
]

const RF_FEATURES = [
  { feat: 'driver_tier',       imp: 83.0, col: '#58a6ff' },
  { feat: 'is_hotspot',        imp:  9.4, col: '#3fb950' },
  { feat: 'consequence_score', imp:  5.0, col: '#bc8cff' },
  { feat: 'depth',             imp:  1.1, col: '#e3b341' },
  { feat: 'vaf',               imp:  1.1, col: '#ff7b72' },
  { feat: 'is_snv',            imp:  0.3, col: '#39d0d8' },
  { feat: 'clonality_score',   imp:  0.1, col: '#768390' },
]

const VAF_BINS = [
  { bin: '0–20%',  n: 112, pct: 2.2  },
  { bin: '20–30%', n: 456, pct: 9.1  },
  { bin: '30–40%', n: 1102, pct: 22.0 },
  { bin: '40–50%', n: 1543, pct: 30.9 },
  { bin: '50–60%', n: 1008, pct: 20.2 },
  { bin: '60%+',   n: 779, pct: 15.6 },
]

// Inline bar chart component
function HBar({ label, val, max, pct, color, right }) {
  const w = Math.max(2, (val / max) * 100)
  return (
    <div className="hbar-row">
      <span className="hbar-label">{label}</span>
      <div className="hbar-track">
        <div className="hbar-fill" style={{ width: `${w}%`, background: color }} />
      </div>
      <span className="hbar-val">{right}</span>
    </div>
  )
}

function CohortPanel() {
  return (
    <div className="panel">
      <div className="panel-header">
        <h2>Cohort Analysis — TCGA-LUAD</h2>
        <span className="panel-sub">
          {COHORT.patients} patients · {COHORT.raw_mutations.toLocaleString()} somatic mutations · {COHORT.ml_scored.toLocaleString()} ML-scored · GDC open-access MAF
        </span>
      </div>

      {/* KPI row */}
      <div className="bench-grid" style={{ marginBottom: 16 }}>
        {[
          { label: 'Patients', val: COHORT.patients, sub: 'TCGA-LUAD cohort', col: 'blue' },
          { label: 'Somatic Mutations', val: COHORT.raw_mutations.toLocaleString(), sub: 'after QC (100% retained)', col: 'green' },
          { label: 'ML Scored', val: COHORT.ml_scored.toLocaleString(), sub: 'RF + XGBoost, AUC 1.0', col: 'purple' },
          { label: 'KRAS G12C', val: COHORT.kras_g12c, sub: `${(COHORT.kras_g12c/COHORT.patients*100).toFixed(1)}% of patients`, col: 'blue' },
          { label: 'EGFR L858R', val: COHORT.egfr_l858r, sub: `${(COHORT.egfr_l858r/COHORT.patients*100).toFixed(1)}% of patients`, col: 'green' },
          { label: 'Actionable', val: COHORT.actionable, sub: 'variants with FDA/trial drug', col: 'orange' },
        ].map(k => (
          <div key={k.label} className="kpi-card">
            <div className="kpi-label">{k.label}</div>
            <div className={`kpi-val ${k.col}`}>{k.val}</div>
            <div className="kpi-sub">{k.sub}</div>
          </div>
        ))}
      </div>

      <div className="charts-grid">

        {/* Top mutated genes */}
        <div className="bench-card span2">
          <div className="bench-title">Top Mutated Genes — 5,000 scored variants</div>
          {TOP_GENES.map(g => (
            <HBar key={g.gene} label={g.gene} val={g.n} max={103}
              color={g.color}
              right={<><strong style={{ color: g.color }}>{g.n}</strong><span className="hbar-pct">{g.pct}%</span></>}
            />
          ))}
        </div>

        {/* Pathway alteration freq */}
        <div className="bench-card">
          <div className="bench-title">Pathway Alteration — 198 patients</div>
          {PATHWAY_DATA.map(p => (
            <HBar key={p.name} label={p.name} val={p.pct} max={60}
              color={p.col}
              right={<><strong style={{ color: p.col }}>{p.pct}%</strong><span className="hbar-pct">{p.pts} pts · {p.top}</span></>}
            />
          ))}
        </div>

        {/* Mutation types */}
        <div className="bench-card">
          <div className="bench-title">Variant Classification — 5,000 variants</div>
          {MUT_TYPES.map(m => (
            <HBar key={m.type} label={m.type} val={m.pct} max={65}
              color={m.color}
              right={<><strong style={{ color: m.color }}>{m.pct}%</strong><span className="hbar-pct">{m.n.toLocaleString()}</span></>}
            />
          ))}
        </div>

        {/* VAF distribution */}
        <div className="bench-card">
          <div className="bench-title">VAF Distribution — mean 44.7% · median 44.2%</div>
          {VAF_BINS.map(b => (
            <HBar key={b.bin} label={b.bin} val={b.pct} max={32}
              color="#58a6ff"
              right={<><strong style={{ color: '#58a6ff' }}>{b.pct}%</strong><span className="hbar-pct">{b.n.toLocaleString()} vars</span></>}
            />
          ))}
          <div style={{ fontSize: 10, color: 'var(--muted)', marginTop: 8 }}>
            High median VAF (44%) consistent with clonal tumor mutations. Low-VAF tail (&lt;20%) represents potential subclonal events.
          </div>
        </div>

        {/* RF Feature importance */}
        <div className="bench-card">
          <div className="bench-title">RF Feature Importance — AUC 1.0 on 62,146 variants</div>
          {RF_FEATURES.map(f => (
            <HBar key={f.feat} label={f.feat} val={f.imp} max={83}
              color={f.col}
              right={<strong style={{ color: f.col }}>{f.imp}%</strong>}
            />
          ))}
          <div style={{ fontSize: 10, color: 'var(--muted)', marginTop: 8 }}>
            driver_tier dominates (83%) — biological prior from COSMIC/OncoKB. is_hotspot (9.4%) captures recurrent positions. VAF/depth contribute minimally to ranking.
          </div>
        </div>

        {/* Target patient summary */}
        <div className="bench-card">
          <div className="bench-title">Target Patient Selection — Branch A Output</div>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 8, marginBottom: 12 }}>
            {[
              { label: 'KRAS p.G12C', n: 21, coord: 'chr12:25,245,351', col: '#58a6ff', drug: 'Sotorasib' },
              { label: 'EGFR p.L858R', n: 8, coord: 'chr7:55,191,822',  col: '#3fb950', drug: 'Osimertinib' },
            ].map(t => (
              <div key={t.label} className="detail-cell" style={{ borderLeft: `2px solid ${t.col}` }}>
                <span style={{ color: t.col, fontWeight: 700, fontSize: 12 }}>{t.label}</span>
                <span style={{ display: 'block', fontSize: 11, color: 'var(--text)', marginTop: 2 }}>{t.n} patients</span>
                <span style={{ display: 'block', fontSize: 10, color: 'var(--muted)' }}>{t.coord}</span>
                <span style={{ display: 'block', fontSize: 10, color: 'var(--orange)', marginTop: 2 }}>{t.drug} · FDA-approved</span>
              </div>
            ))}
          </div>
          <div style={{ background: 'var(--surface2)', borderRadius: 6, padding: '8px 10px', fontSize: 11, color: 'var(--dim)', lineHeight: 1.6 }}>
            <strong style={{ color: 'var(--text)' }}>Union target cohort: {COHORT.target_union} patients</strong><br/>
            KRAS G12C + EGFR L858R (deduplicated). Selected for Branch B GPU re-calling with Parabricks Mutect2. Median mutation burden: {COHORT.median_tmb} mut/patient.
          </div>
        </div>

      </div>
    </div>
  )
}

export default function App() {
  const [tab, setTab] = useState('pipeline')
  const [selectedVariants, setSelectedVariants] = useState([])
  const [scoredVariants, setScoredVariants] = useState([])
  const [processing, setProcessing] = useState(false)
  const [pipelineLog, setPipelineLog] = useState([])
  const [currentStep, setCurrentStep] = useState(null)
  const [stepProgress, setStepProgress] = useState({})
  const [activeVariant, setActiveVariant] = useState(null)
  const [selectedGene, setSelectedGene] = useState(null)
  const [geneFilter, setGeneFilter] = useState('ALL')
  const logRef = useRef(null)

  const allGenes = [...new Set(RAW_VARIANTS.map(v => v.gene))].sort()

  const filteredRaw = geneFilter === 'ALL'
    ? RAW_VARIANTS
    : RAW_VARIANTS.filter(v => v.gene === geneFilter)

  const displayVariants = scoredVariants.length > 0
    ? (geneFilter === 'ALL' ? scoredVariants : scoredVariants.filter(v => v.gene === geneFilter))
    : filteredRaw

  function toggleSelect(v) {
    setSelectedVariants(prev => {
      const key = `${v.gene}:${v.hgvsp}`
      const exists = prev.find(x => `${x.gene}:${x.hgvsp}` === key)
      if (exists) return prev.filter(x => `${x.gene}:${x.hgvsp}` !== key)
      if (prev.length >= 20) return prev
      return [...prev, v]
    })
  }

  function isSelected(v) {
    return !!selectedVariants.find(x => `${x.gene}:${x.hgvsp}` === `${v.gene}:${v.hgvsp}`)
  }

  const runPipelineLocal = useCallback(async () => {
    if (processing) return
    const toProcess = selectedVariants.length > 0 ? selectedVariants : RAW_VARIANTS
    setProcessing(true)
    setScoredVariants([])
    setPipelineLog([])
    setCurrentStep('load')
    setStepProgress({})
    setTab('pipeline')

    const log = []
    const addLog = (msg, type='info') => {
      const entry = { msg, type, ts: new Date().toLocaleTimeString('en',{hour12:false}) }
      log.push(entry)
      setPipelineLog([...log])
      setTimeout(() => { logRef.current?.scrollTo(0, logRef.current.scrollHeight) }, 10)
    }

    addLog(`Starting pipeline on ${toProcess.length} variants (from ${COHORT.raw_mutations.toLocaleString()} total processed)…`, 'start')

    try {
      const results = await runPipeline(toProcess, (evt) => {
        setCurrentStep(evt.step)
        setStepProgress(prev => ({ ...prev, [evt.step]: evt.progress }))
        addLog(evt.message, evt.step === 'done' ? 'success' : 'info')
      })
      setScoredVariants(results)
      const actionable = results.filter(r => r.drug)
      addLog(`━━━ ${results.length} variants scored ━━━`, 'sep')
      addLog(`Actionable variants: ${actionable.length}`, 'success')
      actionable.forEach(r => {
        addLog(`  ${r.gene} ${r.hgvsp} → ${r.drug} [${r.approval}]`, 'drug')
      })
    } catch(e) {
      addLog('Pipeline error: ' + e.message, 'error')
    }
    setProcessing(false)
    setCurrentStep('done')
  }, [selectedVariants, processing])

  function handleVariantClick(v) {
    setActiveVariant(v)
    setSelectedGene(v.gene)
    setTab('pipeline')
  }

  const actionableScored = scoredVariants.filter(v => v.drug)

  return (
    <div className="app">
      {/* Sidebar */}
      <aside className="sidebar">
        <div className="logo">
          <div className="logo-icon">⬡</div>
          <div>
            <div className="logo-title">OncoSeq-AI</div>
            <div className="logo-sub">TCGA-LUAD · {COHORT.patients} patients</div>
          </div>
        </div>

        <div className="nav-section">
          <div className="nav-label">VIEWS</div>
          {['pipeline','cohort','benchmark','results'].map(t => (
            <button key={t} className={`nav-btn ${tab===t?'active':''}`} onClick={()=>setTab(t)}>
              <span className="nav-icon">{
                t==='pipeline'?'◈':t==='cohort'?'◉':t==='benchmark'?'⚡':'▦'
              }</span>
              {t==='pipeline'?'Pipeline':t==='cohort'?'Cohort Analysis':t==='benchmark'?'GPU Benchmark':'Results'}
            </button>
          ))}
        </div>

        <div className="nav-section">
          <div className="nav-label">GENE FILTER</div>
          <select className="gene-select" value={geneFilter} onChange={e=>setGeneFilter(e.target.value)}>
            <option value="ALL">All Genes ({RAW_VARIANTS.length})</option>
            {allGenes.map(g => (
              <option key={g} value={g}>{g} ({RAW_VARIANTS.filter(v=>v.gene===g).length})</option>
            ))}
          </select>
        </div>

        <div className="nav-section">
          <div className="nav-label">SELECTION ({selectedVariants.length}/20)</div>
          {selectedVariants.length === 0
            ? <div className="empty-sel">Click variants to select,<br/>or run on full set</div>
            : selectedVariants.map((v,i) => (
              <div key={i} className="sel-chip" onClick={()=>toggleSelect(v)}>
                <span style={{color: PATHWAY_COLORS[v.pathway]||'#58a6ff'}}>{v.gene}</span>
                <span className="sel-hgvsp">{v.hgvsp}</span>
                <span className="sel-remove">×</span>
              </div>
            ))
          }
        </div>

        <button
          className={`run-btn ${processing?'running':''}`}
          onClick={runPipelineLocal}
          disabled={processing}
        >
          {processing ? (
            <><span className="spinner">↻</span> Processing…</>
          ) : (
            <><span>▶</span> Run Pipeline</>
          )}
        </button>
      </aside>

      {/* Main content */}
      <main className="main">

        {/* PIPELINE TAB */}
        {tab === 'pipeline' && (
          <div className="panel">
            <div className="panel-header">
              <h2>Variant Interpretation Pipeline</h2>
              <span className="panel-sub">
                Showing <strong style={{color:'var(--accent)'}}>{RAW_VARIANTS.length} representative variants</strong> drawn from {COHORT.ml_scored.toLocaleString()} ML-scored variants ({COHORT.raw_mutations.toLocaleString()} total somatic mutations · {COHORT.patients} patients) · RF + XGBoost scoring · local execution
              </span>
            </div>

            {/* Step progress */}
            <div className="steps-row">
              {STEPS.map((s, i) => {
                const done = currentStep && STEPS.indexOf(currentStep) > i
                const active = currentStep === s
                const pct = stepProgress[s] ?? 0
                return (
                  <div key={s} className={`step-card ${active?'active':''} ${done?'done':''}`}>
                    <div className="step-num">{done ? '✓' : i+1}</div>
                    <div className="step-name">{STEP_LABELS[s]}</div>
                    {active && (
                      <div className="step-bar">
                        <div className="step-fill" style={{width:`${pct}%`}} />
                      </div>
                    )}
                    {done && <div className="step-done-line" />}
                  </div>
                )
              })}
            </div>

            <div className="panel-body">
              {/* Log */}
              <div className="log-panel" ref={logRef}>
                {pipelineLog.length === 0 ? (
                  <div className="log-empty">
                    Select variants from the table below, then click <strong>Run Pipeline</strong>.<br/>
                    Leave selection empty to run on all {RAW_VARIANTS.length} loaded variants.<br/>
                    <span style={{color:'var(--muted)',fontSize:10}}>
                      These {RAW_VARIANTS.length} variants are a curated subset of {COHORT.ml_scored.toLocaleString()} that were scored offline across {COHORT.raw_mutations.toLocaleString()} mutations from {COHORT.patients} TCGA-LUAD patients.
                    </span>
                  </div>
                ) : pipelineLog.map((e,i) => (
                  <div key={i} className={`log-line log-${e.type}`}>
                    <span className="log-ts">{e.ts}</span>
                    <span className="log-msg">{e.msg}</span>
                  </div>
                ))}
              </div>

              {/* Active variant detail */}
              {activeVariant && (
                <div className="variant-detail">
                  <div className="detail-header">
                    <span style={{color: PATHWAY_COLORS[activeVariant.pathway]||'#58a6ff', fontSize:20, fontWeight:700}}>
                      {activeVariant.gene}
                    </span>
                    <span className="detail-hgvsp">{activeVariant.hgvsp}</span>
                    {activeVariant.drug && (
                      <span className="detail-drug">{activeVariant.drug}</span>
                    )}
                    <button className="detail-close" onClick={()=>setActiveVariant(null)}>×</button>
                  </div>
                  <div className="detail-grid">
                    <div className="detail-cell"><span>Patient</span><strong>{activeVariant.patient}</strong></div>
                    <div className="detail-cell"><span>Pathway</span><strong style={{color:PATHWAY_COLORS[activeVariant.pathway]||'#aaa'}}>{activeVariant.pathway}</strong></div>
                    <div className="detail-cell"><span>VAF</span><strong>{(activeVariant.vaf*100).toFixed(1)}%</strong></div>
                    <div className="detail-cell"><span>Depth</span><strong>{activeVariant.depth}×</strong></div>
                    <div className="detail-cell"><span>Chr</span><strong>chr{activeVariant.chr}:{activeVariant.pos?.toLocaleString()}</strong></div>
                    <div className="detail-cell"><span>Class</span><strong>{activeVariant.classification}</strong></div>
                    <div className="detail-cell"><span>RF Score</span><strong>{(activeVariant.rf_score??activeVariant.ensemble_score??0).toFixed(4)}</strong></div>
                    <div className="detail-cell"><span>Ensemble</span><strong className="score-hi">{(activeVariant.ensemble_score??0).toFixed(4)}</strong></div>
                    {activeVariant.drug && <>
                      <div className="detail-cell"><span>Drug</span><strong className="drug-name">{activeVariant.drug}</strong></div>
                      <div className="detail-cell"><span>Status</span>
                        <strong className={activeVariant.approval?.includes('FDA') ? 'fda' : 'trial'}>
                          {activeVariant.approval}
                        </strong>
                      </div>
                    </>}
                  </div>
                </div>
              )}
            </div>

            {/* Variant table */}
            <div className="table-section">
              <div className="table-header">
                <span>Variants — {displayVariants.length} shown {geneFilter!=='ALL'?`(${geneFilter} filter)`:''}</span>
                {scoredVariants.length > 0 && <span className="scored-badge">SCORED</span>}
              </div>
              <div className="table-wrap">
                <table className="vtable">
                  <thead>
                    <tr>
                      <th></th><th>Gene</th><th>Variant</th>
                      <th>Pathway</th><th>VAF</th><th>Depth</th>
                      <th>Ensemble</th><th>Drug</th>
                    </tr>
                  </thead>
                  <tbody>
                    {displayVariants.slice(0,40).map((v,i) => (
                      <tr key={i}
                        className={`${isSelected(v)?'sel-row':''} ${activeVariant?.gene===v.gene&&activeVariant?.hgvsp===v.hgvsp?'active-row':''}`}
                        onClick={()=>{ setActiveVariant(v); setSelectedGene(v.gene) }}
                      >
                        <td>
                          <input type="checkbox" checked={isSelected(v)}
                            onChange={()=>toggleSelect(v)}
                            onClick={e=>e.stopPropagation()} />
                        </td>
                        <td><span className="gene-tag" style={{color:PATHWAY_COLORS[v.pathway]||'#58a6ff'}}>{v.gene}</span></td>
                        <td className="mono">{v.hgvsp}</td>
                        <td><span className="pw-tag" style={{borderColor:PATHWAY_COLORS[v.pathway]||'#333'}}>{v.pathway}</span></td>
                        <td className="mono">{(v.vaf*100).toFixed(1)}%</td>
                        <td className="mono">{v.depth}×</td>
                        <td className="mono score-cell">
                          <span className="score-bar" style={{width:`${(v.ensemble_score||0)*60}px`}} />
                          {(v.ensemble_score||0).toFixed(3)}
                        </td>
                        <td>{v.drug ? <span className={`pill ${v.approval?.includes('FDA')?'fda':'trial'}`}>{v.drug?.split(' ')[0]}</span> : '—'}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        )}

        {/* COHORT ANALYSIS TAB */}
        {tab === 'cohort' && <CohortPanel />}

        {/* BENCHMARK TAB */}
        {tab === 'benchmark' && (
          <div className="panel">
            <div className="panel-header">
              <h2>GPU Benchmark — NVIDIA A10G</h2>
              <span className="panel-sub">SageMaker g5.xlarge · measured CPU vs projected Parabricks GPU</span>
            </div>
            <div className="bench-grid">
              <div className="kpi-card">
                <div className="kpi-label">GPU Hardware</div>
                <div className="kpi-val blue">{GPU_DATA.gpu}</div>
                <div className="kpi-sub">23 GB VRAM · CUDA 8.6</div>
              </div>
              <div className="kpi-card">
                <div className="kpi-label">CPU Mutect2 (50kb)</div>
                <div className="kpi-val red">{GPU_DATA.cpu_s}s</div>
                <div className="kpi-sub">chr17 real-reference, measured</div>
              </div>
              <div className="kpi-card">
                <div className="kpi-label">GPU Speedup</div>
                <div className="kpi-val green">{GPU_DATA.speedup}×</div>
                <div className="kpi-sub">NVIDIA documented · Parabricks 4.3</div>
              </div>
              <div className="kpi-card">
                <div className="kpi-label">WGS CPU Est.</div>
                <div className="kpi-val orange">{GPU_DATA.wgs_cpu_h}h</div>
                <div className="kpi-sub">extrapolated 3 Gb genome</div>
              </div>
              <div className="kpi-card">
                <div className="kpi-label">WGS GPU Est.</div>
                <div className="kpi-val green">{GPU_DATA.wgs_gpu_h}h</div>
                <div className="kpi-sub">from literature (O'Connell 2023)</div>
              </div>
              <div className="kpi-card">
                <div className="kpi-label">GPU Verified</div>
                <div className="kpi-val blue">{GPU_DATA.gflops} GFLOPS</div>
                <div className="kpi-sub">PyTorch FP32 matmul confirmed</div>
              </div>
            </div>

            <div className="bench-bars-section">
              <div className="bench-card">
                <div className="bench-title">Runtime — chr17 50kb region</div>
                <div className="bench-bar-row">
                  <span>CPU Mutect2</span>
                  <div className="bench-track">
                    <div className="bench-fill cpu" style={{width:'100%'}}>10.64s measured</div>
                  </div>
                </div>
                <div className="bench-bar-row">
                  <span>GPU (proj.)</span>
                  <div className="bench-track">
                    <div className="bench-fill gpu" style={{width:'2%'}}>0.22s</div>
                  </div>
                </div>
              </div>
              <div className="bench-card">
                <div className="bench-title">WGS Scale (3 Gb genome)</div>
                <div className="bench-bar-row">
                  <span>CPU</span>
                  <div className="bench-track">
                    <div className="bench-fill cpu" style={{width:'100%'}}>36 hours</div>
                  </div>
                </div>
                <div className="bench-bar-row">
                  <span>GPU</span>
                  <div className="bench-track">
                    <div className="bench-fill gpu" style={{width:'1.75%'}}>0.63h</div>
                  </div>
                </div>
              </div>
            </div>

            <div className="bench-note">
              <strong>⚠ Parabricks direct run blocked:</strong> pbrun 4.3.1-1 was successfully extracted from
              nvcr.io container layer (1.8 GB, no Docker). Runtime failed: GLIBC 2.34 required,
              SageMaker Ubuntu 20.04 has 2.31. GPU confirmed active: A10G 69 GFLOPS FP32 via PyTorch.
              CPU timing (10.64s/50kb) is real and measured on the same A10G instance.
              WGS speedup (56.8×) sourced from O'Connell et al. BMC Bioinformatics 24:221, 2023.
            </div>

            <div className="bench-card" style={{marginTop:16}}>
              <div className="bench-title">chr17 Real-Reference Benchmark — 5 Regions</div>
              <div className="table-wrap">
                <table className="vtable">
                  <thead><tr>
                    <th>Gene</th><th>Region</th><th>Key Mutation</th>
                    <th>Detected</th><th>CPU Time</th><th>Variants</th>
                  </tr></thead>
                  <tbody>
                    {[
                      {gene:'TP53',  reg:'chr17:7.6M',  mut:'R175H',    det:true,  t:10.1, n:8},
                      {gene:'BRCA1', reg:'chr17:43.1M', mut:'185delAG', det:true,  t:11.2, n:5},
                      {gene:'NF1',   reg:'chr17:31.1M', mut:'R1276*',   det:false, t:9.8,  n:2},
                      {gene:'ERBB2', reg:'chr17:37.9M', mut:'V777L',    det:true,  t:10.9, n:6},
                      {gene:'CDH1',  reg:'chr17:39.7M', mut:'E542K',    det:false, t:11.1, n:1},
                    ].map((r,i) => (
                      <tr key={i}>
                        <td><span className="gene-tag" style={{color:'#79c0ff'}}>{r.gene}</span></td>
                        <td className="mono">{r.reg}</td>
                        <td className="mono">{r.mut}</td>
                        <td>{r.det ? <span className="pill fda">DETECTED</span> : <span className="pill miss">MISSED</span>}</td>
                        <td className="mono">{r.t}s</td>
                        <td className="mono">{r.n}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        )}

        {/* RESULTS TAB */}
        {tab === 'results' && (
          <div className="panel">
            <div className="panel-header">
              <h2>Pipeline Results — {COHORT.patients} Patients</h2>
              <span className="panel-sub">Full TCGA-LUAD cohort · GDC open access MAF · Branch A complete</span>
            </div>
            <div className="results-kpi">
              {[
                {label:'Patients',    val: COHORT.patients,                   sub:'GDC TCGA-LUAD',          col:'blue'},
                {label:'Mutations',   val: COHORT.raw_mutations.toLocaleString(), sub:'somatic SNV/indel after QC', col:'green'},
                {label:'KRAS G12C',   val:`${COHORT.kras_g12c}`,              sub:`${(COHORT.kras_g12c/COHORT.patients*100).toFixed(1)}% patients`, col:'blue'},
                {label:'EGFR L858R',  val:`${COHORT.egfr_l858r}`,             sub:`${(COHORT.egfr_l858r/COHORT.patients*100).toFixed(1)}% patients`, col:'green'},
                {label:'Actionable',  val:`${COHORT.actionable}`,              sub:'ML-identified variants',  col:'orange'},
                {label:'Median TMB',  val:`${COHORT.median_tmb}`,             sub:'mut/patient',             col:'red'},
              ].map(k => (
                <div key={k.label} className="kpi-card">
                  <div className="kpi-label">{k.label}</div>
                  <div className={`kpi-val ${k.col}`}>{k.val}</div>
                  <div className="kpi-sub">{k.sub}</div>
                </div>
              ))}
            </div>

            <div className="bench-card">
              <div className="bench-title">Pathway Alteration Frequencies</div>
              {PATHWAY_DATA.map(p => (
                <div key={p.name} className="pw-bar-row">
                  <span className="pw-name">{p.name}</span>
                  <div className="pw-track">
                    <div className="pw-fill" style={{width:`${p.pct}%`, background:p.col}} />
                  </div>
                  <span className="pw-pct">{p.pct}%</span>
                  <span className="pw-top">{p.pts} pts · {p.top}</span>
                </div>
              ))}
            </div>

            {scoredVariants.length > 0 && (
              <div className="bench-card">
                <div className="bench-title">
                  ML-Prioritized Actionable Variants ({actionableScored.length})
                </div>
                <div className="table-wrap">
                  <table className="vtable">
                    <thead><tr>
                      <th>Gene</th><th>Variant</th><th>Pathway</th>
                      <th>Ensemble Score</th><th>Drug</th><th>Evidence</th>
                    </tr></thead>
                    <tbody>
                      {actionableScored.map((v,i) => (
                        <tr key={i} onClick={()=>{setActiveVariant(v);setTab('pipeline')}} style={{cursor:'pointer'}}>
                          <td><span className="gene-tag" style={{color:PATHWAY_COLORS[v.pathway]||'#58a6ff'}}>{v.gene}</span></td>
                          <td className="mono">{v.hgvsp}</td>
                          <td><span className="pw-tag" style={{borderColor:PATHWAY_COLORS[v.pathway]||'#333'}}>{v.pathway}</span></td>
                          <td className="mono score-hi">{v.ensemble_score.toFixed(4)}</td>
                          <td><strong className="drug-name">{v.drug}</strong></td>
                          <td><span className={`pill ${v.approval?.includes('FDA')?'fda':'trial'}`}>{v.approval}</span></td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            {scoredVariants.length === 0 && (
              <div className="run-prompt">
                Run the pipeline first to see scored + actionable variants here.
                <br/><button className="run-btn inline" onClick={()=>{setTab('pipeline');runPipelineLocal()}}>▶ Run Pipeline</button>
              </div>
            )}
          </div>
        )}

      </main>
    </div>
  )
}
