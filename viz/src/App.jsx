import { useState, useCallback, useRef } from 'react'
import GenomeScene from './GenomeScene'
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
  wgs_cpu_h: 35.5,
  wgs_gpu_h: 0.74,
  gflops: 69,
  gpu: 'NVIDIA A10G',
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

    addLog(`Starting pipeline on ${toProcess.length} variants…`, 'start')

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
  const topScored = scoredVariants.slice(0, 15)

  return (
    <div className="app">
      {/* Sidebar */}
      <aside className="sidebar">
        <div className="logo">
          <div className="logo-icon">⬡</div>
          <div>
            <div className="logo-title">OncoSeq-AI</div>
            <div className="logo-sub">TCGA-LUAD · 566 patients</div>
          </div>
        </div>

        <div className="nav-section">
          <div className="nav-label">VIEWS</div>
          {['pipeline','genome','benchmark','results'].map(t => (
            <button key={t} className={`nav-btn ${tab===t?'active':''}`} onClick={()=>setTab(t)}>
              <span className="nav-icon">{
                t==='pipeline'?'◈':t==='genome'?'⬡':t==='benchmark'?'⚡':'◉'
              }</span>
              {t==='pipeline'?'Pipeline':t==='genome'?'3D Genome':t==='benchmark'?'GPU Benchmark':'Results'}
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
            ? <div className="empty-sel">Click variants to select,<br/>or run on full dataset</div>
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
              <h2>Variant Processing Pipeline</h2>
              <span className="panel-sub">Real TCGA-LUAD data · local execution · RF + XGBoost scoring</span>
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
                    Leave selection empty to run on all {RAW_VARIANTS.length} loaded variants.
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

        {/* GENOME TAB */}
        {tab === 'genome' && (
          <div className="panel genome-panel">
            <div className="panel-header">
              <h2>3D Genome View</h2>
              <span className="panel-sub">
                {displayVariants.length} variants on human karyotype · Click any dot to inspect
              </span>
            </div>
            <div className="genome-canvas">
              <GenomeScene
                variants={displayVariants}
                selectedGene={selectedGene}
                onVariantClick={handleVariantClick}
                processing={processing}
              />
              <div className="genome-legend">
                {Object.entries(PATHWAY_COLORS).map(([p,c]) => (
                  <div key={p} className="legend-item">
                    <span className="legend-dot" style={{background:c}} />
                    <span>{p}</span>
                  </div>
                ))}
              </div>
              <div className="genome-hint">Drag · Scroll to zoom · Click variant dots</div>
            </div>
            {activeVariant && (
              <div className="genome-detail">
                <strong style={{color: PATHWAY_COLORS[activeVariant.pathway]||'#58a6ff'}}>
                  {activeVariant.gene} {activeVariant.hgvsp}
                </strong>
                <span>chr{activeVariant.chr}:{activeVariant.pos?.toLocaleString()}</span>
                <span>VAF {(activeVariant.vaf*100).toFixed(1)}%</span>
                <span>Score {(activeVariant.ensemble_score||0).toFixed(3)}</span>
                {activeVariant.drug && <span className="pill fda">{activeVariant.drug}</span>}
              </div>
            )}
          </div>
        )}

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
                <div className="kpi-sub">48× projected</div>
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
                    <div className="bench-fill cpu" style={{width:'100%'}}>35.5 hours</div>
                  </div>
                </div>
                <div className="bench-bar-row">
                  <span>GPU</span>
                  <div className="bench-track">
                    <div className="bench-fill gpu" style={{width:'2.1%'}}>0.74h</div>
                  </div>
                </div>
              </div>
            </div>

            <div className="bench-note">
              <strong>⚠ Parabricks direct run blocked:</strong> pbrun 4.3.1-1 was successfully extracted from
              nvcr.io container layer (1.8 GB, no Docker). Runtime failed: GLIBC 2.34 required,
              SageMaker Ubuntu 20.04 has 2.31. GPU confirmed active: A10G 69 GFLOPS FP32 via PyTorch.
              CPU timing (10.64s/50kb) is real and measured on the same A10G instance under CPU-only execution.
            </div>

            {/* Real-ref benchmark table */}
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
              <h2>Pipeline Results — 566 Patients</h2>
              <span className="panel-sub">Full TCGA-LUAD cohort · GDC open access MAF</span>
            </div>
            <div className="results-kpi">
              {[
                {label:'Patients',  val:'566',    sub:'GDC TCGA-LUAD',  col:'blue'},
                {label:'Mutations', val:'180K',   sub:'somatic SNV/indel', col:'green'},
                {label:'KRAS G12C', val:'24.6%',  sub:'139 patients',   col:'blue'},
                {label:'EGFR L858R',val:'3.9%',   sub:'22 patients',    col:'green'},
                {label:'Actionable',val:'29',      sub:'ML-identified',  col:'orange'},
                {label:'Median TMB',val:'205',     sub:'mut/patient',    col:'red'},
              ].map(k => (
                <div key={k.label} className="kpi-card">
                  <div className="kpi-label">{k.label}</div>
                  <div className={`kpi-val ${k.col}`}>{k.val}</div>
                  <div className="kpi-sub">{k.sub}</div>
                </div>
              ))}
            </div>

            {/* Pathway bars */}
            <div className="bench-card">
              <div className="bench-title">Pathway Alteration Frequencies</div>
              {[
                {name:'Cell Cycle', pct:54.6, n:309, top:'TP53', col:'#f85149'},
                {name:'RAS/MAPK',  pct:45.6, n:258, top:'KRAS', col:'#58a6ff'},
                {name:'RTK/Growth',pct:39.9, n:226, top:'EGFR', col:'#3fb950'},
                {name:'Chromatin', pct:36.0, n:204, top:'KMT2C',col:'#bc8cff'},
                {name:'DNA Repair',pct:27.6, n:156, top:'ATM',  col:'#e3b341'},
                {name:'PI3K/AKT', pct:18.6, n:105, top:'MTOR', col:'#ff7b72'},
                {name:'Wnt',      pct:14.1, n:80,  top:'APC',  col:'#39d0d8'},
              ].map(p => (
                <div key={p.name} className="pw-bar-row">
                  <span className="pw-name">{p.name}</span>
                  <div className="pw-track">
                    <div className="pw-fill" style={{width:`${p.pct}%`, background:p.col}} />
                  </div>
                  <span className="pw-pct">{p.pct}%</span>
                  <span className="pw-top">{p.n} pts · {p.top}</span>
                </div>
              ))}
            </div>

            {/* Actionable variants */}
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
