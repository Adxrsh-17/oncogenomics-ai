import { useState, useCallback } from 'react'
import { runPipeline, MODEL_METRICS } from './pipeline'
import Viewer3D from './Viewer3D'
import { kras_pred, egfr_pred, kras_crystal, egfr_crystal } from './pdbData'

const PDB = { kras_pred, egfr_pred, kras_crystal, egfr_crystal }

const PROTEIN_CONFIG = {
  'KRAS_p.G12C': {
    label: 'KRAS G12C', mutResidue: 12, accent: '#ffb347',
    predKey: 'kras_pred', crystalKey: 'kras_crystal',
    drug: 'Sotorasib (AMG-510)', pdbId: '6OIM',
  },
  'EGFR_p.L858R': {
    label: 'EGFR L858R', mutResidue: 858, accent: '#00e5a0',
    predKey: 'egfr_pred', crystalKey: 'egfr_crystal',
    drug: 'Osimertinib (AZD9291)', pdbId: '4ZAU',
  },
}

const SAMPLE_CSV = `patient_id,gene,hgvsp,classification,vaf,depth
TCGA-78-7145,KRAS,p.G12C,Missense_Mutation,0.38,245
TCGA-55-7914,EGFR,p.L858R,Missense_Mutation,0.45,312
TCGA-44-6147,TP53,p.R248W,Missense_Mutation,0.52,410
TCGA-44-6145,BRAF,p.V600E,Missense_Mutation,0.29,198
TCGA-97-7941,KRAS,p.G12D,Missense_Mutation,0.33,267
TCGA-55-8089,STK11,p.Q37*,Nonsense_Mutation,0.41,189
TCGA-50-5935,EGFR,p.T790M,Missense_Mutation,0.18,423
TCGA-86-8073,MET,p.D1010N,Missense_Mutation,0.22,156
TCGA-64-5779,RB1,p.E137*,Nonsense_Mutation,0.47,302
TCGA-44-6778,CDKN2A,p.R80*,Nonsense_Mutation,0.36,221`

const STEPS = [
  { key: 'load',     label: 'Parse Variants',      icon: '⬤' },
  { key: 'features', label: 'Extract Features',     icon: '⬤' },
  { key: 'ml',       label: 'ML Scoring (RF+XGB)',  icon: '⬤' },
  { key: 'action',   label: 'Actionability Lookup', icon: '⬤' },
  { key: 'done',     label: 'Generate Report',      icon: '⬤' },
]

function parseCSV(text) {
  const lines = text.trim().split('\n').filter(l => l.trim())
  if (lines.length < 2) return []
  const header = lines[0].split(',').map(h => h.trim())
  return lines.slice(1).map(line => {
    const vals = line.split(',').map(v => v.trim())
    const obj = {}
    header.forEach((h, i) => { obj[h] = vals[i] ?? '' })
    return {
      patient_id: obj.patient_id || 'UNKNOWN',
      gene: obj.gene || '',
      hgvsp: obj.hgvsp || '',
      classification: obj.classification || 'Missense_Mutation',
      vaf: parseFloat(obj.vaf) || 0.3,
      depth: parseInt(obj.depth) || 200,
    }
  }).filter(v => v.gene)
}

export default function DemoTab() {
  const [csvInput, setCsvInput]         = useState(SAMPLE_CSV)
  const [running, setRunning]           = useState(false)
  const [stepStates, setStepStates]     = useState({})
  const [stepLog, setStepLog]           = useState('')
  const [results, setResults]           = useState(null)
  const [selected3D, setSelected3D]     = useState(null)
  const [viewerMode, setViewerMode]     = useState('crystal')
  const [viewerStyle, setViewerStyle]   = useState('cartoon')

  const handleRun = useCallback(async () => {
    const variants = parseCSV(csvInput)
    if (!variants.length) return

    setRunning(true)
    setResults(null)
    setSelected3D(null)
    setStepStates({})
    setStepLog('')

    const output = await runPipeline(variants, ({ step, progress, message }) => {
      setStepStates(prev => ({ ...prev, [step]: { progress, done: progress === 100 } }))
      setStepLog(message)
    })

    setResults(output)
    setRunning(false)
  }, [csvInput])

  const handleVariantClick = (v) => {
    const key = `${v.gene}_${v.hgvsp}`
    if (PROTEIN_CONFIG[key]) {
      setSelected3D({ key, ...PROTEIN_CONFIG[key] })
      setViewerMode('crystal')
      setViewerStyle('cartoon')
    } else {
      setSelected3D(null)
    }
  }

  const stepStatus = (key) => {
    const s = stepStates[key]
    if (!s) return 'pending'
    if (s.done) return 'done'
    return 'running'
  }

  const actionableCount = results ? results.filter(v => v.drug).length : 0

  return (
    <div className="view-container demo-view">
      <div className="view-intro">
        <div className="intro-badge">Live Pipeline Demo · Input → Report in Seconds</div>
        <h2>Real-time Clinical Mutation Report</h2>
        <p>
          Paste a list of somatic mutations (CSV) and the pipeline parses, scores with purity-corrected ML (PR-AUC&nbsp;
          <strong>{MODEL_METRICS.pooled_pr_auc_mean}</strong> internal&nbsp;/&nbsp;
          <strong>{MODEL_METRICS.external_pr_auc}</strong> CPTAC external), and matches each actionable variant to FDA-approved therapies — automatically, in seconds.
        </p>
      </div>

      <div className="demo-grid">
        {/* Left: Input + Controls */}
        <div className="demo-input-col">
          <div className="dashboard-card demo-input-card">
            <div className="card-title-row">
              <h3 className="card-title">Mutation Input</h3>
              <button className="btn-ghost" onClick={() => setCsvInput(SAMPLE_CSV)}>Load sample</button>
            </div>
            <p className="card-sub">Format: patient_id, gene, hgvsp, classification, vaf, depth</p>
            <textarea
              className="demo-textarea"
              value={csvInput}
              onChange={e => setCsvInput(e.target.value)}
              spellCheck={false}
            />
            <button
              className={`btn-primary demo-run-btn ${running ? 'btn-running' : ''}`}
              onClick={handleRun}
              disabled={running}
            >
              {running ? (
                <><span className="spinner" /> Running pipeline…</>
              ) : (
                '▶  Run Pipeline'
              )}
            </button>
          </div>

          {/* Pipeline Steps */}
          <div className="dashboard-card demo-steps-card">
            <h3 className="card-title">Pipeline Progress</h3>
            <div className="steps-list">
              {STEPS.map(({ key, label }) => {
                const status = stepStatus(key)
                return (
                  <div key={key} className={`step-item step-${status}`}>
                    <div className="step-dot" />
                    <div className="step-body">
                      <span className="step-label">{label}</span>
                      {status === 'running' && (
                        <div className="step-bar-wrap">
                          <div
                            className="step-bar-fill"
                            style={{ width: `${stepStates[key]?.progress ?? 0}%` }}
                          />
                        </div>
                      )}
                    </div>
                    <span className="step-status-icon">
                      {status === 'done' ? '✓' : status === 'running' ? '…' : '○'}
                    </span>
                  </div>
                )
              })}
            </div>
            {stepLog && <div className="step-log">{stepLog}</div>}
          </div>

          {/* Model badge */}
          <div className="demo-model-badge">
            <span className="model-badge-icon">⬛</span>
            <div>
              <div className="model-badge-title">RF + XGBoost · Purity-Corrected CCF</div>
              <div className="model-badge-sub">
                Internal PR-AUC {MODEL_METRICS.pooled_pr_auc_mean} · CPTAC external {MODEL_METRICS.external_pr_auc} · {MODEL_METRICS.enrichment_factor}
              </div>
            </div>
          </div>
        </div>

        {/* Right: Results */}
        <div className="demo-results-col">
          {!results && !running && (
            <div className="demo-empty-state">
              <div className="empty-icon">⚕</div>
              <p>Enter mutations and click <strong>Run Pipeline</strong> to generate a clinical report.</p>
              <p className="empty-sub">The pipeline scores each variant for hotspot-likeness and matches actionable mutations to FDA-approved targeted therapies.</p>
            </div>
          )}

          {(results || running) && (
            <div className="demo-report-card dashboard-card">
              <div className="card-title-row">
                <h3 className="card-title">Clinical Report</h3>
                {results && (
                  <span className="report-summary-badge">
                    {results.length} variants · <span style={{ color: 'var(--emerald, #00e5a0)' }}>{actionableCount} actionable</span>
                  </span>
                )}
              </div>
              <div className="table-responsive">
                <table className="triage-table demo-report-table">
                  <thead>
                    <tr>
                      <th>Patient</th>
                      <th>Gene</th>
                      <th>Variant</th>
                      <th>VAF</th>
                      <th>ML Score</th>
                      <th>Drug / Status</th>
                    </tr>
                  </thead>
                  <tbody>
                    {results ? results.map((v, i) => {
                      const protKey = `${v.gene}_${v.hgvsp}`
                      const has3D = !!PROTEIN_CONFIG[protKey]
                      return (
                        <tr
                          key={i}
                          className={`triage-row ${v.drug ? 'row-actionable' : ''} ${has3D ? 'row-clickable-3d' : ''}`}
                          onClick={() => handleVariantClick(v)}
                          title={has3D ? 'Click to view 3D protein structure' : ''}
                        >
                          <td className="mono-text" style={{ fontSize: 11 }}>{v.patient_id}</td>
                          <td><strong className="gene-name">{v.gene}</strong></td>
                          <td className="mono-text">{v.hgvsp}</td>
                          <td className="mono-text">{(v.vaf * 100).toFixed(0)}%</td>
                          <td>
                            <div className="ml-score-wrap">
                              <div className="score-bar-bg">
                                <div className="score-bar-fill" style={{ width: `${Math.min(100, (v.ensemble_score ?? 0) * 100)}%` }} />
                              </div>
                              <span className="score-num">{(v.ensemble_score ?? 0).toFixed(3)}</span>
                            </div>
                          </td>
                          <td>
                            {v.drug ? (
                              <div className="drug-cell">
                                <span className="drug-badge-sm">{v.drug}</span>
                                <span className={`approval-tag ${v.approval?.includes('FDA') ? 'fda' : 'trial'}`}>
                                  {v.approval?.includes('FDA') ? 'FDA ✓' : 'Trial'}
                                </span>
                                {has3D && <span className="view3d-hint" title="View 3D structure">🔬</span>}
                              </div>
                            ) : (
                              <span className="dim-text">—</span>
                            )}
                          </td>
                        </tr>
                      )
                    }) : (
                      Array.from({ length: 5 }).map((_, i) => (
                        <tr key={i} className="skeleton-row">
                          {Array.from({ length: 6 }).map((_, j) => (
                            <td key={j}><div className="skeleton-cell" /></td>
                          ))}
                        </tr>
                      ))
                    )}
                  </tbody>
                </table>
              </div>

              {results && actionableCount === 0 && (
                <div className="no-action-note">
                  No actionable variants found in this sample. Try including KRAS p.G12C or EGFR p.L858R.
                </div>
              )}
            </div>
          )}
        </div>
      </div>

      {/* 3D Protein Viewer */}
      {selected3D && (
        <div className="demo-viewer-section dashboard-card">
          <div className="viewer-header-row">
            <div>
              <h3 className="card-title" style={{ color: selected3D.accent }}>
                {selected3D.label} — 3D Protein Structure
              </h3>
              <p className="card-sub">
                NVIDIA OpenFold2 GPU-predicted structure · Drug: {selected3D.drug} · PDB {selected3D.pdbId}
              </p>
            </div>
            <div className="viewer-controls">
              <div className="mode-grid" style={{ gap: 6 }}>
                <button
                  className={`mode-btn ${viewerMode === 'gpu' ? 'mode-active' : ''}`}
                  onClick={() => setViewerMode('gpu')}
                >
                  <span style={{ color: '#76b900', fontWeight: 700, fontSize: 12 }}>⬛ GPU Predicted</span>
                </button>
                <button
                  className={`mode-btn ${viewerMode === 'crystal' ? 'mode-active' : ''}`}
                  onClick={() => setViewerMode('crystal')}
                >
                  <span style={{ fontSize: 12 }}>💎 Crystal + Drug</span>
                </button>
              </div>
              <div className="style-row" style={{ marginTop: 6 }}>
                {['cartoon', 'surface', 'ball-stick'].map(s => (
                  <button
                    key={s}
                    className={`style-btn ${viewerStyle === s ? 'style-active' : ''}`}
                    onClick={() => setViewerStyle(s)}
                  >
                    {s.replace('-', ' ')}
                  </button>
                ))}
              </div>
            </div>
            <button className="detail-close-btn" onClick={() => setSelected3D(null)}>×</button>
          </div>
          <div style={{ height: 480 }}>
            <Viewer3D
              key={`${selected3D.key}-${viewerMode}`}
              pdbStr={viewerMode === 'gpu' ? PDB[selected3D.predKey] : PDB[selected3D.crystalKey]}
              accentColor={selected3D.accent}
              mutResidue={selected3D.mutResidue}
              mode={viewerMode}
              style={viewerStyle}
              showDrug={viewerMode === 'crystal'}
              showMutation={true}
              spinning={true}
            />
          </div>
          <div style={{ textAlign: 'center', marginTop: 8, color: 'var(--text-muted)', fontSize: 12 }}>
            Click table rows with 🔬 to switch protein · Powered by NVIDIA BioNeMo + OpenFold2
          </div>
        </div>
      )}
    </div>
  )
}
