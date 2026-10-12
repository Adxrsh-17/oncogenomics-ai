import { useState, Suspense } from 'react';
import Viewer3D from './Viewer3D';
import { kras_pred, egfr_pred, kras_crystal, egfr_crystal } from './pdbData';
import './App.css';

const PDB = { kras_pred, egfr_pred, kras_crystal, egfr_crystal };

const PROTEINS = {
  kras: {
    label: 'KRAS G12C', mutLabel: 'Gly→Cys at position 12', mutResidue: 12,
    drug: 'Sotorasib (AMG-510)', drugShort: 'Sotorasib',
    mechanism: 'Locks KRAS in GDP-bound inactive state via covalent bond to Cys-12, blocking RAS/MAPK oncogenic signalling.',
    approval: 'FDA-approved May 2021',  patients: 21, pdbId: '6OIM',
    accent: '#ffb347', accentDim: 'rgba(255,179,71,0.15)',
    predKey: 'kras_pred', crystalKey: 'kras_crystal',
  },
  egfr: {
    label: 'EGFR L858R', mutLabel: 'Leu→Arg at position 858', mutResidue: 858,
    drug: 'Osimertinib (AZD9291)', drugShort: 'Osimertinib',
    mechanism: 'L858R stabilises the active kinase conformation. Osimertinib binds irreversibly to Cys-797 in the ATP pocket.',
    approval: 'FDA-approved Nov 2015',  patients: 8,  pdbId: '4ZAU',
    accent: '#00e5a0', accentDim: 'rgba(0,229,160,0.15)',
    predKey: 'egfr_pred', crystalKey: 'egfr_crystal',
  },
};

export default function App() {
  const [protein,      setProtein]      = useState('kras');
  const [mode,         setMode]         = useState('gpu');
  const [style,        setStyle]        = useState('cartoon');
  const [showDrug,     setShowDrug]     = useState(true);
  const [showMutation, setShowMutation] = useState(true);
  const [spinning,     setSpinning]     = useState(true);

  const p = PROTEINS[protein];
  const pdbStr = mode === 'gpu' ? PDB[p.predKey] : PDB[p.crystalKey];

  return (
    <div className="app">
      {/* Header */}
      <header className="topbar">
        <div>
          <div className="topbar-title">LUAD MUTATION EXPLORER</div>
          <div className="topbar-sub">KRAS G12C · EGFR L858R · GPU structure prediction + drug-bound crystal</div>
        </div>
        <div className="nvidia-badge">⬛ NVIDIA BioNeMo · OpenFold2</div>
      </header>

      <div className="main">
        {/* Sidebar */}
        <aside className="sidebar">

          {/* Protein selector */}
          <section>
            <div className="section-label">Select Mutation</div>
            <div className="ptab-grid">
              {Object.entries(PROTEINS).map(([key, pr]) => (
                <button
                  key={key}
                  className={`ptab ${protein === key ? 'ptab-active' : ''}`}
                  style={protein === key ? { borderColor: pr.accent, background: pr.accentDim } : {}}
                  onClick={() => setProtein(key)}
                >
                  <div className="ptab-name" style={{ color: pr.accent }}>{pr.label}</div>
                  <div className="ptab-drug">→ {pr.drugShort}</div>
                  <div className="ptab-pts">{pr.patients} patients</div>
                </button>
              ))}
            </div>
          </section>

          {/* Mode */}
          <section>
            <div className="section-label">Structure Source</div>
            <div className="mode-grid">
              <button className={`mode-btn ${mode === 'gpu' ? 'mode-active' : ''}`} onClick={() => setMode('gpu')}>
                <span style={{ color: '#76b900', fontWeight: 700 }}>⬛ GPU Predicted</span>
                <span className="mode-sub">NVIDIA OpenFold2</span>
              </button>
              <button className={`mode-btn ${mode === 'crystal' ? 'mode-active' : ''}`} onClick={() => setMode('crystal')}>
                <span>💎 Crystal + Drug</span>
                <span className="mode-sub">PDB {p.pdbId}</span>
              </button>
            </div>
          </section>

          {/* Display style */}
          <section>
            <div className="section-label">Display Style</div>
            <div className="style-row">
              {['cartoon', 'surface', 'ball-stick'].map(s => (
                <button key={s} className={`style-btn ${style === s ? 'style-active' : ''}`} onClick={() => setStyle(s)}>
                  {s.replace('-', ' ')}
                </button>
              ))}
            </div>
          </section>

          {/* Toggles */}
          <section>
            <div className="section-label">Highlights</div>
            <div className="toggle-stack">
              <Toggle label={`Mutation site (res. ${p.mutResidue})`} checked={showMutation} onChange={setShowMutation} accent={p.accent} />
              {mode === 'crystal' && <Toggle label="Drug molecule" checked={showDrug} onChange={setShowDrug} accent={p.accent} />}
              <Toggle label="Auto-rotate" checked={spinning} onChange={setSpinning} accent={p.accent} />
            </div>
          </section>

          {/* Info card */}
          <section className="info-card">
            <div className="info-name" style={{ color: p.accent }}>{p.label}</div>
            <div className="info-mut">{p.mutLabel}</div>
            <div className="info-drug-row">
              <span className="drug-badge" style={{ background: p.accentDim, color: p.accent, borderColor: p.accent + '55' }}>
                {p.drug}
              </span>
              <span className="fda-badge">FDA ✓</span>
            </div>
            <p className="info-mech">{p.mechanism}</p>
            <div className="info-pts-row">
              <span className="info-pts-num" style={{ color: p.accent }}>{p.patients}</span>
              <span className="info-pts-label"> patients in TCGA-LUAD cohort</span>
            </div>
          </section>

          {/* Source note */}
          <div className="source-note">
            {mode === 'gpu'
              ? <><span style={{ color: '#76b900' }}>◆ GPU:</span> Predicted ab initio with NVIDIA OpenFold2 via BioNeMo hosted inference. Mutant sequence only — no ligand present.</>
              : <><span style={{ color: '#00d4ff' }}>◆ Crystal:</span> X-ray structure from RCSB PDB {p.pdbId}. HETATM records show bound drug molecule.</>
            }
          </div>
        </aside>

        {/* 3D Viewer */}
        <div className="viewer-col">
          <div className="viewer-badges">
            {mode === 'gpu'
              ? <><span className="badge badge-gpu">NVIDIA OPENFOLD2</span><span className="badge badge-pred">GPU PREDICTED</span></>
              : <><span className="badge badge-crystal">PDB: {p.pdbId}</span><span className="badge badge-xray">CRYSTAL + DRUG</span></>
            }
          </div>
          <Viewer3D
            key={`${protein}-${mode}`}
            pdbStr={pdbStr}
            accentColor={p.accent}
            mutResidue={p.mutResidue}
            mode={mode}
            style={style}
            showDrug={showDrug}
            showMutation={showMutation}
            spinning={spinning}
          />
        </div>
      </div>
    </div>
  );
}

function Toggle({ label, checked, onChange, accent }) {
  return (
    <label className="toggle-row">
      <span className="toggle-label">{label}</span>
      <div className="toggle-track" style={checked ? { background: accent + '30', borderColor: accent } : {}}>
        <input type="checkbox" checked={checked} onChange={e => onChange(e.target.checked)} style={{ display: 'none' }} />
        <div className="toggle-thumb" style={checked ? { transform: 'translateX(16px)', background: accent } : {}} />
      </div>
    </label>
  );
}
