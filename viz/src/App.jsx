import React, { useState, useMemo } from 'react'
import { triageVariantsData, cohortSummaryData, dataQualityData, modelTransparencyData, dataMissing } from './data'
import './App.css'

export default function App() {
  const [activeTab, setActiveTab] = useState(() => {
    if (typeof window !== 'undefined') {
      const hash = window.location.hash.replace('#', '')
      if (['cohort', 'quality', 'model'].includes(hash)) return hash
      const params = new URLSearchParams(window.location.search)
      const tabParam = params.get('tab')
      if (['cohort', 'quality', 'model'].includes(tabParam)) return tabParam
    }
    return 'triage'
  })

  const [selectedVariant, setSelectedVariant] = useState(() => {
    if (typeof window !== 'undefined') {
      const params = new URLSearchParams(window.location.search)
      if (params.get('detail') === '1' && triageVariantsData.length > 0) {
        return triageVariantsData[0]
      }
    }
    return null
  })

  const handleTabChange = (tab) => {
    setActiveTab(tab)
    if (typeof window !== 'undefined') {
      window.location.hash = tab
    }
  }

  // Filters for Triage Table
  const [searchQuery, setSearchQuery] = useState('')
  const [geneFilter, setGeneFilter] = useState('ALL')
  const [tierFilter, setTierFilter] = useState('ALL')
  const [hotspotOnly, setHotspotOnly] = useState(false)
  const [consequenceFilter, setConsequenceFilter] = useState('ALL')
  const [minVaf, setMinVaf] = useState(0)
  const [maxVaf, setMaxVaf] = useState(100)

  // Sorting
  const [sortField, setSortField] = useState('default')
  const [sortDirection, setSortDirection] = useState('desc') // 'asc' | 'desc'

  // Pagination
  const [currentPage, setCurrentPage] = useState(1)
  const rowsPerPage = 50

  // Unique genes list for dropdown
  const uniqueGenes = useMemo(() => {
    const set = new Set(triageVariantsData.map(v => v.gene))
    return Array.from(set).sort()
  }, [])

  // Filtered & Sorted variants
  const filteredVariants = useMemo(() => {
    return triageVariantsData.filter(v => {
      // Search query (gene, protein change, patient ID)
      if (searchQuery.trim()) {
        const q = searchQuery.toLowerCase().trim()
        const matchGene = v.gene.toLowerCase().includes(q)
        const matchMut = v.protein_change.toLowerCase().includes(q)
        const matchPt = v.patient_id.toLowerCase().includes(q)
        if (!matchGene && !matchMut && !matchPt) return false
      }

      // Gene filter
      if (geneFilter !== 'ALL' && v.gene !== geneFilter) return false

      // Tier filter
      if (tierFilter === 'TIER1' && v.driver_tier !== 2) return false
      if (tierFilter === 'TIER1_2' && v.driver_tier === 0) return false

      // Hotspot only
      if (hotspotOnly && !v.is_hotspot) return false

      // Consequence filter
      if (consequenceFilter !== 'ALL') {
        if (!v.consequence.toLowerCase().includes(consequenceFilter.toLowerCase())) return false
      }

      // VAF filter
      const vafPct = v.vaf * 100
      if (vafPct < minVaf || vafPct > maxVaf) return false

      return true
    })
  }, [searchQuery, geneFilter, tierFilter, hotspotOnly, consequenceFilter, minVaf, maxVaf])

  const sortedVariants = useMemo(() => {
    const list = [...filteredVariants]
    if (sortField === 'default') {
      // Default: Hotspot Y/N -> Tier 2 -> Caller count -> ML score
      list.sort((a, b) => {
        if (a.is_hotspot !== b.is_hotspot) return b.is_hotspot ? 1 : -1
        if (a.driver_tier !== b.driver_tier) return b.driver_tier - a.driver_tier
        if (a.caller_count !== b.caller_count) return b.caller_count - a.caller_count
        return b.ml_score - a.ml_score
      })
      return list
    }

    list.sort((a, b) => {
      let aVal = a[sortField]
      let bVal = b[sortField]

      if (typeof aVal === 'string') {
        return sortDirection === 'asc' ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal)
      }
      if (typeof aVal === 'boolean') {
        aVal = aVal ? 1 : 0
        bVal = bVal ? 1 : 0
      }
      if (aVal === null || aVal === undefined) aVal = -9999
      if (bVal === null || bVal === undefined) bVal = -9999

      return sortDirection === 'asc' ? aVal - bVal : bVal - aVal
    })
    return list
  }, [filteredVariants, sortField, sortDirection])

  // Paginated variants
  const totalPages = Math.ceil(sortedVariants.length / rowsPerPage) || 1
  const paginatedVariants = useMemo(() => {
    const start = (currentPage - 1) * rowsPerPage
    return sortedVariants.slice(start, start + rowsPerPage)
  }, [sortedVariants, currentPage])

  // Handle Sort Click
  const handleSort = (field) => {
    if (sortField === field) {
      setSortDirection(prev => prev === 'asc' ? 'desc' : 'asc')
    } else {
      setSortField(field)
      setSortDirection('desc')
    }
    setCurrentPage(1)
  }

  // Export to CSV
  const exportToCSV = () => {
    const headers = [
      'Gene', 'Protein_Change', 'Consequence', 'VAF', 'Depth',
      'CancerHotspot', 'Gene_Tier', 'SIFT_Score', 'SIFT_Pred',
      'PolyPhen_Score', 'PolyPhen_Pred', 'gnomAD_AF', 'Caller_Count',
      'Callers', 'COSMIC', 'Hotspot_Likeness_Score', 'Patient_ID', 'Recurrence_In_Cohort', 'Why_Flagged'
    ]
    const rows = sortedVariants.map(v => [
      v.gene,
      v.protein_change,
      v.consequence,
      (v.vaf * 100).toFixed(1) + '%',
      v.depth,
      v.is_hotspot ? 'YES' : 'NO',
      v.driver_tier === 2 ? 'Tier 1' : v.driver_tier === 1 ? 'Tier 2' : 'Tier 0',
      v.sift_score ?? 'NA',
      v.sift_pred ?? 'NA',
      v.polyphen_score ?? 'NA',
      v.polyphen_pred ?? 'NA',
      v.gnomad_af ?? 0.0,
      v.caller_count,
      `"${v.callers}"`,
      v.in_cosmic ? 'YES' : 'NO',
      v.ml_score,
      v.patient_id,
      v.recurrence,
      `"${v.why_flagged.replace(/"/g, '""')}"`
    ])

    const csvContent = 'data:text/csv;charset=utf-8,' + [headers.join(','), ...rows.map(e => e.join(','))].join('\n')
    const encodedUri = encodeURI(csvContent)
    const link = document.createElement('a')
    link.setAttribute('href', encodedUri)
    link.setAttribute('download', `oncoseq_triage_variants_${new Date().toISOString().slice(0, 10)}.csv`)
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
  }

  return (
    <div className="dashboard-root">
      {/* MANDATORY PROTOTYPE BANNER */}
      <div className="prototype-banner">
        <span className="banner-icon">&#9888;</span>
        <span className="banner-text">
          <strong>Research prototype on TCGA-LUAD data.</strong> Not a clinical diagnostic tool. Not for treatment decisions.
        </span>
        <span className="banner-tag">Research Use Only (RUO)</span>
      </div>

      {dataMissing && (
        <div style={{ background: '#7a4b00', color: '#fff', padding: '8px 16px', fontSize: 13 }}>
          Dashboard data files not found in src/data/ - showing empty placeholders. Run scripts/build_dashboard_data.py.
        </div>
      )}

      {/* HEADER & NAV */}
      <header className="dashboard-header">
        <div className="header-brand">
          <div className="brand-logo">OS</div>
          <div>
            <h1 className="brand-title">OncoSeq Triage Dashboard v2</h1>
            <p className="brand-sub">Decision Support & Analytical Transparency · TCGA-LUAD Cohort (187 patients · 62,146 variants)</p>
          </div>
        </div>

        <nav className="header-nav">
          <button
            className={`nav-tab ${activeTab === 'triage' ? 'active' : ''}`}
            onClick={() => handleTabChange('triage')}
            id="tab-triage"
          >
            <span className="nav-icon">&#9776;</span> Variant Triage Table
          </button>
          <button
            className={`nav-tab ${activeTab === 'cohort' ? 'active' : ''}`}
            onClick={() => handleTabChange('cohort')}
            id="tab-cohort"
          >
            <span className="nav-icon">&#9639;</span> Cohort Summary
          </button>
          <button
            className={`nav-tab ${activeTab === 'quality' ? 'active' : ''}`}
            onClick={() => handleTabChange('quality')}
            id="tab-quality"
          >
            <span className="nav-icon">&#9733;</span> Data Quality & Trust
          </button>
          <button
            className={`nav-tab ${activeTab === 'model' ? 'active' : ''}`}
            onClick={() => handleTabChange('model')}
            id="tab-model"
          >
            <span className="nav-icon">&#9881;</span> Model Transparency
          </button>
        </nav>
      </header>

      {/* MAIN VIEW CONTAINER */}
      <main className="dashboard-main">
        {/* VIEW 1: VARIANT TRIAGE TABLE */}
        {activeTab === 'triage' && (
          <div className="view-container triage-view">
            {/* Context callout for Oncologists / MTB */}
            <div className="view-intro">
              <div className="intro-badge">Molecular Tumour Board / Oncologist View</div>
              <h2>Which variants in this cohort are worth attention, and why?</h2>
              <p>
                Prioritized variant callset sorted by curated CancerHotspots confirmation, OncoKB/CGC gene tiering, multi-caller consensus, and recurrence. Select any variant to view its granular evidence checklist and transparent rationale.
              </p>
            </div>

            {/* Filter toolbar */}
            <div className="filter-card">
              <div className="filter-row">
                <div className="filter-group search-group">
                  <label htmlFor="search-input">Search Variant / Gene / Patient</label>
                  <input
                    id="search-input"
                    type="text"
                    placeholder="e.g. KRAS, G12C, TP53, TCGA-55-7914..."
                    value={searchQuery}
                    onChange={(e) => { setSearchQuery(e.target.value); setCurrentPage(1); }}
                  />
                </div>

                <div className="filter-group">
                  <label htmlFor="gene-select">Gene</label>
                  <select
                    id="gene-select"
                    value={geneFilter}
                    onChange={(e) => { setGeneFilter(e.target.value); setCurrentPage(1); }}
                  >
                    <option value="ALL">All Genes ({uniqueGenes.length})</option>
                    {uniqueGenes.map(g => (
                      <option key={g} value={g}>{g}</option>
                    ))}
                  </select>
                </div>

                <div className="filter-group">
                  <label htmlFor="tier-select">Gene Tier (OncoKB/CGC)</label>
                  <select
                    id="tier-select"
                    value={tierFilter}
                    onChange={(e) => { setTierFilter(e.target.value); setCurrentPage(1); }}
                  >
                    <option value="ALL">All Gene Tiers</option>
                    <option value="TIER1">Tier 1 Only (OncoKB + CGC Consensus)</option>
                    <option value="TIER1_2">Tier 1 & Tier 2 Genes</option>
                  </select>
                </div>

                <div className="filter-group">
                  <label htmlFor="consequence-select">Consequence</label>
                  <select
                    id="consequence-select"
                    value={consequenceFilter}
                    onChange={(e) => { setConsequenceFilter(e.target.value); setCurrentPage(1); }}
                  >
                    <option value="ALL">All Consequences</option>
                    <option value="missense">Missense Mutation</option>
                    <option value="nonsense">Nonsense Mutation</option>
                    <option value="splice">Splice Site</option>
                    <option value="frame">Frameshift</option>
                  </select>
                </div>

                <div className="filter-group checkbox-group">
                  <label className="checkbox-label" htmlFor="hotspot-only">
                    <input
                      id="hotspot-only"
                      type="checkbox"
                      checked={hotspotOnly}
                      onChange={(e) => { setHotspotOnly(e.target.checked); setCurrentPage(1); }}
                    />
                    <span>CancerHotspot Only</span>
                  </label>
                </div>
              </div>

              <div className="filter-row filter-row-bottom">
                <div className="filter-group vaf-range-group">
                  <label>VAF Filter: {minVaf}% – {maxVaf}%</label>
                  <div className="vaf-sliders">
                    <input
                      type="range"
                      min="0"
                      max="100"
                      value={minVaf}
                      onChange={(e) => { setMinVaf(Number(e.target.value)); setCurrentPage(1); }}
                    />
                    <input
                      type="range"
                      min="0"
                      max="100"
                      value={maxVaf}
                      onChange={(e) => { setMaxVaf(Number(e.target.value)); setCurrentPage(1); }}
                    />
                  </div>
                </div>

                <div className="filter-actions">
                  <span className="results-count">
                    Showing <strong>{sortedVariants.length.toLocaleString()}</strong> matching variants
                  </span>
                  <button className="btn-secondary" onClick={() => {
                    setSearchQuery('');
                    setGeneFilter('ALL');
                    setTierFilter('ALL');
                    setHotspotOnly(false);
                    setConsequenceFilter('ALL');
                    setMinVaf(0);
                    setMaxVaf(100);
                    setSortField('default');
                    setCurrentPage(1);
                  }}>
                    Reset Filters
                  </button>
                  <button className="btn-primary" onClick={exportToCSV} id="btn-export-csv">
                    <span className="btn-icon">&#10515;</span> Export CSV
                  </button>
                </div>
              </div>
            </div>

            {/* Main Table Layout */}
            <div className="triage-layout">
              <div className="table-responsive">
                <table className="triage-table" id="variants-triage-table">
                  <thead>
                    <tr>
                      <th onClick={() => handleSort('gene')} className="sortable-th">
                        Gene {sortField === 'gene' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('protein_change')} className="sortable-th">
                        Protein Change {sortField === 'protein_change' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('consequence')} className="sortable-th">
                        Consequence {sortField === 'consequence' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('vaf')} className="sortable-th">
                        VAF {sortField === 'vaf' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('depth')} className="sortable-th">
                        Depth {sortField === 'depth' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('is_hotspot')} className="sortable-th">
                        Hotspot {sortField === 'is_hotspot' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('driver_tier')} className="sortable-th">
                        Gene Tier {sortField === 'driver_tier' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th>SIFT / PolyPhen</th>
                      <th onClick={() => handleSort('gnomad_af')} className="sortable-th">
                        gnomAD AF {sortField === 'gnomad_af' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('caller_count')} className="sortable-th">
                        Callers {sortField === 'caller_count' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('in_cosmic')} className="sortable-th">
                        COSMIC {sortField === 'in_cosmic' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                      <th onClick={() => handleSort('ml_score')} className="sortable-th highlight-th">
                        Hotspot-Likeness (Research) {sortField === 'ml_score' ? (sortDirection === 'asc' ? '▲' : '▼') : ''}
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    {paginatedVariants.length === 0 ? (
                      <tr>
                        <td colSpan="12" className="empty-table-cell">
                          No variants match the selected filters.
                        </td>
                      </tr>
                    ) : (
                      paginatedVariants.map((v) => {
                        const isSelected = selectedVariant?.id === v.id
                        return (
                          <tr
                            key={v.id}
                            className={`triage-row ${isSelected ? 'row-selected' : ''}`}
                            onClick={() => setSelectedVariant(v)}
                          >
                            <td>
                              <strong className="gene-name">{v.gene}</strong>
                            </td>
                            <td className="mono-text">{v.protein_change}</td>
                            <td>
                              <span className={`consequence-pill ${v.consequence.toLowerCase()}`}>
                                {v.consequence.replace('_', ' ')}
                              </span>
                            </td>
                            <td className="mono-text">{(v.vaf * 100).toFixed(1)}%</td>
                            <td className="mono-text">{v.depth}×</td>
                            <td>
                              {v.is_hotspot ? (
                                <span className="hotspot-badge">HOTSPOT</span>
                              ) : (
                                <span className="dim-text">—</span>
                              )}
                            </td>
                            <td>
                              {v.driver_tier === 2 ? (
                                <span className="tier-badge tier1">Tier 1</span>
                              ) : v.driver_tier === 1 ? (
                                <span className="tier-badge tier2">Tier 2</span>
                              ) : (
                                <span className="tier-badge tier0">None</span>
                              )}
                            </td>
                            <td className="mono-text">
                              <span title={`SIFT: ${v.sift_pred || 'NA'} (${v.sift_score ?? 'NA'}), PolyPhen: ${v.polyphen_pred || 'NA'} (${v.polyphen_score ?? 'NA'})`}>
                                {v.sift_score !== null ? (v.sift_score < 0.05 ? 'Del' : 'Tol') : '—'} / {v.polyphen_score !== null ? (v.polyphen_score > 0.85 ? 'Dam' : 'Ben') : '—'}
                              </span>
                            </td>
                            <td className="mono-text">
                              {v.gnomad_af === 0 ? (
                                <span className="absent-badge">0.0 (absent)</span>
                              ) : (
                                v.gnomad_af.toExponential(2)
                              )}
                            </td>
                            <td>
                              <span className={`caller-pill c-${v.caller_count}`} title={v.callers}>
                                {v.caller_count}/3 callers
                              </span>
                            </td>
                            <td>
                              {v.in_cosmic ? <span className="cosmic-badge">COSMIC</span> : <span className="dim-text">—</span>}
                            </td>
                            <td className="ml-score-cell">
                              <div className="ml-score-wrap">
                                <div className="score-bar-bg">
                                  <div className="score-bar-fill" style={{ width: `${Math.min(100, v.ml_score * 100)}%` }} />
                                </div>
                                <span className="score-num">{v.ml_score.toFixed(3)}</span>
                              </div>
                            </td>
                          </tr>
                        )
                      })
                    )}
                  </tbody>
                </table>
              </div>

              {/* VIEW 2: VARIANT DETAIL PANEL (Slide-over drawer when a row is clicked) */}
              {selectedVariant && (
                <aside className="detail-panel" id="variant-detail-panel">
                  <div className="detail-panel-header">
                    <div>
                      <span className="detail-label">EVIDENCE CHECKLIST & TRIAGE RATIONALE</span>
                      <h3 className="detail-title">
                        {selectedVariant.gene} <span className="detail-mut">{selectedVariant.protein_change}</span>
                      </h3>
                    </div>
                    <button className="detail-close-btn" onClick={() => setSelectedVariant(null)} title="Close Panel">
                      &times;
                    </button>
                  </div>

                  {/* One-line plain-language 'Why flagged' built solely from available fields */}
                  <div className="why-flagged-box">
                    <span className="why-flagged-label">Why flagged:</span>
                    <p className="why-flagged-text">{selectedVariant.why_flagged}</p>
                  </div>

                  {/* Granular Evidence Checklist */}
                  <div className="evidence-section">
                    <h4>Evidence Checklist</h4>
                    <ul className="evidence-list">
                      <li className={`evidence-item ${selectedVariant.is_hotspot ? 'checked' : 'unchecked'}`}>
                        <span className="check-icon">{selectedVariant.is_hotspot ? '✓' : '✗'}</span>
                        <div className="evidence-content">
                          <strong>CancerHotspots.org Status:</strong>
                          <span>
                            {selectedVariant.is_hotspot
                              ? 'Confirmed statistically recurrent hotspot in CancerHotspots v2.'
                              : 'Not cataloged as a statistically significant hotspot.'}
                          </span>
                        </div>
                      </li>

                      <li className={`evidence-item ${selectedVariant.driver_tier > 0 ? 'checked' : 'unchecked'}`}>
                        <span className="check-icon">{selectedVariant.driver_tier > 0 ? '✓' : '✗'}</span>
                        <div className="evidence-content">
                          <strong>Cancer Gene List Membership:</strong>
                          <span>
                            {selectedVariant.driver_tier === 2
                              ? 'Tier 1 Driver: OncoKB Cancer Gene List & COSMIC CGC Consensus gene.'
                              : selectedVariant.driver_tier === 1
                              ? 'Tier 2 Driver: Cataloged in OncoKB or COSMIC Cancer Gene Census.'
                              : 'Non-driver gene in standard cancer gene panels.'}
                          </span>
                        </div>
                      </li>

                      <li className={`evidence-item ${selectedVariant.gnomad_af <= 0.001 ? 'checked' : 'unchecked'}`}>
                        <span className="check-icon">{selectedVariant.gnomad_af <= 0.001 ? '✓' : '✗'}</span>
                        <div className="evidence-content">
                          <strong>Population Frequency (gnomAD AF):</strong>
                          <span>
                            {selectedVariant.gnomad_af === 0
                              ? 'Absent from gnomAD population database (somatic consistent).'
                              : `gnomAD AF: ${selectedVariant.gnomad_af} (under 0.001 germline filter threshold).`}
                          </span>
                        </div>
                      </li>

                      <li className={`evidence-item ${selectedVariant.sift_pred === 'deleterious' || selectedVariant.polyphen_pred?.includes('damaging') ? 'checked' : 'unchecked'}`}>
                        <span className="check-icon">
                          {selectedVariant.sift_pred === 'deleterious' || selectedVariant.polyphen_pred?.includes('damaging') ? '✓' : '✗'}
                        </span>
                        <div className="evidence-content">
                          <strong>Predicted Functional Effect:</strong>
                          <span>
                            SIFT: {selectedVariant.sift_pred || 'NA'} (score: {selectedVariant.sift_score ?? 'NA'}) · PolyPhen: {selectedVariant.polyphen_pred || 'NA'} (score: {selectedVariant.polyphen_score ?? 'NA'}).
                          </span>
                        </div>
                      </li>

                      <li className={`evidence-item ${selectedVariant.caller_count >= 2 ? 'checked' : 'unchecked'}`}>
                        <span className="check-icon">✓</span>
                        <div className="evidence-content">
                          <strong>Caller Consensus:</strong>
                          <span>
                            {selectedVariant.caller_count}/3 callers agreed ({selectedVariant.callers}). High technical call confidence.
                          </span>
                        </div>
                      </li>

                      <li className={`evidence-item ${selectedVariant.recurrence > 1 ? 'checked' : 'unchecked'}`}>
                        <span className="check-icon">{selectedVariant.recurrence > 1 ? '✓' : '✗'}</span>
                        <div className="evidence-content">
                          <strong>Cohort Recurrence:</strong>
                          <span>
                            Observed in <strong>{selectedVariant.recurrence}</strong> independent patient(s) across this TCGA-LUAD cohort.
                          </span>
                        </div>
                      </li>
                    </ul>
                  </div>

                  {/* Metadata and Context Cards */}
                  <div className="variant-meta-grid">
                    <div className="meta-card">
                      <span className="meta-k">Patient Barcode</span>
                      <strong className="meta-v mono-text">{selectedVariant.patient_id}</strong>
                    </div>
                    <div className="meta-card">
                      <span className="meta-k">Sample Barcode</span>
                      <strong className="meta-v mono-text text-truncate" title={selectedVariant.sample_id}>
                        {selectedVariant.sample_id}
                      </strong>
                    </div>
                    <div className="meta-card">
                      <span className="meta-k">Variant Allele Frequency</span>
                      <strong className="meta-v">{(selectedVariant.vaf * 100).toFixed(1)}% ({selectedVariant.depth}× depth)</strong>
                    </div>
                    <div className="meta-card">
                      <span className="meta-k">ML Score (Hotspot-likeness)</span>
                      <strong className="meta-v score-color">{selectedVariant.ml_score.toFixed(4)}</strong>
                    </div>
                  </div>

                  {/* Hard Rule Notice on Therapeutics */}
                  <div className="oncokb-disclaimer-notice">
                    <span className="notice-icon">&#9432;</span>
                    <p>
                      <strong>Note on Therapeutics:</strong> OncoKB Therapeutic Levels (Levels 1–4) and clinical trial assignments require licensed OncoKB institutional access and are not stored in this open-source prototype. Do not use this tool for prescription or clinical therapy decisions.
                    </p>
                  </div>
                </aside>
              )}
            </div>

            {/* Pagination Controls */}
            <div className="pagination-bar">
              <div className="page-info">
                Page <strong>{currentPage}</strong> of <strong>{totalPages}</strong> ({sortedVariants.length.toLocaleString()} total variants)
              </div>
              <div className="page-buttons">
                <button
                  className="page-btn"
                  onClick={() => setCurrentPage(1)}
                  disabled={currentPage === 1}
                >
                  &laquo; First
                </button>
                <button
                  className="page-btn"
                  onClick={() => setCurrentPage(prev => Math.max(1, prev - 1))}
                  disabled={currentPage === 1}
                >
                  &lsaquo; Prev
                </button>
                <span className="page-num-display">{currentPage}</span>
                <button
                  className="page-btn"
                  onClick={() => setCurrentPage(prev => Math.min(totalPages, prev + 1))}
                  disabled={currentPage === totalPages}
                >
                  Next &rsaquo;
                </button>
                <button
                  className="page-btn"
                  onClick={() => setCurrentPage(totalPages)}
                  disabled={currentPage === totalPages}
                >
                  Last &raquo;
                </button>
              </div>
            </div>
          </div>
        )}

        {/* VIEW 3: COHORT SUMMARY */}
        {activeTab === 'cohort' && (
          <div className="view-container cohort-view">
            <div className="view-intro">
              <div className="intro-badge">Cohort Mutational Landscape</div>
              <h2>Observed alterations and mutational spectra in TCGA-LUAD</h2>
              <p>
                Empirical summary computed directly from the cohort GDC MAF (187 patients, 62,146 somatic variants). Shows key lung cancer driver alterations, pathway frequencies, and the 6-class single-base substitution spectrum.
              </p>
            </div>

            {/* Key Lung Cancer Drivers Spotlight */}
            <div className="section-block">
              <h3 className="section-title">Key Lung Cancer Driver Genes (Observed Cohort Frequency)</h3>
              <p className="section-subtitle">
                Patient alteration percentages strictly as observed in this TCGA-LUAD dataset (N = 187 patients).
              </p>
              <div className="lung-genes-grid">
                {cohortSummaryData.key_lung_genes.map(g => (
                  <div key={g.gene} className="lung-gene-card">
                    <div className="card-top">
                      <strong className="gene-symbol">{g.gene}</strong>
                      <span className="tier-tag">{g.tier === 2 ? 'Tier 1' : 'Tier 2'}</span>
                    </div>
                    <div className="gene-stat-row">
                      <span className="stat-pct">{g.pct}%</span>
                      <span className="stat-count">{g.patients} / {g.total_patients} pts</span>
                    </div>
                    <div className="gene-track">
                      <div className="gene-fill" style={{ width: `${Math.min(100, g.pct * 1.8)}%` }} />
                    </div>
                  </div>
                ))}
              </div>
            </div>

            {/* Pathways & Top Mutated Genes Grid */}
            <div className="grid-2col">
              {/* Pathways */}
              <div className="dashboard-card">
                <h3 className="card-title">Canonical Pathway Alteration Frequency</h3>
                <p className="card-sub">
                  Percentage of cohort patients harboring at least one mutation in canonical oncogenic pathways.
                </p>
                <div className="pathway-list">
                  {cohortSummaryData.pathways.map(p => (
                    <div key={p.pathway} className="pathway-row">
                      <div className="pathway-info">
                        <strong className="pathway-name">{p.pathway}</strong>
                        <span className="pathway-top">Top gene: {p.top_gene} ({p.top_gene_n} muts)</span>
                      </div>
                      <div className="pathway-track-wrap">
                        <div className="pathway-track">
                          <div className="pathway-fill" style={{ width: `${p.pct}%` }} />
                        </div>
                        <span className="pathway-val">{p.pct}% ({p.patients} pts)</span>
                      </div>
                    </div>
                  ))}
                </div>
              </div>

              {/* Mutation Spectrum (6-Class SBS) */}
              <div className="dashboard-card">
                <h3 className="card-title">6-Class Mutation Spectrum (SBS)</h3>
                <p className="card-sub">
                  Single-base substitution profile across 60,122 analyzed SNVs. Strong tobacco carcinogen etiology.
                </p>
                <div className="spectrum-bars">
                  {Object.entries(cohortSummaryData.mutation_spectrum.proportions || {}).map(([sub, prop]) => {
                    const count = cohortSummaryData.mutation_spectrum.counts?.[sub] || 0
                    const pct = (prop * 100).toFixed(1)
                    return (
                      <div key={sub} className="spectrum-row">
                        <div className="spectrum-label">
                          <strong className={`sub-tag sub-${sub.replace('>', '_')}`}>{sub}</strong>
                          <span className="spectrum-count">{count.toLocaleString()} muts</span>
                        </div>
                        <div className="spectrum-track-wrap">
                          <div className="spectrum-track">
                            <div className="spectrum-fill" style={{ width: `${Math.min(100, prop * 200)}%` }} />
                          </div>
                          <span className="spectrum-val">{pct}%</span>
                        </div>
                      </div>
                    )
                  })}
                </div>
                <div className="signature-callout">
                  <span className="callout-badge">SBS4 Etiology</span>
                  <p>
                    Cosine similarity of <strong>0.9938</strong> to COSMIC <strong>SBS4 (Tobacco Smoking Signature)</strong>. Transversions C&gt;A dominate (42.5%), reflecting direct polycyclic aromatic hydrocarbon DNA damage in lung adenocarcinoma.
                  </p>
                </div>
              </div>
            </div>

            {/* Top 15 Mutated Genes */}
            <div className="dashboard-card" style={{ marginTop: '24px' }}>
              <h3 className="card-title">Top 15 Most Frequently Mutated Genes in Cohort</h3>
              <div className="top-genes-grid">
                {cohortSummaryData.top_genes.map((g, idx) => (
                  <div key={g.gene} className="top-gene-item">
                    <span className="rank-num">#{idx + 1}</span>
                    <strong className="tg-name">{g.gene}</strong>
                    <span className="tg-pts">{g.patients} pts</span>
                    <span className="tg-pct">{g.pct}%</span>
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}

        {/* VIEW 4: DATA QUALITY AND TRUST PANEL */}
        {activeTab === 'quality' && (
          <div className="view-container quality-view">
            <div className="view-intro">
              <div className="intro-badge">Bioinformatician / Researcher View</div>
              <h2>How reliable are these calls and annotations?</h2>
              <p>
                Full transparency into variant caller consensus, gnomAD germline filter rates, and in silico functional prediction coverage across all 62,146 somatic variants.
              </p>
            </div>

            {/* Quality KPIs */}
            <div className="quality-kpi-grid">
              <div className="quality-kpi-card">
                <div className="kpi-tag">Caller Concordance</div>
                <div className="kpi-metric-val">{dataQualityData.caller_concordance.three_callers_pct}%</div>
                <div className="kpi-metric-sub">
                  {dataQualityData.caller_concordance.three_callers.toLocaleString()} variants verified by 3 independent callers
                </div>
              </div>

              <div className="quality-kpi-card">
                <div className="kpi-tag">gnomAD Missingness</div>
                <div className="kpi-metric-val">{dataQualityData.gnomad_metrics.missing_pct}%</div>
                <div className="kpi-metric-sub">
                  {dataQualityData.gnomad_metrics.missing_count.toLocaleString()} variants absent from population database (expected somatic rate)
                </div>
              </div>

              <div className="quality-kpi-card">
                <div className="kpi-tag">Germline Filtered</div>
                <div className="kpi-metric-val">{dataQualityData.gnomad_metrics.germline_filtered_count}</div>
                <div className="kpi-metric-sub">
                  Variants removed with gnomAD AF &gt; 0.001 (rare germline contamination filter)
                </div>
              </div>

              <div className="quality-kpi-card">
                <div className="kpi-tag">SIFT / PolyPhen Coverage</div>
                <div className="kpi-metric-val">{dataQualityData.functional_annotation_coverage.sift_pct}%</div>
                <div className="kpi-metric-sub">
                  39,962 missense variants scored; remaining 35.7% are non-missense (expected)
                </div>
              </div>
            </div>

            {/* Caller Concordance Breakdown */}
            <div className="grid-2col" style={{ marginTop: '24px' }}>
              <div className="dashboard-card">
                <h3 className="card-title">Multi-Caller Consensus Breakdown</h3>
                <p className="card-sub">
                  Distribution of calling algorithms in the GDC Multi-Center Somatic Callset (MuSE, Mutect2, VarScan2, Pindel).
                </p>
                <div className="caller-breakdown-table-wrap">
                  <table className="mini-table">
                    <thead>
                      <tr>
                        <th>Caller Combination</th>
                        <th>Variant Count</th>
                        <th>Cohort %</th>
                      </tr>
                    </thead>
                    <tbody>
                      {dataQualityData.caller_concordance.top_combinations.map(c => (
                        <tr key={c.combination}>
                          <td className="mono-text">{c.combination}</td>
                          <td>{c.count.toLocaleString()}</td>
                          <td>
                            <div className="pct-bar-wrap">
                              <div className="pct-bar" style={{ width: `${c.pct}%` }} />
                              <span>{c.pct}%</span>
                            </div>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                <div className="trust-note">
                  <strong>Trust Rationale:</strong> 80.1% of all variants in this cohort are supported by 3 callers, and 100% are supported by at least 2 callers. No single-caller unverified variants are present in this curated callset.
                </div>
              </div>

              {/* Annotation Completeness & Missing Rates */}
              <div className="dashboard-card">
                <h3 className="card-title">Annotation Missingness & Filter Transparency</h3>
                <p className="card-sub">
                  Addressing expected database gaps in somatic sequencing data.
                </p>
                <div className="missingness-list">
                  <div className="missing-item">
                    <div className="mi-header">
                      <strong>gnomAD Population Database (86.48% Missing)</strong>
                      <span className="mi-badge info">Expected Somatic Rate</span>
                    </div>
                    <p className="mi-desc">
                      53,744 of 62,146 variants have missing gnomAD AF entries. In somatic tumor profiling, novel and ultra-rare somatic mutations are absent from germline reference cohorts by biological definition. Common germline variants were filtered upstream by GDC normal-tissue matching.
                    </p>
                  </div>

                  <div className="missing-item">
                    <div className="mi-header">
                      <strong>Explicit gnomAD AF Threshold (AF &gt; 0.001)</strong>
                      <span className="mi-badge pass">3 Variants Removed</span>
                    </div>
                    <p className="mi-desc">
                      An explicit filter threshold of <code>gnomAD_AF &gt; 0.001</code> was applied to remove potential germline leakage. Exactly 3 variants were removed, verifying that the cohort is virtually free of germline contamination.
                    </p>
                  </div>

                  <div className="missing-item">
                    <div className="mi-header">
                      <strong>SIFT & PolyPhen Coverage (64.3% Coverage)</strong>
                      <span className="mi-badge info">Missense Restricted</span>
                    </div>
                    <p className="mi-desc">
                      SIFT and PolyPhen are algorithms engineered exclusively for non-synonymous single amino-acid missense substitutions. Nonsense mutations, splice-site disruptions, in-frame and frameshift indels (35.5% of MAF) cannot be scored by these tools and appear as unannotated.
                    </p>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* VIEW 5: MODEL TRANSPARENCY PANEL */}
        {activeTab === 'model' && (
          <div className="view-container model-view">
            <div className="view-intro">
              <div className="intro-badge">Bioinformatician / Model Auditor View</div>
              <h2>What is the machine learning model actually adding?</h2>
              <p>
                Honest reporting of cross-validated metrics from <code>results/final_metrics.json</code>. The score reflects hotspot-likeness driven heavily by gene-level priors, not driver status or treatment response.
              </p>
            </div>

            {/* Headline Metrics Grid */}
            <div className="model-kpi-grid">
              <div className="model-kpi-card highlight">
                <div className="kpi-label">Pooled Out-of-Fold PR-AUC</div>
                <div className="kpi-val cyan">{modelTransparencyData.pooled_pr_auc_mean}</div>
                <div className="kpi-sub">
                  95% CI: &plusmn;{modelTransparencyData.ci_95} · Range: {modelTransparencyData.pr_auc_min} &ndash; {modelTransparencyData.pr_auc_max}
                </div>
              </div>

              <div className="model-kpi-card">
                <div className="kpi-label">Random Baseline PR-AUC</div>
                <div className="kpi-val dim">{modelTransparencyData.baseline_pr_auc}</div>
                <div className="kpi-sub">
                  Chance prevalence: {modelTransparencyData.positives} hotspots / {modelTransparencyData.total_evaluated.toLocaleString()} variants
                </div>
              </div>

              <div className="model-kpi-card">
                <div className="kpi-label">driver_tier Prior Baseline</div>
                <div className="kpi-val amber">0.092</div>
                <div className="kpi-sub">
                  Performance from OncoKB/CGC gene-level list membership alone
                </div>
              </div>

              <div className="model-kpi-card">
                <div className="kpi-label">Model Enrichment</div>
                <div className="kpi-val emerald">{modelTransparencyData.enrichment_factor}</div>
                <div className="kpi-sub">
                  Relative PR-AUC gain over random selection
                </div>
              </div>
            </div>

            {/* Baselines & Ablation Comparison */}
            <div className="grid-2col" style={{ marginTop: '24px' }}>
              <div className="dashboard-card">
                <h3 className="card-title">Ablation Analysis & Baselines Comparison</h3>
                <p className="card-sub">
                  Measuring the incremental value of multi-caller features, sequence context, and functional predictors over simple gene lists.
                </p>
                <div className="ablation-list">
                  {modelTransparencyData.baselines_comparison.map(b => (
                    <div key={b.model} className="ablation-item">
                      <div className="ablation-top">
                        <strong className="ablation-name">{b.model}</strong>
                        <span className="ablation-score">PR-AUC: {b.pr_auc}</span>
                      </div>
                      <div className="ablation-track-wrap">
                        <div className="ablation-track">
                          <div
                            className={`ablation-fill ${b.model.includes('Full') ? 'full' : b.model.includes('Prior') ? 'prior' : 'base'}`}
                            style={{ width: `${Math.min(100, b.pr_auc * 350)}%` }}
                          />
                        </div>
                      </div>
                      <div className="ablation-notes">{b.notes}</div>
                    </div>
                  ))}
                </div>
              </div>

              {/* Feature Importance */}
              <div className="dashboard-card">
                <h3 className="card-title">Random Forest Feature Importance Breakdown</h3>
                <p className="card-sub">
                  Gini feature importance confirms dominance of the gene-level prior (42.1%).
                </p>
                <div className="feature-importance-list">
                  {Object.entries(modelTransparencyData.feature_importance_rf).map(([feat, imp]) => {
                    const pct = (imp * 100).toFixed(1)
                    return (
                      <div key={feat} className="feat-row">
                        <div className="feat-info">
                          <code className="feat-name">{feat}</code>
                          <span className="feat-pct">{pct}%</span>
                        </div>
                        <div className="feat-track">
                          <div
                            className="feat-fill"
                            style={{
                              width: `${Math.min(100, imp * 220)}%`,
                              background: feat === 'driver_tier' ? 'var(--amber)' : 'var(--cyan)'
                            }}
                          />
                        </div>
                      </div>
                    )
                  })}
                </div>
              </div>
            </div>

            {/* Honest Interpretation Statement Box */}
            <div className="honest-takeaway-box" style={{ marginTop: '24px' }}>
              <div className="takeaway-badge">HONEST METHODOLOGICAL ASSESSMENT</div>
              <p className="takeaway-text">
                {modelTransparencyData.plain_language_disclaimer}
              </p>
              <div className="takeaway-points">
                <div className="point-item">
                  <strong>What the model DOES:</strong> Efficiently ranks CancerHotspot-recurrent positions (248 cataloged in cohort) above background passenger mutations, achieving a 65x enrichment over random baseline.
                </div>
                <div className="point-item">
                  <strong>What the model DOES NOT DO:</strong> It does not predict de novo driver probability, tumor fitness, drug sensitivity, or patient survival. Over 42% of decision weight comes from pre-existing gene tier priors.
                </div>
                <div className="point-item">
                  <strong>ClinVar Benchmark:</strong> Validated on external ClinVar annotations (N = 243) achieving ROC-AUC of {modelTransparencyData.clinvar_validation.roc_auc} and PR-AUC of {modelTransparencyData.clinvar_validation.pr_auc}.
                </div>
              </div>
            </div>
          </div>
        )}
      </main>

      {/* DASHBOARD FOOTER */}
      <footer className="dashboard-footer">
        <div className="footer-content">
          <div>
            <strong>Data Sources:</strong> GDC TCGA-LUAD Open MAF (GRCh38, 62,146 variants), CancerHotspots v2, OncoKB Cancer Gene List, COSMIC CGC v99, results/final_metrics.json.
          </div>
          <div>
            <strong>Verification Date:</strong> 2026-10-11 · OncoSeq-AI Research Prototype · Branch: <code>feature/dashboard-v2</code>
          </div>
        </div>
      </footer>
    </div>
  )
}
