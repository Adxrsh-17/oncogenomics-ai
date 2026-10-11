// Generated JSON (scripts/build_dashboard_data.py -> src/data/) is gitignored.
// Load it if present, otherwise fall back to empty-shaped defaults so the app builds and renders.
const files = import.meta.glob('./data/*.json', { eager: true, import: 'default' })
const load = (name, fallback) => files[`./data/${name}.json`] ?? fallback

export const triageVariantsData = load('triage_variants', [])
export const cohortSummaryData = load('cohort_summary', {
  key_lung_genes: [], pathways: [], top_genes: [],
  mutation_spectrum: { proportions: {}, counts: {} },
})
export const dataQualityData = load('data_quality', {
  caller_concordance: { three_callers_pct: 'n/a', three_callers: 0, top_combinations: [] },
  gnomad_metrics: { missing_pct: 'n/a', missing_count: 0, germline_filtered_count: 'n/a' },
  functional_annotation_coverage: { sift_pct: 'n/a' },
})
export const modelTransparencyData = load('model_transparency', {
  pooled_pr_auc_mean: 'n/a', ci_95: 'n/a', pr_auc_min: 'n/a', pr_auc_max: 'n/a',
  baseline_pr_auc: 'n/a', positives: 0, total_evaluated: 0, enrichment_factor: 'n/a',
  baselines_comparison: [], feature_importance_rf: {},
  plain_language_disclaimer: 'Data files missing - run scripts/build_dashboard_data.py.',
  clinvar_validation: { roc_auc: 'n/a', pr_auc: 'n/a' },
})
export const dataMissing = Object.keys(files).length < 4
