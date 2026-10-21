import { useEffect, useRef, useState } from 'react'
import {
  ArrowRight,
  CheckCircle2,
  Clock3,
  Download,
  ExternalLink,
  MessageSquare,
  Moon,
  Pencil,
  Play,
  Plus,
  Sun,
  Trash2,
  XCircle,
} from 'lucide-react'
import FusionStructureViewer from './FusionStructureViewer'
import {
  buildTemplateCandidates,
  fusionFasta,
  parseProteinSequence,
  summarizePdbConfidence,
  validateFusionPdb,
  type FusionCandidate,
} from './fusion'
import {
  errorMessage,
  foldStartTimeFromRunDir,
  formatElapsedTime,
  getComputeStatus,
  getFoldStatus,
  inDesktop,
  openPdbInPymol,
  startLocalFold,
  type ComputeStatus,
  type FoldStatus,
} from './desktopFold'
import {
  auditFusionCoordinates,
  assessProtectedSiteClearance,
  calculateDomainPreservation,
  parseColabFoldScoreJson,
  parseProtectedSites,
  validateDomainReferencePdb,
  type ColabFoldMetrics,
} from './structureEvidence'
import {
  appendHistoryEntry,
  checksumManifest,
  createStoredZip,
  sha256Hex,
  verifyHistory,
  type HistoryEntry,
  type ResearchBundleFile,
} from './researchBundle'

type RigidityMode = 'flexible' | 'rigid'
type AppPage = 'candidates' | 'structure' | 'results'
type ColorTheme = 'light' | 'dark'

const appPages: Array<{ id: AppPage; label: string; step: string }> = [
  { id: 'candidates', label: 'Candidates', step: '01' },
  { id: 'structure', label: 'Prediction', step: '02' },
  { id: 'results', label: 'Results', step: '03' },
]

function pageFromHash(hash: string): AppPage {
  const route = hash.replace(/^#/, '')
  if (route === 'runs' || route === 'candidates') return 'candidates'
  if (route === 'viewer' || route === 'structure') return 'structure'
  if (route === 'provenance' || route === 'results') return 'results'
  return 'structure'
}

interface RankedCandidateRow {
  rank: number
  candidate_id: string
  method_name: string
  sequence: string
  length: number
  compatibility_score: number
  hard_rejection_pass: boolean
  provenance: string
  rejection_reasons: string[]
  contour_length_A: number
}

interface ApiCandidateRow extends Omit<RankedCandidateRow, 'rejection_reasons' | 'contour_length_A'> {
  scorecard?: {
    rejection_reasons?: string[]
    contour_length_A?: number
  }
}

interface ViewerResult {
  fold: FoldStatus | null
  pdb: string | null
  source: 'colabfold' | 'uploaded' | null
  metrics: ColabFoldMetrics | null
  elapsed: string | null
  running: boolean
}

interface DesignSession {
  id: string
  updatedAt: string
  projectName: string
  domainALabel: string
  domainBLabel: string
  domainAInput: string
  domainBInput: string
  biologicalObjective: string
  expressionSystem: string
  protectedSites: string
  rigidity: RigidityMode
  cleavable: boolean
  targetLength: number
}

interface DomainReference {
  sourceId: string
  sourceKind: 'experimental' | 'predicted'
  filename: string | null
  pdb: string | null
}

interface DomainReferences {
  A: DomainReference
  B: DomainReference
}

interface SavedProjectV1 {
  schema_version: 'linkerforge.project.v1'
  project: {
    name: string
    biological_objective: string
    expression_system: string
    protected_sites: string
  }
  domains: {
    A: { label: string; sequence: string; reference_source_id: string; reference_source_kind: DomainReference['sourceKind'] }
    B: { label: string; sequence: string; reference_source_id: string; reference_source_kind: DomainReference['sourceKind'] }
  }
  design: { rigidity: RigidityMode; cleavable: boolean; target_length: number }
  privacy: { public_msa_consent_saved: false }
}

const PRESET_CANDIDATES: Record<string, RankedCandidateRow[]> = {
  flexible_non_cleavable: [
    {
      rank: 1,
      candidate_id: 'lf_template_1',
      method_name: 'linkerforge_ranker',
      sequence: 'GGGGSGGGGSGGGGS',
      length: 15,
      compatibility_score: 10.3379,
      hard_rejection_pass: true,
      provenance: 'Part 2 Template (Chen2013 Flexible GGGGS)',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 2,
      candidate_id: 'lf_authority_5',
      method_name: 'nearest_neighbour_retrieval',
      sequence: 'ARGGQYSNQGVVAST',
      length: 15,
      compatibility_score: 10.1476,
      hard_rejection_pass: true,
      provenance: 'X2 Authority (mined_swissprot_gap:Q9XTD6:272-286, mined_high_confidence_dfl_call)',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 3,
      candidate_id: 'lf_template_2',
      method_name: 'linkerforge_ranker',
      sequence: 'GSAGSAAGSGEFGSA',
      length: 15,
      compatibility_score: 10.1155,
      hard_rejection_pass: true,
      provenance: 'Part 2 Template (Waldo GFP-Compatible GSAGSAAGSGEF)',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 4,
      candidate_id: 'baseline_random_gs',
      method_name: 'random_length_matched_gs',
      sequence: 'GSGGSGSGGSGSGGS',
      length: 15,
      compatibility_score: 10.0750,
      hard_rejection_pass: true,
      provenance: 'Baseline 3: Random length-matched GS',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 5,
      candidate_id: 'baseline_geometry_only',
      method_name: 'geometry_only',
      sequence: 'AAAGGGAAAGGGAAA',
      length: 15,
      compatibility_score: 9.3547,
      hard_rejection_pass: true,
      provenance: 'Baseline 4: Geometry-only length match',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 6,
      candidate_id: 'synthetic_neg_rigid_helical',
      method_name: 'synthetic_negative:rigid_helical_where_flexible_needed',
      sequence: 'EAAAKEAAAKEAAAK',
      length: 15,
      compatibility_score: 8.1518,
      hard_rejection_pass: true,
      provenance: 'Synthetic Negative (rigid_helical_where_flexible_needed)',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 7,
      candidate_id: 'synthetic_neg_charge_clash',
      method_name: 'synthetic_negative:charge_clash_with_junction_anchors',
      sequence: 'KKKKKKKKKKKKKKK',
      length: 15,
      compatibility_score: 5.5658,
      hard_rejection_pass: true,
      provenance: 'Synthetic Negative (charge_clash_with_junction_anchors)',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 8,
      candidate_id: 'synthetic_neg_too_short',
      method_name: 'synthetic_negative:too_short_by_50pct',
      sequence: 'GGGGSGG',
      length: 7,
      compatibility_score: 5.1972,
      hard_rejection_pass: false,
      provenance: 'Synthetic Negative (too_short_by_50pct)',
      rejection_reasons: ['physical_span_closure_failed(contour=26.6A < required=28.4A, L=7)'],
      contour_length_A: 26.6,
    },
    {
      rank: 9,
      candidate_id: 'synthetic_neg_too_long',
      method_name: 'synthetic_negative:too_long_by_200pct',
      sequence: 'GGGGSGGGGSGGGGSGGGGSGGGGSGGGGSGGGGSGGGGSGGGGS',
      length: 45,
      compatibility_score: 4.5453,
      hard_rejection_pass: false,
      provenance: 'Synthetic Negative (too_long_by_200pct)',
      rejection_reasons: ['domain_clash_or_entropic_excess(L=45, target=15, detour=1.15)'],
      contour_length_A: 171.0,
    },
    {
      rank: 10,
      candidate_id: 'synthetic_neg_buried_protease',
      method_name: 'synthetic_negative:protease_site_buried_in_rigid_helix',
      sequence: 'EAAAKEAAAKENLYFQSEAAAKEAAAK',
      length: 27,
      compatibility_score: 3.6214,
      hard_rejection_pass: false,
      provenance: 'Synthetic Negative (protease_site_buried_in_rigid_helix)',
      rejection_reasons: [
        'unintended_protease_cleavage_motif_present',
        'protease_site_buried_inside_rigid_alpha_helix_without_gs_spacers',
      ],
      contour_length_A: 102.6,
    },
    {
      rank: 11,
      candidate_id: 'synthetic_neg_all_hydrophobic',
      method_name: 'synthetic_negative:all_hydrophobic_core',
      sequence: 'LVLFLILVWMFLVLF',
      length: 15,
      compatibility_score: 3.1204,
      hard_rejection_pass: false,
      provenance: 'Synthetic Negative (all_hydrophobic_core)',
      rejection_reasons: ['hydrophobic_collapse_risk(hydro_frac=1.00, max_run=15)'],
      contour_length_A: 57.0,
    },
    {
      rank: 12,
      candidate_id: 'baseline_direct_0aa',
      method_name: 'direct_0aa_fusion',
      sequence: '',
      length: 0,
      compatibility_score: -2.1164,
      hard_rejection_pass: false,
      provenance: 'Baseline 1: Direct 0-aa fusion',
      rejection_reasons: [
        'physical_span_closure_failed(contour=0.0A < required=28.4A, L=0)',
        'domain_clash_or_entropic_excess(L=0, target=15, detour=1.15)',
      ],
      contour_length_A: 0.0,
    },
  ],
  rigid_non_cleavable: [
    {
      rank: 1,
      candidate_id: 'lf_template_2',
      method_name: 'linkerforge_ranker',
      sequence: 'AEAAAKEAAAKEAAAKA',
      length: 17,
      compatibility_score: 10.4812,
      hard_rejection_pass: true,
      provenance: 'Part 2 Template (A(EAAAK)3A Helical Rod, Chen et al. 2013)',
      rejection_reasons: [],
      contour_length_A: 64.6,
    },
    {
      rank: 2,
      candidate_id: 'lf_template_3',
      method_name: 'linkerforge_ranker',
      sequence: 'PAPAPPAPAPPAPAPPAP',
      length: 18,
      compatibility_score: 10.3190,
      hard_rejection_pass: true,
      provenance: 'Part 2 Template (Chen2013 Proline-Rich PAPAP Rod)',
      rejection_reasons: [],
      contour_length_A: 68.4,
    },
    {
      rank: 3,
      candidate_id: 'lf_authority_1',
      method_name: 'nearest_neighbour_retrieval',
      sequence: 'PEPEPEPEPEPEPEPEPE',
      length: 18,
      compatibility_score: 10.2441,
      hard_rejection_pass: true,
      provenance: 'X2 Authority (mined_uniprot_explicit:P14921:188-205, mined_high_confidence_dfl_call)',
      rejection_reasons: [],
      contour_length_A: 68.4,
    },
    {
      rank: 4,
      candidate_id: 'baseline_fixed_ggggs3',
      method_name: 'fixed_ggggs3',
      sequence: 'GGGGSGGGGSGGGGS',
      length: 15,
      compatibility_score: 8.6120,
      hard_rejection_pass: true,
      provenance: 'Baseline 2: Fixed (GGGGS)3 (mismatched rigidity axis)',
      rejection_reasons: [],
      contour_length_A: 57.0,
    },
    {
      rank: 5,
      candidate_id: 'baseline_direct_0aa',
      method_name: 'direct_0aa_fusion',
      sequence: '',
      length: 0,
      compatibility_score: -2.1164,
      hard_rejection_pass: false,
      provenance: 'Baseline 1: Direct 0-aa fusion',
      rejection_reasons: ['physical_span_closure_failed(contour=0.0A < required=28.4A, L=0)'],
      contour_length_A: 0.0,
    },
  ],
  flexible_cleavable: [
    {
      rank: 1,
      candidate_id: 'lf_template_2',
      method_name: 'linkerforge_ranker',
      sequence: 'GSAGSAENLYFQSGSAGSA',
      length: 19,
      compatibility_score: 10.5240,
      hard_rejection_pass: true,
      provenance: 'Part 2 Cleavable GSA-TEV-GSA Template (Hubbard 1994; Chen 2013)',
      rejection_reasons: [],
      contour_length_A: 72.2,
    },
    {
      rank: 2,
      candidate_id: 'lf_template_1',
      method_name: 'linkerforge_ranker',
      sequence: 'GGGGSENLYFQSGGGGS',
      length: 17,
      compatibility_score: 10.3915,
      hard_rejection_pass: true,
      provenance: 'Part 2 Cleavable GS-TEV-GS Template',
      rejection_reasons: [],
      contour_length_A: 64.6,
    },
    {
      rank: 3,
      candidate_id: 'synthetic_neg_buried_protease',
      method_name: 'synthetic_negative:protease_site_buried_in_rigid_helix',
      sequence: 'EAAAKEAAAKENLYFQSEAAAKEAAAK',
      length: 27,
      compatibility_score: 4.8920,
      hard_rejection_pass: false,
      provenance: 'Synthetic Negative (protease_site_buried_in_rigid_helix)',
      rejection_reasons: ['protease_site_buried_inside_rigid_alpha_helix_without_gs_spacers'],
      contour_length_A: 102.6,
    },
    {
      rank: 4,
      candidate_id: 'baseline_fixed_ggggs3',
      method_name: 'fixed_ggggs3',
      sequence: 'GGGGSGGGGSGGGGS',
      length: 15,
      compatibility_score: 5.1531,
      hard_rejection_pass: false,
      provenance: 'Baseline 2: Fixed (GGGGS)3 (missing TEV cleavage motif)',
      rejection_reasons: ['missing_requested_protease_cleavage_motif'],
      contour_length_A: 57.0,
    },
  ],
}

function downloadText(filename: string, content: string) {
  const url = URL.createObjectURL(new Blob([content], { type: 'text/plain;charset=utf-8' }))
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(url)
}

function downloadJson(filename: string, value: unknown) {
  downloadText(filename, `${JSON.stringify(value, null, 2)}\n`)
}

function downloadBytes(filename: string, bytes: Uint8Array, type: string) {
  const url = URL.createObjectURL(new Blob([bytes as BlobPart], { type }))
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(url)
}

function displayConfidence(value: number | null | undefined) {
  return value == null ? 'Not available' : value.toFixed(1)
}

function ProteinViewer({
  candidate,
  domainALength,
  domainALabel,
  domainBLabel,
  references,
  onReferencesChange,
  msaConsent,
  onResultChange,
}: {
  candidate: FusionCandidate | null
  domainALength: number
  domainALabel: string
  domainBLabel: string
  references: DomainReferences
  onReferencesChange: (references: DomainReferences) => void
  msaConsent: boolean
  onResultChange: (result: ViewerResult) => void
}) {
  const [mode, setMode] = useState<'fusion' | 'A' | 'B'>('fusion')
  const [fullscreen, setFullscreen] = useState(false)
  const [structures, setStructures] = useState<Record<string, { pdb: string; source: 'colabfold' | 'uploaded' }>>({})
  const [uploadError, setUploadError] = useState<string | null>(null)
  const [compute, setCompute] = useState<ComputeStatus | null>(null)
  const [fold, setFold] = useState<FoldStatus | null>(null)
  const [foldSequence, setFoldSequence] = useState<string | null>(null)
  const [foldStarting, setFoldStarting] = useState(false)
  const [foldStartedAt, setFoldStartedAt] = useState<number | null>(null)
  const [clock, setClock] = useState(0)
  const foldLaunchInFlight = useRef(false)
  const [pymolMessage, setPymolMessage] = useState<string | null>(null)
  const [manualMetrics, setManualMetrics] = useState<ColabFoldMetrics | null>(null)
  const selectedSequence = candidate?.fullSequence ?? ''
  const structure = structures[selectedSequence] ?? null
  const pdb = structure?.pdb ?? null
  const linkerLength = candidate?.linker.length ?? 0
  const domainASequence = selectedSequence.slice(0, domainALength)
  const domainBSequence = selectedSequence.slice(domainALength + linkerLength)
  const selectedReference = mode === 'B' ? 'B' : 'A'
  const reference = references[selectedReference]
  const referenceSequence = selectedReference === 'A' ? domainASequence : domainBSequence
  const referenceLabel = selectedReference === 'A' ? domainALabel : domainBLabel

  useEffect(() => {
    if (!inDesktop()) return
    void getComputeStatus().then(setCompute).catch(() => setCompute(null))
  }, [])

  useEffect(() => {
    if (!inDesktop() || !selectedSequence) return
    let cancelled = false
    void getFoldStatus().then(async (status) => {
      const selectedDigest = await sha256Hex(selectedSequence)
      if (cancelled || status.sequenceSha256 !== selectedDigest) return
      setFold(status)
      setFoldStartedAt(foldStartTimeFromRunDir(status.runDir))
      if (status.phase === 'running') {
        setFoldSequence(selectedSequence)
      } else if (status.phase === 'complete' && status.pdb) {
        validateFusionPdb(status.pdb, selectedSequence)
        setStructures((previous) => ({ ...previous, [selectedSequence]: { pdb: status.pdb!, source: 'colabfold' } }))
        setMode('fusion')
        setUploadError(null)
      }
    }).catch((cause) => {
      if (!cancelled) setUploadError(errorMessage(cause, 'Could not restore the local fold status.'))
    })
    return () => { cancelled = true }
  }, [selectedSequence])

  useEffect(() => {
    if (!foldSequence) return
    const timer = window.setInterval(() => {
      void getFoldStatus().then((status) => {
        setFold(status)
        if (status.phase === 'complete' && status.pdb) {
          try {
            validateFusionPdb(status.pdb, foldSequence)
            setStructures((previous) => ({ ...previous, [foldSequence]: { pdb: status.pdb!, source: 'colabfold' } }))
            setMode('fusion')
            setUploadError(null)
          } catch (cause) {
            setUploadError(cause instanceof Error ? cause.message : 'Predicted PDB did not match the fusion.')
          }
          setFoldSequence(null)
        } else if (status.phase === 'failed') {
          setFoldSequence(null)
        }
      }).catch((cause) => {
        setUploadError(cause instanceof Error ? cause.message : 'Could not check the local fold.')
        setFoldSequence(null)
      })
    }, 2000)
    return () => window.clearInterval(timer)
  }, [foldSequence])

  useEffect(() => {
    if (!foldStarting && fold?.phase !== 'running') return
    const timer = window.setInterval(() => setClock(Date.now()), 1000)
    return () => window.clearInterval(timer)
  }, [foldStarting, fold?.phase])

  const effectiveFoldStart = foldStartedAt ?? foldStartTimeFromRunDir(fold?.runDir ?? null)
  const elapsed = effectiveFoldStart === null ? null : formatElapsedTime(clock - effectiveFoldStart)

  useEffect(() => {
    onResultChange({
      fold,
      pdb,
      source: structure?.source ?? null,
      metrics: fold?.metrics ?? manualMetrics,
      elapsed,
      running: foldStarting || fold?.phase === 'running',
    })
  }, [elapsed, fold, foldStarting, manualMetrics, onResultChange, pdb, structure?.source])

  const handleLocalFold = async () => {
    if (!candidate || foldLaunchInFlight.current) return
    foldLaunchInFlight.current = true
    setFoldStarting(true)
    setFoldStartedAt(Date.now())
    setClock(Date.now())
    setUploadError(null)
    try {
      const status = await startLocalFold(candidate.fullSequence, domainALength, candidate.linker.length)
      setFold(status)
      setFoldSequence(candidate.fullSequence)
      setUploadError(null)
    } catch (cause) {
      // A rapid second IPC request can be rejected after the first one has
      // already launched. Reconcile with the backend before reporting failure.
      try {
        const current = await getFoldStatus()
        if (current.phase === 'running') {
          setFold(current)
          setFoldSequence(candidate.fullSequence)
          setFoldStartedAt((previous) => previous ?? foldStartTimeFromRunDir(current.runDir))
          setUploadError(null)
          return
        }
      } catch {
        // Preserve the original launch error below.
      }
      setUploadError(errorMessage(cause, 'Could not start ColabFold.'))
    } finally {
      foldLaunchInFlight.current = false
      setFoldStarting(false)
    }
  }

  useEffect(() => {
    if (!fullscreen) return
    const onKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') setFullscreen(false)
    }
    document.addEventListener('keydown', onKeyDown)
    return () => document.removeEventListener('keydown', onKeyDown)
  }, [fullscreen])

  const handlePdbUpload = async (file: File | undefined) => {
    if (!file || !candidate) return
    try {
      if (file.size > 15_000_000) throw new Error('The local PDB limit is 15 MB.')
      const content = await file.text()
      validateFusionPdb(content, candidate.fullSequence)
      setStructures((previous) => ({ ...previous, [candidate.fullSequence]: { pdb: content, source: 'uploaded' } }))
      setManualMetrics(null)
      setUploadError(null)
      setPymolMessage(null)
      setMode('fusion')
    } catch (cause) {
      setUploadError(cause instanceof Error ? cause.message : 'Could not read the PDB file.')
    }
  }

  const handleScoreUpload = async (file: File | undefined) => {
    if (!file || !candidate) return
    try {
      if (file.size > 20_000_000) throw new Error('The local score JSON limit is 20 MB.')
      const metrics = parseColabFoldScoreJson(await file.text(), candidate.fullSequence.length, domainALength, candidate.linker.length)
      setManualMetrics(metrics)
      setUploadError(null)
    } catch (cause) {
      setManualMetrics(null)
      setUploadError(cause instanceof Error ? cause.message : 'Could not read the ColabFold score JSON.')
    }
  }

  const handleReferenceUpload = async (file: File | undefined) => {
    if (!file || !candidate) return
    try {
      if (file.size > 15_000_000) throw new Error('The local PDB limit is 15 MB.')
      const content = await file.text()
      validateDomainReferencePdb(content, referenceSequence, referenceLabel || `Domain ${selectedReference}`)
      onReferencesChange({
        ...references,
        [selectedReference]: { ...reference, filename: file.name, pdb: content },
      })
      setUploadError(null)
    } catch (cause) {
      setUploadError(cause instanceof Error ? cause.message : 'Could not read the reference PDB.')
    }
  }

  const handleOpenInPymol = async () => {
    if (!candidate || !pdb) return
    try {
      const message = await openPdbInPymol({
        pdb,
        sequence: candidate.fullSequence,
        domainALength,
        linkerLength: candidate.linker.length,
      })
      setPymolMessage(message)
      setUploadError(null)
    } catch (cause) {
      setPymolMessage(null)
      setUploadError(cause instanceof Error ? cause.message : 'Could not open the fusion in PyMOL.')
    }
  }

  const pml = candidate ? [
    '# Run after loading the matching full-fusion PDB in PyMOL.',
    'hide everything, all',
    'show cartoon, all',
    `select domain_a, resi 1-${domainALength}`,
    ...(candidate.linker.length > 0 ? [`select linker, resi ${domainALength + 1}-${domainALength + candidate.linker.length}`] : []),
    `select domain_b, resi ${domainALength + candidate.linker.length + 1}-${candidate.fullSequence.length}`,
    'color green, domain_a',
    ...(candidate.linker.length > 0 ? ['color magenta, linker'] : []),
    'color red, domain_b',
    'zoom all',
  ].join('\n') + '\n' : ''

  return (
    <section className={`section viewer-section ${fullscreen ? 'viewer-fullscreen' : ''}`} id="viewer" aria-label="Structure viewer">
      <div className="viewer-header">
        <div className="section-heading">
          <span>Prediction workspace / full fusion</span>
          <h2>{candidate ? 'Inspect both domains and the connecting linker' : 'Inspect reference structures'}</h2>
        </div>
        <button type="button" className="viewer-expand" onClick={() => setFullscreen((value) => !value)}>
          {fullscreen ? 'Close full screen' : 'Open full screen'}
        </button>
      </div>
      <div className="viewer-mode-row model-access-tabs" role="tablist" aria-label="Prediction models">
        <button type="button" className={mode === 'fusion' ? 'active' : ''} onClick={() => setMode('fusion')}>
          Full fusion model <small>{pdb ? 'PDB loaded' : 'awaits prediction'}</small>
        </button>
        <button type="button" className={mode === 'A' ? 'active' : ''} onClick={() => setMode('A')}>
          {domainALabel} <small>{references.A.pdb ? 'individual model loaded' : 'load individual model'}</small>
        </button>
        <button type="button" className={mode === 'B' ? 'active' : ''} onClick={() => setMode('B')}>
          {domainBLabel} <small>{references.B.pdb ? 'individual model loaded' : 'load individual model'}</small>
        </button>
        {inDesktop() && candidate && <button type="button" className="gpu-fold-button"
          disabled={!compute?.predictorAvailable || !msaConsent || foldStarting || fold?.phase === 'running'}
          onClick={() => { setMode('fusion'); void handleLocalFold() }}>
          <Play size={14} /> {foldStarting ? 'Starting ColabFold…' : fold?.phase === 'running' ? 'Folding locally…' : 'Fold full fusion on this GPU'}
        </button>}
      </div>
      {mode === 'fusion' ? (
        <>
          <div className="fusion-legend">
            <span><i style={{ background: '#3acf78' }} /> Domain A ({domainALength} aa)</span>
            <span><i style={{ background: '#ff4dba' }} /> Linker ({candidate?.linker.length ?? 0} aa)</span>
            <span><i style={{ background: '#f66b62' }} /> Domain B</span>
          </div>
          <div className="viewer-stage">
            <FusionStructureViewer
              pdb={pdb}
              fullSequence={selectedSequence}
              domainALength={domainALength}
              linkerLength={candidate?.linker.length ?? 0}
            />
          </div>
          <div className="viewer-actions">
            {candidate ? (
              <>
                <button type="button" onClick={() => downloadText('linkerforge-fusion.fasta', fusionFasta(candidate.fullSequence))}>
                  <Download size={14} /> Export full-fusion FASTA
                </button>
                <label className="upload-pdb">
                  Load predicted full-fusion PDB
                  <input type="file" accept=".pdb,.ent,text/plain" onChange={(event) => void handlePdbUpload(event.target.files?.[0])} />
                </label>
                {pdb && <label className="upload-pdb">
                  Load matching ColabFold score JSON
                  <input type="file" accept=".json,application/json" onChange={(event) => void handleScoreUpload(event.target.files?.[0])} />
                </label>}
                {pdb && <button type="button" onClick={() => downloadText('linkerforge-fusion.pml', pml)}>
                  <Download size={14} /> PyMOL colouring script
                </button>}
                {inDesktop() && pdb && <button type="button" disabled={!compute?.pymolAvailable}
                  onClick={() => void handleOpenInPymol()}>
                  <ExternalLink size={14} /> Open coloured fusion in PyMOL
                </button>}
              </>
            ) : (
              <p>Enter both protein sequences in the design studio to construct a full fusion first.</p>
            )}
          </div>
          {uploadError && <p className="structure-error" role="alert">{uploadError}</p>}
          {pymolMessage && <p className="structure-success" role="status">{pymolMessage}</p>}
          {inDesktop() && !msaConsent && <p className="structure-warning" role="status">
            Local prediction is locked until the public-MSA data disclosure is acknowledged on the Design input page.
          </p>}
          {inDesktop() && <div className="desktop-compute" aria-live="polite">
            <strong>Local compute:</strong> {compute?.gpu ?? 'Checking GPU…'} · ColabFold {compute?.predictorAvailable ? `available via ${compute.predictorBackend}` : 'not installed'}.
            <p>PyMOL {compute?.pymolAvailable ? `available${compute.pymolPath ? ` at ${compute.pymolPath}` : ''}` : 'not installed or not detected'}.</p>
            {fold && <p>{fold.phase}: {fold.message}</p>}
            {fold?.runDir && <small>Run files: {fold.runDir}</small>}
          </div>}
          <p className="viewer-disclaimer">A manually uploaded PDB is sequence-checked but is not assumed to be an AlphaFold result. Load its matching ColabFold rank-1 score JSON before interpreting PAE or pTM. Linker coordinates appear only when the full fusion has actually been predicted.</p>
        </>
      ) : (
        <>
          <div className="individual-model-heading">
            <div><span>Individual domain model</span><strong>{referenceLabel}</strong></div>
            <small>{reference.pdb ? `${reference.sourceKind} · ${reference.sourceId || reference.filename}` : 'No model loaded for this domain'}</small>
          </div>
          <div className="viewer-stage">
            <FusionStructureViewer
              pdb={reference.pdb}
              fullSequence={referenceSequence}
              domainALength={referenceSequence.length}
              linkerLength={0}
              mode="reference"
              referenceColor={selectedReference === 'A' ? '#3acf78' : '#f66b62'}
              ariaLabel={`${referenceLabel} individual reference structure`}
            />
          </div>
          <div className="reference-controls">
            <label>Reference source identifier
              <input value={reference.sourceId} placeholder="Example: PDB 1UBQ, chain A"
                onChange={(event) => onReferencesChange({ ...references, [selectedReference]: { ...reference, sourceId: event.target.value } })} />
            </label>
            <label>Evidence type
              <select value={reference.sourceKind}
                onChange={(event) => onReferencesChange({ ...references, [selectedReference]: { ...reference, sourceKind: event.target.value as DomainReference['sourceKind'] } })}>
                <option value="experimental">Experimental structure</option>
                <option value="predicted">Predicted structure</option>
              </select>
            </label>
            <label className="upload-pdb">Load {referenceLabel} reference PDB
              <input type="file" accept=".pdb,.ent,text/plain" onChange={(event) => void handleReferenceUpload(event.target.files?.[0])} />
            </label>
          </div>
          <p className="viewer-disclaimer">This reference belongs to the current entered domain only. Its single-chain C-alpha sequence must exactly match the domain sequence. The source identifier and experimental/predicted classification are researcher-supplied provenance; LinkerForge does not infer them.</p>
        </>
      )}
    </section>
  )
}
function newDesignSessionId() {
  return `design-${Date.now()}-${Math.random().toString(16).slice(2, 8)}`
}

function App() {
  const historyStorageKey = 'linkerforge.research-history.v1'
  const sessionStorageKey = 'linkerforge.design-sessions.v1'
  const themeStorageKey = 'linkerforge.color-theme.v1'
  const [colorTheme, setColorTheme] = useState<ColorTheme>(() => {
    const savedTheme = localStorage.getItem(themeStorageKey)
    if (savedTheme === 'light' || savedTheme === 'dark') return savedTheme
    return typeof window.matchMedia === 'function' && window.matchMedia('(prefers-color-scheme: dark)').matches
      ? 'dark'
      : 'light'
  })
  const [page, setPage] = useState<AppPage>(() => pageFromHash(window.location.hash))
  const [composerOpen, setComposerOpen] = useState(true)
  const [activeSessionId, setActiveSessionId] = useState<string | null>(null)
  const [designSessions, setDesignSessions] = useState<DesignSession[]>(() => {
    try {
      const saved = localStorage.getItem(sessionStorageKey)
      return saved ? JSON.parse(saved) as DesignSession[] : []
    } catch {
      return []
    }
  })
  const [projectName, setProjectName] = useState('Untitled linker design')
  const [domainALabel, setDomainALabel] = useState('Domain A')
  const [domainBLabel, setDomainBLabel] = useState('Domain B')
  const [biologicalObjective, setBiologicalObjective] = useState('')
  const [expressionSystem, setExpressionSystem] = useState('Not specified')
  const [protectedSites, setProtectedSites] = useState('')
  const [msaConsent, setMsaConsent] = useState(false)
  const [viewerResult, setViewerResult] = useState<ViewerResult>({ fold: null, pdb: null, source: null, metrics: null, elapsed: null, running: false })
  const [domainReferences, setDomainReferences] = useState<DomainReferences>({
    A: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null },
    B: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null },
  })
  const [domainAInput, setDomainAInput] = useState('')
  const [domainBInput, setDomainBInput] = useState('')
  const [domainALength, setDomainALength] = useState(0)
  const [customCandidates, setCustomCandidates] = useState<FusionCandidate[]>([])
  const [selectedCustomIndex, setSelectedCustomIndex] = useState(0)
  const [customError, setCustomError] = useState<string | null>(null)
  const [rigidity, setRigidity] = useState<RigidityMode>('flexible')
  const [cleavable, setCleavable] = useState<boolean>(false)
  const [targetLength, setTargetLength] = useState<number>(15)
  const [candidates, setCandidates] = useState<RankedCandidateRow[]>([])
  const [selectedCandidateSeq, setSelectedCandidateSeq] = useState<string>('')
  const [apiConnected, setApiConnected] = useState<boolean>(false)
  const [researchHistory, setResearchHistory] = useState<HistoryEntry[]>(() => {
    try {
      const saved = localStorage.getItem(historyStorageKey)
      return saved ? JSON.parse(saved) as HistoryEntry[] : []
    } catch {
      return []
    }
  })
  const [historyVerified, setHistoryVerified] = useState<boolean | null>(null)

  useEffect(() => {
    const syncPage = () => {
      setPage(pageFromHash(window.location.hash))
    }
    window.addEventListener('hashchange', syncPage)
    return () => window.removeEventListener('hashchange', syncPage)
  }, [])

  useEffect(() => {
    if (window.location.hash !== `#${page}`) window.history.replaceState(null, '', `#${page}`)
  }, [page])

  useEffect(() => {
    document.documentElement.dataset.theme = colorTheme
    document.documentElement.style.colorScheme = colorTheme
    localStorage.setItem(themeStorageKey, colorTheme)
  }, [colorTheme])

  useEffect(() => {
    try {
      localStorage.setItem(sessionStorageKey, JSON.stringify(designSessions))
    } catch {
      // Keep the current session usable when browser storage is unavailable.
    }
  }, [designSessions])

  useEffect(() => {
    let cancelled = false
    void verifyHistory(researchHistory).then((verified) => {
      if (!cancelled) setHistoryVerified(verified)
    })
    return () => { cancelled = true }
  }, [researchHistory])

  const goToPage = (nextPage: AppPage) => {
    setPage(nextPage)
  }

  const currentSession = (id: string): DesignSession => ({
    id,
    updatedAt: new Date().toISOString(),
    projectName: projectName.trim() || 'Untitled linker design',
    domainALabel: domainALabel.trim() || 'Domain A',
    domainBLabel: domainBLabel.trim() || 'Domain B',
    domainAInput,
    domainBInput,
    biologicalObjective,
    expressionSystem,
    protectedSites,
    rigidity,
    cleavable,
    targetLength,
  })

  const openSession = (session: DesignSession) => {
    try {
      const domainA = parseProteinSequence(session.domainAInput, 'Domain A')
      const domainB = parseProteinSequence(session.domainBInput, 'Domain B')
      setProjectName(session.projectName)
      setDomainALabel(session.domainALabel)
      setDomainBLabel(session.domainBLabel)
      setDomainAInput(session.domainAInput)
      setDomainBInput(session.domainBInput)
      setBiologicalObjective(session.biologicalObjective)
      setExpressionSystem(session.expressionSystem)
      setProtectedSites(session.protectedSites)
      setMsaConsent(false)
      setRigidity(session.rigidity)
      setCleavable(session.cleavable)
      setTargetLength(session.targetLength)
      setDomainALength(domainA.length)
      setCustomCandidates(buildTemplateCandidates(domainA, domainB, session.targetLength, session.rigidity, session.cleavable))
      setCandidates([])
      setSelectedCustomIndex(0)
      setSelectedCandidateSeq('')
      setViewerResult({ fold: null, pdb: null, source: null, metrics: null, elapsed: null, running: false })
      setDomainReferences({
        A: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null },
        B: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null },
      })
      setActiveSessionId(session.id)
      setComposerOpen(false)
      setCustomError(null)
      goToPage('structure')
    } catch (cause) {
      setComposerOpen(true)
      setCustomError(cause instanceof Error ? cause.message : 'Could not open this design session.')
    }
  }

  const startNewDesign = () => {
    setActiveSessionId(null)
    setComposerOpen(true)
    setProjectName('Untitled linker design')
    setDomainALabel('Domain A')
    setDomainBLabel('Domain B')
    setDomainAInput('')
    setDomainBInput('')
    setBiologicalObjective('')
    setExpressionSystem('Not specified')
    setProtectedSites('')
    setRigidity('flexible')
    setCleavable(false)
    setTargetLength(15)
    setMsaConsent(false)
    setCustomCandidates([])
    setCandidates([])
    setSelectedCustomIndex(0)
    setSelectedCandidateSeq('')
    setViewerResult({ fold: null, pdb: null, source: null, metrics: null, elapsed: null, running: false })
    setDomainReferences({
      A: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null },
      B: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null },
    })
    setCustomError(null)
    goToPage('structure')
  }

  const editSession = (session: DesignSession) => {
    openSession(session)
    setComposerOpen(true)
  }

  const deleteSession = (session: DesignSession) => {
    const confirmed = window.confirm(
      `Delete "${session.projectName}"? This removes the saved design from this device but does not delete any prediction run files.`
    )
    if (!confirmed) return
    setDesignSessions((previous) => previous.filter((item) => item.id !== session.id))
    if (activeSessionId === session.id) startNewDesign()
  }

  const handleRunDesign = async (
    nextRigidity: RigidityMode = rigidity,
    nextCleavable: boolean = cleavable,
    nextLength: number = targetLength
  ) => {
    try {
      const domainA = parseProteinSequence(domainAInput, 'Domain A')
      const domainB = parseProteinSequence(domainBInput, 'Domain B')
      const generated = buildTemplateCandidates(domainA, domainB, nextLength, nextRigidity, nextCleavable)
      const protectedSiteParse = parseProtectedSites(
        protectedSites,
        domainA.length,
        generated[0]?.linker.length ?? nextLength,
        domainB.length
      )
      if (protectedSiteParse.errors.length > 0) {
        throw new Error(`Protected residues must use valid domain-relative positions (for example A:48;B:12-18). Check: ${protectedSiteParse.errors.join(', ')}`)
      }
      setDomainALength(domainA.length)
      setCustomCandidates(generated)
      setCandidates([])
      setSelectedCustomIndex(0)
      setSelectedCandidateSeq('')
      setCustomError(null)
      setApiConnected(false)
      const id = activeSessionId ?? newDesignSessionId()
      const session = currentSession(id)
      setDesignSessions((previous) => [session, ...previous.filter((item) => item.id !== id)])
      setActiveSessionId(id)
      setComposerOpen(false)
      setViewerResult({ fold: null, pdb: null, source: null, metrics: null, elapsed: null, running: false })
      goToPage('structure')
    } catch (cause) {
      setCustomCandidates([])
      setCustomError(cause instanceof Error ? cause.message : 'Invalid sequence input.')
    }
  }

  const handleLoadPresetCandidates = async () => {
    setCustomCandidates([])
    setCustomError(null)
    try {
      const resp = await fetch('http://127.0.0.1:8765/api/design', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          rigidity,
          cleavable,
          protease_family: 'TEV',
          target_length: targetLength,
          required_path_length_A: 28.4,
        }),
      })
      if (resp.ok) {
        const data: { ranked_candidates: ApiCandidateRow[] } = await resp.json()
        setApiConnected(true)
        const mapped: RankedCandidateRow[] = data.ranked_candidates.map((c) => ({
          rank: c.rank,
          candidate_id: c.candidate_id,
          method_name: c.method_name,
          sequence: c.sequence,
          length: c.length,
          compatibility_score: c.compatibility_score,
          hard_rejection_pass: c.hard_rejection_pass,
          provenance: c.provenance,
          rejection_reasons: c.scorecard?.rejection_reasons ?? [],
          contour_length_A: c.scorecard?.contour_length_A ?? c.length * 3.8,
        }))
        setCandidates(mapped)
        if (mapped.length > 0) {
          setSelectedCandidateSeq(mapped[0].sequence)
        }
        goToPage('candidates')
        return
      }
    } catch {
      // The local scientific dataset/API is optional for the labelled UI demonstration.
    }
    setApiConnected(false)
    const key = cleavable
      ? 'flexible_cleavable'
      : rigidity === 'rigid'
        ? 'rigid_non_cleavable'
        : 'flexible_non_cleavable'
    const rows = PRESET_CANDIDATES[key]
    setCandidates(rows)
    if (rows.length > 0) {
      setSelectedCandidateSeq(rows[0].sequence)
    }
  }

  const selectedCustomCandidate = customCandidates[selectedCustomIndex] ?? null
  const selectedDomainBLength = selectedCustomCandidate
    ? selectedCustomCandidate.fullSequence.length - domainALength - selectedCustomCandidate.linker.length
    : 0
  const localComputeRisk = !selectedCustomCandidate
    ? null
    : selectedCustomCandidate.fullSequence.length > 400
      ? 'High risk on this 4 GB GPU; use the college cloud if local execution fails.'
      : selectedCustomCandidate.fullSequence.length > 150
        ? 'Moderate local GPU risk; retain a cloud fallback.'
        : 'Short-sequence smoke-test range; successful execution is still not guaranteed.'
  const confidence = viewerResult.pdb && selectedCustomCandidate
    ? summarizePdbConfidence(viewerResult.pdb, domainALength, selectedCustomCandidate.linker.length)
    : null
  const predictionPhase = viewerResult.pdb ? 'complete' : viewerResult.fold?.phase ?? 'not started'
  const coordinateAudit = viewerResult.pdb ? auditFusionCoordinates(viewerResult.pdb) : null
  const domainPreservation = viewerResult.pdb && selectedCustomCandidate
    ? calculateDomainPreservation(
      viewerResult.pdb,
      domainALength,
      selectedCustomCandidate.linker.length,
      domainReferences.A.pdb,
      domainReferences.B.pdb
    )
    : null
  const protectedSiteAudit = selectedCustomCandidate
    ? parseProtectedSites(protectedSites, domainALength, selectedCustomCandidate.linker.length, selectedDomainBLength)
    : { residueIndices: [], errors: [] }
  const protectedSiteClearance = viewerResult.pdb && selectedCustomCandidate && protectedSiteAudit.errors.length === 0 && protectedSiteAudit.residueIndices.length > 0
    ? assessProtectedSiteClearance(viewerResult.pdb, protectedSiteAudit.residueIndices, domainALength, selectedCustomCandidate.linker.length)
    : null

  const projectRecord = (): SavedProjectV1 => ({
      schema_version: 'linkerforge.project.v1',
      project: {
        name: projectName,
        biological_objective: biologicalObjective,
        expression_system: expressionSystem,
        protected_sites: protectedSites,
      },
      domains: {
        A: { label: domainALabel, sequence: domainAInput, reference_source_id: domainReferences.A.sourceId, reference_source_kind: domainReferences.A.sourceKind },
        B: { label: domainBLabel, sequence: domainBInput, reference_source_id: domainReferences.B.sourceId, reference_source_kind: domainReferences.B.sourceKind },
      },
      design: { rigidity, cleavable, target_length: targetLength },
      privacy: { public_msa_consent_saved: false },
  })

  const exportProject = () => {
    downloadJson('linkerforge-project.v1.json', projectRecord())
  }

  const importProject = async (file: File | undefined) => {
    if (!file) return
    try {
      const saved = JSON.parse(await file.text()) as SavedProjectV1
      if (saved.schema_version !== 'linkerforge.project.v1' || !saved.domains?.A || !saved.domains?.B) {
        throw new Error('This is not a LinkerForge project.v1 file.')
      }
      if (saved.domains.A.sequence.trim() || saved.domains.B.sequence.trim()) {
        parseProteinSequence(saved.domains.A.sequence, 'Domain A')
        parseProteinSequence(saved.domains.B.sequence, 'Domain B')
      }
      if (!['flexible', 'rigid'].includes(saved.design.rigidity)) throw new Error('Saved rigidity setting is invalid.')
      if (!['experimental', 'predicted'].includes(saved.domains.A.reference_source_kind)
        || !['experimental', 'predicted'].includes(saved.domains.B.reference_source_kind)) {
        throw new Error('Saved reference evidence type is invalid.')
      }
      if (!Number.isInteger(saved.design.target_length) || saved.design.target_length < 5 || saved.design.target_length > 30) {
        throw new Error('Saved target linker length must be between 5 and 30 residues.')
      }
      setProjectName(saved.project.name)
      setBiologicalObjective(saved.project.biological_objective)
      setExpressionSystem(saved.project.expression_system)
      setProtectedSites(saved.project.protected_sites)
      setDomainALabel(saved.domains.A.label)
      setDomainBLabel(saved.domains.B.label)
      setDomainAInput(saved.domains.A.sequence)
      setDomainBInput(saved.domains.B.sequence)
      setRigidity(saved.design.rigidity)
      setCleavable(saved.design.cleavable)
      setTargetLength(saved.design.target_length)
      setMsaConsent(false)
      setCustomCandidates([])
      setDomainReferences({
        A: { sourceId: saved.domains.A.reference_source_id, sourceKind: saved.domains.A.reference_source_kind, filename: null, pdb: null },
        B: { sourceId: saved.domains.B.reference_source_id, sourceKind: saved.domains.B.reference_source_kind, filename: null, pdb: null },
      })
      setActiveSessionId(null)
      setComposerOpen(true)
      setCustomError(null)
    } catch (cause) {
      setCustomError(cause instanceof Error ? cause.message : 'Could not open the project file.')
    }
  }

  const resultRecord = async () => {
    if (!selectedCustomCandidate) throw new Error('Select a candidate before creating a result record.')
    const configuration = JSON.stringify({ rigidity, cleavable, targetLength, protectedSites })
    return {
      schema_version: 'linkerforge.result.v1',
      generated_at: new Date().toISOString(),
      status: viewerResult.pdb ? 'unassessed' : viewerResult.fold?.phase === 'failed' ? 'failed' : 'unassessed',
      project: { name: projectName, biological_objective: biologicalObjective, expression_system: expressionSystem },
      candidate: {
        method: selectedCustomCandidate.method,
        evidence_quality: selectedCustomCandidate.evidenceQuality,
        domain_a_label: domainALabel,
        domain_a_length: domainALength,
        linker_sequence: selectedCustomCandidate.linker,
        linker_length: selectedCustomCandidate.linker.length,
        domain_b_label: domainBLabel,
        domain_b_length: selectedDomainBLength,
        full_sequence_length: selectedCustomCandidate.fullSequence.length,
        full_sequence_sha256: await sha256Hex(selectedCustomCandidate.fullSequence),
      },
      structure_prediction: {
        source: viewerResult.source,
        run_directory: viewerResult.fold?.runDir ?? null,
        sequence_sha256: viewerResult.fold?.sequenceSha256 ?? null,
        artifacts: viewerResult.fold?.artifacts ?? null,
      },
      structural_metrics: {
        confidence_source: viewerResult.metrics?.predictor ?? (viewerResult.pdb ? 'PDB B-factor field; predictor unknown' : null),
        plddt: confidence,
        ptm: viewerResult.metrics?.ptm ?? null,
        max_pae: viewerResult.metrics?.maxPae ?? null,
        mean_cross_domain_pae_angstrom: viewerResult.metrics?.meanInterdomainPae ?? null,
        coordinate_pre_screen: coordinateAudit,
        domain_preservation_rmsd_angstrom: domainPreservation,
        protected_site_mapping: { ...protectedSiteAudit, spatial_ca_pre_screen: protectedSiteClearance },
      },
      references: {
        A: {
          source_id: domainReferences.A.sourceId,
          source_kind: domainReferences.A.sourceKind,
          filename: domainReferences.A.filename,
          pdb_sha256: domainReferences.A.pdb ? await sha256Hex(domainReferences.A.pdb) : null,
        },
        B: {
          source_id: domainReferences.B.sourceId,
          source_kind: domainReferences.B.sourceKind,
          filename: domainReferences.B.filename,
          pdb_sha256: domainReferences.B.pdb ? await sha256Hex(domainReferences.B.pdb) : null,
        },
      },
      completed_comparators: [],
      limitations: [
        'No validated LinkerForge compatibility ranker was used.',
        'No matched structural comparison of all required controls was completed.',
        'No Rosetta closure/all-atom clash assessment was completed.',
        'No OpenMM/GROMACS validation or experimental assay was completed.',
      ],
      provenance: {
        application_version: '0.1.0',
        pipeline_commit: import.meta.env.VITE_LINKERFORGE_COMMIT ?? 'not embedded',
        working_tree_state: import.meta.env.VITE_LINKERFORGE_DIRTY === 'true' ? 'modified' : 'clean',
        config_sha256: await sha256Hex(configuration),
        evidence_quality: 'predicted_unvalidated',
      },
      claim_boundary: 'Computational candidate record; not evidence of expression, function, or linker optimality.',
    }
  }

  const exportResultRecord = async () => {
    downloadJson('linkerforge-result.v1.json', await resultRecord())
  }

  const recordResultInHistory = async () => {
    const updated = await appendHistoryEntry(researchHistory, await resultRecord())
    localStorage.setItem(historyStorageKey, JSON.stringify(updated))
    setResearchHistory(updated)
  }

  const exportResearchBundle = async () => {
    if (!selectedCustomCandidate) return
    const files: Record<string, ResearchBundleFile> = {
      'project.json': `${JSON.stringify(projectRecord(), null, 2)}\n`,
      'result.json': `${JSON.stringify(await resultRecord(), null, 2)}\n`,
      'fusion.fasta': fusionFasta(selectedCustomCandidate.fullSequence),
      'history.json': `${JSON.stringify(researchHistory, null, 2)}\n`,
      'README.txt': [
        'LinkerForge computational research bundle',
        '',
        'Verify every file against checksums.sha256 before interpretation.',
        'This bundle records an unvalidated computational candidate, not experimental evidence.',
        'PyMOL is a visualization tool and no molecular simulation is implied.',
        '',
      ].join('\n'),
    }
    if (viewerResult.pdb) files['structures/fusion.pdb'] = viewerResult.pdb
    if (domainReferences.A.pdb) files['references/domain-A.pdb'] = domainReferences.A.pdb
    if (domainReferences.B.pdb) files['references/domain-B.pdb'] = domainReferences.B.pdb
    files['checksums.sha256'] = await checksumManifest(files)
    downloadBytes('linkerforge-research-bundle.zip', createStoredZip(files), 'application/zip')
  }

  const activePage = appPages.find((item) => item.id === page) ?? appPages[0]
  const isDesktop = inDesktop()
  const workspaceState = viewerResult.pdb
    ? 'Structure loaded'
    : customCandidates.length > 0
      ? `${customCandidates.length} candidates ready`
      : candidates.length > 0
        ? 'Labelled example loaded'
        : 'Awaiting input'

  return (
    <div className={`app-shell theme-${colorTheme} ${isDesktop ? 'native-shell' : 'web-shell'}`}>
      <div className="app-workspace">
        <aside className="app-sidebar">
          <div className="sidebar-brand">
            <span className="brand-mark">LF</span>
            <span><strong>LinkerForge</strong><small>Local research workspace</small></span>
          </div>
          <button type="button" className="new-design-button" onClick={startNewDesign}>
            <Plus size={16} /> New linker design
          </button>

          <section className={`sidebar-composer ${composerOpen ? 'open' : ''}`} aria-label="New linker design input">
            <button type="button" className="composer-toggle" aria-expanded={composerOpen}
              onClick={() => setComposerOpen((value) => !value)}>
              <span><strong>{activeSessionId ? 'Edit design inputs' : 'Describe a new fusion'}</strong><small>Two sequences + design constraints</small></span>
              <span>{composerOpen ? '−' : '+'}</span>
            </button>
            {composerOpen && <form className="sidebar-input-form" onSubmit={(event) => { event.preventDefault(); void handleRunDesign(rigidity, cleavable, targetLength) }}>
              <label htmlFor="project-name">Design name
                <input id="project-name" value={projectName} onChange={(event) => setProjectName(event.target.value)} />
              </label>
              <div className="sidebar-domain-grid">
                <label htmlFor="domain-a-label">Domain A name
                  <input id="domain-a-label" value={domainALabel} onChange={(event) => setDomainALabel(event.target.value)} />
                </label>
                <label htmlFor="domain-b-label">Domain B name
                  <input id="domain-b-label" value={domainBLabel} onChange={(event) => setDomainBLabel(event.target.value)} />
                </label>
              </div>
              <label htmlFor="domain-a">{domainALabel || 'Domain A'} amino-acid sequence
                <textarea id="domain-a" value={domainAInput}
                  onChange={(event) => { setDomainAInput(event.target.value); setCustomCandidates([]); setDomainReferences({ A: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null }, B: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null } }) }}
                  placeholder="Paste sequence or FASTA" rows={4} />
              </label>
              <label htmlFor="domain-b">{domainBLabel || 'Domain B'} amino-acid sequence
                <textarea id="domain-b" value={domainBInput}
                  onChange={(event) => { setDomainBInput(event.target.value); setCustomCandidates([]); setDomainReferences({ A: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null }, B: { sourceId: '', sourceKind: 'experimental', filename: null, pdb: null } }) }}
                  placeholder="Paste sequence or FASTA" rows={4} />
              </label>

              <div className="sidebar-constraint-row">
                <label>Linker type
                  <select value={rigidity} onChange={(event) => { const value = event.target.value as RigidityMode; setRigidity(value); if (value === 'rigid') setCleavable(false); setCustomCandidates([]) }}>
                    <option value="flexible">Flexible</option>
                    <option value="rigid">Rigid</option>
                  </select>
                </label>
                <label>Cleavage
                  <select value={cleavable ? 'tev' : 'intact'} onChange={(event) => { const value = event.target.value === 'tev'; setCleavable(value); if (value) setRigidity('flexible'); setCustomCandidates([]) }}>
                    <option value="intact">Non-cleavable</option>
                    <option value="tev">TEV-cleavable</option>
                  </select>
                </label>
              </div>
              <label>Target linker length <strong>{targetLength} aa</strong>
                <input type="range" min={5} max={30} value={targetLength}
                  onChange={(event) => { setTargetLength(Number(event.target.value)); setCustomCandidates([]) }} />
              </label>

              <section className="sidebar-advanced-inputs" aria-labelledby="biological-context-heading">
                <div id="biological-context-heading" className="sidebar-advanced-heading">Biological context and safeguards</div>
                <label htmlFor="biological-objective">Biological objective
                  <textarea id="biological-objective" value={biologicalObjective} rows={3}
                    onChange={(event) => setBiologicalObjective(event.target.value)} placeholder="What must the fusion preserve or enable?" />
                </label>
                <label htmlFor="expression-system">Expression system
                  <select id="expression-system" value={expressionSystem} onChange={(event) => setExpressionSystem(event.target.value)}>
                    <option>Not specified</option><option>E. coli</option><option>Yeast</option><option>Insect cells</option><option>Mammalian cells</option><option>Cell-free</option><option>Other</option>
                  </select>
                </label>
                <label htmlFor="protected-sites">Protected residues
                  <input id="protected-sites" value={protectedSites} onChange={(event) => setProtectedSites(event.target.value)} placeholder="A:48;B:12-18" />
                </label>
                <label className="sidebar-msa-consent">
                  <input type="checkbox" checked={msaConsent} onChange={(event) => setMsaConsent(event.target.checked)} />
                  <span>Allow this complete sequence to use the public ColabFold/MMseqs2 MSA service.</span>
                </label>
                <p>Use exact expressed sequences. Signal peptides, tags and propeptides must be included or removed intentionally. Only canonical amino acids are accepted.</p>
              </section>

              <button type="submit" className="sidebar-submit" disabled={!domainAInput.trim() || !domainBInput.trim()}>
                <Play size={15} /> Create prediction workspace
              </button>
              <div className="sidebar-project-files">
                <button type="button" onClick={exportProject}><Download size={13} /> Save JSON</button>
                <label>Open JSON<input type="file" accept=".json,application/json" onChange={(event) => void importProject(event.target.files?.[0])} /></label>
              </div>
              {customError && <p className="sidebar-error" role="alert">{customError}</p>}
            </form>}
          </section>

          <div className="session-list" aria-label="Saved design sessions">
            <div className="session-list-heading"><span>Designs</span><small>{designSessions.length}</small></div>
            {designSessions.length === 0 ? <p>No designs yet. Your submitted inputs will appear here.</p> : designSessions.map((session) => (
              <div key={session.id} className={`session-item ${activeSessionId === session.id ? 'active' : ''}`}>
                <button type="button" className="session-open" aria-label={`Open ${session.projectName}`}
                  onClick={() => openSession(session)}>
                  <MessageSquare size={15} />
                  <span><strong>{session.projectName}</strong><small>{session.domainALabel} → {session.domainBLabel} · {session.targetLength} aa</small></span>
                </button>
                <div className="session-actions" aria-label={`${session.projectName} actions`}>
                  <button type="button" aria-label={`Edit ${session.projectName}`} title="Edit design"
                    onClick={() => editSession(session)}><Pencil size={13} /></button>
                  <button type="button" aria-label={`Delete ${session.projectName}`} title="Delete design"
                    onClick={() => deleteSession(session)}><Trash2 size={13} /></button>
                </div>
              </div>
            ))}
          </div>
          <div className="sidebar-status"><strong>{workspaceState}</strong><small>Saved locally on this device · computational evidence only</small></div>
        </aside>

        <div className="workspace-main">
          <header className="workspace-toolbar">
            <div className="workspace-title"><span>Current design</span><strong>{activeSessionId ? projectName : 'New prediction workspace'}</strong></div>
            <nav className="workspace-tabs" aria-label="Workspace views">
              {appPages.map((item) => <button type="button" key={item.id} className={page === item.id ? 'active' : ''}
                onClick={() => goToPage(item.id)} aria-current={page === item.id ? 'page' : undefined}>
                {item.label}
              </button>)}
            </nav>
            <button type="button" className="theme-toggle"
              onClick={() => setColorTheme((theme) => theme === 'dark' ? 'light' : 'dark')}
              aria-label={`Switch to ${colorTheme === 'dark' ? 'light' : 'dark'} mode`}
              title={`Switch to ${colorTheme === 'dark' ? 'light' : 'dark'} mode`}>
              {colorTheme === 'dark' ? <Sun size={15} /> : <Moon size={15} />}
              <span>{colorTheme === 'dark' ? 'Light' : 'Dark'}</span>
            </button>
            <div className={`prediction-clock ${viewerResult.running ? 'running' : viewerResult.fold?.phase ?? 'idle'}`} aria-live="polite">
              <Clock3 size={15} />
              <span><small>{viewerResult.running ? 'Prediction running' : viewerResult.pdb ? 'Prediction loaded' : isDesktop ? 'GPU ready' : 'Review mode'}</small>
                <strong>{viewerResult.running ? viewerResult.elapsed ?? '00:00' : viewerResult.pdb ? 'Complete' : activePage.label}</strong></span>
            </div>
          </header>

          {customCandidates.length > 0 && <div className="candidate-switcher" aria-label="Candidate models">
            <span>Candidate</span>
            {customCandidates.map((candidate, index) => <button type="button" key={candidate.id}
              className={selectedCustomIndex === index ? 'active' : ''} onClick={() => setSelectedCustomIndex(index)}>
              <strong>{candidate.method}</strong><small>{candidate.linker.length} aa · {candidate.linker || 'direct'}</small>
            </button>)}
          </div>}

          <main id="app-page" className="page-shell">

        <div className="app-page" hidden={page !== 'structure'}>
          {selectedCustomCandidate ? <>
            <ProteinViewer candidate={selectedCustomCandidate} domainALength={domainALength}
              domainALabel={domainALabel || 'Domain A'} domainBLabel={domainBLabel || 'Domain B'}
              references={domainReferences} onReferencesChange={setDomainReferences}
              msaConsent={msaConsent} onResultChange={setViewerResult} />
            <div className="page-actions section-page-actions">
              <button type="button" onClick={() => goToPage('candidates')}>Review candidates</button>
              <button type="button" className="primary-btn" onClick={() => goToPage('results')}>Open results <ArrowRight size={16} /></button>
            </div>
          </> : <section className="workspace-empty" aria-label="Empty prediction workspace">
            <span className="brand-mark">LF</span>
            <h2>Start a linker prediction</h2>
            <p>Enter both protein-domain sequences and design constraints in the sidebar. LinkerForge will create an auditable candidate set and open the selected full-fusion workspace here.</p>
            <button type="button" onClick={() => setComposerOpen(true)}>Open design inputs</button>
          </section>}
        </div>

        {page === 'candidates' && <div className="app-page"><section className="section candidate-section" id="candidates">
          <div className="candidate-intro-grid">
            <div>
            <div className="section-heading">
                <span>Candidate decision gate</span>
              <h2>Choose which complete fusion is worth predicting</h2>
            </div>
            <p>
              {customCandidates.length > 0
                ? `Generated ${customCandidates.length} auditable design candidates and required controls. Select one exact sequence to send through the same structural-prediction settings.`
                : candidates.length > 0
                  ? apiConnected
                    ? 'The local API returned these candidate rows. They still require the same provenance, validation, and structural evidence gates.'
                    : 'These are labelled preset examples for the fixed sfGFP–mCherry pair. They are not predictions for user-entered proteins.'
                  : 'No candidate set is loaded. Create one from your sequences or open the labelled example set here for interface training.'}
            </p>
            {customCandidates.length === 0 && candidates.length === 0 && <div className="candidate-example-action">
              <div>
                <strong>Need a safe interface walkthrough?</strong>
                <span>The sfGFP–mCherry example is visibly labelled and never treated as a result for your project.</span>
              </div>
              <button type="button" onClick={() => void handleLoadPresetCandidates()}>
                <Play size={15} /> Load labelled preset example
              </button>
            </div>}
            <div className="why-candidates">
              <strong>Why this page comes before Structure</strong>
              <p>ColabFold predicts one complete Domain A–linker–Domain B sequence and consumes GPU time. This gate records the exact sequence, design method, control role, and biological context before compute begins. A professional comparison should later fold the shortlisted design and required controls under identical settings.</p>
            </div>
            </div>
            <aside className="candidate-context" aria-label="Selected construct context">
            <span>Current construct</span>
            <strong>{projectName || 'Untitled linker design'}</strong>
            <dl>
              <div><dt>Domains</dt><dd>{domainALabel || 'Domain A'} → {domainBLabel || 'Domain B'}</dd></div>
              <div><dt>Expression</dt><dd>{expressionSystem}</dd></div>
              <div><dt>Objective</dt><dd>{biologicalObjective || 'Not documented'}</dd></div>
              <div><dt>Protected sites</dt><dd>{protectedSites || 'Not documented'}</dd></div>
            </dl>
            </aside>
          </div>
          {customCandidates.length === 0 && candidates.length > 0 && <p className="notice"><strong>{apiConnected ? 'Unvalidated local API output:' : 'Labelled demonstration data:'}</strong> the ranks and compatibility scores below are not produced by the validated frozen fusion ranker, are not calibrated probabilities, and must not be used to choose a research construct.</p>}
          {customCandidates.length > 0 ? <div className="custom-candidate-list" aria-label="Custom fusion candidates">
            {customCandidates.map((candidate, index) => <button type="button" key={candidate.id}
              className={selectedCustomIndex === index ? 'active' : ''} aria-pressed={selectedCustomIndex === index}
              onClick={() => setSelectedCustomIndex(index)}>
              <span className={`candidate-role ${candidate.role}`}>{candidate.role.replace('_', ' ')}</span>
              <strong>{candidate.method} · {candidate.linker.length} aa</strong>
              <code>{candidate.linker || '(no linker)'}</code>
              <small>Full fusion: {candidate.fullSequence.length} residues · {candidate.role === 'design_template' ? 'unranked design candidate' : 'comparison control'}</small>
              <small>{candidate.rationale}</small>
            </button>)}
          </div> : candidates.length > 0 ? <div className="run-table" role="table" aria-label="Example runs">
            <div className="run-row header" role="row">
              <span>Rank / ID</span>
              <span>Sequence &amp; Provenance</span>
              <span>Method / Gate</span>
              <span>Score</span>
            </div>
            {candidates.map((row) => (
              <div
                className={`run-row clickable-row ${selectedCandidateSeq === row.sequence ? 'selected-row' : ''}`}
                role="row"
                key={row.candidate_id}
                onClick={() => setSelectedCandidateSeq(row.sequence)}
              >
                <strong>
                  #{row.rank} · {row.candidate_id}
                </strong>
                <span>
                  <code>{row.sequence || '(0-aa direct fusion)'}</code> ({row.length} aa,{' '}
                  {row.contour_length_A.toFixed(1)} Å max)
                  <small className="prov-line">{row.provenance}</small>
                  {row.rejection_reasons.length > 0 && (
                    <small className="reject-line">
                      Rejected: {row.rejection_reasons.join('; ')}
                    </small>
                  )}
                </span>
                <span>
                  {row.hard_rejection_pass ? (
                    <span className="badge-pass">
                      <CheckCircle2 size={13} /> PASS ({row.method_name})
                    </span>
                  ) : (
                    <span className="badge-fail">
                      <XCircle size={13} /> REJECT ({row.method_name})
                    </span>
                  )}
                </span>
                <span className="score">{row.compatibility_score.toFixed(2)}</span>
              </div>
            ))}
          </div> : <div className="candidate-empty-state">
            <strong>No candidates yet</strong>
            <p>Use the input composer in the sidebar for a new construct, or load the labelled preset above to inspect the candidate-review interface.</p>
          </div>}
          {selectedCustomCandidate && <div className="selected-candidate-review" aria-label="Selected candidate review">
            <div>
              <span>Selected for the next prediction</span>
              <h3>{selectedCustomCandidate.method}</h3>
              <code>{selectedCustomCandidate.linker || '(direct fusion; no linker)'}</code>
            </div>
            <dl>
              <div><dt>{domainALabel || 'Domain A'}</dt><dd>{domainALength} aa</dd></div>
              <div><dt>Linker</dt><dd>{selectedCustomCandidate.linker.length} aa</dd></div>
              <div><dt>{domainBLabel || 'Domain B'}</dt><dd>{selectedDomainBLength} aa</dd></div>
              <div><dt>Total fusion</dt><dd>{selectedCustomCandidate.fullSequence.length} aa</dd></div>
            </dl>
            <p><strong>Compute preflight:</strong> {localComputeRisk}</p>
            <p><strong>Evidence:</strong> {selectedCustomCandidate.evidenceQuality === 'control' ? 'Required computational control' : 'Deterministic candidate only; no learned compatibility rank'}</p>
            <details>
              <summary>Review complete fusion sequence</summary>
              <code className="full-sequence">{selectedCustomCandidate.fullSequence}</code>
            </details>
          </div>}
          {customCandidates.length > 0 && <div className="baseline-status">
            <h3>Comparison coverage</h3>
            <p><strong>Included now:</strong> direct fusion and fixed (GGGGS)3 controls.</p>
            <p><strong>Still required for a professional recommendation:</strong> random length-matched GS, geometry-only, natural-linker retrieval, and feature-based ranker baselines evaluated under identical structural rules.</p>
          </div>}
          <div className="page-actions candidates-actions">
            <button type="button" onClick={() => setComposerOpen(true)}>Edit inputs in sidebar</button>
            <button type="button" className="primary-btn" onClick={() => goToPage('structure')}>
              {selectedCustomCandidate ? 'Continue with selected candidate' : 'Open reference structure workspace'} <ArrowRight size={16} />
            </button>
          </div>
        </section></div>}

        {page === 'results' && <div className="app-page"><section className="section provenance" id="results">
          <div className="section-heading">
            <span>Selected fusion result</span>
            <h2>Structure prediction, confidence, and validation status</h2>
            <p>Review the measurements available for this candidate, the evidence source for each value, and the analyses that still require additional inputs or future tooling.</p>
          </div>
          {!selectedCustomCandidate ? <div className="results-empty">
            <h3>No user fusion has been selected</h3>
            <p>Use the input composer in the sidebar, enter both domain sequences, and create a prediction workspace before opening a result record.</p>
          </div> : <>
            <div className={`result-status ${predictionPhase.replace(' ', '-')}`}>
              <div><span>Prediction status</span><strong>{predictionPhase}</strong></div>
              <p>{viewerResult.fold?.message ?? (viewerResult.pdb ? 'A sequence-matched full-fusion PDB is loaded.' : 'No structure prediction has been started for this candidate.')}</p>
            </div>

            <div className="result-overview">
              <article>
                <span>Construct</span>
                <h3>{projectName || 'Untitled linker design'}</h3>
                <dl>
                  <div><dt>Domain order</dt><dd>{domainALabel || 'Domain A'} → linker → {domainBLabel || 'Domain B'}</dd></div>
                  <div><dt>Candidate</dt><dd>{selectedCustomCandidate.method}</dd></div>
                  <div><dt>Linker</dt><dd><code>{selectedCustomCandidate.linker || '(none)'}</code> ({selectedCustomCandidate.linker.length} aa)</dd></div>
                  <div><dt>Complete fusion</dt><dd>{selectedCustomCandidate.fullSequence.length} aa</dd></div>
                  <div><dt>Expression context</dt><dd>{expressionSystem}</dd></div>
                </dl>
              </article>
              <article>
                <span>Biological interpretation context</span>
                <dl>
                  <div><dt>Objective</dt><dd>{biologicalObjective || 'Not documented'}</dd></div>
                  <div><dt>Protected sites</dt><dd>{protectedSites || 'Not documented'}</dd></div>
                  <div><dt>Design evidence</dt><dd>{selectedCustomCandidate.evidenceQuality === 'control' ? 'Computational control' : 'Deterministic candidate; unranked'}</dd></div>
                  <div><dt>Structure source</dt><dd>{viewerResult.source === 'colabfold' ? 'Local ColabFold prediction' : viewerResult.source === 'uploaded' ? 'User-provided, sequence-matched PDB' : 'Not available'}</dd></div>
                </dl>
              </article>
            </div>

            <div className="confidence-section">
              <div className="section-heading compact-heading">
                <span>Structure confidence · source identified</span>
                <h3>{viewerResult.source === 'colabfold' || viewerResult.metrics
                  ? 'AlphaFold2 confidence produced through ColabFold'
                  : 'Values read from the uploaded PDB B-factor field'}</h3>
              </div>
              <div className="confidence-grid">
                <article><span>Whole fusion</span><strong>{displayConfidence(confidence?.overall)}</strong><small>mean pLDDT</small></article>
                <article><span>{domainALabel || 'Domain A'}</span><strong>{displayConfidence(confidence?.domainA)}</strong><small>mean pLDDT</small></article>
                <article><span>Linker</span><strong>{displayConfidence(confidence?.linker)}</strong><small>mean pLDDT</small></article>
                <article><span>{domainBLabel || 'Domain B'}</span><strong>{displayConfidence(confidence?.domainB)}</strong><small>mean pLDDT</small></article>
                <article><span>Global fold confidence</span><strong>{viewerResult.metrics?.ptm == null ? 'Not available' : viewerResult.metrics.ptm.toFixed(3)}</strong><small>AlphaFold pTM</small></article>
                <article><span>Relative domain uncertainty</span><strong>{displayConfidence(viewerResult.metrics?.meanInterdomainPae)}</strong><small>mean cross-domain PAE (Å)</small></article>
              </div>
              <p className="result-caveat"><strong>What generated this confidence:</strong> {viewerResult.metrics?.predictor ?? (viewerResult.source === 'colabfold' ? 'AlphaFold2 model executed by ColabFold.' : 'Unknown predictor; only numeric B-factor values were read from the uploaded coordinates.')} ColabFold is the execution/MSA pipeline; AlphaFold2 is the structure model. LinkerForge did not generate or calibrate these confidence values, and its future compatibility ranker is a separate model.</p>
              <p className="result-caveat">pLDDT estimates local-coordinate confidence. PAE estimates uncertainty in the position of one residue relative to another, and pTM summarizes predicted global topology. None proves expression, activity, linker flexibility, or a stable interdomain orientation. A flexible linker may appropriately have lower pLDDT.</p>
            </div>

            <div className="evidence-checks" aria-label="Structural evidence checklist">
              <h3>Evidence status and available alternatives</h3>
              <p className="evidence-legend">Green means measured for this candidate. Amber means an input or optional analysis is still needed. Blue identifies an implemented project-level alternative. Red is reserved for an actual failed or invalid check.</p>
              <div className="evidence-row"><span>Full-chain sequence match</span><strong className={viewerResult.pdb ? 'available' : 'pending'}>{viewerResult.pdb ? 'Passed' : 'Structure required'}</strong><small>The loaded PDB must exactly match the selected fusion.</small></div>
              <div className="evidence-row"><span>Per-region pLDDT</span><strong className={confidence ? 'available' : 'pending'}>{confidence ? 'Calculated' : 'Structure required'}</strong><small>Calculated from C-alpha B-factors in the loaded prediction.</small></div>
              <div className="evidence-row"><span>PAE / interdomain uncertainty</span><strong className={viewerResult.metrics?.meanInterdomainPae == null ? 'pending' : 'available'}>{viewerResult.metrics?.meanInterdomainPae == null ? 'Score JSON required' : `${viewerResult.metrics.meanInterdomainPae.toFixed(2)} Å`}</strong><small>{viewerResult.metrics?.meanInterdomainPae == null ? 'Load the matching ColabFold score JSON or complete a desktop fold.' : 'Mean of both Domain A↔Domain B PAE blocks; lower is more certain. No acceptance threshold is frozen.'}</small></div>
              <div className="evidence-row"><span>Chain continuity / severe Cα overlap pre-screen</span><strong className={!coordinateAudit ? 'pending' : coordinateAudit.chainBreakCount === 0 && coordinateAudit.severeCaOverlapCount === 0 ? 'available' : 'failed'}>{coordinateAudit ? `${coordinateAudit.chainBreakCount} breaks · ${coordinateAudit.severeCaOverlapCount} overlaps` : 'Structure required'}</strong><small>{coordinateAudit ? `Maximum adjacent Cα distance ${coordinateAudit.maxAdjacentCaDistance?.toFixed(2)} Å. This deterministic coordinate check is not Rosetta all-atom closure or clash scoring.` : 'Load a fusion structure to run the coordinate pre-screen.'}</small></div>
              <div className="evidence-row"><span>Reference-domain structural preservation</span><strong className={domainPreservation?.domainA != null && domainPreservation?.domainB != null ? 'available' : 'pending'}>{domainPreservation?.domainA != null && domainPreservation?.domainB != null ? `A ${domainPreservation.domainA.toFixed(2)} Å · B ${domainPreservation.domainB.toFixed(2)} Å` : 'Domain references required'}</strong><small>Sequence-matched Cα RMSD after rigid-body superposition. Load each domain’s own trusted reference on the Structure page; source classification remains part of interpretation.</small></div>
              <div className="evidence-row"><span>Protected-site mapping</span><strong className={protectedSiteAudit.errors.length ? 'failed' : protectedSiteClearance ? 'available' : 'pending'}>{protectedSiteAudit.errors.length ? `Invalid: ${protectedSiteAudit.errors.join(', ')}` : protectedSiteClearance ? `${protectedSiteClearance.sites.length} mapped · ${protectedSiteClearance.minimumDistanceAngstrom?.toFixed(2) ?? 'N/A'} Å minimum linker Cα distance` : protectedSiteAudit.residueIndices.length ? `${protectedSiteAudit.residueIndices.length} mapped · structure required` : 'Optional sites not supplied'}</strong><small>{protectedSiteClearance ? `Coordinate pre-screen at a documented ${protectedSiteClearance.contactCutoffAngstrom.toFixed(1)} Å contact cutoff. This is not solvent accessibility, ligand/interface annotation, or an all-atom obstruction calculation.` : 'Enter domain-relative sites such as A:48 or B:12 and load a sequence-matched full-fusion structure.'}</small></div>
              <div className="evidence-row"><span>Domain surface burial</span><strong className="alternative">CPU report workflow available</strong><small>The report workflow measures residue-level loss of solvent-accessible surface area (SASA) caused by the linker and the other domain at fixed predicted coordinates. This is not calculated automatically for the current candidate. Run the report workflow and supply functional-site annotations to assess whether relevant sites are buried; surface burial alone does not establish loss of function.</small></div>
              <div className="evidence-row"><span>All-atom geometry assessment</span><strong className="alternative">Partial workflow available</strong><small>OpenMM minimization can expose severe force-field strain and relax a finalist, but it has not run for this candidate and is not a loop-reconstruction substitute. The Cα continuity screen runs now; a pinned Reduce/Probe report supplies a complementary hydrogen-aware overlap check, while Rosetta remains future scope.</small></div>
              <div className="evidence-row"><span>Matched control comparison</span><strong className="alternative">Public reference panel available</strong><small>A seven-construct 1UBQ–2GB1 panel has matched ColabFold settings. A durable three-pair public panel is configured to extend this sensitivity check; the current user design still requires its own matched controls.</small></div>
              <div className="evidence-row"><span>Energy minimization / MD</span><strong className="future">Future scope · on hold</strong><small>The CUDA OpenMM runner completed only a three-replica 50 ps execution smoke test. Production 3 × 50 ns MD is not queued to start automatically: a future explicit launch must freeze the finalist, repeat disk/dependency/CUDA/input preflights, use a new durable run directory, and complete replicate-aware convergence and change analysis. The storage-aware protocol retains 10 ps metrics, stores DCD every 50 ps and projects 31.76 GPU hours. Neither the smoke run nor a future simulation proves biological stability.</small></div>
            </div>

            <div className="result-provenance">
              <h3>Run record</h3>
              <p><strong>Run directory:</strong> {viewerResult.fold?.runDir ?? 'No local ColabFold run directory for this result.'}</p>
              <p><strong>Sequence hash:</strong> {viewerResult.fold?.sequenceSha256 ?? 'Not recorded for a manual upload.'}</p>
              <p><strong>Privacy route:</strong> {viewerResult.source === 'colabfold' ? 'Public ColabFold/MMseqs2 MSA service acknowledged for this session.' : 'No local MSA submission recorded in this result view.'}</p>
              <p><strong>Claim boundary:</strong> This is a predicted structure for an unranked candidate or control, not a validated linker recommendation.</p>
            </div>

            <div className="result-provenance" aria-label="Tamper-evident local history">
              <h3>Checksummed run history</h3>
              <p><strong>Local entries:</strong> {researchHistory.length} · <strong>hash-chain verification:</strong> {historyVerified == null ? 'checking' : historyVerified ? 'passed' : 'failed—do not trust this history'}</p>
              <p>This browser history is append-only through the interface and detects later changes through a SHA-256 hash chain. It is not a digital signature and is not a substitute for institutional records.</p>
              <div className="export-actions">
                <button type="button" onClick={() => void recordResultInHistory()}>Append current result to history</button>
                <button type="button" onClick={() => downloadJson('linkerforge-history.v1.json', researchHistory)}>Export history JSON</button>
              </div>
            </div>

            <div className="evidence-checks" aria-label="Comparator audit coverage">
              <h3>Comparator coverage and future integrations</h3>
              <div className="evidence-row"><span>AlphaFold2 ensemble agreement</span><strong className="future">Future scope · step 1</strong><small>Compare all five sequence-matched models after separate domain superpositions, reporting pairwise domain RMSD, linker/junction displacement, interdomain orientation, pLDDT, pTM and cross-domain PAE. This reveals within-family disagreement but is not an independent-predictor comparison.</small></div>
              <div className="evidence-row"><span>Reference-domain lDDT / TM-score</span><strong className="future">Future scope · step 2</strong><small>Map each submitted domain to its own provenance-bearing reference, then run pinned OpenStructure lDDT and US-align TM-score calculations with explicit residue coverage. These complement RMSD; thresholds must be frozen on development data before any pass/fail claim.</small></div>
              <div className="evidence-row"><span>Automatic functional-site mapping</span><strong className="future">Future scope · step 3</strong><small>Resolve dated UniProt/PDB annotations, map source numbering through an auditable sequence alignment, and measure ligand/interface/site distances, contacts and SASA changes. Ambiguous mappings remain unknown; predicted obstruction is a risk flag, not proof of lost function.</small></div>
              <div className="evidence-row"><span>DSSP junction analysis</span><strong className="future">Future scope · step 4</strong><small>Run a pinned mkdssp build on references, controls and every retained ensemble member, then compare per-residue secondary structure and helix/strand persistence across both attachment junctions. Missing residues and model disagreement must remain visible.</small></div>
              <div className="evidence-row"><span>Post-OpenMM change analysis</span><strong className="future">Future scope · step 5</strong><small>After any future replicated production run, compare starting, minimized and trajectory structures using domain-aligned RMSD/RMSF, radius of gyration, linker distance, orientation, contacts, protected-site accessibility, lDDT/TM-score, DSSP and replicate/block convergence. Even 3 × 50 ns remains computational plausibility evidence.</small></div>
              <div className="evidence-row"><span>Published repository datasets</span><strong className="available">Byte-validated</strong><small>The historical DFL split contains 60 train, 28 development, and 9 auxiliary proteins. The separate 27,766-record mined/natural collection is a source inventory—not a validated fusion-compatibility training set; ranker training remains unauthorized.</small></div>
              <div className="evidence-row"><span>Historical hand-feature ranker</span><strong className="future">Future compatibility study</strong><small>Reproduced against rule-generated synthetic corruptions. It is not the frozen ESM-2/geometric-GNN model, is not calibrated to biological success, and cannot supply a candidate score in this workflow.</small></div>
              <div className="evidence-row"><span>DFLpred published detector</span><strong className="available">28 constructs run</strong><small>Four public domain contexts × seven inserted/direct linkers. It detects flexible-linker propensity, not fusion compatibility or function; results are application-level sensitivity evidence, not this candidate’s validation.</small></div>
              <div className="evidence-row"><span>Retrospective CAID3 common-set benchmark</span><strong className="alternative">Spent result available</strong><small>On the same 27 linker-positive proteins, the best historical LinkerForge head scored APS 0.089 / Fmax 0.165, below IPA-AF2-Linker 0.482 / 0.491, LINKER-Pred-Lite 0.429 / 0.473, and DisoFLAG-DFL 0.409 / 0.464. This public retrospective is spent: it cannot be reused for tuning or model selection, does not measure no-linker specificity, and is not evidence of fusion-linker compatibility.</small></div>
              <div className="evidence-row"><span>LINKER-Pred / LINKER-Pred Lite</span><strong className="future">Future reproducibility work</strong><small>The official workflow is available through separate Embedding and predictor packages, with Docker recommended by its authors. It still requires model assets and a firewall-safe common-target evaluation.</small></div>
              <div className="evidence-row"><span>TransDFL / DisoFLAG-DFL</span><strong className="future">Future detector comparison</strong><small>Public implementations exist, but a fair rerun needs pinned environments and a protected-label-safe target set. Spent CAID3 scores cannot guide model selection.</small></div>
              <div className="evidence-row"><span>Alternative structure predictors</span><strong className="future">Future model-family sensitivity</strong><small>ESMFold offers an MSA-free architecture; OpenFold reproduces the AlphaFold2 family; Boltz provides a newer open model family. They require separate validated environments and do not themselves prove linker compatibility.</small></div>
              <div className="evidence-row"><span>Reduce / Probe component screen</span><strong className="alternative">7 public structures screened</strong><small>Reduce 4.15 and Probe 2.18 are locked in a project-owned environment. The GGGGS10 acceptance prediction contained 48 unique severe-overlap pairs (22.22 per 1000 explicit atoms; worst overlap 1.473 Å). No pass threshold is frozen. This transparent pair rate is not the canonical MolProbity clashscore and does not replace Rosetta loop-closure sampling.</small></div>
            </div>
          </>}
          <div className="page-actions results-actions">
            <a href="#structure" onClick={() => goToPage('structure')}>Back to structure</a>
            {selectedCustomCandidate && <button type="button" onClick={() => downloadText('linkerforge-selected-fusion.fasta', fusionFasta(selectedCustomCandidate.fullSequence))}>Export selected FASTA</button>}
            {selectedCustomCandidate && <button type="button" onClick={() => void exportResultRecord()}>Export versioned result JSON</button>}
            {selectedCustomCandidate && <button type="button" onClick={() => void exportResearchBundle()}><Download size={14} /> Export checksummed research ZIP</button>}
            {viewerResult.pdb && <button type="button" onClick={() => downloadText('linkerforge-selected-fusion.pdb', viewerResult.pdb!)}>Export structure PDB</button>}
            <button type="button" onClick={startNewDesign}>Start another design</button>
          </div>
        </section></div>}
      </main>

        </div>
      </div>
    </div>
  )
}

export default App
