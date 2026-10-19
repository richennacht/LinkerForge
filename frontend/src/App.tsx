import { useEffect, useState } from 'react'
import {
  Activity,
  ArrowRight,
  CheckCircle2,
  CircleDot,
  Database,
  Download,
  GitBranch,
  Play,
  ShieldCheck,
  XCircle,
} from 'lucide-react'
import FusionStructureViewer from './FusionStructureViewer'
import { buildTemplateCandidates, fusionFasta, parseProteinSequence, validateFusionPdb, type FusionCandidate } from './fusion'
import { getComputeStatus, getFoldStatus, inDesktop, startLocalFold, type ComputeStatus, type FoldStatus } from './desktopFold'

const stages = [
  ['01', 'Check inputs', 'Confirm both protein files and protected sites are usable.'],
  ['02', 'Measure joining points', 'Describe the shape and surroundings of each possible connection.'],
  ['03', 'Recommend a linker', 'Rank connection direction, linker type, and length.'],
  ['04', 'Predict and compare', 'Predict the combined shape and compare it with the originals.'],
]

const referenceStructures = {
  sfGFP: {
    label: 'sfGFP',
    pdbId: '2B3P',
    description: 'Experimentally determined superfolder green fluorescent protein structure (238 aa beta-barrel; C-terminal port Lys238, rSASA = 0.62).',
  },
  mCherry: {
    label: 'mCherry',
    pdbId: '2H5Q',
    description: 'Experimentally determined mCherry monomeric red fluorescent protein structure (236 aa beta-barrel; N-terminal port Met1, rSASA = 0.58).',
  },
} as const

type ReferenceStructure = keyof typeof referenceStructures
type RigidityMode = 'flexible' | 'rigid'

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

function ProteinViewer({
  candidate,
  domainALength,
}: {
  candidate: FusionCandidate | null
  domainALength: number
}) {
  const [selected, setSelected] = useState<ReferenceStructure>('sfGFP')
  const [mode, setMode] = useState<'fusion' | 'reference'>('fusion')
  const [fullscreen, setFullscreen] = useState(false)
  const [structures, setStructures] = useState<Record<string, string>>({})
  const [uploadError, setUploadError] = useState<string | null>(null)
  const [compute, setCompute] = useState<ComputeStatus | null>(null)
  const [fold, setFold] = useState<FoldStatus | null>(null)
  const [foldSequence, setFoldSequence] = useState<string | null>(null)
  const selectedSequence = candidate?.fullSequence ?? ''
  const pdb = structures[selectedSequence] ?? null
  const reference = referenceStructures[selected]

  useEffect(() => {
    if (!inDesktop()) return
    void getComputeStatus().then(setCompute).catch(() => setCompute(null))
  }, [])

  useEffect(() => {
    if (!foldSequence) return
    const timer = window.setInterval(() => {
      void getFoldStatus().then((status) => {
        setFold(status)
        if (status.phase === 'complete' && status.pdb) {
          try {
            validateFusionPdb(status.pdb, foldSequence)
            setStructures((previous) => ({ ...previous, [foldSequence]: status.pdb! }))
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

  const handleLocalFold = async () => {
    if (!candidate) return
    try {
      const status = await startLocalFold(candidate.fullSequence)
      setFold(status)
      setFoldSequence(candidate.fullSequence)
      setUploadError(null)
    } catch (cause) {
      setUploadError(cause instanceof Error ? cause.message : 'Could not start ColabFold.')
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
      setStructures((previous) => ({ ...previous, [candidate.fullSequence]: content }))
      setUploadError(null)
      setMode('fusion')
    } catch (cause) {
      setUploadError(cause instanceof Error ? cause.message : 'Could not read the PDB file.')
    }
  }

  const pml = candidate ? [
    '# Run after loading the matching full-fusion PDB in PyMOL.',
    'hide everything, all',
    'show cartoon, all',
    `select domain_a, resi 1-${domainALength}`,
    `select linker, resi ${domainALength + 1}-${domainALength + candidate.linker.length}`,
    `select domain_b, resi ${domainALength + candidate.linker.length + 1}-${candidate.fullSequence.length}`,
    'color green, domain_a',
    'color magenta, linker',
    'color red, domain_b',
    'zoom all',
  ].join('\n') + '\n' : ''

  return (
    <section className={`section viewer-section ${fullscreen ? 'viewer-fullscreen' : ''}`} id="viewer" aria-label="Structure viewer">
      <div className="viewer-header">
        <div className="section-heading">
          <span>02 / Full-fusion structure viewer</span>
          <h2>{candidate ? 'Inspect both domains and the connecting linker' : 'Inspect reference structures'}</h2>
        </div>
        <button type="button" className="viewer-expand" onClick={() => setFullscreen((value) => !value)}>
          {fullscreen ? 'Close full screen' : 'Open full screen'}
        </button>
      </div>
      <div className="viewer-mode-row" role="group" aria-label="Structure source">
        <button type="button" className={mode === 'fusion' ? 'active' : ''} onClick={() => setMode('fusion')}>
          Full fusion {pdb ? '· PDB loaded' : '· awaits prediction'}
        </button>
        <button type="button" className={mode === 'reference' ? 'active' : ''} onClick={() => setMode('reference')}>
          Individual experimental references
        </button>
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
                {inDesktop() && <button type="button" disabled={!compute?.predictorAvailable || fold?.phase === 'running'}
                  onClick={() => void handleLocalFold()}>
                  <Play size={14} /> {fold?.phase === 'running' ? 'Folding locally…' : 'Fold full fusion on this GPU'}
                </button>}
                {pdb && <button type="button" onClick={() => downloadText('linkerforge-fusion.pml', pml)}>
                  <Download size={14} /> PyMOL colouring script
                </button>}
              </>
            ) : (
              <p>Enter both protein sequences in the design studio to construct a full fusion first.</p>
            )}
          </div>
          {uploadError && <p className="structure-error" role="alert">{uploadError}</p>}
          {inDesktop() && <div className="desktop-compute" aria-live="polite">
            <strong>Local compute:</strong> {compute?.gpu ?? 'Checking GPU…'} · ColabFold {compute?.predictorAvailable ? `available via ${compute.predictorBackend}` : 'not installed'}.
            {fold && <p>{fold.phase}: {fold.message}</p>}
            {fold?.runDir && <small>Run files: {fold.runDir}</small>}
          </div>}
          <p className="viewer-disclaimer">The uploaded PDB is user-provided and sequence-checked, not independently verified as an AlphaFold result. Linker coordinates appear only when the full fusion has actually been predicted.</p>
        </>
      ) : (
        <>
          <div className="structure-tabs" aria-label="Experimental reference structure">
            {(Object.keys(referenceStructures) as ReferenceStructure[]).map((key) => (
              <button key={key} type="button" className={selected === key ? 'active' : ''}
                aria-pressed={selected === key} onClick={() => setSelected(key)}>
                {referenceStructures[key].label} <small>PDB {referenceStructures[key].pdbId}</small>
              </button>
            ))}
          </div>
          <div className="viewer-stage">
            <iframe key={reference.pdbId} src={`https://molstar.org/viewer/?pdb=${reference.pdbId}`}
              title={`${reference.label} experimental structure in Mol*`} allowFullScreen />
          </div>
          <p className="viewer-disclaimer">{reference.description} This is an individual experimental reference, not a predicted fusion.</p>
        </>
      )}
    </section>
  )
}
function App() {
  const [domainAInput, setDomainAInput] = useState('')
  const [domainBInput, setDomainBInput] = useState('')
  const [domainALength, setDomainALength] = useState(0)
  const [customCandidates, setCustomCandidates] = useState<FusionCandidate[]>([])
  const [selectedCustomIndex, setSelectedCustomIndex] = useState(0)
  const [customError, setCustomError] = useState<string | null>(null)
  const [rigidity, setRigidity] = useState<RigidityMode>('flexible')
  const [cleavable, setCleavable] = useState<boolean>(false)
  const [targetLength, setTargetLength] = useState<number>(15)
  const [candidates, setCandidates] = useState<RankedCandidateRow[]>(
    PRESET_CANDIDATES.flexible_non_cleavable
  )
  const [selectedCandidateSeq, setSelectedCandidateSeq] = useState<string>('GGGGSGGGGSGGGGS')
  const [apiConnected, setApiConnected] = useState<boolean>(false)

  const handleRunDesign = async (
    nextRigidity: RigidityMode = rigidity,
    nextCleavable: boolean = cleavable,
    nextLength: number = targetLength
  ) => {
    if (domainAInput.trim() || domainBInput.trim()) {
      try {
        const domainA = parseProteinSequence(domainAInput, 'Domain A')
        const domainB = parseProteinSequence(domainBInput, 'Domain B')
        const generated = buildTemplateCandidates(domainA, domainB, nextLength, nextRigidity, nextCleavable)
        setDomainALength(domainA.length)
        setCustomCandidates(generated)
        setSelectedCustomIndex(0)
        setCustomError(null)
        setApiConnected(false)
      } catch (cause) {
        setCustomCandidates([])
        setCustomError(cause instanceof Error ? cause.message : 'Invalid sequence input.')
      }
      return
    }
    setCustomCandidates([])
    setCustomError(null)
    try {
      const resp = await fetch('http://127.0.0.1:8765/api/design', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          rigidity: nextRigidity,
          cleavable: nextCleavable,
          protease_family: 'TEV',
          target_length: nextLength,
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
        return
      }
    } catch {
      // The local scientific dataset/API is optional for the labelled UI demonstration.
    }
    setApiConnected(false)
    const key = nextCleavable
      ? 'flexible_cleavable'
      : nextRigidity === 'rigid'
        ? 'rigid_non_cleavable'
        : 'flexible_non_cleavable'
    const rows = PRESET_CANDIDATES[key]
    setCandidates(rows)
    if (rows.length > 0) {
      setSelectedCandidateSeq(rows[0].sequence)
    }
  }

  return (
    <div className="app-shell">
      <header className="topbar">
        <a className="brand" href="#top" aria-label="LinkerForge home">
          <span className="brand-mark">LF</span>
          <span>LinkerForge</span>
        </a>
        <nav aria-label="Main navigation">
          <a href="#workflow">Workflow</a>
          <a href="#studio">Design studio</a>
          <a href="#viewer">3D viewer</a>
          <a href="#runs">Ranked candidates</a>
          <a href="#provenance">Provenance</a>
        </nav>
        <span className="release">{inDesktop() ? 'Local desktop' : 'Computational prototype'} · {customCandidates.length ? 'template mode' : 'example mode'}</span>
      </header>

      <main id="top">
        <section className="hero">
          <div className="eyebrow">
            <CircleDot size={14} /> Computational research workspace · local demonstration
          </div>
          <h1>Design protein connections that can be inspected, repeated, and challenged.</h1>
          <p className="hero-copy">
            Enter two protein sequences to build complete fusion candidates, then inspect a
            sequence-matched predicted structure. Sequence-only templates are available locally;
            {inDesktop()
              ? ' this desktop app can launch a full-fusion prediction when LocalColabFold is separately installed.'
              : ' learned ranking and AlphaFold/ColabFold prediction require separate scientific compute.'}
          </p>
          <div className="hero-actions">
            <a className="primary-btn" href="#studio">
              Start a design <ArrowRight size={17} />
            </a>
            <a href="#viewer">Open full-fusion viewer</a>
          </div>
          <p className="notice">
            <ShieldCheck size={17} /> Demonstration only. The current site uses labelled example
            results and makes no claim of biological function.
          </p>
        </section>

        <section className="metrics" aria-label="Project status">
          <article>
            <span>Reported candidate inventory</span>
            <strong>27,766 records</strong>
            <small>27,595 model-scored gaps; source data not bundled here</small>
          </article>
          <article>
            <span>Attachment-site diagnostic</span>
            <strong>Constructed task</strong>
            <small>Terminal-vs-interior scores are not independent site validation</small>
          </article>
          <article>
            <span>Ranking diagnostic</span>
            <strong>Synthetic controls</strong>
            <small>No matched published-method or experimental comparison yet</small>
          </article>
          <article>
            <span>Traceability</span>
            <strong>Deterministic</strong>
            <small>Hard-rejection scorecard (no fake probabilities)</small>
          </article>
        </section>

        <section className="section" id="workflow">
          <div className="section-heading">
            <span>01 / Workflow</span>
            <h2>One traceable path from input to result</h2>
          </div>
          <div className="stage-grid">
            {stages.map(([number, title, body]) => (
              <article className="stage" key={number}>
                <span>{number}</span>
                <h3>{title}</h3>
                <p>{body}</p>
              </article>
            ))}
          </div>
        </section>

        <section className="section studio-section" id="studio">
          <div className="section-heading">
            <span>01B / Interactive End-to-End Linker Design Studio</span>
            <h2>Enter Domain A and Domain B to construct a full fusion</h2>
          </div>
          <div className="studio-controls">
            <div className="sequence-inputs">
              <label htmlFor="domain-a">Domain A amino-acid sequence (FASTA accepted)</label>
              <textarea id="domain-a" value={domainAInput} onChange={(event) => { setDomainAInput(event.target.value); setCustomCandidates([]) }}
                placeholder="Paste the first protein sequence" rows={5} />
              <label htmlFor="domain-b">Domain B amino-acid sequence (FASTA accepted)</label>
              <textarea id="domain-b" value={domainBInput} onChange={(event) => { setDomainBInput(event.target.value); setCustomCandidates([]) }}
                placeholder="Paste the second protein sequence" rows={5} />
              <small>Leave both blank to explore the labelled sfGFP–mCherry preset. Custom inputs use template enumeration only; no compatibility ranking or structural prediction is implied.</small>
            </div>
            <div className="control-group">
              <label>Rigidity Axis (Component 3)</label>
              <div className="pill-row">
                <button
                  type="button"
                  className={rigidity === 'flexible' ? 'active' : ''}
                  onClick={() => {
                    setRigidity('flexible')
                    handleRunDesign('flexible', cleavable, targetLength)
                  }}
                >
                  Flexible Coil (GGGGS / Waldo / Mined DFL)
                </button>
                <button
                  type="button"
                  className={rigidity === 'rigid' ? 'active' : ''}
                  onClick={() => {
                    setRigidity('rigid')
                    setCleavable(false)
                    setTargetLength(18)
                    handleRunDesign('rigid', false, 18)
                  }}
                >
                  Rigid Spacer (EAAAK Helical / PAPAP Rod)
                </button>
              </div>
            </div>

            <div className="control-group">
              <label>Protease Cleavability (Hubbard 1994 Spacer Rule)</label>
              <div className="pill-row">
                <button
                  type="button"
                  className={!cleavable ? 'active' : ''}
                  onClick={() => {
                    setCleavable(false)
                    handleRunDesign(rigidity, false, targetLength)
                  }}
                >
                  Intact (Non-cleavable)
                </button>
                <button
                  type="button"
                  className={cleavable ? 'active' : ''}
                  onClick={() => {
                    setCleavable(true)
                    setRigidity('flexible')
                    setTargetLength(19)
                    handleRunDesign('flexible', true, 19)
                  }}
                >
                  TEV Cleavable (ENLYFQS with GS spacers)
                </button>
              </div>
            </div>

            <div className="control-group">
              <label>Target Linker Length: {targetLength} aa (3D Path = 28.4 Å)</label>
              <input
                type="range"
                min={8}
                max={30}
                value={targetLength}
                onChange={(e) => setTargetLength(Number(e.target.value))}
              />
            </div>

            <button
              type="button"
              className="run-pipeline-btn"
              onClick={() => handleRunDesign(rigidity, cleavable, targetLength)}
            >
              <Play size={15} /> {domainAInput.trim() || domainBInput.trim()
                ? 'Generate full-fusion template candidates'
                : apiConnected ? 'Re-Rank Candidates (Live Local API)' : 'Load Preset Example Candidates'}
            </button>
            {customError && <p className="structure-error" role="alert">{customError}</p>}
            {!apiConnected && customCandidates.length === 0 && !domainAInput.trim() && !domainBInput.trim() && (
              <p className="notice">
                Example mode: controls switch preset results; the length slider does not recalculate
                scores until the dataset-backed local API is running.
              </p>
            )}
          </div>
        </section>

        <ProteinViewer candidate={customCandidates[selectedCustomIndex] ?? null} domainALength={domainALength} />

        <section className="section split" id="runs">
          <div>
            <div className="section-heading">
              <span>03 / Ranked Fusion Candidates &amp; Baselines</span>
              <h2>Results visible to the whole team</h2>
            </div>
            <p>
              {customCandidates.length > 0
                ? 'These are unranked sequence templates, not predicted successes. Choose one, export its full-fusion FASTA, and load a matching predicted PDB to inspect the linker.'
                : 'These are labelled preset examples for the fixed sfGFP–mCherry pair. They are not predictions for user-entered proteins.'}
            </p>
          </div>
          {customCandidates.length > 0 ? <div className="custom-candidate-list" aria-label="Custom fusion candidates">
            {customCandidates.map((candidate, index) => <button type="button" key={candidate.id}
              className={selectedCustomIndex === index ? 'active' : ''} onClick={() => setSelectedCustomIndex(index)}>
              <strong>{candidate.method} · {candidate.linker.length} aa</strong>
              <code>{candidate.linker}</code>
              <small>Full fusion: {candidate.fullSequence.length} residues · unranked template</small>
            </button>)}
          </div> : <div className="run-table" role="table" aria-label="Example runs">
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
          </div>}
        </section>

        <section className="section provenance" id="provenance">
          <div className="section-heading">
            <span>04 / Provenance</span>
            <h2>Every answer carries its history</h2>
          </div>
          <div className="provenance-grid">
            <article>
              <GitBranch />
              <div>
                <strong>Code version</strong>
                <p>Exact source-code revision and working configuration (`codex/classifier-foundation`).</p>
              </div>
            </article>
            <article>
              <Database />
              <div>
                <strong>Reported input inventory</strong>
                <p>
                  The source report lists 27,595 model-scored candidate gaps and 171
                  natural/template records. The dataset is not bundled with this web demo.
                </p>
              </div>
            </article>
            <article>
              <Activity />
              <div>
                <strong>Calculation record</strong>
                <p>
                  Preset example ranking unless the local API is connected. The earlier
                  conformal precision claim is not certified for proteome-wide use.
                </p>
              </div>
            </article>
            <article>
              <CheckCircle2 />
              <div>
                <strong>Decision record</strong>
                <p>Every passed rule, failed rule, score, and stated limitation (`SPEC.md` §5.8).</p>
              </div>
            </article>
          </div>
        </section>
      </main>

      <footer>
        <strong>LinkerForge</strong>
        <span>Encoding-guided design of structurally plausible fusion proteins</span>
        <span>Computational predictions require future experimental validation.</span>
      </footer>
    </div>
  )
}

export default App
