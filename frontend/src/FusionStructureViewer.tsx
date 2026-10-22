import { useEffect, useMemo, useRef, useState } from 'react'
import {
  Box,
  Camera,
  ChevronDown,
  Crosshair,
  Download,
  Eye,
  EyeOff,
  Focus,
  Info,
  Layers3,
  MousePointer2,
  Pause,
  Play,
  RotateCcw,
  Ruler,
  Search,
} from 'lucide-react'
import type {
  AtomSelectionSpec,
  AtomSpec,
  AtomStyleSpec,
  ColorschemeSpec,
  GLViewer,
  SurfaceStyleSpec,
} from '3dmol'
import { validateFusionPdb } from './fusion'
import { validateDomainReferencePdb } from './structureEvidence'

interface Props {
  pdb: string | null
  fullSequence: string
  domainALength: number
  linkerLength: number
  mode?: 'fusion' | 'reference'
  referenceColor?: string
  ariaLabel?: string
  structureId?: string
  structureTitle?: string
  structureSource?: string
  sourceUrl?: string
}

type Representation = 'cartoon' | 'sticks' | 'ball-stick' | 'spheres' | 'surface'
type ColourMode = 'domains' | 'confidence' | 'spectrum'
type CanvasTheme = 'dark' | 'light'
type ComponentKey = 'domainA' | 'linker' | 'domainB' | 'reference'

interface ComponentDefinition {
  key: ComponentKey
  label: string
  range: string
  length: number
  selection: AtomSelectionSpec
  color: string
}

interface ViewerStyleOptions {
  representation: Representation
  colourMode: ColourMode
  mode: 'fusion' | 'reference'
  referenceColor: string
  domainALength: number
  linkerLength: number
  sequenceLength: number
  hiddenComponents: Set<ComponentKey>
  selectedResidue: number | null
}

interface Measurement {
  start: number
  end: number
  distance: number
}

const CONFIDENCE_COLOURS = {
  prop: 'b',
  gradient: 'linear',
  min: 0,
  max: 100,
  colors: ['#ff7d45', '#ffdb13', '#65cbf3', '#0053d6'],
}

function representationStyle(representation: Representation, colour: string | Exclude<ColorschemeSpec, string>): AtomStyleSpec {
  const colourSpec = typeof colour === 'string' ? { color: colour } : { colorscheme: colour }
  if (representation === 'sticks') return { stick: { ...colourSpec, radius: 0.22 } }
  if (representation === 'ball-stick') return {
    stick: { ...colourSpec, radius: 0.18 },
    sphere: { ...colourSpec, scale: 0.28 },
  }
  if (representation === 'spheres') return { sphere: { ...colourSpec, scale: 0.72 } }
  return { cartoon: colourSpec }
}

function componentDefinitions(options: Pick<ViewerStyleOptions, 'mode' | 'referenceColor' | 'domainALength' | 'linkerLength' | 'sequenceLength'>): ComponentDefinition[] {
  if (options.mode === 'reference') return [{
    key: 'reference',
    label: 'Sequence-matched polymer',
    range: `1–${options.sequenceLength}`,
    length: options.sequenceLength,
    selection: {},
    color: options.referenceColor,
  }]

  const linkerStart = options.domainALength + 1
  const linkerEnd = options.domainALength + options.linkerLength
  const definitions: ComponentDefinition[] = [{
    key: 'domainA',
    label: 'Domain A',
    range: `1–${options.domainALength}`,
    length: options.domainALength,
    selection: { resi: `1-${options.domainALength}` },
    color: '#3acf78',
  }]
  if (options.linkerLength > 0) definitions.push({
    key: 'linker',
    label: 'Linker',
    range: `${linkerStart}–${linkerEnd}`,
    length: options.linkerLength,
    selection: { resi: `${linkerStart}-${linkerEnd}` },
    color: '#ff4dba',
  })
  definitions.push({
    key: 'domainB',
    label: 'Domain B',
    range: `${linkerEnd + 1}–${options.sequenceLength}`,
    length: Math.max(0, options.sequenceLength - linkerEnd),
    selection: { resi: `${linkerEnd + 1}-${options.sequenceLength}` },
    color: '#f66b62',
  })
  return definitions
}

function componentForResidue(residue: number, mode: 'fusion' | 'reference', domainALength: number, linkerLength: number) {
  if (mode === 'reference') return 'reference'
  if (residue <= domainALength) return 'domainA'
  if (residue <= domainALength + linkerLength) return 'linker'
  return 'domainB'
}

function applyViewerStyle(viewer: GLViewer, options: ViewerStyleOptions) {
  viewer.removeAllSurfaces()
  viewer.setStyle({}, {})
  const components = componentDefinitions(options).filter((item) => !options.hiddenComponents.has(item.key))
  const mappedColours = options.colourMode === 'confidence'
    ? CONFIDENCE_COLOURS
    : { prop: 'resi', gradient: 'roygb', min: 1, max: options.sequenceLength }

  for (const item of components) {
    if (options.representation === 'surface') {
      const style: SurfaceStyleSpec = options.colourMode === 'domains'
        ? { opacity: 0.84, color: item.color }
        : { opacity: 0.84, colorscheme: mappedColours }
      void viewer.addSurface('VDW', style, item.selection)
      viewer.setStyle(item.selection, { cartoon: { color: '#d5ded8', opacity: 0.24 } })
    } else {
      const colour = options.colourMode === 'domains' ? item.color : mappedColours
      viewer.setStyle(item.selection, representationStyle(options.representation, colour))
    }
  }

  if (options.selectedResidue !== null) {
    viewer.addStyle({ resi: options.selectedResidue }, {
      stick: { color: '#ffd94a', radius: 0.25 },
      sphere: { color: '#ffd94a', scale: 0.38 },
    })
  }
  viewer.render()
}

function downloadViewerImage(viewer: GLViewer, filename: string) {
  const anchor = document.createElement('a')
  anchor.href = viewer.pngURI()
  anchor.download = filename
  document.body.appendChild(anchor)
  anchor.click()
  anchor.remove()
}

function downloadCoordinates(pdb: string, filename: string) {
  const url = URL.createObjectURL(new Blob([pdb], { type: 'chemical/x-pdb' }))
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(url)
}

function atomPoint(atom: AtomSpec) {
  if (atom.x == null || atom.y == null || atom.z == null) return null
  return { x: atom.x, y: atom.y, z: atom.z }
}

function distanceBetween(first: AtomSpec, second: AtomSpec) {
  const a = atomPoint(first)
  const b = atomPoint(second)
  if (!a || !b) return null
  return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)
}

export default function FusionStructureViewer({
  pdb,
  fullSequence,
  domainALength,
  linkerLength,
  mode = 'fusion',
  referenceColor = '#55a7ff',
  ariaLabel = 'Predicted full fusion structure',
  structureId,
  structureTitle,
  structureSource,
  sourceUrl,
}: Props) {
  const containerRef = useRef<HTMLDivElement | null>(null)
  const viewerRef = useRef<GLViewer | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [representation, setRepresentation] = useState<Representation>('cartoon')
  const [colourMode, setColourMode] = useState<ColourMode>(mode === 'fusion' ? 'domains' : 'spectrum')
  const [canvasTheme, setCanvasTheme] = useState<CanvasTheme>('dark')
  const [spinning, setSpinning] = useState(false)
  const [viewerReady, setViewerReady] = useState(0)
  const [selectedResidue, setSelectedResidue] = useState<number | null>(null)
  const [hoveredResidue, setHoveredResidue] = useState<string>('Hover or select a residue')
  const [hiddenComponents, setHiddenComponents] = useState<Set<ComponentKey>>(new Set())
  const [selectionStart, setSelectionStart] = useState(1)
  const [selectionEnd, setSelectionEnd] = useState(Math.max(1, fullSequence.length))
  const [measureStart, setMeasureStart] = useState(1)
  const [measureEnd, setMeasureEnd] = useState(Math.max(1, fullSequence.length))
  const [measurement, setMeasurement] = useState<Measurement | null>(null)

  const components = useMemo(() => componentDefinitions({
    mode,
    referenceColor,
    domainALength,
    linkerLength,
    sequenceLength: fullSequence.length,
  }), [domainALength, fullSequence.length, linkerLength, mode, referenceColor])

  useEffect(() => {
    const container = containerRef.current
    if (!container || !pdb) return
    let disposed = false
    let observer: ResizeObserver | null = null
    const load = async () => {
      try {
        if (mode === 'fusion') validateFusionPdb(pdb, fullSequence)
        else validateDomainReferencePdb(pdb, fullSequence, 'Domain')
        const mol = await import('3dmol')
        if (disposed) return
        container.replaceChildren()
        const viewer = mol.createViewer(container, { backgroundColor: canvasTheme === 'dark' ? '#08120f' : '#f4f7f5' })
        viewer.addModel(pdb, 'pdb')
        viewer.setClickable({ atom: 'CA' }, true, (atom: AtomSpec) => {
          if (atom.resi == null) return
          setSelectedResidue(atom.resi)
          setHoveredResidue(`${atom.resn ?? 'Residue'} ${atom.resi} · chain ${atom.chain || 'A'}`)
        })
        viewer.setHoverable({ atom: 'CA' }, true, (atom: AtomSpec) => {
          if (atom.resi == null) return
          setHoveredResidue(`${atom.resn ?? 'Residue'} ${atom.resi} · chain ${atom.chain || 'A'} · B ${atom.b?.toFixed(1) ?? 'n/a'}`)
        }, () => setHoveredResidue(selectedResidue ? `Residue ${selectedResidue} selected` : 'Hover or select a residue'))
        viewerRef.current = viewer
        viewer.zoomTo()
        applyViewerStyle(viewer, {
          representation,
          colourMode,
          mode,
          referenceColor,
          domainALength,
          linkerLength,
          sequenceLength: fullSequence.length,
          hiddenComponents,
          selectedResidue,
        })
        observer = new ResizeObserver(() => viewer.resize())
        observer.observe(container)
        setError(null)
        setViewerReady((value) => value + 1)
      } catch (cause) {
        setError(cause instanceof Error ? cause.message : 'Could not display this structure.')
      }
    }
    void load()
    return () => {
      disposed = true
      observer?.disconnect()
      viewerRef.current?.spin(false)
      viewerRef.current = null
      container.replaceChildren()
    }
    // Viewer creation is tied to the coordinate record. Display controls are applied below without rebuilding it.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [ariaLabel, domainALength, fullSequence, linkerLength, mode, pdb, referenceColor])

  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer) return
    viewer.setBackgroundColor(canvasTheme === 'dark' ? '#08120f' : '#f4f7f5', 1)
    applyViewerStyle(viewer, {
      representation,
      colourMode,
      mode,
      referenceColor,
      domainALength,
      linkerLength,
      sequenceLength: fullSequence.length,
      hiddenComponents,
      selectedResidue,
    })
  }, [canvasTheme, colourMode, domainALength, fullSequence.length, hiddenComponents, linkerLength, mode, referenceColor, representation, selectedResidue, viewerReady])

  useEffect(() => {
    viewerRef.current?.spin(spinning ? 'y' : false, 0.7)
    return () => viewerRef.current?.spin(false)
  }, [spinning, viewerReady])

  const focusSelection = (start: number, end: number) => {
    const safeStart = Math.min(Math.max(1, start), fullSequence.length)
    const safeEnd = Math.min(Math.max(safeStart, end), fullSequence.length)
    setSelectionStart(safeStart)
    setSelectionEnd(safeEnd)
    setSelectedResidue(safeStart === safeEnd ? safeStart : null)
    viewerRef.current?.zoomTo({ resi: `${safeStart}-${safeEnd}` }, 450)
    viewerRef.current?.render()
  }

  const toggleComponent = (key: ComponentKey) => {
    setHiddenComponents((previous) => {
      const next = new Set(previous)
      if (next.has(key)) next.delete(key)
      else next.add(key)
      return next
    })
  }

  const addMeasurement = () => {
    const viewer = viewerRef.current
    if (!viewer) return
    const first = viewer.selectedAtoms({ atom: 'CA', resi: measureStart })[0]
    const second = viewer.selectedAtoms({ atom: 'CA', resi: measureEnd })[0]
    const distance = first && second ? distanceBetween(first, second) : null
    const start = first && atomPoint(first)
    const end = second && atomPoint(second)
    if (distance == null || !start || !end) {
      setError('The selected residues do not both contain C-alpha coordinates.')
      return
    }
    viewer.removeAllShapes()
    viewer.removeAllLabels()
    viewer.addLine({ start, end, color: '#ffd94a', dashed: true })
    viewer.addLabel(`${distance.toFixed(2)} Å`, {
      position: { x: (start.x + end.x) / 2, y: (start.y + end.y) / 2, z: (start.z + end.z) / 2 },
      backgroundColor: '#111916',
      backgroundOpacity: 0.88,
      fontColor: '#ffd94a',
      fontSize: 14,
      borderThickness: 1,
      borderColor: '#ffd94a',
      inFront: true,
    })
    viewer.render()
    setMeasurement({ start: measureStart, end: measureEnd, distance })
    setError(null)
  }

  const clearMeasurement = () => {
    viewerRef.current?.removeAllShapes()
    viewerRef.current?.removeAllLabels()
    viewerRef.current?.render()
    setMeasurement(null)
  }

  if (!pdb) {
    return (
      <div className="structure-empty">
        <Box size={34} aria-hidden="true" />
        <strong>{mode === 'fusion' ? 'No full-fusion structure loaded' : 'No reference structure loaded'}</strong>
        <p>{mode === 'fusion'
          ? 'Run the complete fusion through local ColabFold in the desktop app, or load a matching single-chain PDB. The viewer never invents coordinates.'
          : 'Enter a four-character PDB accession below or load a matching local structure. LinkerForge verifies the domain sequence before displaying it.'}</p>
      </div>
    )
  }

  return (
    <div className="molecular-workspace">
      <header className="molecular-topbar">
        <div className="structure-identity">
          <span className="structure-id">{structureId || (mode === 'fusion' ? 'FULL FUSION' : 'REFERENCE')}</span>
          <div>
            <strong>{structureTitle || ariaLabel}</strong>
            <small>{structureSource || `${fullSequence.length} residues · sequence verified`}</small>
          </div>
        </div>
        <div className="molecular-top-actions">
          <button type="button" aria-pressed={spinning} onClick={() => setSpinning((value) => !value)}>
            {spinning ? <Pause size={15} /> : <Play size={15} />} {spinning ? 'Stop' : 'Spin'}
          </button>
          <button type="button" onClick={() => { viewerRef.current?.zoomTo(); viewerRef.current?.render() }}>
            <RotateCcw size={15} /> Reset
          </button>
          <button type="button" onClick={() => viewerRef.current && downloadViewerImage(viewerRef.current, mode === 'fusion' ? 'linkerforge-fusion.png' : 'linkerforge-reference.png')}>
            <Camera size={15} /> Image
          </button>
        </div>
      </header>

      <section className="sequence-browser" aria-label="Interactive protein sequence">
        <div className="sequence-browser-heading">
          <span><MousePointer2 size={14} /> Sequence</span>
          <small>{selectedResidue ? `Residue ${selectedResidue} selected` : `${fullSequence.length} amino acids`}</small>
        </div>
        <div className="sequence-track">
          {Array.from(fullSequence).map((residue, index) => {
            const residueNumber = index + 1
            const component = componentForResidue(residueNumber, mode, domainALength, linkerLength)
            return <button key={residueNumber} type="button"
              className={`sequence-residue residue-${component} ${selectedResidue === residueNumber ? 'selected' : ''}`}
              title={`Residue ${residueNumber}: ${residue}`}
              aria-label={`Select residue ${residueNumber} ${residue}`}
              onClick={() => focusSelection(residueNumber, residueNumber)}>
              {residueNumber === 1 || residueNumber % 10 === 0 ? <small>{residueNumber}</small> : null}
              <span>{residue}</span>
            </button>
          })}
        </div>
      </section>

      <div className="molecular-body">
        <nav className="canvas-tools" aria-label="Canvas tools">
          <button type="button" title="Reset camera" aria-label="Reset molecular camera" onClick={() => { viewerRef.current?.zoomTo(); viewerRef.current?.render() }}><RotateCcw size={18} /></button>
          <button type="button" title="Focus selected residue" aria-label="Focus selected residue" disabled={!selectedResidue} onClick={() => selectedResidue && focusSelection(selectedResidue, selectedResidue)}><Crosshair size={18} /></button>
          <button type="button" title="Toggle rotation" aria-label="Toggle molecular rotation" aria-pressed={spinning} onClick={() => setSpinning((value) => !value)}>{spinning ? <Pause size={18} /> : <Play size={18} />}</button>
          <button type="button" title="Save image" aria-label="Save molecular-view image" onClick={() => viewerRef.current && downloadViewerImage(viewerRef.current, 'linkerforge-structure.png')}><Camera size={18} /></button>
        </nav>

        <main className="molecular-canvas-shell">
          <div ref={containerRef} className={`fusion-3d-canvas canvas-${canvasTheme}`} aria-label={ariaLabel} />
          <div className="canvas-readout" aria-live="polite"><Info size={13} /> {hoveredResidue}</div>
          <div className="canvas-help">Drag: rotate · wheel: zoom · right-drag: pan · click: select residue</div>
          {error && <p role="alert" className="structure-error canvas-error">{error}</p>}
        </main>

        <aside className="structure-control-panel" aria-label="Structure controls">
          <details open>
            <summary><Search size={15} /> Selection <ChevronDown size={14} /></summary>
            <div className="control-panel-content">
              <div className="range-fields">
                <label>From<input type="number" min={1} max={fullSequence.length} value={selectionStart} onChange={(event) => setSelectionStart(Number(event.target.value))} /></label>
                <label>To<input type="number" min={1} max={fullSequence.length} value={selectionEnd} onChange={(event) => setSelectionEnd(Number(event.target.value))} /></label>
              </div>
              <button type="button" className="panel-action" onClick={() => focusSelection(selectionStart, selectionEnd)}><Focus size={14} /> Focus range</button>
              <p>{selectedResidue ? `Selected residue ${selectedResidue} (${fullSequence[selectedResidue - 1]}).` : 'Select from the sequence strip or the 3D canvas.'}</p>
            </div>
          </details>

          <details open>
            <summary><Layers3 size={15} /> Components <ChevronDown size={14} /></summary>
            <div className="control-panel-content component-list">
              {components.map((component) => <div className="component-row" key={component.key}>
                <button type="button" className="visibility-button" aria-label={`${hiddenComponents.has(component.key) ? 'Show' : 'Hide'} ${component.label}`} onClick={() => toggleComponent(component.key)}>
                  {hiddenComponents.has(component.key) ? <EyeOff size={15} /> : <Eye size={15} />}
                </button>
                <button type="button" className="component-focus" onClick={() => {
                  const [start, end = start] = component.range.split('–').map(Number)
                  focusSelection(start, end)
                }}>
                  <i style={{ background: component.color }} />
                  <span><strong>{component.label}</strong><small>{component.range} · {component.length} aa</small></span>
                </button>
              </div>)}
            </div>
          </details>

          <details open>
            <summary><Box size={15} /> Representation <ChevronDown size={14} /></summary>
            <div className="control-panel-content control-fields">
              <label>Style<select aria-label="Structure representation" value={representation} onChange={(event) => setRepresentation(event.target.value as Representation)}>
                <option value="cartoon">Cartoon</option>
                <option value="sticks">Sticks</option>
                <option value="ball-stick">Ball and stick</option>
                <option value="spheres">Space-filling spheres</option>
                <option value="surface">Molecular surface</option>
              </select></label>
              <label>Colour<select aria-label="Structure colour scheme" value={colourMode} onChange={(event) => setColourMode(event.target.value as ColourMode)}>
                <option value="domains">{mode === 'fusion' ? 'Domains and linker' : 'Reference component'}</option>
                <option value="confidence">B-factor / pLDDT</option>
                <option value="spectrum">Residue spectrum</option>
              </select></label>
              <label>Canvas<select aria-label="Viewer background" value={canvasTheme} onChange={(event) => setCanvasTheme(event.target.value as CanvasTheme)}>
                <option value="dark">Dark</option>
                <option value="light">Light</option>
              </select></label>
            </div>
          </details>

          <details open>
            <summary><Ruler size={15} /> Measurements <ChevronDown size={14} /></summary>
            <div className="control-panel-content">
              <div className="range-fields">
                <label>Residue 1<input aria-label="First measurement residue" type="number" min={1} max={fullSequence.length} value={measureStart} onChange={(event) => setMeasureStart(Number(event.target.value))} /></label>
                <label>Residue 2<input aria-label="Second measurement residue" type="number" min={1} max={fullSequence.length} value={measureEnd} onChange={(event) => setMeasureEnd(Number(event.target.value))} /></label>
              </div>
              <div className="measurement-actions">
                <button type="button" className="panel-action" onClick={addMeasurement}>Measure Cα</button>
                {measurement && <button type="button" onClick={clearMeasurement}>Clear</button>}
              </div>
              <p>{measurement ? `${measurement.start} ↔ ${measurement.end}: ${measurement.distance.toFixed(2)} Å` : 'Adds a Cα–Cα distance to the canvas.'}</p>
            </div>
          </details>

          <details open>
            <summary><Download size={15} /> Export & source <ChevronDown size={14} /></summary>
            <div className="control-panel-content export-panel-actions">
              <button type="button" onClick={() => viewerRef.current && downloadViewerImage(viewerRef.current, 'linkerforge-structure.png')}><Camera size={14} /> PNG image</button>
              <button type="button" onClick={() => downloadCoordinates(pdb, `${(structureId || 'linkerforge-structure').toLowerCase()}.pdb`)}><Download size={14} /> PDB coordinates</button>
              {sourceUrl && <a href={sourceUrl} target="_blank" rel="noreferrer">View source record</a>}
              <small>Coordinates are displayed only after sequence validation. RCSB accessions link back to their authoritative public record.</small>
            </div>
          </details>
        </aside>
      </div>
    </div>
  )
}
