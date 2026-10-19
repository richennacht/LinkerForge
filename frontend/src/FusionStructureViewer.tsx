import { useEffect, useRef, useState } from 'react'
import { validateFusionPdb } from './fusion'

interface Props {
  pdb: string | null
  fullSequence: string
  domainALength: number
  linkerLength: number
}

export default function FusionStructureViewer({ pdb, fullSequence, domainALength, linkerLength }: Props) {
  const containerRef = useRef<HTMLDivElement | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const container = containerRef.current
    if (!container || !pdb) return
    let disposed = false
    let observer: ResizeObserver | null = null
    const load = async () => {
      try {
        validateFusionPdb(pdb, fullSequence)
        const mol = await import('3dmol')
        if (disposed) return
        container.replaceChildren()
        const viewer = mol.createViewer(container, { backgroundColor: '#0b1712' })
        viewer.addModel(pdb, 'pdb')
        viewer.setStyle({}, { cartoon: { color: '#9aa49b' } })
        viewer.setStyle({ resi: `1-${domainALength}` }, { cartoon: { color: '#3acf78' } })
        const linkerStart = domainALength + 1
        const linkerEnd = domainALength + linkerLength
        if (linkerLength > 0) {
          viewer.setStyle(
            { resi: `${linkerStart}-${linkerEnd}` },
            { cartoon: { color: '#ff4dba' }, stick: { color: '#ff4dba', radius: 0.22 } }
          )
        }
        viewer.setStyle(
          { resi: `${linkerEnd + 1}-${fullSequence.length}` },
          { cartoon: { color: '#f66b62' } }
        )
        viewer.zoomTo()
        viewer.render()
        observer = new ResizeObserver(() => viewer.resize())
        observer.observe(container)
        setError(null)
      } catch (cause) {
        setError(cause instanceof Error ? cause.message : 'Could not display this structure.')
      }
    }
    void load()
    return () => {
      disposed = true
      observer?.disconnect()
      container.replaceChildren()
    }
  }, [pdb, fullSequence, domainALength, linkerLength])

  if (!pdb) {
    return (
      <div className="structure-empty">
        <strong>No full-fusion structure loaded</strong>
        <p>Export the selected full-fusion FASTA, predict that whole sequence with AlphaFold/ColabFold, then load its single-chain PDB here. No coordinates are invented by this page.</p>
      </div>
    )
  }
  return <>{error && <p role="alert" className="structure-error">{error}</p>}<div ref={containerRef} className="fusion-3d-canvas" aria-label="Predicted full fusion structure" /></>
}
