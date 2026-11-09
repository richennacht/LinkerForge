import { fireEvent, render, screen, waitFor } from '@testing-library/react'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import FusionStructureViewer from './FusionStructureViewer'

const viewer = vi.hoisted(() => ({
  addModel: vi.fn(),
  setStyle: vi.fn(),
  addStyle: vi.fn(),
  setClickable: vi.fn(),
  setHoverable: vi.fn(),
  selectedAtoms: vi.fn(({ resi }: { resi?: number }) => [{ resi, atom: 'CA', x: Number(resi) * 3.8, y: 0, z: 0 }]),
  removeAllSurfaces: vi.fn(),
  addSurface: vi.fn(),
  removeAllShapes: vi.fn(),
  removeAllLabels: vi.fn(),
  addLine: vi.fn(),
  addLabel: vi.fn(),
  setBackgroundColor: vi.fn(),
  zoomTo: vi.fn(),
  render: vi.fn(),
  resize: vi.fn(),
  spin: vi.fn(),
  pngURI: vi.fn(() => 'data:image/png;base64,AAAA'),
}))

vi.mock('3dmol', () => ({ createViewer: () => viewer }))

const pdb = [
  'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C',
  'ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 90.00           C',
  'ATOM      3  CA  CYS A   3       7.600   0.000   0.000  1.00 90.00           C',
].join('\n')

describe('full-fusion structure viewer', () => {
  beforeEach(() => {
    vi.clearAllMocks()
    vi.stubGlobal('ResizeObserver', class {
      observe() {}
      disconnect() {}
    })
  })

  it('does not create a molecular model without coordinates', () => {
    render(<FusionStructureViewer pdb={null} fullSequence="AGC" domainALength={1} linkerLength={1} />)
    expect(screen.getByText(/no full-fusion structure loaded/i)).toBeInTheDocument()
    expect(viewer.addModel).not.toHaveBeenCalled()
  })

  it('colours the verified linker separately from both domains', async () => {
    render(<FusionStructureViewer pdb={pdb} fullSequence="AGC" domainALength={1} linkerLength={1} />)
    await waitFor(() => expect(viewer.render).toHaveBeenCalled())
    expect(viewer.addModel).toHaveBeenCalledWith(pdb, 'pdb')
    expect(viewer.setStyle).toHaveBeenCalledWith(
      { resi: '2-2' },
      { cartoon: { color: '#ff4dba' } }
    )
  })

  it('changes the embedded representation and controls spin without opening another app', async () => {
    render(<FusionStructureViewer pdb={pdb} fullSequence="AGC" domainALength={1} linkerLength={1} />)
    await waitFor(() => expect(viewer.render).toHaveBeenCalled())

    fireEvent.change(screen.getByLabelText(/structure representation/i), { target: { value: 'ball-stick' } })
    await waitFor(() => expect(viewer.setStyle).toHaveBeenCalledWith(
      { resi: '2-2' },
      {
        stick: { color: '#ff4dba', radius: 0.18 },
        sphere: { color: '#ff4dba', scale: 0.28 },
      }
    ))

    fireEvent.click(screen.getByRole('button', { name: /^spin$/i }))
    expect(viewer.spin).toHaveBeenCalledWith('y', 0.7)
  })

  it('provides sequence selection, components, and a C-alpha measurement panel', async () => {
    render(<FusionStructureViewer pdb={pdb} fullSequence="AGC" domainALength={1} linkerLength={1} />)
    await waitFor(() => expect(viewer.render).toHaveBeenCalled())

    fireEvent.click(screen.getByRole('button', { name: /select residue 2 g/i }))
    expect(viewer.zoomTo).toHaveBeenCalledWith({ resi: '2-2' }, 450)
    expect(screen.getByText(/selected residue 2/i)).toBeInTheDocument()

    fireEvent.change(screen.getByLabelText(/first measurement residue/i), { target: { value: '1' } })
    fireEvent.change(screen.getByLabelText(/second measurement residue/i), { target: { value: '3' } })
    fireEvent.click(screen.getByRole('button', { name: /measure cα/i }))
    expect(viewer.addLine).toHaveBeenCalled()
    expect(screen.getByText(/1 ↔ 3: 7.60 Å/i)).toBeInTheDocument()
  })
})
