import { render, screen, waitFor } from '@testing-library/react'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import FusionStructureViewer from './FusionStructureViewer'

const viewer = vi.hoisted(() => ({
  addModel: vi.fn(),
  setStyle: vi.fn(),
  zoomTo: vi.fn(),
  render: vi.fn(),
  resize: vi.fn(),
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
      { cartoon: { color: '#ff4dba' }, stick: { color: '#ff4dba', radius: 0.22 } }
    )
  })
})
