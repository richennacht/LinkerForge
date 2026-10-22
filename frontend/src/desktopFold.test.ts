import { beforeEach, describe, expect, it, vi } from 'vitest'

const tauri = vi.hoisted(() => ({
  invoke: vi.fn(),
  isTauri: vi.fn(() => true),
}))

vi.mock('@tauri-apps/api/core', () => tauri)

import { errorMessage, foldStartTimeFromRunDir, formatElapsedTime, inDesktop, openPdbInPymol, openReferenceInPymol } from './desktopFold'

describe('desktop folding bridge', () => {
  it('formats elapsed prediction time and recovers the run start timestamp', () => {
    expect(formatElapsedTime(3_723_900)).toBe('01:02:03')
    expect(foldStartTimeFromRunDir('M:\\runs\\folds\\1792481307-abc123')).toBe(1_792_481_307_000)
    expect(foldStartTimeFromRunDir(null)).toBeNull()
  })

  it('preserves string errors returned by Tauri commands', () => {
    expect(errorMessage('A local fold is already running.', 'fallback')).toBe('A local fold is already running.')
    expect(errorMessage(null, 'fallback')).toBe('fallback')
  })

  beforeEach(() => {
    tauri.invoke.mockReset()
  })

  it('forwards exact sequence and segment boundaries to the native PyMOL gate', async () => {
    tauri.invoke.mockResolvedValue('opened')
    const request = {
      pdb: 'ATOM',
      sequence: 'AAAGGBBB',
      domainALength: 3,
      linkerLength: 2,
    }

    await expect(openPdbInPymol(request)).resolves.toBe('opened')
    expect(tauri.invoke).toHaveBeenCalledWith('open_pdb_in_pymol', { request })
    expect(inDesktop()).toBe(true)
  })

  it('forwards a sequence-checked domain reference to the native PyMOL gate', async () => {
    tauri.invoke.mockResolvedValue('opened')
    const request = { pdb: 'ATOM', sequence: 'ACD', label: 'Ubiquitin', color: '#3acf78' }

    await expect(openReferenceInPymol(request)).resolves.toBe('opened')
    expect(tauri.invoke).toHaveBeenCalledWith('open_reference_in_pymol', { request })
  })
})
