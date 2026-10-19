import { describe, expect, it } from 'vitest'
import { buildTemplateCandidates, parseProteinSequence, validateFusionPdb } from './fusion'

describe('custom fusion workflow', () => {
  it('accepts FASTA and rejects ambiguous residues', () => {
    expect(parseProteinSequence('>domain A\nACD\nEFG', 'A')).toBe('ACDEFG')
    expect(() => parseProteinSequence('AXD', 'A')).toThrow(/20 standard/)
  })

  it('assembles both supplied domains around each actual linker sequence', () => {
    const candidates = buildTemplateCandidates('ACD', 'EFG', 15, 'flexible', false)
    expect(candidates[0].linker).toBe('GGGGSGGGGSGGGGS')
    expect(candidates[0].fullSequence).toBe('ACDGGGGSGGGGSGGGGSEFG')
  })

  it('requires the full predicted chain to match the chosen fusion', () => {
    const pdb = [
      'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C',
      'ATOM      2  CA  CYS A   2       3.800   0.000   0.000  1.00 90.00           C',
      'END',
    ].join('\n')
    expect(() => validateFusionPdb(pdb, 'AC')).not.toThrow()
    expect(() => validateFusionPdb(pdb, 'AD')).toThrow(/does not match/)
    expect(() => validateFusionPdb(pdb.replace('A   2', 'A   9'), 'AC')).toThrow(/consecutively/)
  })
})
