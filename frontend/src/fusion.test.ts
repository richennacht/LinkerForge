import { describe, expect, it } from 'vitest'
import { buildTemplateCandidates, parseProteinSequence, summarizePdbConfidence, validateFusionPdb } from './fusion'

describe('custom fusion workflow', () => {
  it('accepts FASTA and rejects ambiguous residues', () => {
    expect(parseProteinSequence('>domain A\nACD\nEFG', 'A')).toBe('ACDEFG')
    expect(() => parseProteinSequence('AXD', 'A')).toThrow(/20 standard/)
  })

  it('assembles both supplied domains around each actual linker sequence', () => {
    const candidates = buildTemplateCandidates('ACD', 'EFG', 15, 'flexible', false)
    expect(candidates[0].linker).toBe('GGGGSGGGGSGGGGS')
    expect(candidates[0].fullSequence).toBe('ACDGGGGSGGGGSGGGGSEFG')
    expect(candidates.some((candidate) => candidate.role === 'fixed_control')).toBe(true)
    expect(candidates.some((candidate) => candidate.role === 'direct_control' && candidate.linker === '')).toBe(true)
  })

  it('labels the requested length separately from shorter sensitivity candidates', () => {
    const candidates = buildTemplateCandidates('ACD', 'EFG', 8, 'flexible', false)
    expect(candidates[0].linker).toHaveLength(8)
    expect(candidates[0].method).toMatch(/target-length/i)
    expect(candidates.some((candidate) => candidate.linker.length === 5 && candidate.method.match(/shorter/i))).toBe(true)
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

  it('summarizes AlphaFold confidence stored in PDB B-factors by region', () => {
    const pdb = [
      'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C',
      'ATOM      2  CA  CYS A   2       3.800   0.000   0.000  1.00 80.00           C',
      'ATOM      3  CA  GLY A   3       7.600   0.000   0.000  1.00 50.00           C',
      'ATOM      4  CA  ASP A   4      11.400   0.000   0.000  1.00 70.00           C',
      'END',
    ].join('\n')
    const summary = summarizePdbConfidence(pdb, 2, 1)
    expect(summary.overall).toBe(72.5)
    expect(summary.domainA).toBe(85)
    expect(summary.linker).toBe(50)
    expect(summary.domainB).toBe(70)
  })
})
