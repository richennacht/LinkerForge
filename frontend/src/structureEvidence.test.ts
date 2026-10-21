import { describe, expect, it } from 'vitest'
import {
  auditFusionCoordinates,
  assessProtectedSiteClearance,
  calculateDomainPreservation,
  parseColabFoldScoreJson,
  parseProtectedSites,
  validateDomainReferencePdb,
} from './structureEvidence'

const pdb = [
  'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C',
  'ATOM      2  CA  CYS A   2       3.800   0.000   0.000  1.00 80.00           C',
  'ATOM      3  CA  GLY A   3       7.600   0.000   0.000  1.00 50.00           C',
  'ATOM      4  CA  ASP A   4      11.400   0.000   0.000  1.00 70.00           C',
  'ATOM      5  CA  GLU A   5      15.200   0.000   0.000  1.00 75.00           C',
  'END',
].join('\n')

describe('structural evidence calculations', () => {
  it('requires each domain reference to match its entered sequence', () => {
    expect(validateDomainReferencePdb(pdb, 'ACGDE', 'A')).toHaveLength(5)
    expect(() => validateDomainReferencePdb(pdb, 'ACGDF', 'A')).toThrow(/does not exactly match/i)
  })

  it('audits C-alpha continuity and overlaps', () => {
    const result = auditFusionCoordinates(pdb)
    expect(result.chainBreakCount).toBe(0)
    expect(result.severeCaOverlapCount).toBe(0)
    expect(result.maxAdjacentCaDistance).toBeCloseTo(3.8)
  })

  it('calculates rigid-body invariant domain RMSD', () => {
    const translatedDomainA = [
      'ATOM      1  CA  ALA A   1       5.000   5.000   5.000  1.00 90.00           C',
      'ATOM      2  CA  CYS A   2       8.800   5.000   5.000  1.00 80.00           C',
      'ATOM      3  CA  GLY A   3      12.600   5.000   5.000  1.00 50.00           C',
    ].join('\n')
    const result = calculateDomainPreservation(pdb, 3, 0, translatedDomainA, pdb.split('\n').slice(3, 5).join('\n'))
    expect(result.domainA).toBeCloseTo(0)
    expect(result.domainB).toBeNull()
  })

  it('parses protected sites with domain-relative numbering', () => {
    expect(parseProtectedSites('A:2; B:1-2', 3, 2, 4)).toEqual({ residueIndices: [2, 6, 7], errors: [] })
    expect(parseProtectedSites('A:9', 3, 2, 4).errors).toEqual(['A:9'])
  })

  it('measures protected-site C-alpha clearance from the linker', () => {
    const result = assessProtectedSiteClearance(pdb, [1, 5], 2, 1)
    expect(result.minimumDistanceAngstrom).toBeCloseTo(7.6)
    expect(result.sites.map((site) => site.withinContactCutoff)).toEqual([true, true])
  })

  it('loads pTM and cross-domain PAE from ColabFold score JSON', () => {
    const metrics = parseColabFoldScoreJson(JSON.stringify({
      plddt: [90, 80, 50, 70],
      pae: [[0, 1, 2, 10], [1, 0, 2, 20], [2, 2, 0, 3], [30, 40, 3, 0]],
      ptm: 0.6,
      max_pae: 40,
    }), 4, 2, 1)
    expect(metrics.meanPlddt).toBe(72.5)
    expect(metrics.meanInterdomainPae).toBe(25)
    expect(metrics.ptm).toBe(0.6)
  })
})
