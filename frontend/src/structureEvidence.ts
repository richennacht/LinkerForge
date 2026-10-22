const THREE_TO_ONE: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLU: 'E', GLN: 'Q', GLY: 'G',
  HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
  THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
}

export interface CaRecord {
  chain: string
  residueNumber: string
  aminoAcid: string
  x: number
  y: number
  z: number
}

export interface CoordinateAudit {
  residueCount: number
  chainBreakCount: number
  maxAdjacentCaDistance: number | null
  severeCaOverlapCount: number
}

export interface DomainPreservationResult {
  domainA: number | null
  domainB: number | null
}

export interface ProtectedSiteParseResult {
  residueIndices: number[]
  errors: string[]
}

export interface ProtectedSiteClearanceResult {
  method: 'C-alpha linker-distance pre-screen'
  contactCutoffAngstrom: number
  minimumDistanceAngstrom: number | null
  sites: Array<{
    fusionResidueIndex: number
    minimumLinkerCaDistanceAngstrom: number | null
    withinContactCutoff: boolean | null
  }>
}

export interface ColabFoldMetrics {
  predictor: string
  modelName: string | null
  meanPlddt: number | null
  ptm: number | null
  maxPae: number | null
  meanInterdomainPae: number | null
}

export interface RcsbReferenceStructure {
  pdbId: string
  chainId: string
  pdb: string
  sourceUrl: string
}

export function normalizePdbId(value: string): string {
  const pdbId = value.trim().toUpperCase()
  if (!/^[0-9][A-Z0-9]{3}$/.test(pdbId)) {
    throw new Error('Enter a four-character PDB accession such as 1UBQ.')
  }
  return pdbId
}

function matchingChainFragment(pdb: string, expectedSequence: string, label: string) {
  const chains = new Map<string, Array<{ key: string; aminoAcid: string }>>()
  const seen = new Set<string>()
  for (const line of pdb.split(/\r?\n/)) {
    if (!line.startsWith('ATOM  ') || line.slice(12, 16).trim() !== 'CA') continue
    const chain = line.slice(21, 22)
    const residueNumber = line.slice(22, 27).trim()
    const key = `${chain}:${residueNumber}`
    if (seen.has(key)) continue
    const aminoAcid = THREE_TO_ONE[line.slice(17, 20).trim()]
    if (!aminoAcid) continue
    seen.add(key)
    const records = chains.get(chain) ?? []
    records.push({ key, aminoAcid })
    chains.set(chain, records)
  }

  const matches: Array<{ chain: string; keys: Set<string> }> = []
  for (const [chain, records] of chains) {
    const observed = records.map((record) => record.aminoAcid).join('')
    let start = observed.indexOf(expectedSequence)
    while (start >= 0) {
      matches.push({
        chain,
        keys: new Set(records.slice(start, start + expectedSequence.length).map((record) => record.key)),
      })
      start = observed.indexOf(expectedSequence, start + 1)
    }
  }
  if (!matches.length) {
    throw new Error(`${label} sequence was not found as an exact chain or chain segment in this PDB entry.`)
  }

  const match = matches[0]
  const coordinateLines = pdb.split(/\r?\n/).filter((line) => {
    if (!line.startsWith('ATOM  ') && !line.startsWith('HETATM')) return false
    if (line.slice(21, 22) !== match.chain) return false
    return match.keys.has(`${match.chain}:${line.slice(22, 27).trim()}`)
  })
  const chainId = match.chain.trim() || '(blank)'
  return {
    chainId,
    pdb: [`REMARK 950 LINKERFORGE EXACT SEQUENCE-MATCHED CHAIN ${chainId}`, ...coordinateLines, 'TER', 'END', ''].join('\n'),
  }
}

export async function fetchRcsbReferenceStructure(
  value: string,
  expectedSequence: string,
  label: string,
  fetcher: typeof fetch = fetch
): Promise<RcsbReferenceStructure> {
  const pdbId = normalizePdbId(value)
  const sourceUrl = `https://files.rcsb.org/download/${pdbId}.pdb`
  const response = await fetcher(sourceUrl, { headers: { Accept: 'chemical/x-pdb,text/plain' } })
  if (!response.ok) {
    throw new Error(`RCSB could not provide PDB ${pdbId} (HTTP ${response.status}). Check the accession and network connection.`)
  }
  const bytes = await response.arrayBuffer()
  if (bytes.byteLength > 15_000_000) throw new Error('The downloaded PDB exceeds the 15 MB safety limit.')
  const downloaded = new TextDecoder().decode(bytes)
  const matched = matchingChainFragment(downloaded, expectedSequence, label)
  validateDomainReferencePdb(matched.pdb, expectedSequence, label)
  return { pdbId, chainId: matched.chainId, pdb: matched.pdb, sourceUrl }
}

export function parseColabFoldScoreJson(
  json: string,
  sequenceLength: number,
  domainALength: number,
  linkerLength: number
): ColabFoldMetrics {
  const value = JSON.parse(json) as { plddt?: unknown; pae?: unknown; ptm?: unknown; max_pae?: unknown }
  if (!Array.isArray(value.plddt) || value.plddt.length !== sequenceLength || !value.plddt.every((item) => typeof item === 'number')) {
    throw new Error(`ColabFold score JSON must contain ${sequenceLength} numeric pLDDT values.`)
  }
  if (!Array.isArray(value.pae) || value.pae.length !== sequenceLength || !value.pae.every((row) => Array.isArray(row) && row.length === sequenceLength)) {
    throw new Error(`ColabFold score JSON must contain a ${sequenceLength} × ${sequenceLength} PAE matrix.`)
  }
  const pae = value.pae as number[][]
  const domainBStart = domainALength + linkerLength
  const crossDomain: number[] = []
  for (let left = 0; left < domainALength; left += 1) {
    for (let right = domainBStart; right < sequenceLength; right += 1) crossDomain.push(pae[left][right])
  }
  for (let left = domainBStart; left < sequenceLength; left += 1) {
    for (let right = 0; right < domainALength; right += 1) crossDomain.push(pae[left][right])
  }
  return {
    predictor: 'AlphaFold2 model executed by ColabFold',
    modelName: null,
    meanPlddt: value.plddt.reduce((sum, item) => sum + item, 0) / value.plddt.length,
    ptm: typeof value.ptm === 'number' ? value.ptm : null,
    maxPae: typeof value.max_pae === 'number' ? value.max_pae : null,
    meanInterdomainPae: crossDomain.length ? crossDomain.reduce((sum, item) => sum + item, 0) / crossDomain.length : null,
  }
}

export function readCaRecords(pdb: string): CaRecord[] {
  const records: CaRecord[] = []
  const seen = new Set<string>()
  for (const line of pdb.split(/\r?\n/)) {
    if (!line.startsWith('ATOM  ') || line.slice(12, 16).trim() !== 'CA') continue
    const chain = line.slice(21, 22)
    const residueNumber = line.slice(22, 27).trim()
    const key = `${chain}:${residueNumber}`
    if (seen.has(key)) continue
    const aminoAcid = THREE_TO_ONE[line.slice(17, 20).trim()]
    const x = Number(line.slice(30, 38).trim())
    const y = Number(line.slice(38, 46).trim())
    const z = Number(line.slice(46, 54).trim())
    if (!aminoAcid || ![x, y, z].every(Number.isFinite)) {
      throw new Error('The PDB contains an unsupported residue or malformed C-alpha coordinate.')
    }
    seen.add(key)
    records.push({ chain, residueNumber, aminoAcid, x, y, z })
  }
  if (!records.length) throw new Error('The PDB contains no readable C-alpha trace.')
  if (new Set(records.map((record) => record.chain)).size !== 1) {
    throw new Error('Provide one protein chain per reference PDB.')
  }
  return records
}

export function validateDomainReferencePdb(pdb: string, expectedSequence: string, label: string): CaRecord[] {
  const records = readCaRecords(pdb)
  const observed = records.map((record) => record.aminoAcid).join('')
  if (observed !== expectedSequence) {
    throw new Error(`${label} reference sequence does not exactly match the entered domain (${observed.length} of ${expectedSequence.length} residues).`)
  }
  return records
}

function distance(a: CaRecord, b: CaRecord) {
  return Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z)
}

export function auditFusionCoordinates(pdb: string): CoordinateAudit {
  const records = readCaRecords(pdb)
  const adjacent = records.slice(1).map((record, index) => distance(records[index], record))
  let severeCaOverlapCount = 0
  for (let left = 0; left < records.length; left += 1) {
    for (let right = left + 3; right < records.length; right += 1) {
      if (distance(records[left], records[right]) < 3.0) severeCaOverlapCount += 1
    }
  }
  return {
    residueCount: records.length,
    chainBreakCount: adjacent.filter((value) => value > 4.5).length,
    maxAdjacentCaDistance: adjacent.length ? Math.max(...adjacent) : null,
    severeCaOverlapCount,
  }
}

type Vector = [number, number, number]

function alignedRmsd(mobile: CaRecord[], target: CaRecord[]): number | null {
  if (mobile.length !== target.length || mobile.length < 3) return null
  const centroid = (records: CaRecord[]): Vector => [
    records.reduce((sum, item) => sum + item.x, 0) / records.length,
    records.reduce((sum, item) => sum + item.y, 0) / records.length,
    records.reduce((sum, item) => sum + item.z, 0) / records.length,
  ]
  const cm = centroid(mobile)
  const ct = centroid(target)
  const m = mobile.map((item): Vector => [item.x - cm[0], item.y - cm[1], item.z - cm[2]])
  const t = target.map((item): Vector => [item.x - ct[0], item.y - ct[1], item.z - ct[2]])
  const s = Array.from({ length: 3 }, () => [0, 0, 0])
  for (let index = 0; index < m.length; index += 1) {
    for (let row = 0; row < 3; row += 1) {
      for (let column = 0; column < 3; column += 1) s[row][column] += m[index][row] * t[index][column]
    }
  }
  const [sxx, sxy, sxz] = s[0]
  const [syx, syy, syz] = s[1]
  const [szx, szy, szz] = s[2]
  const n = [
    [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
    [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
    [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
    [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
  ]
  let q = [1, 0, 0, 0]
  for (let iteration = 0; iteration < 80; iteration += 1) {
    const next = n.map((row) => row.reduce((sum, value, index) => sum + value * q[index], 0))
    const norm = Math.hypot(...next)
    if (!norm) break
    q = next.map((value) => value / norm)
  }
  const [w, x, y, z] = q
  const rotation = [
    [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
    [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
    [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
  ]
  const squared = m.reduce((sum, vector, index) => {
    const rotated = rotation.map((row) => row.reduce((value, coefficient, axis) => value + coefficient * vector[axis], 0))
    return sum + rotated.reduce((value, coordinate, axis) => value + (coordinate - t[index][axis]) ** 2, 0)
  }, 0)
  return Math.sqrt(squared / mobile.length)
}

export function calculateDomainPreservation(
  fusionPdb: string,
  domainALength: number,
  linkerLength: number,
  domainAReference: string | null,
  domainBReference: string | null
): DomainPreservationResult {
  const fusion = readCaRecords(fusionPdb)
  const domainBStart = domainALength + linkerLength
  return {
    domainA: domainAReference ? alignedRmsd(fusion.slice(0, domainALength), readCaRecords(domainAReference)) : null,
    domainB: domainBReference ? alignedRmsd(fusion.slice(domainBStart), readCaRecords(domainBReference)) : null,
  }
}

export function parseProtectedSites(input: string, domainALength: number, linkerLength: number, domainBLength: number): ProtectedSiteParseResult {
  if (!input.trim()) return { residueIndices: [], errors: [] }
  const residueIndices = new Set<number>()
  const errors: string[] = []
  for (const token of input.split(/[;,]+/).map((value) => value.trim()).filter(Boolean)) {
    const match = token.match(/^(?:domain\s*)?([AB])\s*:\s*(\d+)(?:\s*[-–]\s*(\d+))?$/i)
    if (!match) {
      errors.push(token)
      continue
    }
    const domain = match[1].toUpperCase()
    const start = Number(match[2])
    const end = Number(match[3] ?? match[2])
    const limit = domain === 'A' ? domainALength : domainBLength
    if (start < 1 || end < start || end > limit) {
      errors.push(token)
      continue
    }
    const offset = domain === 'A' ? 0 : domainALength + linkerLength
    for (let value = start; value <= end; value += 1) residueIndices.add(offset + value)
  }
  return { residueIndices: [...residueIndices].sort((a, b) => a - b), errors }
}

export function assessProtectedSiteClearance(
  pdb: string,
  residueIndices: number[],
  domainALength: number,
  linkerLength: number,
  contactCutoffAngstrom = 8.0
): ProtectedSiteClearanceResult {
  const records = readCaRecords(pdb)
  const linker = records.slice(domainALength, domainALength + linkerLength)
  const sites = residueIndices.map((fusionResidueIndex) => {
    const site = records[fusionResidueIndex - 1]
    if (!site) throw new Error(`Protected-site residue ${fusionResidueIndex} is outside the fusion structure.`)
    const minimum = linker.length ? Math.min(...linker.map((record) => distance(site, record))) : null
    return {
      fusionResidueIndex,
      minimumLinkerCaDistanceAngstrom: minimum,
      withinContactCutoff: minimum == null ? null : minimum <= contactCutoffAngstrom,
    }
  })
  const measured = sites
    .map((site) => site.minimumLinkerCaDistanceAngstrom)
    .filter((value): value is number => value != null)
  return {
    method: 'C-alpha linker-distance pre-screen',
    contactCutoffAngstrom,
    minimumDistanceAngstrom: measured.length ? Math.min(...measured) : null,
    sites,
  }
}
