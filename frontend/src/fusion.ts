export const MAX_DOMAIN_LENGTH = 1000

const CANONICAL_AA = /^[ACDEFGHIKLMNPQRSTVWY]+$/

export function parseProteinSequence(input: string, name: string): string {
  const sequence = input
    .split(/\r?\n/)
    .filter((line) => !line.trim().startsWith('>'))
    .join('')
    .replace(/\s+/g, '')
    .toUpperCase()
  if (!sequence) throw new Error(`${name} needs an amino-acid sequence.`)
  if (!CANONICAL_AA.test(sequence)) {
    throw new Error(`${name} must contain only the 20 standard amino-acid letters.`)
  }
  if (sequence.length > MAX_DOMAIN_LENGTH) {
    throw new Error(`${name} exceeds the ${MAX_DOMAIN_LENGTH}-residue local input limit.`)
  }
  return sequence
}

export interface FusionCandidate {
  id: string
  linker: string
  fullSequence: string
  method: string
}

export function buildTemplateCandidates(
  domainA: string,
  domainB: string,
  targetLength: number,
  rigidity: 'flexible' | 'rigid',
  cleavable: boolean
): FusionCandidate[] {
  const requested = Math.max(5, Math.min(30, Math.round(targetLength)))
  const lengths = [...new Set([requested, Math.max(5, requested - 5), Math.min(30, requested + 5)])]
  const linkers = lengths.map((length) => {
    if (cleavable) {
      const motif = 'ENLYFQS'
      const flank = Math.max(0, length - motif.length)
      const left = 'G'.repeat(Math.floor(flank / 2))
      const right = 'S'.repeat(Math.ceil(flank / 2))
      return `${left}${motif}${right}`
    }
    const unit = rigidity === 'rigid' ? 'EAAAK' : 'GGGGS'
    return unit.repeat(Math.ceil(length / unit.length)).slice(0, length)
  })
  return [...new Set(linkers)].map((linker, index) => ({
    id: `template-${index + 1}`,
    linker,
    fullSequence: `${domainA}${linker}${domainB}`,
    method: cleavable ? 'TEV motif template' : rigidity === 'rigid' ? 'EAAAK template' : 'GGGGS template',
  }))
}

export function fusionFasta(sequence: string): string {
  return `>LinkerForge_candidate_full_fusion\n${sequence.match(/.{1,80}/g)?.join('\n') ?? ''}\n`
}

const THREE_TO_ONE: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLU: 'E', GLN: 'Q', GLY: 'G',
  HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
  THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
}

export function validateFusionPdb(pdb: string, expectedSequence: string): void {
  const residues = new Map<string, string>()
  const chains = new Set<string>()
  for (const line of pdb.split(/\r?\n/)) {
    if (!line.startsWith('ATOM  ') || line.slice(12, 16).trim() !== 'CA') continue
    const chain = line.slice(21, 22)
    const residue = line.slice(17, 20).trim()
    const key = `${chain}:${line.slice(22, 27)}`
    const aa = THREE_TO_ONE[residue]
    if (!aa) throw new Error(`Unrecognized PDB residue ${residue || '(blank)'}.`)
    chains.add(chain)
    residues.set(key, aa)
  }
  if (chains.size !== 1 || residues.size === 0) {
    throw new Error('Provide a single-chain full-fusion PDB with C-alpha atoms.')
  }
  const observed = [...residues.values()].join('')
  const numbering = [...residues.keys()].map((key) => Number(key.split(':')[1].trim()))
  if (numbering.some((number, index) => number !== index + 1)) {
    throw new Error('PDB residues must be numbered consecutively from 1 for linker colouring.')
  }
  if (observed !== expectedSequence) {
    throw new Error(`PDB sequence does not match the selected full fusion (${observed.length} of ${expectedSequence.length} residues).`)
  }
}
