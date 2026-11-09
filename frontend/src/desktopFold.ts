import { invoke, isTauri } from '@tauri-apps/api/core'

export interface ComputeStatus {
  predictorAvailable: boolean
  predictorBackend: string
  gpu: string
  pymolAvailable: boolean
  pymolPath: string | null
  note: string
}

export interface FoldStatus {
  phase: 'idle' | 'running' | 'complete' | 'failed'
  message: string
  runDir: string | null
  sequenceSha256: string | null
  pdb: string | null
  metrics: FoldMetrics | null
  artifacts: FoldArtifacts | null
}

export interface FoldMetrics {
  predictor: string
  modelName: string | null
  meanPlddt: number | null
  ptm: number | null
  maxPae: number | null
  meanInterdomainPae: number | null
}

export interface FoldArtifacts {
  pdbPath: string | null
  scoreJsonPath: string | null
  paePlotPath: string | null
  plddtPlotPath: string | null
  logPath: string | null
}

export interface PymolLaunchRequest {
  pdb: string
  sequence: string
  domainALength: number
  linkerLength: number
}

export interface ReferencePymolLaunchRequest {
  pdb: string
  sequence: string
  label: string
  color: string
}

export const foldStartTimeFromRunDir = (runDir: string | null) => {
  const match = runDir?.match(/(?:^|[\\/])(\d{10,13})-[^\\/]+$/)
  if (!match) return null
  const timestamp = Number(match[1])
  return Number.isFinite(timestamp) ? timestamp * (timestamp < 1_000_000_000_000 ? 1000 : 1) : null
}

export const formatElapsedTime = (elapsedMs: number) => {
  const totalSeconds = Math.max(0, Math.floor(elapsedMs / 1000))
  const hours = Math.floor(totalSeconds / 3600)
  const minutes = Math.floor((totalSeconds % 3600) / 60)
  const seconds = totalSeconds % 60
  return [hours, minutes, seconds].map((value) => value.toString().padStart(2, '0')).join(':')
}

export const errorMessage = (cause: unknown, fallback: string) => {
  if (cause instanceof Error) return cause.message
  if (typeof cause === 'string' && cause.trim()) return cause
  return fallback
}

export const inDesktop = () => isTauri()

export const getComputeStatus = () => invoke<ComputeStatus>('compute_status')

export const startLocalFold = (sequence: string, domainALength: number, linkerLength: number) =>
  invoke<FoldStatus>('start_fold', { sequence, domainALength, linkerLength })

export const getFoldStatus = () => invoke<FoldStatus>('fold_status')

export const openPdbInPymol = (request: PymolLaunchRequest) =>
  invoke<string>('open_pdb_in_pymol', { request })

export const openReferenceInPymol = (request: ReferencePymolLaunchRequest) =>
  invoke<string>('open_reference_in_pymol', { request })
