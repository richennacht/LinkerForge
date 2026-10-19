import { invoke, isTauri } from '@tauri-apps/api/core'

export interface ComputeStatus {
  predictorAvailable: boolean
  predictorBackend: string
  gpu: string
  note: string
}

export interface FoldStatus {
  phase: 'idle' | 'running' | 'complete' | 'failed'
  message: string
  runDir: string | null
  sequenceSha256: string | null
  pdb: string | null
}

export const inDesktop = () => isTauri()

export const getComputeStatus = () => invoke<ComputeStatus>('compute_status')

export const startLocalFold = (sequence: string) =>
  invoke<FoldStatus>('start_fold', { sequence })

export const getFoldStatus = () => invoke<FoldStatus>('fold_status')
