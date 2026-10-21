const encoder = new TextEncoder()

export type ResearchBundleFile = string | Uint8Array

export async function sha256Hex(value: string | Uint8Array): Promise<string> {
  const bytes = typeof value === 'string' ? encoder.encode(value) : value
  const copy = new Uint8Array(bytes.byteLength)
  copy.set(bytes)
  const digest = await crypto.subtle.digest('SHA-256', copy.buffer)
  return Array.from(new Uint8Array(digest), (byte) => byte.toString(16).padStart(2, '0')).join('')
}

function crc32(bytes: Uint8Array): number {
  let crc = 0xffffffff
  for (const byte of bytes) {
    crc ^= byte
    for (let bit = 0; bit < 8; bit += 1) crc = (crc >>> 1) ^ (0xedb88320 & -(crc & 1))
  }
  return (crc ^ 0xffffffff) >>> 0
}

function write16(view: DataView, offset: number, value: number) {
  view.setUint16(offset, value, true)
}

function write32(view: DataView, offset: number, value: number) {
  view.setUint32(offset, value, true)
}

function joinBytes(parts: Uint8Array[]): Uint8Array {
  const output = new Uint8Array(parts.reduce((sum, part) => sum + part.length, 0))
  let offset = 0
  for (const part of parts) {
    output.set(part, offset)
    offset += part.length
  }
  return output
}

export async function checksumManifest(files: Record<string, ResearchBundleFile>): Promise<string> {
  const lines = await Promise.all(Object.keys(files).sort().map(async (name) => {
    const content = files[name]
    return `${await sha256Hex(content)}  ${name}`
  }))
  return `${lines.join('\n')}\n`
}

/** Build a deterministic ZIP archive using the portable uncompressed method. */
export function createStoredZip(files: Record<string, ResearchBundleFile>): Uint8Array {
  const localParts: Uint8Array[] = []
  const centralParts: Uint8Array[] = []
  let localOffset = 0
  for (const name of Object.keys(files).sort()) {
    if (!name || name.includes('..') || name.startsWith('/') || name.includes('\\')) {
      throw new Error(`Unsafe research-bundle path: ${name}`)
    }
    const nameBytes = encoder.encode(name)
    const content = typeof files[name] === 'string' ? encoder.encode(files[name] as string) : files[name] as Uint8Array
    const checksum = crc32(content)
    const local = new Uint8Array(30 + nameBytes.length)
    const localView = new DataView(local.buffer)
    write32(localView, 0, 0x04034b50)
    write16(localView, 4, 20)
    write16(localView, 6, 0x0800)
    write16(localView, 8, 0)
    write32(localView, 14, checksum)
    write32(localView, 18, content.length)
    write32(localView, 22, content.length)
    write16(localView, 26, nameBytes.length)
    local.set(nameBytes, 30)
    localParts.push(local, content)

    const central = new Uint8Array(46 + nameBytes.length)
    const centralView = new DataView(central.buffer)
    write32(centralView, 0, 0x02014b50)
    write16(centralView, 4, 20)
    write16(centralView, 6, 20)
    write16(centralView, 8, 0x0800)
    write16(centralView, 10, 0)
    write32(centralView, 16, checksum)
    write32(centralView, 20, content.length)
    write32(centralView, 24, content.length)
    write16(centralView, 28, nameBytes.length)
    write32(centralView, 42, localOffset)
    central.set(nameBytes, 46)
    centralParts.push(central)
    localOffset += local.length + content.length
  }
  const central = joinBytes(centralParts)
  const end = new Uint8Array(22)
  const endView = new DataView(end.buffer)
  write32(endView, 0, 0x06054b50)
  write16(endView, 8, centralParts.length)
  write16(endView, 10, centralParts.length)
  write32(endView, 12, central.length)
  write32(endView, 16, localOffset)
  return joinBytes([...localParts, central, end])
}

export interface HistoryEntry {
  schema_version: 'linkerforge.history-entry.v1'
  entry_id: string
  created_at: string
  previous_entry_sha256: string | null
  record_sha256: string
  entry_sha256: string
  record: unknown
}

const historyPayload = (entry: Omit<HistoryEntry, 'entry_sha256'>) => JSON.stringify(entry)

export async function appendHistoryEntry(history: HistoryEntry[], record: unknown, createdAt = new Date().toISOString()): Promise<HistoryEntry[]> {
  const recordText = JSON.stringify(record)
  const partial: Omit<HistoryEntry, 'entry_sha256'> = {
    schema_version: 'linkerforge.history-entry.v1',
    entry_id: `${createdAt}-${(await sha256Hex(recordText)).slice(0, 12)}`,
    created_at: createdAt,
    previous_entry_sha256: history.at(-1)?.entry_sha256 ?? null,
    record_sha256: await sha256Hex(recordText),
    record,
  }
  return [...history, { ...partial, entry_sha256: await sha256Hex(historyPayload(partial)) }]
}

export async function verifyHistory(history: HistoryEntry[]): Promise<boolean> {
  for (let index = 0; index < history.length; index += 1) {
    const entry = history[index]
    const { entry_sha256: entrySha256, ...partial } = entry
    if (entry.previous_entry_sha256 !== (history[index - 1]?.entry_sha256 ?? null)) return false
    if (entry.record_sha256 !== await sha256Hex(JSON.stringify(entry.record))) return false
    if (entrySha256 !== await sha256Hex(historyPayload(partial))) return false
  }
  return true
}
