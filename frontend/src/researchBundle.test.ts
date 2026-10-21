import { describe, expect, it } from 'vitest'
import { appendHistoryEntry, checksumManifest, createStoredZip, verifyHistory } from './researchBundle'

describe('research bundle and history', () => {
  it('creates a deterministic ZIP containing portable paths', () => {
    const first = createStoredZip({ 'fusion.fasta': '>x\nACD\n', 'result.json': '{}\n' })
    const second = createStoredZip({ 'result.json': '{}\n', 'fusion.fasta': '>x\nACD\n' })
    expect(Array.from(first)).toEqual(Array.from(second))
    expect(new DataView(first.buffer).getUint32(0, true)).toBe(0x04034b50)
    expect(() => createStoredZip({ '../escape': 'x' })).toThrow(/unsafe/i)
  })

  it('writes sorted SHA-256 checksums', async () => {
    const manifest = await checksumManifest({ 'b.txt': 'b', 'a.txt': 'a' })
    expect(manifest.split('\n')[0]).toMatch(/^[0-9a-f]{64} {2}a\.txt$/)
  })

  it('detects changes to an append-only hash chain', async () => {
    let history = await appendHistoryEntry([], { run: 1 }, '2026-10-20T00:00:00.000Z')
    history = await appendHistoryEntry(history, { run: 2 }, '2026-10-20T00:01:00.000Z')
    expect(await verifyHistory(history)).toBe(true)
    const changed = history.map((entry, index) => index === 0 ? { ...entry, record: { run: 9 } } : entry)
    expect(await verifyHistory(changed)).toBe(false)
  })
})
