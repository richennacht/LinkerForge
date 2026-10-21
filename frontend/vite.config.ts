import { defineConfig } from 'vitest/config'
import react from '@vitejs/plugin-react'
import { execFileSync } from 'node:child_process'

function gitValue(args: string[], fallback: string) {
  try {
    return execFileSync('git', args, { cwd: '..', encoding: 'utf8' }).trim()
  } catch {
    return fallback
  }
}

const commit = gitValue(['rev-parse', 'HEAD'], 'not-available')
const dirty = gitValue(['status', '--porcelain'], '') ? 'true' : 'false'

export default defineConfig({
  plugins: [react()],
  define: {
    'import.meta.env.VITE_LINKERFORGE_COMMIT': JSON.stringify(commit),
    'import.meta.env.VITE_LINKERFORGE_DIRTY': JSON.stringify(dirty),
  },
  // Relative assets allow the same verified build to run on Vercel previews,
  // the Vercel production domain, and the existing GitHub Pages fallback.
  base: './',
  server: {
    strictPort: true,
    watch: { ignored: ['**/src-tauri/**'] },
  },
  test: {
    environment: 'jsdom',
    setupFiles: './src/test/setup.ts',
    exclude: ['src-tauri/**', 'node_modules/**'],
  },
})
