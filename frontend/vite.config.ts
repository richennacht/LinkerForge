import { defineConfig } from 'vitest/config'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
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
