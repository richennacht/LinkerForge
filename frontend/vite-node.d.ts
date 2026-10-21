declare module 'node:child_process' {
  export function execFileSync(
    file: string,
    args: string[],
    options: { cwd: string; encoding: 'utf8' }
  ): string
}
