import { fireEvent, render, screen, within } from '@testing-library/react'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import App from './App'

function submitTestDesign(name = 'Public test fusion') {
  fireEvent.change(screen.getByLabelText(/design name/i), { target: { value: name } })
  fireEvent.change(screen.getByLabelText(/domain a amino-acid sequence/i), { target: { value: 'ACDEFG' } })
  fireEvent.change(screen.getByLabelText(/domain b amino-acid sequence/i), { target: { value: 'HIKLMN' } })
  fireEvent.click(screen.getByRole('button', { name: /create prediction workspace/i }))
}

describe('LinkerForge interface', () => {
  beforeEach(() => {
    localStorage.clear()
    localStorage.setItem('linkerforge.color-theme.v1', 'light')
    document.documentElement.removeAttribute('data-theme')
    document.documentElement.style.colorScheme = ''
    window.location.hash = '#structure'
  })

  it('switches and persists the color theme shared by web and desktop', () => {
    render(<App />)
    expect(document.documentElement).toHaveAttribute('data-theme', 'light')

    fireEvent.click(screen.getByRole('button', { name: /switch to dark mode/i }))

    expect(document.documentElement).toHaveAttribute('data-theme', 'dark')
    expect(document.documentElement.style.colorScheme).toBe('dark')
    expect(localStorage.getItem('linkerforge.color-theme.v1')).toBe('dark')
    expect(screen.getByRole('button', { name: /switch to light mode/i })).toBeInTheDocument()
  })

  it('opens as a single-window application with its input composer in the sidebar', () => {
    render(<App />)
    expect(screen.getByLabelText(/new linker design input/i)).toBeInTheDocument()
    expect(screen.getByRole('navigation', { name: /workspace views/i })).toBeInTheDocument()
    expect(screen.getByLabelText(/saved design sessions/i)).toBeInTheDocument()
    expect(screen.getByText(/computational evidence only/i)).toBeInTheDocument()
    expect(screen.getByLabelText(/domain a amino-acid sequence/i)).toBeInTheDocument()
    expect(screen.getByRole('region', { name: /biological context and safeguards/i })).toBeInTheDocument()
    expect(screen.getByLabelText(/biological objective/i)).toBeVisible()
    expect(screen.getByRole('heading', { name: /start a linker prediction/i })).toBeInTheDocument()
    expect(screen.queryByText(/^overview$/i)).not.toBeInTheDocument()
  })

  it('keeps the labelled preset action in the Candidates view', () => {
    render(<App />)
    expect(screen.queryByRole('button', { name: /load labelled preset example/i })).not.toBeInTheDocument()
    fireEvent.click(screen.getByRole('button', { name: /^candidates$/i }))
    expect(screen.getByRole('button', { name: /load labelled preset example/i })).toBeInTheDocument()
  })

  it('turns a submitted input into a sidebar session and opens prediction models', () => {
    render(<App />)
    expect(screen.getByRole('slider')).toHaveAttribute('min', '5')
    submitTestDesign()

    const sessions = screen.getByLabelText(/saved design sessions/i)
    expect(within(sessions).getByRole('button', { name: /open public test fusion/i })).toBeInTheDocument()
    expect(within(sessions).getByRole('button', { name: /edit public test fusion/i })).toBeInTheDocument()
    expect(within(sessions).getByRole('button', { name: /delete public test fusion/i })).toBeInTheDocument()
    expect(screen.getByRole('tablist', { name: /prediction models/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /full fusion model/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /domain a.*individual model/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /domain b.*individual model/i })).toBeInTheDocument()
    expect(screen.getByLabelText(/candidate models/i)).toBeInTheDocument()
    expect(screen.getByText(/no full-fusion structure loaded/i)).toBeInTheDocument()
  })

  it('reopens a saved design from the chat-style session list', () => {
    render(<App />)
    submitTestDesign('Reusable fusion')
    fireEvent.click(screen.getByRole('button', { name: /new linker design/i }))
    expect(screen.getByRole('heading', { name: /start a linker prediction/i })).toBeInTheDocument()

    fireEvent.click(within(screen.getByLabelText(/saved design sessions/i)).getByRole('button', { name: /open reusable fusion/i }))
    fireEvent.click(screen.getByRole('button', { name: /edit design inputs/i }))
    expect(screen.getByLabelText(/domain a amino-acid sequence/i)).toHaveValue('ACDEFG')
    expect(screen.getByLabelText(/domain b amino-acid sequence/i)).toHaveValue('HIKLMN')
  })

  it('edits a saved design in place and deletes it with confirmation', () => {
    const confirm = vi.spyOn(window, 'confirm').mockReturnValue(true)
    render(<App />)
    submitTestDesign('Duplicate-prone design')

    const sessions = screen.getByLabelText(/saved design sessions/i)
    fireEvent.click(within(sessions).getByRole('button', { name: /edit duplicate-prone design/i }))
    fireEvent.change(screen.getByLabelText(/design name/i), { target: { value: 'Edited design' } })
    fireEvent.click(screen.getByRole('button', { name: /create prediction workspace/i }))

    expect(within(sessions).queryByRole('button', { name: /open duplicate-prone design/i })).not.toBeInTheDocument()
    expect(within(sessions).getByRole('button', { name: /open edited design/i })).toBeInTheDocument()
    fireEvent.click(within(sessions).getByRole('button', { name: /delete edited design/i }))
    expect(confirm).toHaveBeenCalledOnce()
    expect(within(sessions).getByText(/no designs yet/i)).toBeInTheDocument()
    confirm.mockRestore()
  })

  it('uses per-input domain models instead of hard-coded references', () => {
    render(<App />)
    submitTestDesign()
    fireEvent.click(screen.getByRole('button', { name: /domain a.*individual model/i }))
    expect(screen.getByText(/no reference structure loaded/i)).toBeInTheDocument()
    expect(screen.queryByTitle(/sfGFP experimental structure/i)).not.toBeInTheDocument()
    expect(screen.getByText(/does not guess a reference structure/i)).toBeInTheDocument()
  })

  it('keeps candidate review and evidence results accessible in the same window', () => {
    render(<App />)
    submitTestDesign()
    fireEvent.click(screen.getByRole('button', { name: /^candidates$/i }))
    expect(screen.getByText(/generated 4 auditable design candidates and required controls/i)).toBeInTheDocument()
    expect(within(screen.getByLabelText(/custom fusion candidates/i)).getAllByRole('button')).toHaveLength(4)

    fireEvent.click(screen.getByRole('button', { name: /^results$/i }))
    expect(screen.getByRole('heading', { name: /structure prediction, confidence, and validation status/i })).toBeInTheDocument()
    expect(screen.getByText(/evidence status and available alternatives/i)).toBeInTheDocument()
    expect(screen.getByText(/cpu report workflow available/i)).toBeInTheDocument()
    expect(screen.getByText(/partial workflow available/i)).toBeInTheDocument()
    expect(screen.getByText(/has not run for this candidate/i)).toBeInTheDocument()
    expect(screen.getByText(/7 public structures screened/i)).toBeInTheDocument()
    expect(screen.getByText(/48 unique severe-overlap pairs/i)).toBeInTheDocument()
    expect(screen.getAllByText(/future scope · on hold/i)).not.toHaveLength(0)
    expect(screen.getByText(/31\.76 GPU hours/i)).toBeInTheDocument()
    expect(screen.getByText(/AlphaFold2 ensemble agreement/i)).toBeInTheDocument()
    expect(screen.getByText(/Reference-domain lDDT \/ TM-score/i)).toBeInTheDocument()
    expect(screen.getByText(/Automatic functional-site mapping/i)).toBeInTheDocument()
    expect(screen.getByText(/DSSP junction analysis/i)).toBeInTheDocument()
    expect(screen.getByText(/Post-OpenMM change analysis/i)).toBeInTheDocument()
    expect(screen.getByText(/not the canonical MolProbity clashscore/i)).toBeInTheDocument()
    expect(screen.getByText(/spent result available/i)).toBeInTheDocument()
    expect(screen.getByText(/best historical LinkerForge head scored APS 0\.089 \/ Fmax 0\.165/i)).toBeInTheDocument()
    expect(screen.getByText(/cannot be reused for tuning or model selection/i)).toBeInTheDocument()
    expect(screen.getByText(/no structure prediction has been started/i)).toBeInTheDocument()
    expect(screen.queryByText(/^experimental validation$/i)).not.toBeInTheDocument()
  })
})
