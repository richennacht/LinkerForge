import { fireEvent, render, screen } from '@testing-library/react'
import { describe, expect, it } from 'vitest'
import App from './App'

describe('LinkerForge interface', () => {
  it('states the scientific limitation', () => {
    render(<App />)
    expect(screen.getByText(/makes no claim of biological function/i)).toBeInTheDocument()
    expect(screen.getByText(/source data not bundled here/i)).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /load preset example candidates/i })).toBeInTheDocument()
  })

  it('shows the complete four-stage workflow', () => {
    render(<App />)
    expect(screen.getByText('Check inputs')).toBeInTheDocument()
    expect(screen.getByText('Predict and compare')).toBeInTheDocument()
  })

  it('does not present an individual reference as a predicted fusion', () => {
    render(<App />)
    expect(screen.getByText(/no full-fusion structure loaded/i)).toBeInTheDocument()
    fireEvent.click(screen.getByRole('button', { name: /individual experimental references/i }))
    expect(screen.getByTitle(/sfGFP experimental structure/i)).toHaveAttribute('src', 'https://molstar.org/viewer/?pdb=2B3P')
  })

  it('accepts two custom domains and offers full-fusion templates without invented coordinates', () => {
    render(<App />)
    fireEvent.change(screen.getByLabelText(/domain a amino-acid sequence/i), { target: { value: 'ACDEFG' } })
    fireEvent.change(screen.getByLabelText(/domain b amino-acid sequence/i), { target: { value: 'HIKLMN' } })
    fireEvent.click(screen.getByRole('button', { name: /generate full-fusion template candidates/i }))
    expect(screen.getByText(/full fusion: 27 residues/i)).toBeInTheDocument()
    expect(screen.getByText(/no full-fusion structure loaded/i)).toBeInTheDocument()
    fireEvent.click(screen.getByRole('button', { name: /open full screen/i }))
    expect(screen.getByRole('button', { name: /close full screen/i })).toBeInTheDocument()
  })
})
