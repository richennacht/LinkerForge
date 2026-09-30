import { render, screen } from '@testing-library/react'
import { describe, expect, it } from 'vitest'
import App from './App'

describe('LinkerForge interface', () => {
  it('states the scientific limitation', () => {
    render(<App />)
    expect(screen.getByText(/makes no claim of biological function/i)).toBeInTheDocument()
  })

  it('shows the complete four-stage workflow', () => {
    render(<App />)
    expect(screen.getByText('Check inputs')).toBeInTheDocument()
    expect(screen.getByText('Predict and compare')).toBeInTheDocument()
  })

  it('labels reference structures separately from predictions', () => {
    render(<App />)
    expect(screen.getByText(/does not run AlphaFold/i)).toBeInTheDocument()
    expect(screen.getAllByText(/PDB 2B3P/i)).toHaveLength(2)
    expect(screen.getByTitle(/sfGFP experimental structure/i)).toHaveAttribute('src', 'https://molstar.org/viewer/?pdb=2B3P')
  })
})
