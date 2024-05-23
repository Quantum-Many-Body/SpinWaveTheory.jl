```@meta
CurrentModule = SpinWaveTheory
```

# Square lattice antiferromagnet

## Magnon bands by linear spin wave theory

The following codes could compute the spin wave dispersions of the antiferromagnetic Heisenberg model on the square lattice.

```@example AFM
using QuantumLattices
using TightBindingApproximation
using SpinWaveTheory
using Plots

lattice = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])
cell = Lattice(
    [0.0, 0.0], [1.0, 0.0];
    name=:MagneticCell,
    vectors=[[1.0, 1.0], [1.0, -1.0]]
)
hilbert = Hilbert(site=>Spin{1//2}() for site=1:length(cell))
J = Heisenberg(:J, 1.0, 1)
magneticstructure = MagneticStructure(
    cell,
    Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))
)
antiferromagnet = Algorithm(:SquareAFM, LSWT(lattice, hilbert, J, magneticstructure))

path = ReciprocalPath(reciprocals(lattice), rectangle"Γ-X-M-Γ", length=100)
spectra = antiferromagnet(
    :INSS,
    InelasticNeutronScatteringSpectra(
        path, range(0.0, 2.5, length=251);
        fwhm=0.05, scale=x->log(1+log(1+log(1+x)))
        )
)
energybands = antiferromagnet(:EB, EnergyBands(path))

plt = plot()
plot!(plt, spectra)
plot!(plt, energybands, color=:white, linestyle=:dash)
```

## Auto-generation of the analytical expression of the Hamiltonian matrix by linear spin wave theory

Combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), it is also possible to get the analytical expression of the Hamiltonian in the matrix form obtained by linear wave theory:

```@example AFM-analytical
using SymPy: Sym, symbols
using QuantumLattices
using SpinWaveTheory

lattice = Lattice(
    [zero(Sym), zero(Sym)];
    name=:Square,
    vectors=[[one(Sym), zero(Sym)], [zero(Sym), one(Sym)]]
)
cell = Lattice(
    [zero(Sym), zero(Sym)], [one(Sym), zero(Sym)];
    name=:MagneticCell,
    vectors=[[one(Sym), one(Sym)], [one(Sym), -one(Sym)]]
)
hilbert = Hilbert(site=>Spin{1//2}() for site=1:length(cell))
J = Heisenberg(:J, symbols("J", real=true), 1)
magneticstructure = MagneticStructure(
    cell,
    Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))
)
antiferromagnet = LSWT(lattice, hilbert, J, magneticstructure)

k₁ = symbols("k₁", real=true)
k₂ = symbols("k₂", real=true)
m = matrix(antiferromagnet; k=[k₁, k₂], infinitesimal=0)
```