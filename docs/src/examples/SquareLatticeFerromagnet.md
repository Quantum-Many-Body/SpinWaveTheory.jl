```@meta
CurrentModule = SpinWaveTheory
```

# Square lattice ferromagnet

The following codes could compute the spin wave dispersions of the ferromagnetic Heisenberg model on the square lattice.

```@example FM
using QuantumLattices
using TightBindingApproximation
using SpinWaveTheory
using Plots; pyplot()

lattice = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])
hilbert = Hilbert(site=>Spin{1//2}() for site=1:length(lattice))
J = SpinTerm(:J, -1.0, 1, MatrixCoupling(:, SID, Heisenberg""))
magneticstructure = MagneticStructure(lattice, Dict(site=>[0, 0, 1] for site=1:length(lattice)))
ferromagnet = Algorithm(:SquareFM, LSWT(lattice, hilbert, (J,), magneticstructure))

path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", length=100)
ins = ferromagnet(
    :INS,
    InelasticNeutronScatteringSpectra(path, range(0.0, 5.0, length=501); η=0.3, log=true)
)
energybands = ferromagnet(:dispersion, EnergyBands(path))

plt = plot()
plot!(plt, ins)
plot!(plt, energybands, color=:white, linestyle=:dash)
```
