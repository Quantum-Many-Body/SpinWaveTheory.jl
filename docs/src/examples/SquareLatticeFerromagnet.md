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

lattice = Lattice(:Square,
    [Point(PID(1), [0.0, 0.0])],
    vectors=[[1.0, 0.0], [0.0, 1.0]],
    neighbors=1
    )

hilbert = Hilbert(pid=>Spin{1//2}(1) for pid in lattice.pids)

J = SpinTerm(:J, -1.0, 1, heisenberg"+-z")

magneticstructure = MagneticStructure(lattice, Dict(pid=>[0, 0, 1] for pid in lattice.pids))

ferromagnet = Algorithm(:SquareFM, LSWT(lattice, hilbert, (J,), magneticstructure))

path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", length=100)

ins = ferromagnet(:INS,
    InelasticNeutronSpectra(path, range(0.0, 5.0, length=501); η=0.3, log=true)
    )
energybands = ferromagnet(:dispersion, EnergyBands(path))

plt = plot()
plot!(plt, ins)
plot!(plt, energybands, color=:white, linestyle=:dash)
```
