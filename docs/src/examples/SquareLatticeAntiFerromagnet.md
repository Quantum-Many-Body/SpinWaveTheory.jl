```@meta
CurrentModule = SpinWaveTheory
```

# Square lattice antiferromagnet

The following codes could compute the spin wave dispersions of the antiferromagnetic Heisenberg model on the square lattice.

```@example AFM
using QuantumLattices
using TightBindingApproximation
using SpinWaveTheory
using Plots

lattice = Lattice(:Square,
    [Point(PID(1), [0.0, 0.0])],
    vectors=[[1.0, 0.0], [0.0, 1.0]],
    neighbors=1
    )

cell = Lattice(:MagneticCell,
    [Point(PID(1), [0.0, 0.0]), Point(PID(2), [1.0, 0.0])],
    vectors=[[1.0, 1.0], [1.0, -1.0]],
    neighbors=1
    )

hilbert = Hilbert(pid=>Spin{1//2}(1) for pid in cell.pids)

J = SpinTerm(:J, 1.0, 1, heisenberg"+-z")
h = SpinTerm(:h, -0.01, 0, sᶻ""*sᶻ"")

magneticstructure = MagneticStructure(cell,
    Dict(pid=>(iseven(pid.site) ? [0, 0, 1] : [0, 0, -1]) for pid in cell.pids)
    )

antiferromagnet = Algorithm(:SquareAFM, LSWT(lattice, hilbert, (J, h), magneticstructure))

path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", length=100)
ins = antiferromagnet(:INS,
    InelasticNeutronSpectra(path, range(0.0, 2.5, length=251); η=0.1, log=true)
    )
energybands = antiferromagnet(:EB, EnergyBands(path))

plt = plot()
plot!(plt, ins)
plot!(plt, energybands, color=:white, linestyle=:dash)
```