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

lattice = Lattice("Square", [Point(PID(1), (0.0, 0.0), (0.0, 0.0))],
    vectors=[[1.0, 0.0], [0.0, 1.0]],
    neighbors=1
    )

cell = Lattice("MagneticCell", [
        Point(PID(1), (0.0, 0.0), (0.0, 0.0)),
        Point(PID(2), (1.0, 0.0), (0.0, 0.0))
        ],
    vectors=[[1.0, 1.0], [1.0, -1.0]],
    neighbors=1
    )

hilbert = Hilbert(pid=>Spin{1//2}(1) for pid in cell.pids)

J = SpinTerm(:J, 1.0, 1, heisenberg"+-z")

magneticstructure = MagneticStructure(cell,
    Dict(pid=>(iseven(pid.site) ? [0, 0, 1] : [0, 0, -1]) for pid in cell.pids)
    )

antiferromagnet = Algorithm("SquareAFM", LSWT(lattice, hilbert, (J,), magneticstructure))

path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", len=100)
energybands = register!(antiferromagnet, :EB, TBAEB(path))
plot(energybands)
```