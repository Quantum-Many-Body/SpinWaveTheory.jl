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
using Plots; pyplot()

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

## Auto-generation of the analytical expression of the Hamiltonian matrix by linear spin wave theory

Combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), it is also possible to get the analytical expression of the Hamiltonian in the matrix form obtained by linear wave theory:

```@example AFM-analytical
using SymPy: Sym, symbols
using QuantumLattices
using SpinWaveTheory

lattice = Lattice(:Square,
    [Point(PID(1), [zero(Sym), zero(Sym)])],
    vectors=[[one(Sym), zero(Sym)], [zero(Sym), one(Sym)]],
    neighbors=1
    )

cell = Lattice(:MagneticCell,
    [Point(PID(1), [zero(Sym), zero(Sym)]), Point(PID(2), [one(Sym), zero(Sym)])],
    vectors=[[one(Sym), one(Sym)], [one(Sym), -one(Sym)]],
    neighbors=1
    )

hilbert = Hilbert(pid=>Spin{1//2}(1) for pid in cell.pids)

J = SpinTerm(:J, symbols("J", real=true), 1, heisenberg"+-z")

magneticstructure = MagneticStructure(cell,
    Dict(pid=>(iseven(pid.site) ? [0, 0, 1] : [0, 0, -1]) for pid in cell.pids)
    )

antiferromagnet = LSWT(lattice, hilbert, (J,), magneticstructure)
k₁ = symbols("k₁", real=true)
k₂ = symbols("k₂", real=true)
m = matrix(antiferromagnet; k=[k₁, k₂], atol=0)
```