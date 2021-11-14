```@meta
CurrentModule = SpinWaveTheory
```

# Ferromagnetic spin waves in the monolayer CrBr₃

The following codes could compute the ferromagnetic spin wave dispersions of the monolayer CrBr₃ by the K-Γ model:
!!! note
    K-Γ model is not the correct model to describe CrBr₃. See our [paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.104.L020402).


```@example CrBr3
using QuantumLattices
using TightBindingApproximation
using SpinWaveTheory
using Plots

function kitaev(bond)
    theta = bond|>rcoord|>azimuthd
    (theta≈90  || theta≈270) && return ising"z"
    (theta≈210 || theta≈30 ) && return ising"x"
    (theta≈330 || theta≈150) && return ising"y"
    error("kitaev error: wrong $theta.")
end
function gamma(bond)
    theta = bond|>rcoord|>azimuthd
    (theta≈90  || theta≈270) && return gamma"z"
    (theta≈210 || theta≈30 ) && return gamma"x"
    (theta≈330 || theta≈150) && return gamma"y"
    error("kitaev error: wrong $theta.")
end

lattice = Lattice( "H2", [
                Point(PID(1), (0.0, 0.0), (0.0, 0.0)),
                Point(PID(2), (0.0, √3/3), (0.0, 0.0))
                ],
            vectors=[[1.0, 0.0], [0.5, √3/2]],
            neighbors=3
            )
hilbert = Hilbert(pid=>Spin{3//2}(1) for pid in lattice.pids)

J₁ = SpinTerm(:J₁, 0.0, 1, couplings=heisenberg"xyz", modulate=true)
J₂ = SpinTerm(:J₂, -0.178, 2, couplings=heisenberg"xyz", modulate=true)
J₃ = SpinTerm(:J₃, 0.051, 3, couplings=heisenberg"xyz", modulate=true)
K = SpinTerm(:K, -4.288, 1, couplings=kitaev, modulate=true)
Γ = SpinTerm(:Γ, -0.044, 1, couplings=gamma, modulate=true)

magneticstructure = MagneticStructure(lattice, Dict(pid=>[1, 1, 1] for pid in lattice.pids))
CrBr3 = Algorithm("CrBr3", LSWT(lattice, hilbert, (J₁, J₂, J₃, K, Γ), magneticstructure))
path = ReciprocalPath(lattice.reciprocals, (-2, -1)=>(2, 1), len=400)

ins = register!(CrBr3, :INS,
    InelasticNeutronSpectra(path, range(0.0, 15.0, length=301); eta=0.5, log=true)
    )
eb = register!(CrBr3, :EB, EnergyBands(path))
plt = plot()
plot!(plt, ins)
plot!(plt, eb, color=:white, linestyle=:dash)
yticks!(plt, range(0.0, 15.0, length=16))
```