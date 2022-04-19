```@meta
CurrentModule = SpinWaveTheory
```

# Ferromagnetic spin waves in the monolayer CrBr₃

## Magnon bands and inelastic neutron spectra by linear spin wave theory

The following codes could compute the ferromagnetic spin wave dispersions of the monolayer CrBr₃ by the K-Γ model:
!!! note
    K-Γ model is not the correct model to describe CrBr₃. See our [paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.104.L020402). Here it is just for illustration on how to use this package.


```@example CrBr3
using QuantumLattices
using TightBindingApproximation
using SpinWaveTheory
using Plots; pyplot()

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

lattice = Lattice(:H2,
    [Point(PID(1), (0.0, 0.0)), Point(PID(2), (0.0, √3/3))],
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
CrBr3 = Algorithm(:CrBr3, LSWT(lattice, hilbert, (J₁, J₂, J₃, K, Γ), magneticstructure))
path = ReciprocalPath(lattice.reciprocals, (-2, -1)=>(2, 1), length=400)

ins = CrBr3(:INS,
    InelasticNeutronSpectra(path, range(0.0, 15.0, length=301); η=0.5, log=true)
    )
eb = CrBr3(:EB, EnergyBands(path))
plt = plot()
plot!(plt, ins)
plot!(plt, eb, color=:white, linestyle=:dash)
yticks!(plt, range(0.0, 15.0, length=16))
```

## Berry curvature and Chern number of the magnon bands
The Berry curvatures and the Chern numbers of the magnon bands could also be computed:
```@example CrBr3
brillouin = BrillouinZone(lattice.reciprocals, 90)
berry = CrBr3(:BerryCurvature, BerryCurvature(brillouin, [1, 2]));
plot(berry)
```