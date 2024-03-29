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
using Plots

function kitaev(bond)
    ϕ = bond|>rcoordinate|>azimuthd
    any(≈(ϕ), ( 90, 270)) && return Coupling(:, SID, ('z', 'z'))
    any(≈(ϕ), ( 30, 210)) && return Coupling(:, SID, ('x', 'x'))
    any(≈(ϕ), (150, 330)) && return Coupling(:, SID, ('y', 'y'))
    error("kitaev error: wrong azimuth angle($ϕ) bond.")
end
function gamma(bond)
    ϕ = bond|>rcoordinate|>azimuthd
    any(≈(ϕ), ( 90, 270)) && return MatrixCoupling(:, SID, Γ"z")
    any(≈(ϕ), ( 30, 210)) && return MatrixCoupling(:, SID, Γ"x")
    any(≈(ϕ), (150, 330)) && return MatrixCoupling(:, SID, Γ"y")
    error("kitaev error: wrong azimuth angle($ϕ) of bond.")
end

lattice = Lattice(
    (0.0, 0.0), (0.0, √3/3);
    name=:H2,
    vectors=[[1.0, 0.0], [0.5, √3/2]]
    )
hilbert = Hilbert(site=>Spin{3//2}() for site=1:length(lattice))

J₁ = SpinTerm(:J₁, 0.0, 1, MatrixCoupling(:, SID, Heisenberg""))
J₂ = SpinTerm(:J₂, -0.178, 2, MatrixCoupling(:, SID, Heisenberg""))
J₃ = SpinTerm(:J₃, 0.051, 3, MatrixCoupling(:, SID, Heisenberg""))
K = SpinTerm(:K, -4.288, 1, kitaev)
Γ = SpinTerm(:Γ, -0.044, 1, gamma)

magneticstructure = MagneticStructure(lattice, Dict(site=>[1, 1, 1] for site=1:length(lattice)))
CrBr₃ = Algorithm(:CrBr₃, LSWT(lattice, hilbert, (J₁, J₂, J₃, K, Γ), magneticstructure))
path = ReciprocalPath(reciprocals(lattice), (-2, -1)=>(2, 1), length=400)

spectra = CrBr₃(
    :INSS,
    InelasticNeutronScatteringSpectra(path, range(0.0, 15.0, length=301); fwhm=1.0, scale=log)
    )
energybands = CrBr₃(:EB, EnergyBands(path))
plt = plot()
plot!(plt, spectra)
plot!(plt, energybands, color=:white, linestyle=:dash)
yticks!(plt, range(0.0, 15.0, length=16))
```

## Berry curvature and Chern number of the magnon bands
The Berry curvatures and the Chern numbers of the magnon bands could be computed in the reciprocal unitcell:
```@example CrBr3
brillouin = BrillouinZone(reciprocals(lattice), 90)
berry = CrBr₃(:BerryCurvature, BerryCurvature(brillouin, [1, 2]));
plot(berry)
```
Here, k₁ and k₂ denote the coordinates in the reciprocal space along the two reciprocal vectors.

The Berry curvatures can also be computed in a reciprocal zone beyond the reciprocal unitcell:
```@example CrBr3
b₁, b₂ = 4*pi/√3*[1.0, 0.0], 4*pi/√3*[0.0, 1.0]
reciprocalzone = ReciprocalZone(
    [b₁, b₂], [-1.0=>1.0, -1.0=>1.0];
    length=301, ends=(true, true)
)
berry = CrBr₃(:BerryCurvatureEx, BerryCurvature(reciprocalzone, [1, 2]))
plot(berry)
```