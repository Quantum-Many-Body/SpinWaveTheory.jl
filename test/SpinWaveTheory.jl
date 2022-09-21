using Plots: plot, plot!, savefig
using QuantumLattices: atol, Algorithm, Hilbert, Lattice, MatrixCoupling, ReciprocalPath, SID, Spin, SpinTerm, @Heisenberg_str, @rectangle_str
using SpinWaveTheory
using TightBindingApproximation: EnergyBands, InelasticNeutronScatteringSpectra

@testset "SquareFM" begin
    lattice = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])
    hilbert = Hilbert(site=>Spin{1//2}() for site=1:length(lattice))
    J = SpinTerm(:J, -1.0, 1, MatrixCoupling(:, SID, Heisenberg""))
    ms₁ = MagneticStructure(lattice, Dict(site=>[0, 0, 1] for site=1:length(lattice)))
    ms₂ = MagneticStructure(lattice, Dict(site=>[0, 0] for site=1:length(lattice)))
    @test ms₁.rotations == ms₂.rotations
    lswt = Algorithm(:FM, LSWT(lattice, hilbert, (J,), ms₂))

    path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", length=8)
    data = lswt(:EBS, EnergyBands(path))[2].data[2]
    A(; k) = 2-cos(k[1])-cos(k[2])
    for (i, params) in enumerate(pairs(path))
        @test isapprox(A(; params...), data[i, 1], atol=atol)
        @test isapprox(A(; params...), data[i, 2], atol=atol)
    end

    path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", length=100)
    ebs = lswt(:EBS, EnergyBands(path))
    ins = lswt(:INS, InelasticNeutronScatteringSpectra(path, range(0.0, 5.0, length=501); η=0.3, log=true))
    plt = plot()
    plot!(plt, ins)
    plot!(plt, ebs, color=:white, linestyle=:dash)
    display(plt)
    savefig("spectra.png")
end
