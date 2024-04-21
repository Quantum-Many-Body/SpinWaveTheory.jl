using Plots: plot, plot!, savefig
using QuantumLattices: atol, Algorithm, Heisenberg, Hilbert, Lattice, ReciprocalPath, Spin, Zeeman, reciprocals, @rectangle_str
using SpinWaveTheory
using TightBindingApproximation: EnergyBands, InelasticNeutronScatteringSpectra

@testset "SquareFM" begin
    lattice = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])
    hilbert = Hilbert(site=>Spin{1//2}() for site=1:length(lattice))
    J = Heisenberg(:J, -1.0, 1)
    h = Zeeman(:h, -0.5, 'z')
    ms = MagneticStructure(lattice, Dict(site=>[0, 0, 1] for site=1:length(lattice)))
    @test ms.rotations == MagneticStructure(lattice, Dict(site=>(0, 0) for site=1:length(lattice))).rotations
    lswt = Algorithm(:FM, LSWT(lattice, hilbert, (J, h), ms))

    path = ReciprocalPath(reciprocals(lattice), rectangle"Γ-X-M-Γ", length=8)
    data = lswt(:EBS, EnergyBands(path))[2].data[2]
    A(; k) = 2.5-cos(k[1])-cos(k[2])
    for (i, params) in enumerate(pairs(path))
        @test isapprox(A(; params...), data[i, 1], atol=atol)
        @test isapprox(A(; params...), data[i, 2], atol=atol)
    end

    path = ReciprocalPath(reciprocals(lattice), rectangle"Γ-X-M-Γ", length=100)
    eb = lswt(:EB, EnergyBands(path))
    spectra = lswt(:INSS, InelasticNeutronScatteringSpectra(path, range(0.0, 5.0, length=501); fwhm=0.1, scale=log))
    plt = plot()
    plot!(plt, spectra)
    plot!(plt, eb, color=:white, linestyle=:dash)
    display(plt)
    savefig("inelastic.png")
end
