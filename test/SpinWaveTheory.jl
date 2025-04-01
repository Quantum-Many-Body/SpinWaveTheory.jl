using LinearAlgebra: norm
using Plots: plot, plot!, savefig
using QuantumLattices: atol, Algorithm, Generator, Heisenberg, Hilbert, Lattice, Operator, Operators, ReciprocalPath, Spin, Zeeman, 𝕓, azimuth, azimuthd, bonds, expand, polar, polard, reciprocals, update!, @rectangle_str
using SpinWaveTheory
using TightBindingApproximation: EnergyBands, InelasticNeutronScatteringSpectra

@time @testset "rotation" begin
    input = rand(3)
    dest = input/norm(input)
    @test rotation(input)*[0, 0, 1] ≈ dest
    @test rotation((polar(input), azimuth(input)))*[0, 0, 1] ≈ rotation((polar(input), azimuth(input)); unit=:radian)*[0, 0, 1] ≈ dest
    @test rotation((polard(input), azimuthd(input)); unit=:degree)*[0, 0, 1] ≈ dest
end

@time @testset "MagneticStructure" begin
    cell = Lattice([0.0, 0.0], [1.0, 0.0])
    moments = Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))
    magneticstructure = MagneticStructure(cell, moments)
    @test magneticstructure.rotations[1] == [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]
    @test magneticstructure.rotations[2] == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
end

@time @testset "HolsteinPrimakoff" begin
    lattice = Lattice([0.0, 0.0], [1.0, 0.0])
    hilbert = Hilbert(Spin{1//2}(), length(lattice))
    J = Heisenberg(:J, -1.0, 1)
    ms = MagneticStructure(lattice, Dict(site=>iseven(site) ? [0, 0, 1] : [0, 0, -1] for site=1:length(lattice)))
    spins = Generator(bonds(lattice, 1), hilbert, J; half=false)
    hp = HolsteinPrimakoff{valtype(spins)}(ms)
    @test valtype(hp) == valtype(typeof(hp)) == valtype(typeof(hp), valtype(spins)) == valtype(typeof(hp), eltype(spins))
    bosons = expand(hp(spins))
    @test bosons == Operators(
        Operator(0.5, 𝕓(2, 1, 0, 1, [1.0, 0.0], [0.0, 0.0]), 𝕓(1, 1, 0, 1, [0.0, 0.0], [0.0, 0.0])),
        Operator(0.5, 𝕓(2, 1, 0, 2, [1.0, 0.0], [0.0, 0.0]), 𝕓(1, 1, 0, 2, [0.0, 0.0], [0.0, 0.0])),
        Operator(0.25),
        Operator(-0.5, 𝕓(1, 1, 0, 2, [0.0, 0.0], [0.0, 0.0]), 𝕓(1, 1, 0, 1, [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.5, 𝕓(2, 1, 0, 2, [1.0, 0.0], [0.0, 0.0]), 𝕓(2, 1, 0, 1, [1.0, 0.0], [0.0, 0.0])),
        Operator(1.0, 𝕓(2, 1, 0, 2, [1.0, 0.0], [0.0, 0.0]), 𝕓(2, 1, 0, 1, [1.0, 0.0], [0.0, 0.0]), 𝕓(1, 1, 0, 2, [0.0, 0.0], [0.0, 0.0]), 𝕓(1, 1, 0, 1, [0.0, 0.0], [0.0, 0.0]))
    )
    @test hp(bosons) == bosons
end

@time @testset "SquareFM" begin
    lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]])
    hilbert = Hilbert(Spin{1//2}(), length(lattice))
    J = Heisenberg(:J, -1.0, 1)
    h = Zeeman(:h, 0.0, 'z')
    ms = MagneticStructure(lattice, Dict(site=>[0, 0, 1] for site=1:length(lattice)))
    @test ms.rotations == MagneticStructure(lattice, Dict(site=>(0, 0) for site=1:length(lattice))).rotations
    lswt = Algorithm(:FM, LSWT(lattice, hilbert, (J, h), ms))

    update!(lswt; h=-0.5)
    path = ReciprocalPath(reciprocals(lattice), rectangle"Γ-X-M-Γ", length=8)
    data = lswt(:EBS, EnergyBands(path)).data.values
    A(k) = 2.5-cos(k[1])-cos(k[2])
    for (i, k) in enumerate(path)
        @test isapprox(A(k), data[i, 1], atol=10*atol)
        @test isapprox(A(k), data[i, 2], atol=10*atol)
    end

    path = ReciprocalPath(reciprocals(lattice), rectangle"Γ-X-M-Γ", length=100)
    eb = lswt(:EB, EnergyBands(path))
    spectra = lswt(:INSS, InelasticNeutronScatteringSpectra(path, range(0.0, 5.0, length=501); fwhm=0.1, rescale=x->log(1+x)))
    plt = plot()
    plot!(plt, spectra)
    plot!(plt, eb, color=:white, linestyle=:dash)
    display(plt)
    savefig("inelastic.png")
end
