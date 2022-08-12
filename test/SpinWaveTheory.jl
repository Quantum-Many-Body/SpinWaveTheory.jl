using SpinWaveTheory
using Plots: plot, plot!, savefig
using QuantumLattices: Lattice, Point, PID, Hilbert, Spin, SpinTerm, @heisenberg_str, Algorithm, ReciprocalPath, @rectangle_str, atol
using TightBindingApproximation: EnergyBands, InelasticNeutronScatteringSpectra

@testset "SquareFM" begin
    lattice = Lattice(:Square,
        [Point(PID(1), [0.0, 0.0])],
        vectors=[[1.0, 0.0], [0.0, 1.0]],
        neighbors=1
        )
    hilbert = Hilbert(pid=>Spin{1//2}(1) for pid in lattice.pids)
    J = SpinTerm(:J, -1.0, 1, heisenberg"xyz")
    ms₁ = MagneticStructure(lattice, Dict(pid=>[0, 0, 1] for pid in lattice.pids))
    ms₂ = MagneticStructure(lattice, Dict(pid=>[0, 0] for pid in lattice.pids))
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
