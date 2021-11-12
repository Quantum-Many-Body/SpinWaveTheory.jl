using SpinWaveTheory
using QuantumLattices: Lattice, Point, PID, Hilbert, Spin, SpinTerm, @heisenberg_str, Algorithm, ReciprocalPath, @rectangle_str, register!, atol
using TightBindingApproximation: TBAEB

@testset "SquareFM" begin
    lattice = Lattice("S₁", [Point(PID(1), (0.0, 0.0), (0.0, 0.0))],
        vectors=[[1.0, 0.0], [0.0, 1.0]],
        neighbors=1
        )
    hilbert = Hilbert(pid=>Spin{1//2}(1) for pid in lattice.pids)
    J = SpinTerm(:J, -1.0, 1, heisenberg"xyz")
    ms = MagneticStructure(lattice, Dict(pid=>[0, 0, 1] for pid in lattice.pids))
    lsw = Algorithm("SquareFM", LSWT(lattice, hilbert, (J,), ms))
    path = ReciprocalPath(lattice.reciprocals, rectangle"Γ-X-M-Γ", len=8)
    data = register!(lsw, :EB, TBAEB(path))[2].data[2]

    A(; k) = 2-cos(k[1])-cos(k[2])
    for (i, params) in enumerate(path)
        @test isapprox(A(; params...), data[i, 1], atol=atol)
        @test isapprox(A(; params...), data[i, 2], atol=atol)
    end
end
