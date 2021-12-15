var documenterSearchIndex = {"docs":
[{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/#Square-lattice-antiferromagnet","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"The following codes could compute the spin wave dispersions of the antiferromagnetic Heisenberg model on the square lattice.","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice(\"Square\", [Point(PID(1), (0.0, 0.0), (0.0, 0.0))],\n    vectors=[[1.0, 0.0], [0.0, 1.0]],\n    neighbors=1\n    )\n\ncell = Lattice(\"MagneticCell\", [\n        Point(PID(1), (0.0, 0.0), (0.0, 0.0)),\n        Point(PID(2), (1.0, 0.0), (0.0, 0.0))\n        ],\n    vectors=[[1.0, 1.0], [1.0, -1.0]],\n    neighbors=1\n    )\n\nhilbert = Hilbert(pid=>Spin{1//2}(1) for pid in cell.pids)\n\nJ = SpinTerm(:J, 1.0, 1, heisenberg\"+-z\")\n\nmagneticstructure = MagneticStructure(cell,\n    Dict(pid=>(iseven(pid.site) ? [0, 0, 1] : [0, 0, -1]) for pid in cell.pids)\n    )\n\nantiferromagnet = Algorithm(\"SquareAFM\", LSWT(lattice, hilbert, (J,), magneticstructure))\n\npath = ReciprocalPath(lattice.reciprocals, rectangle\"Γ-X-M-Γ\", len=100)\nenergybands = register!(antiferromagnet, :EB, TBAEB(path))\nplot(energybands)","category":"page"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/SquareLatticeFerromagnet/#Square-lattice-ferromagnet","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"","category":"section"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"The following codes could compute the spin wave dispersions of the ferromagnetic Heisenberg model on the square lattice.","category":"page"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice(\"Square\", [Point(PID(1), (0.0, 0.0), (0.0, 0.0))],\n    vectors=[[1.0, 0.0], [0.0, 1.0]],\n    neighbors=1\n    )\n\nhilbert = Hilbert(pid=>Spin{1//2}(1) for pid in lattice.pids)\n\nJ = SpinTerm(:J, -1.0, 1, heisenberg\"+-z\")\n\nmagneticstructure = MagneticStructure(lattice, Dict(pid=>[0, 0, 1] for pid in lattice.pids))\n\nferromagnet = Algorithm(\"SquareFM\", LSWT(lattice, hilbert, (J,), magneticstructure))\n\npath = ReciprocalPath(lattice.reciprocals, rectangle\"Γ-X-M-Γ\", len=100)\nenergybands = register!(ferromagnet, :dispersion, TBAEB(path))\nplot(energybands)","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/Introduction/#examples","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Here are some examples to illustrate how this package could be used.","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Pages = [\n        \"SquareLatticeFerromagnet.md\",\n        \"SquareLatticeAntiFerromagnet.md\",\n        ]\nDepth = 1","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"#SpinWaveTheory","page":"Home","title":"SpinWaveTheory","text":"","category":"section"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Spin wave theory for magnetically ordered quantum lattice systems based on the QuantumLattices pack and the TightBindingApproximation pack.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In Julia v1.6+, please type ] in the REPL to use the package mode, then type this command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add SpinWaveTheory","category":"page"},{"location":"#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Examples of spin wave theory for magnetically ordered quantum lattice systems","category":"page"},{"location":"#Manuals","page":"Home","title":"Manuals","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [SpinWaveTheory]","category":"page"},{"location":"#QuantumLattices.Essentials.DegreesOfFreedom.Hilbert-Tuple{QuantumLattices.Essentials.DegreesOfFreedom.Hilbert{var\"#s4\", P, M} where {var\"#s4\"<:QuantumLattices.Essentials.QuantumSystems.Spin, P<:QuantumLattices.Essentials.Spatials.AbstractPID, M<:Function}, MagneticStructure}","page":"Home","title":"QuantumLattices.Essentials.DegreesOfFreedom.Hilbert","text":"Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure)\n\nGet the corresponding Hilbert space of the original one after the Holstein-Primakoff transformation. \n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.HPTransformation","page":"Home","title":"SpinWaveTheory.HPTransformation","text":"HPTransformation{S<:Operators, U<:OID{<:Index{<:AbstractPID, <:SimpleIID}}, M<:MagneticStructure} <: AbstractUnitSubstitution{U, S}\n\nHolstein-Primakoff transformation.\n\n\n\n\n\n","category":"type"},{"location":"#SpinWaveTheory.LSWT","page":"Home","title":"SpinWaveTheory.LSWT","text":"LSWT{L<:Lattice, H<:Generator, HP<:HPTransformation, E<:SimplifiedGenerator, F<:SimplifiedGenerator, G<:AbstractMatrix} <: AbstractTBA{TBAKind(:BdG), H, G}\n\nLinear spin wave theory for magnetically ordered quantum lattice systems.\n\n\n\n\n\n","category":"type"},{"location":"#SpinWaveTheory.LSWT-Tuple{QuantumLattices.Essentials.Spatials.Lattice, QuantumLattices.Essentials.DegreesOfFreedom.Hilbert, Tuple{Vararg{QuantumLattices.Essentials.DegreesOfFreedom.Term, N} where N}, MagneticStructure}","page":"Home","title":"SpinWaveTheory.LSWT","text":"LSWT(lattice::Lattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; boundary::Boundary=plain)\n\nConstruct a LSWT.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.MagneticStructure","page":"Home","title":"SpinWaveTheory.MagneticStructure","text":"MagneticStructure{L<:Lattice, P<:AbstractPID, D<:Number}\n\nThe magnetic structure of an ordered quantum lattice system.\n\n\n\n\n\n","category":"type"},{"location":"#SpinWaveTheory.MagneticStructure-Tuple{QuantumLattices.Essentials.Spatials.Lattice, Dict{var\"#s1\", var\"#s6\"} where {var\"#s1\"<:QuantumLattices.Essentials.Spatials.AbstractPID, var\"#s6\"<:(AbstractVector{T} where T)}}","page":"Home","title":"SpinWaveTheory.MagneticStructure","text":"MagneticStructure(cell::Lattice, moments::Dict{<:AbstractPID, <:AbstractVector})\n\nConstruct the magnetic structure on a given lattice with the given moments.\n\n\n\n\n\n","category":"method"},{"location":"#QuantumLattices.Essentials.QuantumOperators.matrix!-Tuple{LSWT}","page":"Home","title":"QuantumLattices.Essentials.QuantumOperators.matrix!","text":"matrix!(lswt::LSWT; k=nothing, atol=atol/5, kwargs...) -> TBAMatrix\n\nGet the tight-binding matrix representation of the linear spin waves.\n\nHere, the atol parameter is used to ensure that the matrix is positive-definite so that the Cholesky decomposition can be performed numerically.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.rotation-Tuple{AbstractVector{var\"#s3\"} where var\"#s3\"<:Number}","page":"Home","title":"SpinWaveTheory.rotation","text":"rotation(destination::AbstractVector{<:Number}) -> SMatrix{3, 3}\n\nGet the rotation matrix which rotates the [0, 0, 1] vector to the direction of the destination vector.\n\n\n\n\n\n","category":"method"}]
}