var documenterSearchIndex = {"docs":
[{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/#Square-lattice-antiferromagnet","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/#Magnon-bands-by-linear-spin-wave-theory","page":"Square lattice antiferromagnet","title":"Magnon bands by linear spin wave theory","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"The following codes could compute the spin wave dispersions of the antiferromagnetic Heisenberg model on the square lattice.","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]])\ncell = Lattice([0.0, 0.0], [1.0, 0.0]; vectors=[[1.0, 1.0], [1.0, -1.0]])\nhilbert = Hilbert(site=>Spin{1//2}() for site=1:length(cell))\nJ = Heisenberg(:J, 1.0, 1)\nmagneticstructure = MagneticStructure(\n    cell,\n    Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))\n)\nantiferromagnet = Algorithm(:SquareAFM, LSWT(lattice, hilbert, J, magneticstructure))\n\npath = ReciprocalPath(reciprocals(lattice), rectangle\"Γ-X-M-Γ\", length=100)\nspectra = antiferromagnet(\n    :INSS,\n    InelasticNeutronScatteringSpectra(\n        path, range(0.0, 2.5, length=251);\n        fwhm=0.05, rescale=x->log(1+log(1+log(1+x)))\n        )\n)\nenergybands = antiferromagnet(:EB, EnergyBands(path))\n\nplt = plot()\nplot!(plt, spectra)\nplot!(plt, energybands, color=:white, linestyle=:dash)","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/#Auto-generation-of-the-analytical-expression-of-the-Hamiltonian-matrix-by-linear-spin-wave-theory","page":"Square lattice antiferromagnet","title":"Auto-generation of the analytical expression of the Hamiltonian matrix by linear spin wave theory","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"Combined with SymPy, it is also possible to get the analytical expression of the Hamiltonian in the matrix form obtained by linear wave theory:","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"using SymPy: Sym, symbols\nusing QuantumLattices\nusing SpinWaveTheory\n\nlattice = Lattice(\n    [zero(Sym), zero(Sym)];\n    name=:Square,\n    vectors=[[one(Sym), zero(Sym)], [zero(Sym), one(Sym)]]\n)\ncell = Lattice(\n    [zero(Sym), zero(Sym)], [one(Sym), zero(Sym)];\n    name=:MagneticCell,\n    vectors=[[one(Sym), one(Sym)], [one(Sym), -one(Sym)]]\n)\nhilbert = Hilbert(site=>Spin{1//2}() for site=1:length(cell))\nJ = Heisenberg(:J, symbols(\"J\", real=true), 1)\nmagneticstructure = MagneticStructure(\n    cell,\n    Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))\n)\nantiferromagnet = LSWT(lattice, hilbert, J, magneticstructure)\n\nk₁ = symbols(\"k₁\", real=true)\nk₂ = symbols(\"k₂\", real=true)\nm = matrix(antiferromagnet, [k₁, k₂]; infinitesimal=0)","category":"page"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/SquareLatticeFerromagnet/#Square-lattice-ferromagnet","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"","category":"section"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"The following codes could compute the spin wave dispersions of the ferromagnetic Heisenberg model on the square lattice.","category":"page"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]])\nhilbert = Hilbert(site=>Spin{1//2}() for site=1:length(lattice))\nJ = Heisenberg(:J, -1.0, 1)\nmagneticstructure = MagneticStructure(lattice, Dict(site=>[0, 0, 1] for site=1:length(lattice)))\nferromagnet = Algorithm(:SquareFM, LSWT(lattice, hilbert, J, magneticstructure))\n\npath = ReciprocalPath(reciprocals(lattice), rectangle\"Γ-X-M-Γ\", length=100)\nspectra = ferromagnet(\n    :INSS,\n    InelasticNeutronScatteringSpectra(path, range(0.0, 5.0, length=501); fwhm=0.2, rescale=x->log(1+x))\n)\nenergybands = ferromagnet(:dispersion, EnergyBands(path))\n\nplt = plot()\nplot!(plt, spectra)\nplot!(plt, energybands, color=:white, linestyle=:dash)","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/Introduction/#examples","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Here are some examples to illustrate how this package could be used.","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Pages = [\n        \"SquareLatticeFerromagnet.md\",\n        \"SquareLatticeAntiFerromagnet.md\",\n        \"CrBr3.md\",\n        ]\nDepth = 2","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"Modules = [SpinWaveTheory]","category":"page"},{"location":"manual/#QuantumLattices.DegreesOfFreedom.Hilbert-Tuple{QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Spin}, MagneticStructure}","page":"Manual","title":"QuantumLattices.DegreesOfFreedom.Hilbert","text":"Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure)\n\nGet the corresponding Hilbert space of the original one after the Holstein-Primakoff transformation. \n\n\n\n\n\n","category":"method"},{"location":"manual/#QuantumLattices.DegreesOfFreedom.Metric-Tuple{Magnonic, QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Fock{:b}}}","page":"Manual","title":"QuantumLattices.DegreesOfFreedom.Metric","text":"Metric(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) -> OperatorIndexToTuple\n\nGet the index-to-tuple metric for a quantum spin system after the Holstein-Primakoff transformation.\n\n\n\n\n\n","category":"method"},{"location":"manual/#SpinWaveTheory.HolsteinPrimakoff","page":"Manual","title":"SpinWaveTheory.HolsteinPrimakoff","text":"HolsteinPrimakoff{S<:Operators, U<:CoordinatedIndex, M<:MagneticStructure} <: UnitSubstitution{U, S}\n\nHolstein-Primakoff transformation.\n\n\n\n\n\n","category":"type"},{"location":"manual/#SpinWaveTheory.LSWT","page":"Manual","title":"SpinWaveTheory.LSWT","text":"LSWT{\n    K<:TBAKind{:BdG},\n    L<:AbstractLattice,\n    S<:OperatorGenerator,\n    HP<:HolsteinPrimakoff,\n    H₀<:CategorizedGenerator,\n    H₂<:CategorizedGenerator,\n    H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}},\n    C<:AbstractMatrix\n} <: TBA{K, H, C}\n\nLinear spin wave theory for magnetically ordered quantum lattice systems.\n\n\n\n\n\n","category":"type"},{"location":"manual/#SpinWaveTheory.LSWT-Tuple{QuantumLattices.Spatials.AbstractLattice, QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Spin}, Union{Tuple{QuantumLattices.DegreesOfFreedom.Term, Vararg{QuantumLattices.DegreesOfFreedom.Term, N}} where N, QuantumLattices.DegreesOfFreedom.Term}, MagneticStructure}","page":"Manual","title":"SpinWaveTheory.LSWT","text":"LSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, terms::OneOrMore{Term}, magneticstructure::MagneticStructure; neighbors::Union{Int, Neighbors}=nneighbor(terms))\n\nConstruct a LSWT.\n\n\n\n\n\n","category":"method"},{"location":"manual/#SpinWaveTheory.MagneticStructure","page":"Manual","title":"SpinWaveTheory.MagneticStructure","text":"MagneticStructure{L<:AbstractLattice, D<:Number}\n\nThe magnetic structure of an ordered quantum lattice system.\n\n\n\n\n\n","category":"type"},{"location":"manual/#SpinWaveTheory.MagneticStructure-Tuple{QuantumLattices.Spatials.AbstractLattice, Dict{Int64, <:Union{Tuple{Number, Number}, AbstractVector{<:Number}}}}","page":"Manual","title":"SpinWaveTheory.MagneticStructure","text":"MagneticStructure(cell::AbstractLattice, moments::Dict{Int, <:Union{AbstractVector{<:Number}, NTuple{2, Number}}}; unit::Symbol=:radian)\n\nConstruct the magnetic structure on a given lattice with the given moments.\n\n\n\n\n\n","category":"method"},{"location":"manual/#SpinWaveTheory.Magnonic","page":"Manual","title":"SpinWaveTheory.Magnonic","text":"Magnonic <: TBAKind{:BdG}\n\nMagnonic quantum lattice system.\n\n\n\n\n\n","category":"type"},{"location":"manual/#QuantumLattices.add!-Tuple{QuantumLattices.QuantumOperators.OperatorSum, TightBindingApproximation.Quadraticization{Magnonic}, QuantumLattices.QuantumOperators.Operator{<:Number, <:Tuple{QuantumLattices.DegreesOfFreedom.CoordinatedIndex{<:QuantumLattices.DegreesOfFreedom.Index{<:QuantumLattices.QuantumSystems.FockIndex{:b}}}, QuantumLattices.DegreesOfFreedom.CoordinatedIndex{<:QuantumLattices.DegreesOfFreedom.Index{<:QuantumLattices.QuantumSystems.FockIndex{:b}}}}}}","page":"Manual","title":"QuantumLattices.add!","text":"add!(dest::OperatorSum, qf::Quadraticization{Magnonic}, m::Operator{<:Number, <:ID{CoordinatedIndex{<:Index{<:FockIndex{:b}}}, 2}}; kwargs...) -> typeof(dest)\n\nGet the unified quadratic form of a rank-2 operator and add it to destination.\n\n\n\n\n\n","category":"method"},{"location":"manual/#SpinWaveTheory.rotation-Tuple{AbstractVector{<:Number}}","page":"Manual","title":"SpinWaveTheory.rotation","text":"rotation(destination::AbstractVector{<:Number}; kwargs...) -> SMatrix{3, 3}\nrotation(destination::Tuple{Number, Number}; unit::Symbol=:radian) -> SMatrix{3, 3}\n\nGet the rotation matrix which rotates [0, 0, 1] to the direction of the destination vector.\n\n\n\n\n\n","category":"method"},{"location":"manual/#TightBindingApproximation.commutator-Tuple{Magnonic, QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Fock{:b}}}","page":"Manual","title":"TightBindingApproximation.commutator","text":"commutator(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) -> Diagonal\n\nGet the commutation relation of the Holstein-Primakoff bosons.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory","page":"Home","title":"SpinWaveTheory","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: CI) (Image: codecov) (Image: ) (Image: ) (Image: 996.icu) (Image: LICENSE) (Image: LICENSE) (Image: Code Style: Blue) (Image: ColPrac: Contributor's Guide on Collaborative Practices for Community Packages)","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Spin wave theory for magnetically ordered quantum lattice systems based on the QuantumLattices pack and the TightBindingApproximation pack.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In Julia v1.8+, please type ] in the REPL to use the package mode, then type this command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add SpinWaveTheory","category":"page"},{"location":"#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Examples of spin wave theory for magnetically ordered quantum lattice systems","category":"page"},{"location":"#Note","page":"Home","title":"Note","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Due to the fast development of this package, releases with different minor version numbers are not guaranteed to be compatible with previous ones before the release of v1.0.0. Comments are welcomed in the issues.","category":"page"},{"location":"#Contact","page":"Home","title":"Contact","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"waltergu1989@gmail.com","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/CrBr3/#Ferromagnetic-spin-waves-in-the-monolayer-CrBr","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"","category":"section"},{"location":"examples/CrBr3/#Magnon-bands-and-inelastic-neutron-spectra-by-linear-spin-wave-theory","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Magnon bands and inelastic neutron spectra by linear spin wave theory","text":"","category":"section"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"The following codes could compute the ferromagnetic spin wave dispersions of the monolayer CrBr₃ by the K-Γ model:","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"note: Note\nK-Γ model is not the correct model to describe CrBr₃. See our paper. Here it is just for illustration on how to use this package.","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])\nhilbert = Hilbert(site=>Spin{3//2}() for site=1:length(lattice))\n\nJ₁ = Heisenberg(:J₁, 0.0, 1)\nJ₂ = Heisenberg(:J₂, -0.178, 2)\nJ₃ = Heisenberg(:J₃, 0.051, 3)\nK = Kitaev(:K, -4.288, 1; x=[30], y=[150], z=[270], unit=:degree)\nG = Γ(:Γ, -0.044, 1; x=[30], y=[150], z=[270], unit=:degree)\n\nmagneticstructure = MagneticStructure(lattice, Dict(site=>[1, 1, 1] for site=1:length(lattice)))\nCrBr₃ = Algorithm(:CrBr₃, LSWT(lattice, hilbert, (J₁, J₂, J₃, K, G), magneticstructure))\npath = ReciprocalPath(reciprocals(lattice), (-2, -1)=>(2, 1), length=400)\n\nspectra = CrBr₃(\n    :INSS,\n    InelasticNeutronScatteringSpectra(path, range(0.0, 15.0, length=301); fwhm=1.0, rescale=x->log(1+x))\n    )\nenergybands = CrBr₃(:EB, EnergyBands(path))\nplt = plot()\nplot!(plt, spectra)\nplot!(plt, energybands, color=:white, linestyle=:dash)\nyticks!(plt, range(0.0, 15.0, length=16))","category":"page"},{"location":"examples/CrBr3/#Berry-curvature-and-Chern-number-of-the-magnon-bands","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Berry curvature and Chern number of the magnon bands","text":"","category":"section"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"The Berry curvatures and the Chern numbers of the magnon bands could be computed in the reciprocal unitcell:","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"brillouin = BrillouinZone(reciprocals(lattice), 90)\nberry = CrBr₃(:BerryCurvature, BerryCurvature(brillouin, [1, 2]));\nplot(berry)","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"Here, k₁ and k₂ denote the coordinates in the reciprocal space along the two reciprocal vectors.","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"The Berry curvatures can also be computed in a reciprocal zone beyond the reciprocal unitcell:","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"b₁, b₂ = 4*pi/√3*[1.0, 0.0], 4*pi/√3*[0.0, 1.0]\nreciprocalzone = ReciprocalZone(\n    [b₁, b₂], [-1.0=>1.0, -1.0=>1.0];\n    length=301, ends=(true, true)\n)\nberry = CrBr₃(:BerryCurvatureEx, BerryCurvature(reciprocalzone, [1, 2]))\nplot(berry)","category":"page"}]
}
