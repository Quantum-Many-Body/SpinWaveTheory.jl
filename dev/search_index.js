var documenterSearchIndex = {"docs":
[{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/#Square-lattice-antiferromagnet","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/#Magnon-bands-by-linear-spin-wave-theory","page":"Square lattice antiferromagnet","title":"Magnon bands by linear spin wave theory","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"The following codes could compute the spin wave dispersions of the antiferromagnetic Heisenberg model on the square lattice.","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])\ncell = Lattice(\n    [0.0, 0.0], [1.0, 0.0];\n    name=:MagneticCell,\n    vectors=[[1.0, 1.0], [1.0, -1.0]]\n)\nhilbert = Hilbert(site=>Spin{1//2}() for site=1:length(cell))\nJ = SpinTerm(:J, 1.0, 1, MatrixCoupling(:, SID, Heisenberg\"\"))\nmagneticstructure = MagneticStructure(\n    cell,\n    Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))\n)\nantiferromagnet = Algorithm(:SquareAFM, LSWT(lattice, hilbert, J, magneticstructure))\n\npath = ReciprocalPath(reciprocals(lattice), rectangle\"Γ-X-M-Γ\", length=100)\nspectra = antiferromagnet(\n    :INSS,\n    InelasticNeutronScatteringSpectra(\n        path, range(0.0, 2.5, length=251);\n        fwhm=0.05, scale=x->log(1+log(1+log(1+x)))\n        )\n)\nenergybands = antiferromagnet(:EB, EnergyBands(path))\n\nplt = plot()\nplot!(plt, spectra)\nplot!(plt, energybands, color=:white, linestyle=:dash)","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/#Auto-generation-of-the-analytical-expression-of-the-Hamiltonian-matrix-by-linear-spin-wave-theory","page":"Square lattice antiferromagnet","title":"Auto-generation of the analytical expression of the Hamiltonian matrix by linear spin wave theory","text":"","category":"section"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"Combined with SymPy, it is also possible to get the analytical expression of the Hamiltonian in the matrix form obtained by linear wave theory:","category":"page"},{"location":"examples/SquareLatticeAntiFerromagnet/","page":"Square lattice antiferromagnet","title":"Square lattice antiferromagnet","text":"using SymPy: Sym, symbols\nusing QuantumLattices\nusing SpinWaveTheory\n\nlattice = Lattice(\n    [zero(Sym), zero(Sym)];\n    name=:Square,\n    vectors=[[one(Sym), zero(Sym)], [zero(Sym), one(Sym)]]\n)\ncell = Lattice(\n    [zero(Sym), zero(Sym)], [one(Sym), zero(Sym)];\n    name=:MagneticCell,\n    vectors=[[one(Sym), one(Sym)], [one(Sym), -one(Sym)]]\n)\nhilbert = Hilbert(site=>Spin{1//2}() for site=1:length(cell))\nJ = SpinTerm(:J, symbols(\"J\", real=true), 1, MatrixCoupling(:, SID, Heisenberg\"\"))\nmagneticstructure = MagneticStructure(\n    cell,\n    Dict(site=>(iseven(site) ? [0, 0, 1] : [0, 0, -1]) for site=1:length(cell))\n)\nantiferromagnet = LSWT(lattice, hilbert, J, magneticstructure)\n\nk₁ = symbols(\"k₁\", real=true)\nk₂ = symbols(\"k₂\", real=true)\nm = matrix(antiferromagnet; k=[k₁, k₂], atol=0)","category":"page"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/SquareLatticeFerromagnet/#Square-lattice-ferromagnet","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"","category":"section"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"The following codes could compute the spin wave dispersions of the ferromagnetic Heisenberg model on the square lattice.","category":"page"},{"location":"examples/SquareLatticeFerromagnet/","page":"Square lattice ferromagnet","title":"Square lattice ferromagnet","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nlattice = Lattice([0.0, 0.0]; name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])\nhilbert = Hilbert(site=>Spin{1//2}() for site=1:length(lattice))\nJ = SpinTerm(:J, -1.0, 1, MatrixCoupling(:, SID, Heisenberg\"\"))\nmagneticstructure = MagneticStructure(lattice, Dict(site=>[0, 0, 1] for site=1:length(lattice)))\nferromagnet = Algorithm(:SquareFM, LSWT(lattice, hilbert, J, magneticstructure))\n\npath = ReciprocalPath(reciprocals(lattice), rectangle\"Γ-X-M-Γ\", length=100)\nspectra = ferromagnet(\n    :INSS,\n    InelasticNeutronScatteringSpectra(path, range(0.0, 5.0, length=501); fwhm=0.2, scale=log)\n)\nenergybands = ferromagnet(:dispersion, EnergyBands(path))\n\nplt = plot()\nplot!(plt, spectra)\nplot!(plt, energybands, color=:white, linestyle=:dash)","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/Introduction/#examples","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Here are some examples to illustrate how this package could be used.","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Pages = [\n        \"SquareLatticeFerromagnet.md\",\n        \"SquareLatticeAntiFerromagnet.md\",\n        \"CrBr3.md\",\n        ]\nDepth = 2","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"#SpinWaveTheory","page":"Home","title":"SpinWaveTheory","text":"","category":"section"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Spin wave theory for magnetically ordered quantum lattice systems based on the QuantumLattices pack and the TightBindingApproximation pack.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In Julia v1.8+, please type ] in the REPL to use the package mode, then type this command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add SpinWaveTheory","category":"page"},{"location":"#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Examples of spin wave theory for magnetically ordered quantum lattice systems","category":"page"},{"location":"#Manuals","page":"Home","title":"Manuals","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [SpinWaveTheory]","category":"page"},{"location":"#QuantumLattices.DegreesOfFreedom.Hilbert-Tuple{QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Spin}, MagneticStructure}","page":"Home","title":"QuantumLattices.DegreesOfFreedom.Hilbert","text":"Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure)\n\nGet the corresponding Hilbert space of the original one after the Holstein-Primakoff transformation. \n\n\n\n\n\n","category":"method"},{"location":"#QuantumLattices.DegreesOfFreedom.Metric-Tuple{Magnonic, QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Fock{:b}}}","page":"Home","title":"QuantumLattices.DegreesOfFreedom.Metric","text":"Metric(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) -> OperatorUnitToTuple\n\nGet the index-to-tuple metric for a quantum spin system after the Holstein-Primakoff transformation.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.HPTransformation","page":"Home","title":"SpinWaveTheory.HPTransformation","text":"HPTransformation{S<:Operators, U<:CompositeIndex, M<:MagneticStructure} <: UnitSubstitution{U, S}\n\nHolstein-Primakoff transformation.\n\n\n\n\n\n","category":"type"},{"location":"#SpinWaveTheory.LSWT","page":"Home","title":"SpinWaveTheory.LSWT","text":"LSWT{K<:TBAKind{:BdG}, L<:AbstractLattice, Hₛ<:OperatorGenerator, HP<:HPTransformation, Ω<:Image, H<:Image} <: AbstractTBA{K, H, AbstractMatrix}\n\nLinear spin wave theory for magnetically ordered quantum lattice systems.\n\n\n\n\n\n","category":"type"},{"location":"#SpinWaveTheory.LSWT-Tuple{QuantumLattices.Spatials.AbstractLattice, QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Spin}, QuantumLattices.DegreesOfFreedom.Term, MagneticStructure}","page":"Home","title":"SpinWaveTheory.LSWT","text":"LSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, term::Term, magneticstructure::MagneticStructure; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)\nLSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)\n\nConstruct a LSWT.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.MagneticStructure","page":"Home","title":"SpinWaveTheory.MagneticStructure","text":"MagneticStructure{L<:AbstractLattice, D<:Number}\n\nThe magnetic structure of an ordered quantum lattice system.\n\n\n\n\n\n","category":"type"},{"location":"#SpinWaveTheory.MagneticStructure-Tuple{QuantumLattices.Spatials.AbstractLattice, Dict{Int64, <:AbstractVector}}","page":"Home","title":"SpinWaveTheory.MagneticStructure","text":"MagneticStructure(cell::AbstractLattice, moments::Dict{Int, <:AbstractVector})\n\nConstruct the magnetic structure on a given lattice with the given moments.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.Magnonic","page":"Home","title":"SpinWaveTheory.Magnonic","text":"Magnonic <: TBAKind{:BdG}\n\nMagnonic quantum lattice system.\n\n\n\n\n\n","category":"type"},{"location":"#QuantumLattices.add!-Tuple{Matrix, TightBindingApproximation.TBAMatrixRepresentation{Magnonic}, QuantumLattices.QuantumOperators.Operator{<:Number, <:Tuple{QuantumLattices.DegreesOfFreedom.CompositeIndex{<:QuantumLattices.DegreesOfFreedom.Index{Int64, <:QuantumLattices.QuantumSystems.FID{:b}}}, QuantumLattices.DegreesOfFreedom.CompositeIndex{<:QuantumLattices.DegreesOfFreedom.Index{Int64, <:QuantumLattices.QuantumSystems.FID{:b}}}}}}","page":"Home","title":"QuantumLattices.add!","text":"add!(dest::Matrix, mr::TBAMatrixRepresentation{Magnonic}, m::Operator{<:Number, <:ID{CompositeIndex{<:Index{Int, <:FID{:b}}}, 2}}; atol=atol/5, kwargs...) -> typeof(dest)\n\nGet the matrix representation of an operator and add it to destination.\n\n\n\n\n\n","category":"method"},{"location":"#SpinWaveTheory.rotation-Tuple{AbstractVector{<:Number}}","page":"Home","title":"SpinWaveTheory.rotation","text":"rotation(destination::AbstractVector{<:Number}) -> SMatrix{3, 3}\n\nGet the rotation matrix which rotates the [0, 0, 1] or [θ, ϕ] vector to the direction of the destination vector.\n\n\n\n\n\n","category":"method"},{"location":"#TightBindingApproximation.commutator-Tuple{Magnonic, QuantumLattices.DegreesOfFreedom.Hilbert{<:QuantumLattices.QuantumSystems.Fock{:b}}}","page":"Home","title":"TightBindingApproximation.commutator","text":"commutator(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) -> Diagonal\n\nGet the commutation relation of the Holstein-Primakoff bosons.\n\n\n\n\n\n","category":"method"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"CurrentModule = SpinWaveTheory","category":"page"},{"location":"examples/CrBr3/#Ferromagnetic-spin-waves-in-the-monolayer-CrBr","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"","category":"section"},{"location":"examples/CrBr3/#Magnon-bands-and-inelastic-neutron-spectra-by-linear-spin-wave-theory","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Magnon bands and inelastic neutron spectra by linear spin wave theory","text":"","category":"section"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"The following codes could compute the ferromagnetic spin wave dispersions of the monolayer CrBr₃ by the K-Γ model:","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"note: Note\nK-Γ model is not the correct model to describe CrBr₃. See our paper. Here it is just for illustration on how to use this package.","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"using QuantumLattices\nusing TightBindingApproximation\nusing SpinWaveTheory\nusing Plots\n\nfunction kitaev(bond)\n    ϕ = bond|>rcoordinate|>azimuthd\n    any(≈(ϕ), ( 90, 270)) && return Coupling(:, SID, ('z', 'z'))\n    any(≈(ϕ), ( 30, 210)) && return Coupling(:, SID, ('x', 'x'))\n    any(≈(ϕ), (150, 330)) && return Coupling(:, SID, ('y', 'y'))\n    error(\"kitaev error: wrong azimuth angle($ϕ) bond.\")\nend\nfunction gamma(bond)\n    ϕ = bond|>rcoordinate|>azimuthd\n    any(≈(ϕ), ( 90, 270)) && return MatrixCoupling(:, SID, Γ\"z\")\n    any(≈(ϕ), ( 30, 210)) && return MatrixCoupling(:, SID, Γ\"x\")\n    any(≈(ϕ), (150, 330)) && return MatrixCoupling(:, SID, Γ\"y\")\n    error(\"kitaev error: wrong azimuth angle($ϕ) of bond.\")\nend\n\nlattice = Lattice(\n    (0.0, 0.0), (0.0, √3/3);\n    name=:H2,\n    vectors=[[1.0, 0.0], [0.5, √3/2]]\n    )\nhilbert = Hilbert(site=>Spin{3//2}() for site=1:length(lattice))\n\nJ₁ = SpinTerm(:J₁, 0.0, 1, MatrixCoupling(:, SID, Heisenberg\"\"))\nJ₂ = SpinTerm(:J₂, -0.178, 2, MatrixCoupling(:, SID, Heisenberg\"\"))\nJ₃ = SpinTerm(:J₃, 0.051, 3, MatrixCoupling(:, SID, Heisenberg\"\"))\nK = SpinTerm(:K, -4.288, 1, kitaev)\nΓ = SpinTerm(:Γ, -0.044, 1, gamma)\n\nmagneticstructure = MagneticStructure(lattice, Dict(site=>[1, 1, 1] for site=1:length(lattice)))\nCrBr₃ = Algorithm(:CrBr₃, LSWT(lattice, hilbert, (J₁, J₂, J₃, K, Γ), magneticstructure))\npath = ReciprocalPath(reciprocals(lattice), (-2, -1)=>(2, 1), length=400)\n\nspectra = CrBr₃(\n    :INSS,\n    InelasticNeutronScatteringSpectra(path, range(0.0, 15.0, length=301); fwhm=1.0, scale=log)\n    )\nenergybands = CrBr₃(:EB, EnergyBands(path))\nplt = plot()\nplot!(plt, spectra)\nplot!(plt, energybands, color=:white, linestyle=:dash)\nyticks!(plt, range(0.0, 15.0, length=16))","category":"page"},{"location":"examples/CrBr3/#Berry-curvature-and-Chern-number-of-the-magnon-bands","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Berry curvature and Chern number of the magnon bands","text":"","category":"section"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"The Berry curvatures and the Chern numbers of the magnon bands could be computed in the reciprocal unitcell:","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"brillouin = BrillouinZone(reciprocals(lattice), 90)\nberry = CrBr₃(:BerryCurvature, BerryCurvature(brillouin, [1, 2]));\nplot(berry)","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"Here, k₁ and k₂ denote the coordinates in the reciprocal space along the two reciprocal vectors.","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"The Berry curvatures can also be computed in a reciprocal zone beyond the reciprocal unitcell:","category":"page"},{"location":"examples/CrBr3/","page":"Ferromagnetic spin waves in the monolayer CrBr₃","title":"Ferromagnetic spin waves in the monolayer CrBr₃","text":"b₁, b₂ = 4*pi/√3*[1.0, 0.0], 4*pi/√3*[0.0, 1.0]\nreciprocalzone = ReciprocalZone(\n    [b₁, b₂], [-1.0=>1.0, -1.0=>1.0];\n    length=301, ends=(true, true)\n)\nberry = CrBr₃(:BerryCurvatureEx, BerryCurvature(reciprocalzone, [1, 2]))\nplot(berry)","category":"page"}]
}
