module SpinWaveTheory

using LinearAlgebra: Diagonal, Hermitian, dot, eigen, norm
using QuantumLattices: atol, lazy, plain, rtol
using QuantumLattices: AbstractLattice, Action, Algorithm, Assignment, CompositeIndex, FID, Fock, Hilbert, ID, Image, Index, Neighbors, Operator, OperatorGenerator, Operators, OperatorSum, OperatorUnitToTuple, RankFilter, SID, Spin, Table, Term, UnitSubstitution
using QuantumLattices: bonds, delta, dimension, direction, dtype, expand, fulltype, icoordinate, idtype, indextype, kind, mul!, rcoordinate, reparameter, sub!
using StaticArrays: SVector, SMatrix, @SMatrix
using TightBindingApproximation: AbstractTBA, InelasticNeutronScatteringSpectra, Quadratic, Quadraticization, TBAKind
using TimerOutputs: @timeit_debug

import QuantumLattices: add!, matrix, update!
import QuantumLattices: Metric
import QuantumLattices: run!
import QuantumLattices: optype
import QuantumLattices: contentnames
import TightBindingApproximation: commutator

export HPTransformation, LSWT, MagneticStructure, Magnonic, rotation

"""
    rotation(destination::AbstractVector{<:Number}; kwargs...) -> SMatrix{3, 3}
    rotation(destination::Tuple{Number, Number}; unit::Symbol=:radian) -> SMatrix{3, 3}

Get the rotation matrix which rotates `[0, 0, 1]` to the direction of the `destination` vector.
"""
function rotation(destination::AbstractVector{<:Number}; kwargs...)
    @assert length(destination)==3 "rotation error: destination vector must be 3-dimensional."
    x, y, z = destination[1], destination[2], destination[3]
    lenxy, lenxyz = √(x^2+y^2), √(x^2+y^2+z^2)
    cosθ, sinθ = z/lenxyz, lenxy/lenxyz
    cosφ, sinφ = if isapprox(lenxy, 0, atol=atol, rtol=rtol)
        one(eltype(destination)), zero(eltype(destination))
    else
        x/lenxy, y/lenxy
    end
    return @SMatrix [cosθ*cosφ -sinφ sinθ*cosφ;
                     cosθ*sinφ  cosφ sinθ*sinφ;
                         -sinθ     0      cosθ]
end
function rotation(destination::Tuple{Number, Number}; unit::Symbol=:radian)
    @assert unit∈(:degree, :radian) "rotation error: unit must be either `:degree` or `:radian`."
    θ₁, ϕ₁ = destination[1], destination[2]
    cosθ, sinθ = unit==:radian ? (cos(θ₁), sin(θ₁)) : (cosd(θ₁), sind(θ₁))
    cosφ, sinφ = unit==:radian ? (cos(ϕ₁), sin(ϕ₁)) : (cosd(ϕ₁), sind(ϕ₁))
    return @SMatrix [cosθ*cosφ -sinφ sinθ*cosφ;
                     cosθ*sinφ  cosφ sinθ*sinφ;
                         -sinθ     0      cosθ]
end

"""
    MagneticStructure{L<:AbstractLattice, D<:Number}

The magnetic structure of an ordered quantum lattice system.
"""
struct MagneticStructure{L<:AbstractLattice, D<:Number}
    cell::L
    moments::Dict{Int, Vector{D}}
    rotations::Dict{Int, SMatrix{3, 3, D, 9}}
end

"""
    MagneticStructure(cell::AbstractLattice, moments::Dict{Int, <:Union{AbstractVector, NTuple{2, Number}}}; unit::Symbol=:radian)

Construct the magnetic structure on a given lattice with the given moments.
"""
function MagneticStructure(cell::AbstractLattice, moments::Dict{Int, <:Union{AbstractVector, NTuple{2, Number}}}; unit::Symbol=:radian)
    @assert length(cell)==length(moments) "MagneticStructure error: mismatched magnetic cell and moments."
    datatype = promote_type(dtype(cell), eltype(valtype(moments)))
    new = Dict{Int, Vector{datatype}}()
    rotations = Dict{Int, SMatrix{3, 3, datatype, 9}}()
    for site=1:length(cell)
        destination = map(i->convert(datatype, i), moments[site])
        new[site] = direction(destination, unit)
        rotations[site] = rotation(destination; unit=unit)
    end
    return MagneticStructure(cell, new, rotations)
end

"""
    HPTransformation{S<:Operators, U<:CompositeIndex, M<:MagneticStructure} <: UnitSubstitution{U, S}

Holstein-Primakoff transformation.
"""
struct HPTransformation{S<:Operators, U<:CompositeIndex, M<:MagneticStructure} <: UnitSubstitution{U, S}
    magneticstructure::M
    function HPTransformation{S}(magneticstructure::MagneticStructure) where {S<:Operators}
        O = optype(HPTransformation, S)
        new{Operators{O, idtype(O)}, eltype(eltype(S)), typeof(magneticstructure)}(magneticstructure)
    end
end
@inline Base.valtype(hp::HPTransformation) = valtype(typeof(hp))
@inline Base.valtype(::Type{<:HPTransformation{S}}) where {S<:Operators} = S
@inline Base.valtype(::Type{<:HPTransformation{S}}, ::Type{<:Operator}) where {S<:Operators} = S
@inline Base.valtype(::Type{<:HPTransformation{S}}, ::Type{<:Operators}) where {S<:Operators} = S
@inline function optype(::Type{<:HPTransformation}, ::Type{S}) where {S<:Operators}
    V = promote_type(valtype(eltype(S)), Complex{Int})
    Iₒ = indextype(eltype(eltype(S)))
    Iₜ = reparameter(Iₒ, :iid, FID{:b, Int, Rational{Int}, Int})
    I = Iₜ<:Iₒ ? Iₒ : Iₜ
    U = reparameter(eltype(eltype(S)), :index, I)
    return fulltype(eltype(S), NamedTuple{(:value, :id), Tuple{V, ID{U}}})
end
@inline (hp::HPTransformation)(index::CompositeIndex; kwargs...) = Operator(1, index)
function (hp::HPTransformation)(index::CompositeIndex{<:Index{Int, <:SID{S}}}; zoff::Bool=false) where S
    datatype = valtype(eltype(valtype(hp)))
    factor = √(2*S*one(datatype))/2
    op₁ = Operator(1, replace(index, index=replace(index.index, iid=FID{:b}(1, 0, 1))))
    op₂ = Operator(1, replace(index, index=replace(index.index, iid=FID{:b}(1, 0, 2))))
    xₗ = add!(add!(zero(valtype(hp)), factor*op₁), factor*op₂)
    yₗ = sub!(add!(zero(valtype(hp)), factor*op₁/1im), factor*op₂/1im)
    zₗ = zero(valtype(hp))
    zoff || sub!(add!(zₗ, S), op₂*op₁)
    x, y, z = hp.magneticstructure.rotations[index.index.site]*SVector(xₗ, yₗ, zₗ)
    index.index.iid.tag=='x' && return x
    index.index.iid.tag=='y' && return y
    index.index.iid.tag=='z' && return z
    index.index.iid.tag=='+' && return add!(x, mul!(y, 1im))
    return sub!(x, mul!(y, 1im))
end

"""
    Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure)

Get the corresponding Hilbert space of the original one after the Holstein-Primakoff transformation. 
"""
@inline Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure) = Hilbert(site=>Fock{:b}(1, 1) for site=1:length(magneticstructure.cell))

"""
    Magnonic <: TBAKind{:BdG}

Magnonic quantum lattice system.
"""
struct Magnonic <: TBAKind{:BdG} end

"""
    Metric(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) -> OperatorUnitToTuple

Get the index-to-tuple metric for a quantum spin system after the Holstein-Primakoff transformation.
"""
@inline @generated Metric(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) = OperatorUnitToTuple(:nambu, :site)

"""
    commutator(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) -> Diagonal

Get the commutation relation of the Holstein-Primakoff bosons.
"""
@inline commutator(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) = Diagonal(kron([1, -1], ones(Int64, sum(length, values(hilbert))÷2)))

"""
    add!(dest::OperatorSum, qf::Quadraticization{Magnonic}, m::Operator{<:Number, <:ID{CompositeIndex{<:Index{Int, <:FID{:b}}}, 2}}; kwargs...) -> typeof(dest)

Get the unified quadratic form of a rank-2 operator and add it to `destination`.
"""
function add!(dest::OperatorSum, qf::Quadraticization{Magnonic}, m::Operator{<:Number, <:ID{CompositeIndex{<:Index{Int, <:FID{:b}}}, 2}}; kwargs...)
    rcoord, icoord = rcoordinate(m), icoordinate(m)
    if m[1]==m[2]'
        seq₁, seq₂ = qf.table[m[1]], qf.table[m[2]]
        add!(dest, Quadratic(m.value, (seq₁, seq₁), rcoord, icoord))
        add!(dest, Quadratic(m.value, (seq₂, seq₂), -rcoord, -icoord))
    else
        seq₁, seq₂ = qf.table[m[1]'], qf.table[m[2]]
        add!(dest, Quadratic(m.value, (seq₁, seq₂), rcoord, icoord))
        add!(dest, Quadratic(m.value', (seq₂, seq₁), -rcoord, -icoord))
    end
    return dest
end

"""
    LSWT{K<:TBAKind{:BdG}, L<:AbstractLattice, Hₛ<:OperatorGenerator, HP<:HPTransformation, Ω<:Image, H<:Image, Hₘ<:Image, C<:AbstractMatrix} <: AbstractTBA{K, Hₘ, C}

Linear spin wave theory for magnetically ordered quantum lattice systems.
"""
struct LSWT{K<:TBAKind{:BdG}, L<:AbstractLattice, Hₛ<:OperatorGenerator, HP<:HPTransformation, Ω<:Image, H<:Image, Hₘ<:Image, C<:AbstractMatrix} <: AbstractTBA{K, Hₘ, C}
    lattice::L
    Hₛ::Hₛ
    hp::HP
    Ω::Ω
    H::H
    Hₘ::Hₘ
    commutator::C
    function LSWT{K}(lattice::AbstractLattice, Hₛ::OperatorGenerator, hp::HPTransformation) where {K<:TBAKind{:BdG}}
        temp = hp(Hₛ)
        Ω = RankFilter(0)(temp)
        H = RankFilter(2)(temp)
        hilbert = Hilbert(Hₛ.hilbert, hp.magneticstructure)
        Hₘ = Quadraticization{K}(Table(hilbert, Metric(K(), hilbert)))(H)
        commt = commutator(K(), hilbert)
        new{K, typeof(lattice), typeof(Hₛ), typeof(hp), typeof(Ω), typeof(H), typeof(Hₘ), typeof(commt)}(lattice, Hₛ, hp, Ω, H, Hₘ, commt)
    end
end
@inline contentnames(::Type{<:LSWT}) = (:lattice, :Hₛ, :hp, :Ω, :H, :commutator)
@inline function update!(lswt::LSWT; k=nothing, kwargs...)
    if length(kwargs)>0
        update!(lswt.Hₛ; kwargs...)
        update!(lswt.Ω; kwargs...)
        update!(lswt.H; kwargs...)
        update!(lswt.Hₘ; kwargs...)
    end
    return lswt
end

"""
    LSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, terms::Union{Term, Tuple{Term, Vararg{Term}}}, magneticstructure::MagneticStructure; neighbors::Union{Nothing, Int, Neighbors}=nothing)

Construct a LSWT.
"""
@inline function LSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, terms::Union{Term, Tuple{Term, Vararg{Term}}}, magneticstructure::MagneticStructure; neighbors::Union{Nothing, Int, Neighbors}=nothing)
    terms = wrapper(terms)
    isnothing(neighbors) && (neighbors=maximum(term->term.bondkind, terms))
    Hₛ = OperatorGenerator(terms, bonds(magneticstructure.cell, neighbors), hilbert, plain, lazy; half=false)
    hp = HPTransformation{valtype(Hₛ)}(magneticstructure)
    return LSWT{Magnonic}(lattice, Hₛ, hp)
end
@inline wrapper(x) = (x,)
@inline wrapper(xs::Tuple) = xs

# Inelastic neutron scattering spectra of magnetically ordered local spin systems.
function run!(lswt::Algorithm{<:LSWT{Magnonic}}, inss::Assignment{<:InelasticNeutronScatteringSpectra})
    operators = spinoperators(lswt.frontend.Hₛ.hilbert, lswt.frontend.hp)
    m = zeros(promote_type(valtype(lswt.frontend), Complex{Int}), dimension(lswt.frontend), dimension(lswt.frontend))
    σ = get(inss.action.options, :fwhm, 0.1)/2/√(2*log(2))
    data = zeros(Complex{Float64}, size(inss.data[3]))
    for (i, momentum) in enumerate(inss.action.reciprocalspace)
        eigenvalues, eigenvectors = eigen(lswt; k=momentum, inss.action.options...)
        @timeit_debug lswt.timer "spectra" for α=1:3, β=1:3
            factor = (delta(α, β) - ((norm(momentum)==0 || α>length(momentum) || β>length(momentum)) ? 0 : momentum[α]*momentum[β]/dot(momentum, momentum)))/√(2pi)/σ
            if !isapprox(abs(factor), 0, atol=atol, rtol=rtol)
                matrix!(m, operators, α, β, lswt.frontend.Hₘ.transformation.table, momentum)
                diag = Diagonal(eigenvectors'*m*eigenvectors)
                for (nₑ, e) in enumerate(inss.action.energies)
                    for j = (dimension(lswt.frontend)÷2+1):dimension(lswt.frontend)
                        data[nₑ, i] += factor*diag[j, j]*exp(-(e-eigenvalues[j])^2/2/σ^2)
                    end
                end
            end
        end
    end
    isapprox(norm(imag(data)), 0, atol=atol, rtol=rtol) || @warn "run! warning: not negligible imaginary part ($(norm(imag(data))))."
    inss.data[3][:, :] .= real(data)[:, :]
    inss.data[3][:, :] = get(inss.action.options, :scale, identity).(inss.data[3].+1)
end
function spinoperators(hilbert::Hilbert{<:Spin}, hp::HPTransformation{S, U}) where {S<:Operators, U<:CompositeIndex{<:Index{Int, <:SID}}}
    M = fulltype(Operator, NamedTuple{(:value, :id), Tuple{valtype(eltype(S)), Tuple{U}}})
    x, y, z = Operators{M}(), Operators{M}(), Operators{M}()
    for (site, spin) in pairs(hilbert)
        rcoordinate = hp.magneticstructure.cell[site]
        icoordinate = zero(rcoordinate)
        for sid in spin
            sid.tag=='x' && add!(x, Operator(1, CompositeIndex(Index(site, sid), rcoordinate, icoordinate)))
            sid.tag=='y' && add!(y, Operator(1, CompositeIndex(Index(site, sid), rcoordinate, icoordinate)))
            sid.tag=='z' && add!(z, Operator(1, CompositeIndex(Index(site, sid), rcoordinate, icoordinate)))
        end
    end
    x = hp(x, zoff=true)|>RankFilter(1)
    y = hp(y, zoff=true)|>RankFilter(1)
    z = hp(z, zoff=true)|>RankFilter(1)
    return @SMatrix [x*x x*y x*z; y*x y*y y*z; z*x z*y z*z]
end
function matrix!(m::Matrix{<:Number}, operators::SMatrix{3, 3, <:Operators, 9}, i::Int, j::Int, table::Table, k)
    m[:, :] .= zero(eltype(m))
    for op in operators[i, j]
        phase = convert(eltype(m), exp(1im*dot(k, rcoordinate(op))))
        seq₁ = table[op[1].index']
        seq₂ = table[op[2].index]
        m[seq₁, seq₂] += op.value*phase
    end
    return m
end

end
