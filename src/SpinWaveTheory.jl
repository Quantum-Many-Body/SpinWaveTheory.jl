module SpinWaveTheory

using StaticArrays: SVector, SMatrix, @SMatrix
using LinearAlgebra: norm, dot, cross, Hermitian
using QuantumLattices: ID, Operator, Operators, AbstractUnitSubstitution, RankFilter,  Generator, SimplifiedGenerator
using QuantumLattices: AbstractPID, Lattice, Bonds, OID, Index, SimpleIID, SID, FID, Metric, OIDToTuple, Table, Hilbert, Spin, Fock, Term, Boundary
using QuantumLattices: atol, rtol, dtype, pidtype, indextype, fulltype, idtype, reparameter, add!, sub!, mul!, plain, expand, rcoord
using TightBindingApproximation: TBAKind, AbstractTBA, TBAMatrix

import QuantumLattices: optype, contentnames, update!, matrix!, statistics, dimension
import TightBindingApproximation: commutator

export rotation, MagneticStructure, HPTransformation, LSWT

"""
    rotation(destination::AbstractVector{<:Number}) -> SMatrix{3, 3}

Get the rotation matrix which rotates the `[0, 0, 1]` vector to the direction of the `destination` vector.
"""
function rotation(destination::AbstractVector{<:Number})
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

"""
    MagneticStructure{L<:Lattice, P<:AbstractPID, D<:Number}

The magnetic structure of an ordered quantum lattice system.
"""
struct MagneticStructure{L<:Lattice, P<:AbstractPID, D<:Number}
    cell::L
    moments::Dict{P, SVector{3, D}}
    rotations::Dict{P, SMatrix{3, 3, D, 9}}
end

"""
    MagneticStructure(cell::Lattice, moments::Dict{<:AbstractPID, <:AbstractVector})

Construct the magnetic structure on a given lattice with the given moments.
"""
function MagneticStructure(cell::Lattice, moments::Dict{<:AbstractPID, <:AbstractVector})
    @assert length(cell)==length(moments) "MagneticStructure error: mismatched magnetic cell and moments."
    datatype = promote_type(dtype(cell), eltype(valtype(moments)))
    rotations = Dict{keytype(moments), SMatrix{3, 3, datatype, 9}}()
    for pid in cell.pids
        rotations[pid] = rotation(moments[pid])
    end
    return MagneticStructure(cell, convert(Dict{keytype(moments), SVector{3, datatype}}, moments), rotations)
end

"""
    HPTransformation{S<:Operators, U<:OID{<:Index{<:AbstractPID, <:SimpleIID}}, M<:MagneticStructure} <: AbstractUnitSubstitution{U, S}

Holstein-Primakoff transformation.
"""
struct HPTransformation{S<:Operators, U<:OID{<:Index{<:AbstractPID, <:SimpleIID}}, M<:MagneticStructure} <: AbstractUnitSubstitution{U, S}
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
@inline @generated function optype(::Type{<:HPTransformation}, ::Type{S}) where {S<:Operators}
    V = promote_type(valtype(eltype(S)), Complex{Int})
    Iₒ = indextype(eltype(eltype(S)))
    Iₜ = reparameter(Iₒ, :iid, FID{:b, Int, Int, Int})
    I = Iₜ<:Iₒ ? Iₒ : Iₜ
    U = reparameter(eltype(eltype(S)), :index, I)
    M = fulltype(eltype(S), NamedTuple{(:value, :id), Tuple{V, ID{U}}})
end
function (hp::HPTransformation)(oid::OID{<:Index{<:AbstractPID, <:SID{S}}}) where S
    datatype = valtype(eltype(valtype(hp)))
    factor = √(2*S*one(datatype))/2
    op₁ = Operator(1, replace(oid, index=replace(oid.index, iid=FID{:b}(oid.index.iid.orbital, 1, 1))))
    op₂ = Operator(1, replace(oid, index=replace(oid.index, iid=FID{:b}(oid.index.iid.orbital, 1, 2))))
    xₗ = add!(add!(zero(valtype(hp)), factor*op₁), factor*op₂)
    yₗ = sub!(add!(zero(valtype(hp)), factor*op₁/1im), factor*op₂/1im)
    zₗ = sub!(add!(zero(valtype(hp)), S), op₂*op₁)
    x, y, z = hp.magneticstructure.rotations[oid.index.pid]*SVector(xₗ, yₗ, zₗ)
    oid.index.iid.tag=='x' && return x
    oid.index.iid.tag=='y' && return y
    oid.index.iid.tag=='z' && return z
    oid.index.iid.tag=='+' && return add!(x, mul!(y, 1im))
    return sub!(x, mul!(y, 1im))
end

"""
    Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure)

Get the corresponding Hilbert space of the original one after the Holstein-Primakoff transformation. 
"""
@inline function Hilbert(hilbert::Hilbert{<:Spin}, magneticstructure::MagneticStructure)
    return Hilbert(pid=>Fock{:b}(norbital=hilbert[pid].norbital, nspin=1, nnambu=2) for pid in magneticstructure.cell.pids)
end

"""
    LSWT{L<:Lattice, H<:Generator, HP<:HPTransformation, E<:SimplifiedGenerator, F<:SimplifiedGenerator, G<:AbstractMatrix} <: AbstractTBA{TBAKind(:BdG), H, G}

Linear spin wave theory for magnetically ordered quantum lattice systems.
"""
struct LSWT{L<:Lattice, H<:Generator, HP<:HPTransformation, E<:SimplifiedGenerator, F<:SimplifiedGenerator, G<:AbstractMatrix} <: AbstractTBA{TBAKind(:BdG), H, G}
    lattice::L
    H::H
    hp::HP
    H₀::E
    H₂::F
    commutator::G
    function LSWT(lattice::Lattice, H::Generator, hp::HPTransformation)
        temp = hp(H)
        hilbert = Hilbert(H.hilbert, hp.magneticstructure)
        table = Table(hilbert, Metric(TBAKind(:BdG), hilbert))
        H₀ = RankFilter(0)(temp, table=table)
        H₂ = RankFilter(2)(temp, table=table)
        commt = commutator(TBAKind(:BdG), hilbert)
        new{typeof(lattice), typeof(H), typeof(hp), typeof(H₀), typeof(H₂), typeof(commt)}(lattice, H, hp, H₀, H₂, commt)
    end
end
@inline contentnames(::Type{<:LSWT}) = (:lattice, :H, :hp, :H₀, :H₂, :commutator)
@inline Base.eltype(T::Type{<:LSWT}) = eltype(fieldtype(T, :H₂))
@inline dimension(lswt::LSWT) = length(lswt.H₂.table)
@inline statistics(::Type{<:LSWT}) = :b
@inline function update!(lswt::LSWT; kwargs...)
    update!(lswt.H; kwargs...)
    update!(lswt.H₀; kwargs...)
    update!(lswt.H₂; kwargs...)
    return lswt
end

"""
    LSWT(lattice::Lattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; boundary::Boundary=plain)

Construct a LSWT.
"""
@inline function LSWT(lattice::Lattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; boundary::Boundary=plain)
    H = Generator(terms, Bonds(magneticstructure.cell), hilbert; half=false, boundary=boundary)
    hp = HPTransformation{valtype(H)}(magneticstructure)
    return LSWT(lattice, H, hp)
end

"""
    matrix!(lswt::LSWT; k=nothing, atol=atol/5, kwargs...) -> TBAMatrix

Get the tight-binding matrix representation of the linear spin waves.

Here, the `atol` parameter is used to ensure that the matrix is positive-definite so that the Cholesky decomposition can be performed numerically.
"""
function matrix!(lswt::LSWT; k=nothing, atol=atol/5, kwargs...)
    length(kwargs)>0 && update!(lswt; kwargs...)
    table = lswt.H₂.table
    result = zeros(valtype(lswt, k), dimension(lswt), dimension(lswt))
    for op in expand(lswt.H₂)
        if op[1]==op[2]'
            seq₁ = table[op[1].index]
            seq₂ = table[op[2].index]
            result[seq₁, seq₁] += op.value+atol
            result[seq₂, seq₂] += op.value+atol
        else
            phase = isnothing(k) ? one(valtype(lswt, k)) : convert(valtype(lswt, k), exp(-1im*dot(k, rcoord(op))))
            seq₁ = table[op[1].index']
            seq₂ = table[op[2].index]
            result[seq₁, seq₂] += op.value*phase
            result[seq₂, seq₁] += op.value'*phase'
        end
    end
    return TBAMatrix(Hermitian(result), lswt.commutator)
end

end