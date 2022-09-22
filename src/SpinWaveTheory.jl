module SpinWaveTheory

using LinearAlgebra: Diagonal, Hermitian, dot, eigen, norm
using QuantumLattices: dtype, expand, mul!, sub!
using QuantumLattices: plain, Boundary, CompositeIndex, Hilbert, Index, OperatorUnitToTuple, Table, Term, indextype
using QuantumLattices: Action, Algorithm, Assignment, Image, OperatorGenerator
using QuantumLattices: AbstractUnitSubstitution, ID, Operator, Operators, RankFilter, idtype
using QuantumLattices: FID, Fock, SID, Spin
using QuantumLattices: AbstractLattice, Neighbors, bonds, icoordinate, rcoordinate
using QuantumLattices: atol, rtol, delta, fulltype, reparameter
using StaticArrays: SVector, SMatrix, @SMatrix
using TightBindingApproximation: AbstractTBA, InelasticNeutronScatteringSpectra, TBAKind, TBAMatrix, TBAMatrixRepresentation
using TimerOutputs: @timeit

import QuantumLattices: add!, dimension, matrix, matrix!, update!
import QuantumLattices: Metric
import QuantumLattices: run!
import QuantumLattices: optype
import QuantumLattices: contentnames
import TightBindingApproximation: commutator

export HPTransformation, LSWT, MagneticStructure, Magnonic, rotation

"""
    rotation(destination::AbstractVector{<:Number}) -> SMatrix{3, 3}

Get the rotation matrix which rotates the `[0, 0, 1]` or `[θ, ϕ]` vector to the direction of the `destination` vector.
"""
function rotation(destination::AbstractVector{<:Number})
    @assert length(destination)==3 || length(destination)==2 "rotation error: destination vector must be 3(or 2)-dimensional."
    if length(destination) == 3
        x, y, z = destination[1], destination[2], destination[3]
        lenxy, lenxyz = √(x^2+y^2), √(x^2+y^2+z^2)
        cosθ, sinθ = z/lenxyz, lenxy/lenxyz
        cosφ, sinφ = if isapprox(lenxy, 0, atol=atol, rtol=rtol)
            one(eltype(destination)), zero(eltype(destination))
        else
            x/lenxy, y/lenxy
        end
    else
        θ₁, ϕ₁ = destination[1], destination[2]
        cosθ, sinθ = cos(θ₁), sin(θ₁)
        cosφ, sinφ = cos(ϕ₁), sin(ϕ₁)
    end
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
    MagneticStructure(cell::AbstractLattice, moments::Dict{Int, <:AbstractVector})

Construct the magnetic structure on a given lattice with the given moments.
"""
function MagneticStructure(cell::AbstractLattice, moments::Dict{Int, <:AbstractVector})
    @assert length(cell)==length(moments) "MagneticStructure error: mismatched magnetic cell and moments."
    datatype = promote_type(dtype(cell), eltype(valtype(moments)))
    moments = convert(Dict{Int, Vector{datatype}}, moments)
    rotations = Dict{Int, SMatrix{3, 3, datatype, 9}}()
    for site=1:length(cell)
        rotations[site] = rotation(moments[site])
    end
    return MagneticStructure(cell, moments, rotations)
end

"""
    HPTransformation{S<:Operators, U<:CompositeIndex, M<:MagneticStructure} <: AbstractUnitSubstitution{U, S}

Holstein-Primakoff transformation.
"""
struct HPTransformation{S<:Operators, U<:CompositeIndex, M<:MagneticStructure} <: AbstractUnitSubstitution{U, S}
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
    Iₜ = reparameter(Iₒ, :iid, FID{:b, Int, Int, Int})
    I = Iₜ<:Iₒ ? Iₒ : Iₜ
    U = reparameter(eltype(eltype(S)), :index, I)
    return fulltype(eltype(S), NamedTuple{(:value, :id), Tuple{V, ID{U}}})
end
@inline (hp::HPTransformation)(index::CompositeIndex; kwargs...) = Operator(1, index)
function (hp::HPTransformation)(index::CompositeIndex{<:Index{Int, <:SID{S}}}; zoff::Bool=false) where S
    datatype = valtype(eltype(valtype(hp)))
    factor = √(2*S*one(datatype))/2
    op₁ = Operator(1, replace(index, index=replace(index.index, iid=FID{:b}(1, 1, 1))))
    op₂ = Operator(1, replace(index, index=replace(index.index, iid=FID{:b}(1, 1, 2))))
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
@inline commutator(::Magnonic, hilbert::Hilbert{<:Fock{:b}}) = Diagonal(kron([1, -1], ones(Int64, sum(dimension, values(hilbert))÷2)))


# When the type of the field `commutator` is a type parameter of `LSWT`, a strange bug would occur when the terms of a quantum
# spin system include both exchange interactions (e.g. the Heisenberg term) and linear onsite interactions (e.g. the potential
# energy under an external magnetic field). This bug does not exist for QuantumLattices@v0.8.6 and SpinWaveTheory@v0.1.1, but
# does for QuantumLattices@v0.8.9 and SpinWaveTheory@v0.1.2. This bug can be fixed by removing the type of the field `commutator`
# as a type parameter of `LSWT` (After several tries, I found this solution although I still don't know why. I guess it may be
# caused by a Julia bug.). This may cause a little bit of type instability but it turns out to be unimportant because the total
# time of the code execution is at most of several seconds, which differs little compared to previous versions.
"""
    LSWT{K<:TBAKind{:BdG}, L<:AbstractLattice, H<:OperatorGenerator, HP<:HPTransformation, E<:Image, F<:Image} <: AbstractTBA{K, H, AbstractMatrix}

Linear spin wave theory for magnetically ordered quantum lattice systems.
"""
struct LSWT{K<:TBAKind{:BdG}, L<:AbstractLattice, H<:OperatorGenerator, HP<:HPTransformation, E<:Image, F<:Image} <: AbstractTBA{K, H, AbstractMatrix}
    lattice::L
    H::H
    hp::HP
    H₀::E
    H₂::F
    commutator::AbstractMatrix
    function LSWT{K}(lattice::AbstractLattice, H::OperatorGenerator, hp::HPTransformation) where {K<:TBAKind{:BdG}}
        temp = hp(H)
        hilbert = Hilbert(H.hilbert, hp.magneticstructure)
        table = Table(hilbert, Metric(K(), hilbert))
        H₀ = RankFilter(0)(temp, table=table)
        H₂ = RankFilter(2)(temp, table=table)
        commt = commutator(K(), hilbert)
        new{K, typeof(lattice), typeof(H), typeof(hp), typeof(H₀), typeof(H₂)}(lattice, H, hp, H₀, H₂, commt)
    end
end
@inline contentnames(::Type{<:LSWT}) = (:lattice, :H, :hp, :H₀, :H₂, :commutator)
@inline Base.valtype(T::Type{<:LSWT}) = valtype(eltype(fieldtype(T, :H₂)))
@inline dimension(lswt::LSWT) = length(lswt.H₂.table)
@inline function update!(lswt::LSWT; k=nothing, kwargs...)
    if length(kwargs)>0
        update!(lswt.H; kwargs...)
        update!(lswt.H₀; kwargs...)
        update!(lswt.H₂; kwargs...)
    end
    return lswt
end

"""
    LSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)

Construct a LSWT.
"""
@inline function LSWT(lattice::AbstractLattice, hilbert::Hilbert{<:Spin}, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)
    isnothing(neighbors) && (neighbors=maximum(term->term.bondkind, terms))
    H = OperatorGenerator(terms, bonds(magneticstructure.cell, neighbors), hilbert; half=false, boundary=boundary)
    hp = HPTransformation{valtype(H)}(magneticstructure)
    return LSWT{Magnonic}(lattice, H, hp)
end

"""
    matrix(lswt::LSWT; k=nothing, gauge=:icoordinate, atol=atol/5, kwargs...) -> TBAMatrix

Get the tight-binding matrix representation of the linear spin waves.

Here, the `atol` parameter is used to ensure that the matrix is positive-definite so that the Cholesky decomposition can be performed numerically.
"""
@inline function matrix(lswt::LSWT; k=nothing, gauge=:icoordinate, atol=atol/5, kwargs...)
    return TBAMatrix(Hermitian(TBAMatrixRepresentation(lswt, k, gauge)(expand(lswt.H₂); atol=atol, kwargs...)), lswt.commutator)
end

"""
    TBAMatrixRepresentation(lswt::LSWT, k=nothing, gauge::Symbol=:icoordinate)

Construct the matrix representation transformation of a quantum spin system using the linear spin wave theory.
"""
@inline function TBAMatrixRepresentation(lswt::LSWT, k=nothing, gauge::Symbol=:icoordinate)
    return TBAMatrixRepresentation{typeof(lswt)}(k, lswt.H₂.table, gauge)
end

"""
    add!(dest::Matrix, mr::TBAMatrixRepresentation{<:LSWT}, m::Operator{<:Number, <:ID{CompositeIndex{<:Index{Int, <:FID{:b}}}}}; atol=atol/5, kwargs...) -> typeof(dest)

Get the matrix representation of an operator and add it to destination.
"""
function add!(dest::Matrix, mr::TBAMatrixRepresentation{<:LSWT}, m::Operator{<:Number, <:ID{CompositeIndex{<:Index{Int, <:FID{:b}}}}}; atol=atol/5, kwargs...)
    if m[1]==m[2]'
        seq₁ = mr.table[m[1].index]
        seq₂ = mr.table[m[2].index]
        dest[seq₁, seq₁] += m.value+atol
        dest[seq₂, seq₂] += m.value+atol
    else
        coordinate = mr.gauge==:rcoordinate ? rcoordinate(m) : icoordinate(m)
        phase = isnothing(mr.k) ? one(eltype(dest)) : convert(eltype(dest), exp(1im*dot(mr.k, coordinate)))
        seq₁ = mr.table[m[1].index']
        seq₂ = mr.table[m[2].index]
        dest[seq₁, seq₂] += m.value*phase
        dest[seq₂, seq₁] += m.value'*phase'
    end
    return dest
end

# Inelastic neutron scattering spectra of magnetically ordered local spin systems.
function run!(lswt::Algorithm{<:LSWT{Magnonic}}, ins::Assignment{<:InelasticNeutronScatteringSpectra})
    operators = spinoperators(lswt.frontend.H.hilbert, lswt.frontend.hp)
    m = zeros(promote_type(valtype(lswt.frontend), Complex{Int}), dimension(lswt.frontend), dimension(lswt.frontend))
    data = zeros(Complex{Float64}, size(ins.data[3]))
    η = get(ins.action.options, :η, 0.1)
    for (i, params) in enumerate(pairs(ins.action.path))
        update!(lswt; params...)
        @timeit lswt.timer "matrix" (mr = matrix(lswt.frontend; gauge=get(ins.action.options, :gauge, :icoordinate), params...))
        @timeit lswt.timer "eigen" ((eigenvalues, eigenvectors) = eigen(mr))
        @timeit lswt.timer "spectra" for α=1:3, β=1:3
            factor = delta(α, β) - ((norm(params.k)==0 || α>length(params.k) || β>length(params.k)) ? 0 : params.k[α]*params.k[β]/dot(params.k, params.k))
            if !isapprox(abs(factor), 0, atol=atol, rtol=rtol)
                matrix!(m, operators, α, β, lswt.frontend.H₂.table, params.k)
                diag = Diagonal(eigenvectors'*m*eigenvectors)
                for (nₑ, e) in enumerate(ins.action.energies)
                    for j = (dimension(lswt.frontend)÷2+1):dimension(lswt.frontend)
                        data[nₑ, i] += factor*diag[j, j]*η^2/(η^2+(e-eigenvalues[j])^2)/pi
                    end
                end
            end
        end
    end
    isapprox(norm(imag(data)), 0, atol=atol, rtol=rtol) || @warn "run! warning: not negligible imaginary part ($(norm(imag(data))))."
    ins.data[3][:, :] .= real(data)[:, :]
    get(ins.action.options, :log, true) && (ins.data[3][:, :] = log.(ins.data[3].+1))
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
