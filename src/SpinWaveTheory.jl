module SpinWaveTheory

using TimerOutputs: @timeit
using StaticArrays: SVector, SMatrix, @SMatrix
using LinearAlgebra: Hermitian, Diagonal, dot, eigen, norm
using QuantumLattices: ID, CompositeOID, Operator, Operators, AbstractUnitSubstitution, RankFilter,  Generator, Image, Action, Algorithm, Assignment
using QuantumLattices: AbstractPID, Lattice, Bonds, OID, Index, SimpleIID, SID, FID, Metric, Table, Hilbert, Spin, Fock, Term, Boundary, ReciprocalPath, LatticeIndex
using QuantumLattices: atol, rtol, dtype, indextype, fulltype, idtype, reparameter, sub!, mul!, expand, plain, rcoord, icoord, delta, fulltype
using TightBindingApproximation: TBAKind, AbstractTBA, TBAMatrix, TBAMatrixRepresentation

import QuantumLattices: optype, contentnames, update!, matrix, matrix!, add!, dimension, prepare!, run!
import TightBindingApproximation: commutator

export rotation, MagneticStructure, HPTransformation
export LSWT, InelasticNeutronSpectra

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
    moments = convert(Dict{keytype(moments), SVector{3, datatype}}, moments)
    rotations = Dict{keytype(moments), SMatrix{3, 3, datatype, 9}}()
    for pid in cell.pids
        rotations[pid] = rotation(moments[pid])
    end
    return MagneticStructure(cell, moments, rotations)
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
@inline function optype(::Type{<:HPTransformation}, ::Type{S}) where {S<:Operators}
    V = promote_type(valtype(eltype(S)), Complex{Int})
    Iₒ = indextype(eltype(eltype(S)))
    Iₜ = reparameter(Iₒ, :iid, FID{:b, Int, Int, Int})
    I = Iₜ<:Iₒ ? Iₒ : Iₜ
    U = reparameter(eltype(eltype(S)), :index, I)
    M = fulltype(eltype(S), NamedTuple{(:value, :id), Tuple{V, ID{U}}})
end
@inline (hp::HPTransformation)(oid::OID; kwargs...) = Operator(1, oid)
function (hp::HPTransformation)(oid::OID{<:Index{<:AbstractPID, <:SID{S}}}; zoff::Bool=false) where S
    datatype = valtype(eltype(valtype(hp)))
    factor = √(2*S*one(datatype))/2
    op₁ = Operator(1, replace(oid, index=replace(oid.index, iid=FID{:b}(oid.index.iid.orbital, 1, 1))))
    op₂ = Operator(1, replace(oid, index=replace(oid.index, iid=FID{:b}(oid.index.iid.orbital, 1, 2))))
    xₗ = add!(add!(zero(valtype(hp)), factor*op₁), factor*op₂)
    yₗ = sub!(add!(zero(valtype(hp)), factor*op₁/1im), factor*op₂/1im)
    zₗ = zero(valtype(hp))
    zoff || sub!(add!(zₗ, S), op₂*op₁)
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

# When the type of the field `commutator` is a type parameter of `LSWT`, a strange bug would occur when the terms of a quantum
# spin system include both exchange interactions (e.g. the Heisenberg term) and linear onsite interactions (e.g. the potential
# energy under an external magnetic field). This bug does not exist for QuantumLattices@v0.8.6 and SpinWaveTheory@v0.1.1, but
# does for QuantumLattices@v0.8.9 and SpinWaveTheory@v0.1.2. This bug can be fixed by removing the type of the field `commutator`
# as a type parameter of `LSWT` (After several tries, I found this solution although I still don't know why. I guess it may be
# caused by a Julia bug.). This may cause a little bit of type instability but it turns out to be unimportant because the total
# time of the code execution is at most of several seconds, which differs little compared to previous versions.
"""
    LSWT{L<:Lattice, H<:Generator, HP<:HPTransformation, E<:Image, F<:Image} <: AbstractTBA{TBAKind(:BdG), H, AbstractMatrix}

Linear spin wave theory for magnetically ordered quantum lattice systems.
"""
struct LSWT{L<:Lattice, H<:Generator, HP<:HPTransformation, E<:Image, F<:Image} <: AbstractTBA{TBAKind(:BdG), H, AbstractMatrix}
    lattice::L
    H::H
    hp::HP
    H₀::E
    H₂::F
    commutator::AbstractMatrix
    function LSWT(lattice::Lattice, H::Generator, hp::HPTransformation)
        temp = hp(H)
        hilbert = Hilbert(H.hilbert, hp.magneticstructure)
        table = Table(hilbert, Metric(TBAKind(:BdG), hilbert))
        H₀ = RankFilter(0)(temp, table=table)
        H₂ = RankFilter(2)(temp, table=table)
        commt = commutator(TBAKind(:BdG), hilbert)
        new{typeof(lattice), typeof(H), typeof(hp), typeof(H₀), typeof(H₂)}(lattice, H, hp, H₀, H₂, commt)
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
    LSWT(lattice::Lattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; boundary::Boundary=plain)

Construct a LSWT.
"""
@inline function LSWT(lattice::Lattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, magneticstructure::MagneticStructure; boundary::Boundary=plain)
    H = Generator(terms, Bonds(magneticstructure.cell), hilbert; half=false, boundary=boundary)
    hp = HPTransformation{valtype(H)}(magneticstructure)
    return LSWT(lattice, H, hp)
end

"""
    matrix(lswt::LSWT; k=nothing, gauge=:icoord, atol=atol/5, kwargs...) -> TBAMatrix

Get the tight-binding matrix representation of the linear spin waves.

Here, the `atol` parameter is used to ensure that the matrix is positive-definite so that the Cholesky decomposition can be performed numerically.
"""
@inline function matrix(lswt::LSWT; k=nothing, gauge=:icoord, atol=atol/5, kwargs...)
    return TBAMatrix(Hermitian(TBAMatrixRepresentation(lswt, k, gauge)(expand(lswt.H₂); atol=atol, kwargs...)), lswt.commutator)
end

"""
    TBAMatrixRepresentation(lswt::LSWT, k=nothing, gauge::Symbol=:icoord)

Construct the matrix representation transformation of a quantum spin system using the linear spin wave theory.
"""
@inline function TBAMatrixRepresentation(lswt::LSWT, k=nothing, gauge::Symbol=:icoord)
    return TBAMatrixRepresentation{typeof(lswt)}(k, lswt.H₂.table, gauge)
end

"""
    add!(dest::Matrix,
        mr::TBAMatrixRepresentation{<:LSWT},
        m::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID{:b}}}}};
        atol=atol/5,
        kwargs...
        ) -> typeof(dest)

Get the matrix representation of an operator and add it to destination.
"""
function add!(dest::Matrix, mr::TBAMatrixRepresentation{<:LSWT}, m::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID{:b}}}}}; atol=atol/5, kwargs...)
    if m[1]==m[2]'
        seq₁ = mr.table[m[1].index]
        seq₂ = mr.table[m[2].index]
        dest[seq₁, seq₁] += m.value+atol
        dest[seq₂, seq₂] += m.value+atol
    else
        coord = mr.gauge==:rcoord ? rcoord(m) : icoord(m)
        phase = isnothing(mr.k) ? one(eltype(dest)) : convert(eltype(dest), exp(-1im*dot(mr.k, coord)))
        seq₁ = mr.table[m[1].index']
        seq₂ = mr.table[m[2].index]
        dest[seq₁, seq₂] += m.value*phase
        dest[seq₂, seq₁] += m.value'*phase'
    end
    return dest
end

"""
    InelasticNeutronSpectra{P<:ReciprocalPath, E<:AbstractVector} <: Action

Inelastic neutron spectra of magnetically ordered quantum lattice systems by linear spin wave theory.
"""
struct InelasticNeutronSpectra{P<:ReciprocalPath, E<:AbstractVector} <: Action
    path::P
    energies::E
    η::Float64
    log::Bool
    function InelasticNeutronSpectra(path::ReciprocalPath, energies::AbstractVector, η::Float64, log::Bool)
        @assert keys(path)==(:k,) "InelasticNeutronSpectra error: the name of the momenta in the path must be :k."
        new{typeof(path), typeof(energies)}(path, energies, η, log)
    end
end
@inline InelasticNeutronSpectra(path::ReciprocalPath, energies::AbstractVector; η::Float64=0.1, log::Bool=true) = InelasticNeutronSpectra(path, energies, η, log)
@inline function prepare!(ins::InelasticNeutronSpectra, lswt::LSWT)
    x = collect(Float64, 0:(length(ins.path)-1))
    y = collect(Float64, ins.energies)
    z = zeros(Float64, length(y), length(x))
    return (x, y, z)
end
function run!(lswt::Algorithm{<:LSWT}, ins::Assignment{<:InelasticNeutronSpectra})
    operators = spinoperators(lswt.engine.H.hilbert, lswt.engine.hp)
    m = zeros(promote_type(valtype(lswt.engine), Complex{Int}), dimension(lswt.engine), dimension(lswt.engine))
    data = zeros(Complex{Float64}, size(ins.data[3]))
    for (i, params) in enumerate(ins.action.path)
        update!(lswt; params...)
        @timeit lswt.timer "matrix" (mr = matrix(lswt.engine; params...))
        @timeit lswt.timer "eigen" ((eigenvalues, eigenvectors) = eigen(mr))
        @timeit lswt.timer "spectra" for α=1:3, β=1:3
            factor = delta(α, β) - ((norm(params.k)==0 || α>length(params.k) || β>length(params.k)) ? 0 : params.k[α]*params.k[β]/dot(params.k, params.k))
            if !isapprox(abs(factor), 0, atol=atol, rtol=rtol)
                matrix!(m, operators, α, β, lswt.engine.H₂.table, params.k)
                diag = Diagonal(eigenvectors'*m*eigenvectors)
                for (nₑ, e) in enumerate(ins.action.energies)
                    for j = (dimension(lswt.engine)÷2+1):dimension(lswt.engine)
                        data[nₑ, i] += factor*diag[j, j]*ins.action.η^2/(ins.action.η^2+(e-eigenvalues[j])^2)/pi
                    end
                end
            end
        end
    end
    isapprox(norm(imag(data)), 0, atol=atol, rtol=rtol) || @warn "run! warning: not negligible imaginary part ($(norm(imag(data))))."
    ins.data[3][:, :] .= real(data)[:, :]
    ins.action.log && (ins.data[3][:, :] = log.(ins.data[3].+1))
end
function spinoperators(hilbert::Hilbert{<:Spin}, hp::HPTransformation{S, U}) where {S<:Operators, U<:OID{<:Index{<:AbstractPID, <:SID}}}
    M = fulltype(Operator, NamedTuple{(:value, :id), Tuple{valtype(eltype(S)), Tuple{U}}})
    x, y, z = Operators{M}(), Operators{M}(), Operators{M}()
    for (pid, spin) in pairs(hilbert)
        rcoord = hp.magneticstructure.cell[LatticeIndex{'R'}(pid)] 
        icoord = hp.magneticstructure.cell[LatticeIndex{'I'}(pid)]
        for sid in spin
            sid.tag=='x' && add!(x, Operator(1, OID(Index(pid, sid), rcoord, icoord)))
            sid.tag=='y' && add!(y, Operator(1, OID(Index(pid, sid), rcoord, icoord)))
            sid.tag=='z' && add!(z, Operator(1, OID(Index(pid, sid), rcoord, icoord)))
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
        phase = convert(eltype(m), exp(1im*dot(k, rcoord(op))))
        seq₁ = table[op[1].index']
        seq₂ = table[op[2].index]
        m[seq₁, seq₂] += op.value*phase
    end
    return m
end

end
