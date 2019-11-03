module PiecewiseIncreasingRanges
export PiecewiseIncreasingRange, NoNearestSampleError, findnearest, resample,
PiecewiseUnitRange,PiecewiseIncreasingStepRange,IncreasingStepRange

using Compat, Base.Order
import Base: @propagate_inbounds

struct IncreasingStepRange{T, S} <: OrdinalRange{T, S}
    start :: T
    step :: S
    stop :: T

    function IncreasingStepRange{T,S}(start,step,stop) where {T,S,R<:OrdinalRange}
        step <= 0 && throw(DomainError(step,"step has to be positive"))
        new{T,S}(start,step,stop)
    end
end

IncreasingStepRange(start::T,step::S,stop::T) where {T<:Integer,S<:Integer} = IncreasingStepRange{T,S}(start,step,stop)
IncreasingStepRange(r::OrdinalRange{T,S}) where {T<:Integer,S<:Integer} = IncreasingStepRange{T,S}(r.start,r.step,r.stop)

Base.step(r::IncreasingStepRange) = r.step
Base.isempty(r::IncreasingStepRange) =
    (r.start != r.stop) & ((r.step > zero(r.step)) != (r.stop > r.start))

function Base.length(r::IncreasingStepRange{T}) where T<:Union{Int,UInt,Int64,UInt64,Int128,UInt128}
    isempty(r) && return zero(T)
    if r.step > 1
        return Base.checked_add(convert(T, div(unsigned(r.stop - r.start), r.step)), one(T))
    else
        return Base.checked_add(div(Base.checked_sub(r.stop, r.start), r.step), one(T))
    end
end
Base.firstindex(::IncreasingStepRange) = 1
@inline @propagate_inbounds function Base.getindex(r::IncreasingStepRange, 
    s::AbstractRange{<:Integer})

    @boundscheck checkbounds(r, s)
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    range(st, step=step(r)*step(s), length=length(s))
end

abstract type AbstractPiecewiseRange{T,R} <: AbstractVector{T} end
abstract type AbstractPiecewiseOrdinalRange{T,R} <: AbstractPiecewiseRange{T,R} end

constructrange(::Type{T}, start, step, stop) where {T<:UnitRange} = T(start, stop)
constructrange(::Type{T}, start, step, stop) where {T<:StepRange} = T(start, step, stop)

function combine_ranges(ranges::Vector{R}, firstrg::R, firstrgidx::Int) where {R<:AbstractRange}
    newranges = R[]
    offsets = Int[1]

    step(firstrg) < zero(eltype(R)) && throw(ArgumentError("ranges must be strictly monotonically increasing"))
    curstart = first(firstrg)
    curstop = last(firstrg)
    curstep = step(firstrg)

    for i = firstrgidx:length(ranges)
        newrg = ranges[i]
        isempty(newrg) && continue

        if step(newrg) == curstep && first(newrg) == curstop+curstep
            # Can extend current range
            curstop = last(newrg)
        else
            if first(newrg) <= curstop || step(newrg) < zero(eltype(R))
                throw(ArgumentError("ranges must be strictly monotonically increasing"))
            end

            # Need to make a new range
            currg = constructrange(R, curstart, curstep, curstop)
            push!(newranges, currg)
            push!(offsets, offsets[end]+length(currg))

            curstart = first(newrg)
            curstop = last(newrg)
            curstep = step(newrg)
        end
    end

    push!(newranges, constructrange(R, curstart, curstep, curstop))
    (newranges, offsets)
end

function combine_ranges(ranges::Vector{R}, firstrg::R, firstrgidx::Int) where {R<:StepRangeLen}
    newranges = R[]
    offsets = Int[1]

    if step(firstrg) <= zero(eltype(R))
        throw(ArgumentError("ranges must be strictly monotonically increasing"))
    end
    currg = firstrg

    for i = firstrgidx:length(ranges)
        newrg = ranges[i]
        isempty(newrg) && continue

        if step(newrg) == step(currg) && first(newrg) == last(currg) + step(currg)
            currg = R(currg.ref, step(currg),
                        length(currg) + length(newrg),currg.offset)
        else
            if first(newrg) <= last(currg) || step(newrg) <= zero(eltype(R))
                throw(ArgumentError("ranges must be strictly monotonically increasing"))
            end

            # Need to make a new range
            push!(newranges, currg)
            push!(offsets, offsets[end]+length(currg))
            currg = newrg
        end
    end

    push!(newranges, currg)
    (newranges, offsets)
end 

struct PiecewiseIncreasingRange{T,R<:AbstractRange,S} <: AbstractPiecewiseRange{T,R}
    ranges::Vector{R}
    offsets::Vector{Int}
    divisor::S

    function PiecewiseIncreasingRange{T,R,S}(ranges::AbstractVector{R}, divisor::S) where {T,R<:AbstractRange,S}
        ranges = convert(Vector{R}, ranges)
        isempty(ranges) && return new{T,R,S}(ranges, Int[])

        # Find first non-empty range
        j = findfirst(!isempty,ranges)
        isnothing(j) && return new{T,R,S}(R[], Int[], divisor)

        firstrg = ranges[j]
        newranges, offsets = combine_ranges(ranges, firstrg, j+1)
        new{T,R,S}(newranges, offsets, divisor)
    end
end
PiecewiseIncreasingRange(ranges::Vector{R}, divisor) where {R<:AbstractRange} = 
    PiecewiseIncreasingRange{typeof(inv(one(eltype(R)))),R,typeof(divisor)}(ranges, divisor)
PiecewiseIncreasingRange(ranges::Vector{R}) where {R<:AbstractRange} = 
    PiecewiseIncreasingRange{eltype(R),R,@compat(Nothing)}(ranges, nothing)
PiecewiseIncreasingRange{T,R,@compat(Nothing)}(ranges::Vector{R}) where {T,R<:AbstractRange} = 
    PiecewiseIncreasingRange{T,R,@compat(Nothing)}(ranges, nothing)

# const PiecewiseUnitRange{T,R<:AbstractUnitRange} = PiecewiseIncreasingRange{T,R,@compat(Nothing)}
struct PiecewiseUnitRange{T,R<:AbstractUnitRange} <: AbstractPiecewiseOrdinalRange{T,R}
    ranges::Vector{R}
    offsets::Vector{Int}

    function PiecewiseUnitRange{T,R}(ranges::AbstractVector) where {T,R<:AbstractUnitRange}
        ranges = convert(Vector{R}, ranges)
        isempty(ranges) && return new{T,R}(ranges, Int[])

        # Find first non-empty range
        j = findfirst(!isempty,ranges)
        isnothing(j) && return new{T,R}(R[], Int[])

        firstrg = ranges[j]
        newranges, offsets = combine_ranges(ranges, firstrg, j+1)
        new{T,R}(newranges, offsets)
    end
end
PiecewiseUnitRange(ranges::Vector{R}) where {R<:AbstractUnitRange} = 
    PiecewiseUnitRange{eltype(R),R}(ranges)

# const PiecewiseIncreasingStepRange{T,R<:StepRange} = PiecewiseIncreasingRange{T,R,@compat(Nothing)}
struct PiecewiseIncreasingStepRange{T,R<:StepRange} <: AbstractPiecewiseOrdinalRange{T,R}
    ranges::Vector{R}
    offsets::Vector{Int}

    function PiecewiseIncreasingStepRange{T,R}(ranges::AbstractVector) where {T,R<:StepRange}
        ranges = convert(Vector{R}, ranges)
        isempty(ranges) && return new{T,R}(ranges, Int[])

        # Find first non-empty range
        j = findfirst(!isempty,ranges)
        isnothing(j) && return new{T,R}(R[], Int[])

        firstrg = ranges[j]
        newranges, offsets = combine_ranges(ranges, firstrg, j+1)
        new{T,R}(newranges, offsets)
    end
end
PiecewiseIncreasingStepRange(ranges::Vector{R}) where {R<:StepRange} = 
    PiecewiseIncreasingStepRange{eltype(R),R}(ranges)
function PiecewiseIncreasingStepRange(ranges::Vector{R}) where {R<:AbstractUnitRange}
    T = eltype(R)
    PiecewiseIncreasingStepRange{T,StepRange{T,T}}(ranges)
end
PiecewiseIncreasingStepRange(x::PiecewiseUnitRange) = PiecewiseIncreasingStepRange(x.ranges)

Base.convert(::Type{PiecewiseIncreasingRange{T,R,S}}, x::PiecewiseIncreasingRange) where {T,R,S} =
    PiecewiseIncreasingRange{T,R,S}(x.ranges, x.divisor)

function Base.convert(::Type{PiecewiseIncreasingStepRange{T,R}},x::PiecewiseUnitRange) where {T,R}
    PiecewiseIncreasingStepRange{T,R}(x.ranges)
end

# Avoid applying the divisor if it is one, to get types right
divide_divisor(::PiecewiseIncreasingRange{T,R,@compat(Nothing)}, x) where {T,R} = x
multiply_divisor(::PiecewiseIncreasingRange{T,R,@compat(Nothing)}, x) where {T,R} = x
divide_divisor(r::PiecewiseIncreasingRange{T,R,S}, x) where {T,R,S} = x/r.divisor
multiply_divisor(r::PiecewiseIncreasingRange{T,R,S}, x) where {T,R,S} = x*r.divisor

divide_divisor(::AbstractPiecewiseOrdinalRange, x) = x
multiply_divisor(::AbstractPiecewiseOrdinalRange, x) = x

function Base.size(r::AbstractPiecewiseRange)
    isempty(r.ranges) && return (0,)
    return (r.offsets[end]+length(r.ranges[end])-1,)
end

function Base.getindex(r::AbstractPiecewiseRange, i::Integer)
    rgidx = searchsortedlast(r.offsets, i, Forward)
    divide_divisor(r, r.ranges[rgidx][i-r.offsets[rgidx]+1])
end

@inline _inds(r,x::OrdinalRange,ind,elmax) = first(x) - r.offsets[ind]+1:step(x):elmax
@inline _inds(r,x::AbstractUnitRange,ind,elmax) = first(x)-r.offsets[ind]+1:elmax

@inline function _newrgx(r,x,ind)
    elmax = min(last(x)-r.offsets[ind]+1, length(r.ranges[ind]))
    newrg = r.ranges[ind][_inds(r,x,ind,elmax)]
    x = x[length(newrg)+1:end]
    return newrg,x
end

@inline function _fill_newrgs!(newrgs,r,x,firstrgidx,lastrgidx)
    for irange in firstrgidx:lastrgidx
        newrg,x = _newrgx(r,x,irange)
        newrgs[irange-firstrgidx+1] = newrg
    end
    return newrgs
end

@inline function _compute_newrgs(Tret::Type,r::AbstractPiecewiseRange,x)
    firstrgidx = searchsortedlast(r.offsets, first(x), Forward)
    lastrgidx = searchsortedlast(r.offsets, last(x), Forward)
    newrgs = Vector{Tret}(undef, lastrgidx-firstrgidx+1)
    _fill_newrgs!(newrgs,r,x,firstrgidx,lastrgidx)
end

@inline function _getindex(r::AbstractPiecewiseRange{T,R},
    x::AbstractUnitRange{<:Integer}) where {T,R}
    _compute_newrgs(R,r,x)
end

# These methods are needed to avoid ambiguities
@inline function _getindex(r::AbstractPiecewiseRange{T,R},
    x::AbstractUnitRange{<:Integer}) where {T,R<:Union{OrdinalRange,StepRangeLen}}
    _compute_newrgs(R,r,x)
end

@inline function _getindex(r::AbstractPiecewiseRange{T,R},
    x::AbstractUnitRange{<:Integer}) where {T,R<:AbstractUnitRange}
    _compute_newrgs(R,r,x)
end

# Slicing a Rational StepRange produces a StepRangeLen, so need to handle this separately
@inline function _getindex(r::AbstractPiecewiseRange{Rational{T},StepRange{Rational{T},Rational{T}}},
    x::AbstractUnitRange{<:Integer}) where {T<:Integer}

    Tret = StepRangeLen{Rational{T},Rational{T}}
    _compute_newrgs(Tret,r,x)
end

@inline function _getindex(r::AbstractPiecewiseRange{Rational{T},StepRange{Rational{T},Rational{T}}},
    x::OrdinalRange{<:Integer}) where {T<:Integer}

    Tret = StepRangeLen{Rational{T},Rational{T}}
    _compute_newrgs(Tret,r,x)
end

# Slicing an AbsractUnitRange using a StepRange produces another StepRange
@inline function _getindex(r::AbstractPiecewiseRange{T,R},
    x::OrdinalRange{<:Integer}) where {T,R<:AbstractUnitRange}

    Tret = StepRange{eltype(R),eltype(R)}
    _compute_newrgs(Tret,r,x)
end

# Slicing an integral StepRange or a StepRangeLen produces the same type
@inline function _getindex(r::AbstractPiecewiseRange{T,R},
    x::OrdinalRange{<:Integer}) where {T,R<:Union{OrdinalRange,StepRangeLen}}

    _compute_newrgs(R,r,x)
end

function Base.getindex(r::PiecewiseIncreasingRange{T,R,S}, x::AbstractRange{Int}) where {T,R,S}
    isempty(x) && return PiecewiseIncreasingRange{T,R,S}(R[], Int[], r.divisor)
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    # Not type-stable!
    if step(x) < 0
        argtype = Tuple{PiecewiseIncreasingRange{T,R,S}, AbstractVector{Int}}
        return invoke(getindex, argtype, r, x)
    end

    newrgs = _getindex(r,x)
    PiecewiseIncreasingRange{T,eltype(newrgs),S}(newrgs, r.divisor)
end

function Base.getindex(r::PiecewiseUnitRange{T,R}, x::AbstractUnitRange{Int}) where {T,R}
    isempty(x) && return PiecewiseUnitRange{T,R}(R[], Int[])
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    PiecewiseUnitRange(_getindex(r,x))
end

@inline @propagate_inbounds function Base.getindex(r::PiecewiseUnitRange{T,R}, x::IncreasingStepRange{Int,Int}) where {T,R}
    isempty(x) && return PiecewiseIncreasingStepRange{T,R}(R[], Int[])
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    PiecewiseIncreasingStepRange(_getindex(r,x))
end

@inline @propagate_inbounds function Base.getindex(r::PiecewiseUnitRange{T,R}, x::StepRange{Int,Int}) where {T,R}
    isempty(x) && return PiecewiseIncreasingStepRange{T,R}(R[], Int[])
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    # Not type-stable!
    if step(x) < 0 
        argtype = Tuple{PiecewiseUnitRange{T,R}, AbstractVector{Int}}
        return invoke(getindex, argtype , r, x)
    end

    PiecewiseIncreasingStepRange(_getindex(r,x))
end

@inline @propagate_inbounds function Base.getindex(r::PiecewiseIncreasingStepRange{T,R}, x::StepRange{Int,Int}) where {T,R}
    isempty(x) && return PiecewiseIncreasingStepRange{T,R}(R[], Int[])
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    # Not type-stable!
    if step(x) < 0 
        argtype = Tuple{PiecewiseIncreasingStepRange{T,R}, AbstractVector{Int}}
        return invoke(getindex, argtype , r, x)
    end

    PiecewiseIncreasingStepRange(_getindex(r,x))
end

@inline @propagate_inbounds function Base.getindex(r::PiecewiseIncreasingStepRange{T,R}, 
    x::Union{AbstractUnitRange{Int},IncreasingStepRange{Int,Int}}) where {T,R}
    
    isempty(x) && return PiecewiseIncreasingStepRange{T,R}(R[], Int[])
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    PiecewiseIncreasingStepRange(_getindex(r,x))
end

@inline @propagate_inbounds Base.getindex(r::AbstractPiecewiseRange,::Colon) = r

function resample(r::PiecewiseIncreasingRange{T,R,S}, ratio::Rational{Int}) where {T,R,S}
    excess = zero(eltype(R))
    curpt = 0
    newrgs = R[]
    for i = 1:length(r.ranges)
        endpt = last(r.ranges[i])*numerator(ratio)
        newrg = first(r.ranges[i])*numerator(ratio)+excess*step(r.ranges[i]):step(r.ranges[i])*denominator(ratio):endpt
        curpt += length(newrg)
        push!(newrgs, newrg)
        if i != length(r.ranges)
            # Need to linearly interpolate between endpoint and next range
            # nextind is the next index in r (= 1 + curpt/ratio) times numerator(ratio)
            nextind = numerator(ratio) + curpt*denominator(ratio)
            if nextind < r.offsets[i+1]*numerator(ratio)
                # Need to interpolate between this point and the next
                weight = mod(nextind, numerator(ratio))
                interpt = last(r.ranges[i])*(numerator(ratio)-weight) + first(r.ranges[i+1])*weight
                interpstep = (first(r.ranges[i+1]) - last(r.ranges[i]))*denominator(ratio)
                newrg = interpt:interpstep:first(r.ranges[i+1])*numerator(ratio)
                push!(newrgs, newrg)
                curpt += length(newrg)
                nextind += length(newrg)*denominator(ratio)
            end
            excess = nextind - r.offsets[i+1]*numerator(ratio)
        end
    end
    PiecewiseIncreasingRange(newrgs, multiply_divisor(r, numerator(ratio)))
end

# searchsortedfirst, searchsortedlast
struct PiecewiseIncreasingRangeFirstOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeFirstOrdering, a, b) = isless(first(a), first(b))
struct PiecewiseIncreasingRangeLastOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeLastOrdering, a, b) = isless(last(a), last(b))

function Base.searchsortedfirst(r::AbstractPiecewiseRange, x)
    isempty(r.ranges) && return 1
    xd = multiply_divisor(r, x)

    rgidx = searchsortedfirst(r.ranges, xd, PiecewiseIncreasingRangeLastOrdering())
    rgidx > length(r.ranges) && return length(r) + 1
    searchsortedfirst(r.ranges[rgidx], xd, Forward) + r.offsets[rgidx] - 1
end

function Base.searchsortedlast(r::AbstractPiecewiseRange, x)
    isempty(r.ranges) && return 1
    xd = multiply_divisor(r, x)

    rgidx = searchsortedlast(r.ranges, xd, PiecewiseIncreasingRangeFirstOrdering())
    rgidx == 0 && return 0
    searchsortedlast(r.ranges[rgidx], xd, Forward) + r.offsets[rgidx] - 1
end

struct NoNearestSampleError <: Exception end

function findnearest(r::AbstractPiecewiseRange, x, within_half_step::Bool=false)
    isempty(r.ranges) && throw(NoNearestSampleError())
    xd = multiply_divisor(r, x)

    rgidx = searchsortedfirst(r.ranges, xd, PiecewiseIncreasingRangeLastOrdering())
    if rgidx > length(r.ranges)
        rgend = r.ranges[end]
        within_half_step && xd > rgend[end]+step(rgend)/2 && throw(NoNearestSampleError())
        return length(r)
    end

    rg = r.ranges[rgidx]
    idxinrg = searchsortedfirst(rg, xd, Forward)
    idx = idxinrg + r.offsets[rgidx] - 1
    d = rg[idxinrg] - xd

    if idxinrg == 1
        if rgidx == 1
            # First element of all elements
            within_half_step && d > step(rg)/2 && throw(NoNearestSampleError())
            return 1
        end

        # Could be closer to last element of preceding range
        rgprev = r.ranges[rgidx-1]
        dprev = xd - rgprev[end]
        within_half_step && d > step(rg)/2 && dprev > step(rgprev)/2 && throw(NoNearestSampleError())
        return idx - (dprev <= d)
    end

    ifelse(d >= step(rg)/2, idx-1, idx)
end

# These functions let us use PiecewiseIncreasingRange as an array axis
Base.checkindex(::Type{Bool}, inds::AbstractPiecewiseRange, i) =
    throw(ArgumentError("unable to check bounds for indices of type $(typeof(i))"))

function Base.checkindex(::Type{Bool},inds::AbstractPiecewiseRange,i::Real)
    inbound=false
    for r in inds.ranges
        inbound |= divide_divisor(inds,first(r)) <= i <= divide_divisor(inds,last(r))
        inbound && break
    end
    inbound
end

Base.checkindex(::Type{Bool},inds::AbstractPiecewiseRange,::Colon) = true
Base.checkindex(::Type{Bool},inds::AbstractPiecewiseRange,::Base.Slice) = true
function Base.checkindex(::Type{Bool}, indx::AbstractPiecewiseRange, I::AbstractVector{Bool}) 
    indx == Base.axes1(parent(I))
end
Base.checkindex(::Type{Bool}, indx::AbstractPiecewiseRange, I::AbstractArray{Bool}) = false
function Base.checkindex(::Type{Bool}, inds::AbstractPiecewiseRange, I::AbstractVector)
    b = true
    for i in I
        b &= checkindex(Bool, inds, i)
    end
    b
end

end # module
