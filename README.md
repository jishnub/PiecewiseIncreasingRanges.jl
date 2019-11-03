# PiecewiseIncreasingRanges

[![Build Status](https://travis-ci.com/jishnub/PiecewiseIncreasingRanges.jl.svg?branch=dev)](https://travis-ci.com/jishnub/PiecewiseIncreasingRanges.jl)
[![Coverage Status](https://coveralls.io/repos/simonster/PiecewiseIncreasingRanges.jl/badge.svg?branch=master)](https://coveralls.io/r/simonster/PiecewiseIncreasingRanges.jl?branch=master)

PiecewiseIncreasingRanges.jl provides a `PiecewiseIncreasingRange` type that corresponds to a set of concatenated, monotonically increasing ranges. It implements indexing as well as `searchsortedfirst`, `searchsortedlast`, and a `findnearest` function. This type is useful for holding the potentially discontinuous timestamps corresponding to large quantities of sampled data, where creating and indexing a vector would be undesirably slow.

There are two other types that are provided for convenience: `PiecewiseUnitRange` and `PiecewiseIncreasingStepRange`, which are concatenated `Unitrange`s and `StepRange`s respectively. These behave like aliases of `PiecewiseIncreasingRange{T,UnitRange{T},Nothing}` and `PiecewiseIncreasingRange{T,StepRange{T},Nothing}`, but have been implemented as separate types.


An extra type `IncreasingStepRange` is introduced for type stability. This can be used to slice any of the piecewise ranges. This behaves exactly like a `StepRange` with the extra constraint that the `step` is guaranteed to be positive. The instability is only in the final return type in `getindex` which returns a `Union` instead of a concrete type, so the performance is almost identical because of [union-splitting](https://julialang.org/blog/2018/08/union-splitting). However `IncreasingStepRange` is included just in case it leads to more efficient code.

## Construction

You can construct a `PiecewiseIncreasingRange` as

```julia
julia> PiecewiseIncreasingRange([0:5:15, 18:2:20])
6-element PiecewiseIncreasingRanges.PiecewiseIncreasingRange{Int64,StepRange{Int64,Int64},Nothing}:
  0
  5
 10
 15
 18
 20
```

`PiecewiseIncreasingRange`s also accept an optional divisor, which is divided out of all elements:

```julia
julia> PiecewiseIncreasingRange(UnitRange{Int}[0:3, 15:16], 4)
6-element PiecewiseIncreasingRanges.PiecewiseIncreasingRange{Float64,UnitRange{Int64},Int64}:
 0.0 
 0.25
 0.5 
 0.75
 3.75
 4.0
```

The types `PiecewiseUnitRange` and `PiecewiseIncreasingStepRange` can be constructed as 

```julia
julia> PiecewiseUnitRange([0:2, 5:7])
6-element PiecewiseUnitRange{Int64,UnitRange{Int64}}:
 0
 1
 2
 5
 6
 7

julia> PiecewiseIncreasingStepRange([0:2:2,5:3:11])
5-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
  0
  2
  5
  8
 11
```

## Slicing

Slicing a `PiecewiseIncreasingRange` produces another `PiecewiseIncreasingRange` for ranges with positive steps. Otherwise it defaults to an array.

```julia
julia> p=PiecewiseIncreasingRange([0:4:20, 61:15:80])
8-element PiecewiseIncreasingRange{Int64,StepRange{Int64,Int64},Nothing}:
  0
  4
  8
 12
 16
 20
 61
 76

julia> p[1:4]
4-element PiecewiseIncreasingRange{Int64,StepRange{Int64,Int64},Nothing}:
  0
  4
  8
 12

julia> p[1:3:7] # not type stable as return type depends on the sign of the step
3-element PiecewiseIncreasingRange{Int64,StepRange{Int64,Int64},Nothing}:
  0
 12
 61

julia> p[IncreasingStepRange(1:3:7)]  # type stable
3-element PiecewiseIncreasingRange{Int64,StepRange{Int64,Int64},Nothing}:
  0
 12
 61

julia> p[4:-1:1]
4-element Array{Int64,1}:
 12
  8
  4
  0
```

Slicing a `PiecewiseUnitRange` produces a `PiecewiseUnitRange` if sliced with an `UnitRange`, a `PiecewiseIncreasingStepRange` if sliced with a `StepRange` or an `IncreasingStepRange`, and an array otherwise.

```julia
julia> p=PiecewiseUnitRange([1:3,6:7,10:12])
8-element PiecewiseUnitRange{Int64,UnitRange{Int64}}:
  1
  2
  3
  6
  7
 10
 11
 12

julia> p[1:4]
4-element PiecewiseUnitRange{Int64,UnitRange{Int64}}:
 1
 2
 3
 6

julia> p[1:2:5] # not type stable as return type depends on the sign of the step
3-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
 1
 3
 7

julia> p[IncreasingStepRange(1:2:5)] # type stable
3-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
 1
 3
 7

julia> p[5:-2:1]
3-element Array{Int64,1}:
 7
 3
 1
```

Slicing a `PiecewiseIncreasingStepRange` produces another `PiecewiseIncreasingStepRange` if the step size is positive, or an array otherwise. Type stability is achieved using `IncreasingStepRange`.

```julia
julia> p=PiecewiseIncreasingStepRange([1:2:5,10:4:18])
6-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
  1
  3
  5
 10
 14
 18

julia> p[1:4]
4-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
  1
  3
  5
 10

julia> p[1:2:5] # not type stable as return type depends on the sign of the step
3-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
  1
  5
 14

julia> p[IncreasingStepRange(1:2:5)] # type stable
3-element PiecewiseIncreasingStepRange{Int64,StepRange{Int64,Int64}}:
  1
  5
 14

julia> p[5:-2:1]
3-element Array{Int64,1}:
 14
  5
  1
```

## Searchsortedfirst and searchsortedlast

These functions are implemented for all the types.

```julia
julia> p=PiecewiseUnitRange([1:3,6:7,10:12])
8-element PiecewiseUnitRange{Int64,UnitRange{Int64}}:
  1
  2
  3
  6
  7
 10
 11
 12

julia> searchsortedfirst(p,4)
4

julia> searchsortedlast(p,4)
3
```

## findnearest

`findnearest(rg, x, within_half_sample)` finds the index of the element of `rg` closest to `x`. If `x` is equidistant between two samples, it chooses the earlier one. If `within_half_sample` is true, `findnearest` throws a `NoNearestSampleError` if there is no sample within half of a step of `x`.