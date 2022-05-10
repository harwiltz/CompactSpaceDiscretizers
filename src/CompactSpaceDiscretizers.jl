module CompactSpaceDiscretizers

using Discretizers
using IntervalSets
using PartialFunctions

export AbstractCompactSpaceDiscretizer
export UniformLinearDiscretizer
export UniformLatticeDiscretizer
export discrete_coords

abstract type AbstractCompactSpaceDiscretizer end
Discretizers.encode(d::AbstractCompactSpaceDiscretizer, s) = s
Discretizers.decode(d::AbstractCompactSpaceDiscretizer, k) = k

struct UniformLinearDiscretizer{S, I, R} <: AbstractCompactSpaceDiscretizer
    discretizers::Vector{S}
    bins_per_dim::I
    ϕ::Vector{Int64}
    n::Int64
    intervals::Vector{R}
end

intervals(disc::AbstractCompactSpaceDiscretizer) = disc.intervals

Base.in(x, disc::AbstractCompactSpaceDiscretizer) = (x .∈ intervals(disc)) |> all

function UniformLinearDiscretizer(
    space::S,
    bins_per_dim::Int;
    I = UInt8
) where {S <: Vector{<: AbstractInterval}}
    lows = flip(getproperty) $ :left .<| space
    highs = flip(getproperty) $ :right .<| space
    UniformLinearDiscretizer(lows, highs, bins_per_dim; I = I, intervals = space)
end

function UniformLinearDiscretizer(
    low::Vector{F},
    high::Vector{F},
    bins_per_dim::Int;
    I = UInt8,
    intervals::Union{Nothing, S} = nothing
) where {F <: Number, S <: Vector{<: AbstractInterval}}
    @assert (bins_per_dim > 0) "bins_per_dim must be greater than 0"
    @assert (bins_per_dim <= typemax(I)) "Too many bins per dimension, try I > $(I)"
    ϕ = bins_per_dim .^ (Base.OneTo(length(low)) .- 1)
    n = bins_per_dim ^ length(low)
    i = isnothing(intervals) ? zip(low, high) .|> Base.splat(ClosedInterval) : intervals
    UniformLinearDiscretizer(
        LinearDiscretizer.(LinRange.(low, high, bins_per_dim)),
        bins_per_dim,
        ϕ,
        n,
        i
    )
end

function Discretizers.encode(d::UniformLinearDiscretizer, x)
    x̅ = clamp.(x, d.intervals)
    coords = Discretizers.encode.(d.discretizers, x̅)
    (coords .- 1)' * d.ϕ
end

"""
Let n = 5, d = 3
[a, b, c] -> (a-1) + 5(b-1) + 25(c-1) = q
q -> {c = q ÷ 25 + 1, b = (q % 25) ÷ 5 + 1, a = ((q % 25) % 5) ÷ 1 + 1 = q % 5 + 1}
Ex: [2, 4, 3] -> 1 + 15 + 50 = 66
           66 -> {c = 3, b = 4, a = 2}
"""
function Discretizers.decode(d::AbstractCompactSpaceDiscretizer, n::Int)
    @assert (n <= d.n) "Value $(n) is not in the encoding space of the discretizer"
    coords = discrete_coords(d, n)
    Discretizers.decode.(d.discretizers, coords)
end

discrete_coords(d::UniformLinearDiscretizer, n::Int) = (n .% (d.bins_per_dim .* d.ϕ)) .÷ d.ϕ .+ 1 

struct UniformLatticeDiscretizer{S, F, R} <: AbstractCompactSpaceDiscretizer
    discretizers::Vector{S}
    ϵ::F
    dim_nbins::Vector{Int}
    ϕ::Vector{Int64}
    n::Int64
    intervals::Vector{R}
end

function UniformLatticeDiscretizer(
    space::S,
    ϵ::Real,
    R = Float64,
) where {S <: Vector{<: AbstractInterval}}
    lows = flip(getproperty) $ :left .<| space
    highs = flip(getproperty) $ :right .<| space
    UniformLatticeDiscretizer(lows, highs, ϵ; R = R, intervals = space)
end

function UniformLatticeDiscretizer(
    low::Vector{F},
    high::Vector{F},
    ϵ::Real;
    R = Float64,
    intervals::Union{Nothing, S} = nothing
) where {F <: Number, S <: Vector{<: AbstractInterval}}
    @assert (ϵ > 0)
    dim_nbins = floor.((high - low) ./ ϵ .+ 1)
    ϕ = fill(R(1), length(low))
    ϕ[2:end] = cumprod(dim_nbins[1:end - 1])
    n = prod(dim_nbins)
    i = isnothing(intervals) ? zip(low, high) .|> Base.splat(ClosedInterval) : intervals
    UniformLatticeDiscretizer(
        LinearDiscretizer.(LinRange.(low, high, Int.(dim_nbins))),
        ϵ,
        Int.(dim_nbins),
        Int.(ϕ),
        Int(n),
        i
    )
end

function Discretizers.encode(d::UniformLatticeDiscretizer, x)
    x̅ = clamp.(x, d.intervals)
    coords = Discretizers.encode.(d.discretizers, x̅)
    (coords .- 1)' * d.ϕ
end

discrete_coords(d::UniformLatticeDiscretizer, n::Int) = (n .% (d.dim_nbins .* d.ϕ)) .÷ d.ϕ .+ 1 

end # module
