module CompactSpaceDiscretizers

using Discretizers
using IntervalSets
using PartialFunctions

export AbstractCompactSpaceDiscretizer
export UniformLinearDiscretizer

abstract type AbstractCompactSpaceDiscretizer end
Discretizers.encode(d::AbstractCompactSpaceDiscretizer, s) = s
Discretizers.decode(d::AbstractCompactSpaceDiscretizer, k) = k

struct UniformLinearDiscretizer{S} <: AbstractCompactSpaceDiscretizer
    discretizers::Vector{S}
    n::Int64
    dotter::Vector{Int64}
    undotter::Vector{Int64}
end

function UniformLinearDiscretizer(space::S, n::Int) where {S <: Vector{<: AbstractInterval}}
    lows = flip(getproperty) $ :left .<| space
    highs = flip(getproperty) $ :right .<| space
    dotter = n .^ (Base.OneTo(length(space)) .- 1)
    UniformLinearDiscretizer(
        LinearDiscretizer.(LinRange.(lows, highs, n)),
        n,
        dotter,
        reverse(dotter)
    )
end

function Discretizers.encode(d::UniformLinearDiscretizer, x)
    grid = Discretizers.encode.(d.discretizers, x)
    (grid .- 1)' * d.dotter
end

function Discretizers.decode(d::UniformLinearDiscretizer, n)
    grid = (n ./ d.undotter) .% d.undotter
    Discretizers.decode.(d.discretizers, grid .+ 1)
end

end # module
