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
    ϕ::Vector{Int64}
end

function UniformLinearDiscretizer(space::S, n::Int) where {S <: Vector{<: AbstractInterval}}
    lows = flip(getproperty) $ :left .<| space
    highs = flip(getproperty) $ :right .<| space
    UniformLinearDiscretizer(lows, highs, n)
end

function UniformLinearDiscretizer(low::Vector{F}, high::Vector{F}, n::Int) where {F <: Number}
    ϕ = n .^ (Base.OneTo(length(low)) .- 1)
    UniformLinearDiscretizer(
        LinearDiscretizer.(LinRange.(low, high, n)),
        n,
        ϕ
    )
end

function Discretizers.encode(d::UniformLinearDiscretizer, x)
    coords = Discretizers.encode.(d.discretizers, x)
    (coords .- 1)' * d.ϕ
end

"""
Let n = 5, d = 3
[a, b, c] -> (a-1) + 5(b-1) + 25(c-1) = q
q -> {c = q ÷ 25 + 1, b = (q % 25) ÷ 5 + 1, a = ((q % 25) % 5) ÷ 1 + 1 = q % 5 + 1}
Ex: [2, 4, 3] -> 1 + 15 + 50 = 66
           66 -> {c = 3, b = 4, a = 2}
"""
function Discretizers.decode(d::UniformLinearDiscretizer, n::Int)
    coords = (n .% (d.n .* d.ϕ)) .÷ d.ϕ .+ 1
    Discretizers.decode.(d.discretizers, coords)
end

end # module
