struct FlexibleIntervals{T} <: Intervals{T}
    n::Int
    flexibility::T
    t_start::Vector{T}
end