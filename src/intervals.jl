# Abstraction 

"""
    AbstractIntervals

Supertype for.
"""
abstract type AbstractIntervals end

"""
    AbstractIntervalsMesh{
        PM<:AbstractPointsMesh,
        MM<:AbstractMethodMesh,
        BM<:AbstractBoundsMesh
    }
    
    Supertype for.
"""
abstract type AbstractIntervalsMesh{
    PM<:AbstractPointsMesh,
    MM<:AbstractMethodMesh,
    BM<:AbstractBoundsMesh
} end

function get_intervals_length end
function get_points_meshes end
function get_bounds_mesh end
function build_intervals_mesh end


# Fixed

struct FixedIntervals <: AbstractIntervals
    number::Int64
    points::Vector{Float64}

    function FixedIntervals(
        number::Integer;
        points::Vector{<:Real}=collect(range(0.0, 1.0, number + 1)),
    )
        if !(number ≥ 1)
            throw(DomainError(number, "Please ensure number ≥ 1."))
        end

        _throw_points_mismatch(number, points)

        return new(number, points)
    end
end

function _throw_points_mismatch(number::Integer, points::Vector{<:Real})

    if !(length(points) == number + 1)
        throw(DomainError(points, "Please ensure length(points) == number + 1."))

    elseif !(points[1] == 0.0) && !(points[end] == 1.0)
        throw(DomainError(points, "Please ensure points[1] == 0.0 and points[end] == 1.0."))

    elseif !issorted(points)
        throw(DomainError(points, "Please ensure issorted(points) == true."))

    elseif !(all(0.0 .≤ points .≤ 1.0))
        throw(DomainError(points, "Please ensure all(0 ≤ points ≤ 1) == true."))
    end

    return nothing
end

struct FixedIntervalsMesh{PM,MM,BM} <: AbstractIntervalsMesh{PM,MM,BM}
    points_meshes::Vector{PM}
    method_meshes::Vector{MM}
    bounds_mesh::BM
end

function mesh_type(
    ::Type{FixedIntervals},
    PM::Type{<:AbstractPointsMesh},
    MM::Type{<:AbstractMethodMesh},
    BM::Type{<:AbstractBoundsMesh},
)
    return FixedIntervalsMesh{PM,MM,BM}
end

get_intervals_length(mesh::FixedIntervalsMesh) = length(mesh.points_meshes)
get_points_meshes(mesh::FixedIntervalsMesh) = mesh.points_meshes
get_bounds_mesh(mesh::FixedIntervalsMesh) = mesh.bounds_mesh

get_points_dif_length(mesh::FixedIntervalsMesh) = get_points_dif_length(mesh.points_meshes[1])
get_points_alg_length(mesh::FixedIntervalsMesh) = get_points_alg_length(mesh.points_meshes[1])
get_points_quad_length(mesh::FixedIntervalsMesh) = get_points_quad_length(mesh.method_meshes[1].quad_points_mesh)


function build_intervals_mesh(
    intervals::FixedIntervals,
    points::AbstractPoints,
    method::AbstractMethod,
    bounds::AbstractBounds,
    t_0::Real,
    t_f::Real,
)
    _throw_if_invalid_bounds(t_0, t_f)

    Δt = t_f - t_0

    points_meshes = [build_points_mesh(
        points,
        t_0 + Δt * intervals.points[i],
        t_0 + Δt * intervals.points[i+1], 
    ) for i in 1:intervals.number]

    method_meshes = [build_method_mesh(
        method,
        points_meshes[i],
    ) for i in 1:intervals.number]

    bounds_mesh = build_bounds_mesh(bounds, points)

    return FixedIntervalsMesh(points_meshes, method_meshes, bounds_mesh)
end


# Flexible

struct FlexibleIntervals <: AbstractIntervals
    number::Int64
    flexibility::Float64
    points::Vector{Float64}

    function FlexibleIntervals(
        number::Integer,
        flexibility::Real;
        points::Vector{Float64}=collect(range(0.0, 1.0, number + 1)),
    )
        if !(number ≥ 2)
            throw(DomainError(number, "Please ensure number ≥ 2."))

        elseif !(0 ≤ flexibility ≤ 1)
            throw(DomainError(flexibility, "Please ensure 0 ≤ flexibility ≤ 1."))
        end

        _throw_points_mismatch(number, points)

        return new(number, flexibility, points)
    end
end

struct FlexibleIntervalsMesh{PM,MM,BM} <: AbstractIntervalsMesh{PM,MM,BM}
    fixed::FixedIntervalsMesh{PM,MM,BM}
    points_mesh::PM
    method_mesh::MM
    Δt_min::Float64
    Δt_max::Float64
end

function mesh_type(
    ::Type{FlexibleIntervals},
    PM::Type{<:AbstractPointsMesh},
    MM::Type{<:AbstractMethodMesh},
    BM::Type{<:AbstractBoundsMesh},
)
    return FlexibleIntervalsMesh{PM,MM,BM}
end

get_intervals_length(mesh::FlexibleIntervalsMesh) = get_intervals_length(mesh.fixed)
get_points_meshes(mesh::FlexibleIntervalsMesh) = get_points_meshes(mesh.fixed)
get_bounds_mesh(mesh::FlexibleIntervalsMesh) = get_bounds_mesh(mesh.fixed)

get_points_dif_length(mesh::FlexibleIntervalsMesh) = get_points_dif_length(mesh.fixed)
get_points_alg_length(mesh::FlexibleIntervalsMesh) = get_points_alg_length(mesh.fixed)
get_points_quad_length(mesh::FlexibleIntervalsMesh) = get_points_quad_length(mesh.fixed)

function build_intervals_mesh(
    intervals::FlexibleIntervals,
    points::AbstractPoints,
    method::AbstractMethod,
    bounds::AbstractBounds,
    t_0::Real,
    t_f::Real,
)
    fixed = build_intervals_mesh(
        FixedIntervals(intervals.number; points=intervals.points),
        points,
        method,
        bounds,
        t_0,
        t_f,
    )
    
    points_mesh = build_points_mesh(points, -1.0, 1.0)
    method_mesh = build_method_mesh(method, points_mesh)

    Δt = t_f - t_0
    Δt_min = (1 - intervals.flexibility) * Δt / intervals.number
    Δt_max = Δt_min + intervals.flexibility * Δt

    return FlexibleIntervalsMesh(fixed, points_mesh, method_mesh, Δt_min, Δt_max)
end