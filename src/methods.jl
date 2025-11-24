# Abstraction

abstract type AbstractMethod end
abstract type AbstractMethodMesh end

function build_method_mesh end


# Quadrature Interpolation

struct PM_MM_Interpolation <: AbstractInterpolant
    interpolant_dif::Matrix{Float64}
    interpolant_alg::Matrix{Float64}
    
    function PM_MM_Interpolation(
        mesh::AbstractPointsMesh,
        quad_mesh::AbstractPointsMesh,
    )

        interpolant_dif = _interpolation_matrix(
            mesh.points_dif,
            mesh.bary_weights_dif,
            quad_mesh.points_alg,
        )

        interpolant_alg = _interpolation_matrix(
            mesh.points_alg,
            mesh.bary_weights_alg,
            quad_mesh.points_alg,
        )

        return new(interpolant_dif, interpolant_alg)
    end

    function PM_MM_Interpolation(
        interpolant_dif::Matrix{Float64},
        interpolant_alg::Matrix{Float64}
    )
        return new(interpolant_dif, interpolant_alg)
    end
end

function _identity_matrix(n::Integer)
    I = zeros(Float64, n, n)
    @inbounds for i in 1:n
        I[i, i] = 1.0
    end
    return I
end

# Collocation

struct Collocation <: AbstractMethod end

struct CollocationMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractMethodMesh
    quad_points_mesh::P
    interpolant::I
end

function CollocationMesh(mesh::AbstractPointsMesh)
    interpolant_dif = _identity_matrix(length(mesh.points_dif))
    interpolant_alg = _identity_matrix(length(mesh.points_alg))
    return CollocationMesh(mesh, PM_MM_Interpolation(interpolant_dif, interpolant_alg))
end

mesh_type(::Type{Collocation}) = CollocationMesh

build_method_mesh(::Collocation, mesh::AbstractPointsMesh) = CollocationMesh(mesh)


# Integrated Residual

abstract type AbstractIntRes <: AbstractMethod end

struct DAIR{T<:AbstractPoints} <: AbstractIntRes
    quad_points::T
end

DAIR(number::Integer) = DAIR(GLPoints(number))

struct QPM{T<:AbstractPoints} <: AbstractIntRes
    quad_points::T
    pen_param::Real
end

QPM(number::Integer; pen_param::Real=1.0) = QPM(GLPoints(number), pen_param)

struct SAIR{T<:AbstractPoints} <: AbstractIntRes
    quad_points::T
end

SAIR(number::Integer) = SAIR(GLPoints(number))


"""
    quad_var_vars[dyn_var][i] should be Vector{MathOptInterface.ScalarAffineFunction{Float64}}
    quad_var_vars[dyn_var][i][q] should replace dyn_var_vars[dyn_var][i][q] in the original transcribe_dyn_fun

    dyn_var_vars::DYN_VAR_VARS
"""

abstract type AbstractIntResMesh <: AbstractMethodMesh end

# DAIR
struct DAIRMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractIntResMesh
    quad_points_mesh::P
    interpolant::I
end

function DAIRMesh(points::DAIR, mesh::AbstractPointsMesh)

    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)

    return DAIRMesh(quad_mesh, interpolant)
end

mesh_type(::Type{DAIR}) = DAIRMesh

build_method_mesh(points::DAIR, mesh::AbstractPointsMesh) = DAIRMesh(points, mesh)

# QPM
struct QPMMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractIntResMesh
    quad_points_mesh::P
    interpolant::I
    pen_param::Real
end

QPMMesh(quad_points_mesh::P, interpolant::I) where {P<:AbstractPointsMesh,I<:AbstractInterpolant} =
    QPMMesh(quad_points_mesh, interpolant, 1.0)

function QPMMesh(points::QPM, mesh::AbstractPointsMesh)

    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)

    return QPMMesh(quad_mesh, interpolant, points.pen_param)
end

mesh_type(::Type{QPM}) = QPMMesh

build_method_mesh(points::QPM, mesh::AbstractPointsMesh) = QPMMesh(points, mesh)

# SAIR
struct SAIRMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractIntResMesh
    quad_points_mesh::P
    interpolant::I
end

function SAIRMesh(points::SAIR, mesh::AbstractPointsMesh)

    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)

    return SAIRMesh(quad_mesh, interpolant)
end

mesh_type(::Type{SAIR}) = SAIRMesh

build_method_mesh(points::SAIR, mesh::AbstractPointsMesh) = SAIRMesh(points, mesh)
