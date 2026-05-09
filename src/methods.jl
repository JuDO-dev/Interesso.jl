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

abstract type AbstractIntResMesh <: AbstractMethodMesh end

"""
    quad_var_vars[dyn_var][i] should be Vector{MathOptInterface.ScalarAffineFunction{Float64}}
    quad_var_vars[dyn_var][i][q] should replace dyn_var_vars[dyn_var][i][q] in the original transcribe_dyn_fun

    dyn_var_vars::DYN_VAR_VARS
"""


# Petrov-Galerkin weighted residual method
struct Galerkin <: AbstractIntRes
    quad_points::GLPoints
end

Galerkin(number::Integer) = Galerkin(GLPoints(number))

struct GalerkinMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractIntResMesh
    quad_points_mesh::P
    interpolant::I
    test_values_dif::Matrix{Float64}
    test_values_alg::Matrix{Float64}
end

function GalerkinMesh(method::Galerkin, mesh::AbstractPointsMesh)
    quad_mesh = GLPointsMesh(method.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)

    n_quad = length(method.quad_points.points_alg_τ)
    moments_dif = get_points_dif_length(mesh) - 1
    moments_alg = get_points_alg_length(mesh)
    max_moments = max(moments_dif, moments_alg)

    if moments_dif < 1
        throw(DomainError(moments_dif, "Please ensure get_points_dif_length(mesh) - 1 ≥ 1."))
    elseif moments_alg < 1
        throw(DomainError(moments_alg, "Please ensure get_points_alg_length(mesh) ≥ 1."))
    elseif max_moments > n_quad
        throw(DomainError(n_quad, "Please ensure the number of quadrature points is at least max(get_points_dif_length(mesh) - 1, get_points_alg_length(mesh))."))
    end

    test_values = Matrix{Float64}(undef, max_moments, n_quad)

    # ϕ_n(τ) = sqrt((2n+1)/2) P_n(τ)
    for (j, τ_j) in enumerate(method.quad_points.points_alg_τ)
        test_values[1, j] = sqrt(0.5)
        max_moments == 1 && continue
        # degree 1
        p_nm2 = 1.0
        p_nm1 = τ_j
        test_values[2, j] = sqrt(1.5) * p_nm1
        # degree >= 2
        for n in 2:(max_moments - 1)
            p_n = ((2n - 1) * τ_j * p_nm1 - (n - 1) * p_nm2) / n
            test_values[n + 1, j] = sqrt(n + 0.5) * p_n
            p_nm2 = p_nm1
            p_nm1 = p_n
        end
    end

    return GalerkinMesh(
        quad_mesh,
        interpolant,
        test_values[1:moments_dif, :],
        test_values[1:moments_alg, :],
    )
end

mesh_type(::Type{Galerkin}) = GalerkinMesh

build_method_mesh(method::Galerkin, mesh::AbstractPointsMesh) = GalerkinMesh(method, mesh)


# DAIR family
abstract type AbstractDAIR <: AbstractIntRes end

struct DAIR{T<:AbstractPoints} <: AbstractDAIR
    quad_points::T
end

DAIR(number::Integer) = DAIR(GLPoints(number))

abstract type AbstractDAIRMesh <: AbstractIntResMesh end

# DAIRFeas
struct DAIRFeas{T<:AbstractPoints} <: AbstractDAIR
    quad_points::T
end

DAIRFeas(number::Integer) = DAIRFeas(GLPoints(number))
DAIRFeas(d::AbstractDAIR) = DAIRFeas(d.quad_points)
struct DAIRFeasMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractDAIRMesh
    quad_points_mesh::P
    interpolant::I
end

function DAIRFeasMesh(points::DAIRFeas, mesh::AbstractPointsMesh)
    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)
    return DAIRFeasMesh(quad_mesh, interpolant)
end

mesh_type(::Type{DAIRFeas}) = DAIRFeasMesh

build_method_mesh(points::DAIRFeas, mesh::AbstractPointsMesh) = DAIRFeasMesh(points, mesh)

# DAIROpti
struct DAIROpti{T<:AbstractPoints} <: AbstractDAIR
    quad_points::T
    tolerance::Float64
end

DAIROpti(number::Integer; tolerance::Float64=1e-8) = DAIROpti(GLPoints(number), tolerance)
DAIROpti(d::AbstractDAIR; tolerance::Float64=1e-8) = DAIROpti(d.quad_points, tolerance)

struct DAIROptiMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant} <: AbstractDAIRMesh
    quad_points_mesh::P
    interpolant::I
    tolerance::Float64
end

function DAIROptiMesh(points::DAIROpti, mesh::AbstractPointsMesh)
    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)
    return DAIROptiMesh(quad_mesh, interpolant, points.tolerance)
end

mesh_type(::Type{DAIROpti}) = DAIROptiMesh

build_method_mesh(points::DAIROpti, mesh::AbstractPointsMesh) = DAIROptiMesh(points, mesh)


# QPM
struct QPM{T<:AbstractPoints, R<:Real} <: AbstractIntRes
    quad_points::T
    penalty::R
end

QPM(number::Integer; penalty::Real=1.0) = QPM(GLPoints(number), penalty)

struct QPMMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant, R<:Real} <: AbstractIntResMesh
    quad_points_mesh::P
    interpolant::I
    penalty::R
end

QPMMesh(quad_points_mesh::P, interpolant::I) where {P<:AbstractPointsMesh,I<:AbstractInterpolant} =
    QPMMesh(quad_points_mesh, interpolant, 1.0)

function QPMMesh(points::QPM, mesh::AbstractPointsMesh)
    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)
    return QPMMesh(quad_mesh, interpolant, points.penalty)
end

mesh_type(::Type{QPM}) = QPMMesh

build_method_mesh(points::QPM, mesh::AbstractPointsMesh) = QPMMesh(points, mesh)


# SAIR
struct SAIR{T<:AbstractPoints} <: AbstractIntRes
    quad_points::T
end

SAIR(number::Integer) = SAIR(GLPoints(number))

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


# SAPM
struct SAPM{T<:AbstractPoints, R<:Real} <: AbstractIntRes
    quad_points::T
    penalty::R
end

SAPM(number::Integer; penalty::Real=1.0) = SAPM(GLPoints(number), penalty)

struct SAPMMesh{P<:AbstractPointsMesh, I<:AbstractInterpolant, R<:Real} <: AbstractIntResMesh
    quad_points_mesh::P
    interpolant::I
    penalty::R
end

SAPMMesh(quad_points_mesh::P, interpolant::I) where {P<:AbstractPointsMesh,I<:AbstractInterpolant} =
    SAPMMesh(quad_points_mesh, interpolant, 1.0)

function SAPMMesh(points::SAPM, mesh::AbstractPointsMesh)

    quad_mesh = GLPointsMesh(points.quad_points, mesh.t_a, mesh.t_b)
    interpolant = PM_MM_Interpolation(mesh, quad_mesh)

    return SAPMMesh(quad_mesh, interpolant, points.penalty)
end

mesh_type(::Type{SAPM}) = SAPMMesh

build_method_mesh(points::SAPM, mesh::AbstractPointsMesh) = SAPMMesh(points, mesh)
