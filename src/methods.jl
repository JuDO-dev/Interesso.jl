# Abstraction

abstract type AbstractMethod end
abstract type AbstractMethodMesh end

function build_method_mesh end


# Collocation

struct Collocation <: AbstractMethod end

struct CollocationMesh <: AbstractMethodMesh end

mesh_type(::Type{Collocation}) = CollocationMesh

build_method_mesh(::Collocation, ::AbstractPointsMesh) = CollocationMesh()