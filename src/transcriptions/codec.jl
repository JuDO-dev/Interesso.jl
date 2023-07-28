struct Metadata
    n_x::Int
    n_u::Int
    n_p::Int
    n_τ_x::Int
    n_τ_u::Int
    x_offset::Int
    u_offset::Int
    z_offsets::Vector{Int}
    
    
    p_offset::Int
    


    p_indices::UnitRange{Int}
end

function Metadata(dop::DOProblem, transcription::Transcription{T, I}) where {T, I<:RigidIntervals}
    
    n_x = dop.x isa Nothing ? 0 : length(dop.x);
    n_u = dop.u isa Nothing ? 0 : length(dop.u);
    n_p = dop.p isa Nothing ? 0 : length(dop.p);
    n_τ_x = transcription.x_rule_rule.n_p;
    n_τ_u = transcription.u.n_p;
    u_offset = n_x * (n_τ_x - 1);
    x_offset = u_offset + n_u * n_τ_u;
    z_offsets = [(i-1)*x_offset+1 for i in 1:transcription.intervals.n];


    
    #indices = [(i-1)*x_offset+1 : i*x_offset+n_x for i in 1:transcription.intervals.n];
    p_offset = indices[end][end];
    p_indices = p_offset + 1 : p_offset + n_p;

    return Metadata(n_x, n_u, n_p, n_τ_x, n_τ_u, x_offset, u_offset, p_offset, indices, p_indices)
end

#function Metadata(dop::DOProblem, transcription::Transcription{T, I}) where {T, I<:FlexibleIntervals} end

function encode_warmstart(warmstart::Warmstart, metadata::Metadata, transcription::Transcription{T, I}) where {T<:Real, I<:RigidIntervals}
    
    n_x = metadata.n_x;
    n_u = metadata.n_u;
    n_p = metadata.n_p;
    n_τ_x = metadata.n_τ_x;
    n_τ_u = metadata.n_τ_u;
    x_offset = metadata.x_offset;
    u_offset = metadata.u_offset;
    p_offset = metadata.p_offset;
    indices = metadata.indices;
    
    z_encoded = Vector{T}(undef, p_offset + n_p);

    τ_x = Vector{T}(undef, n_τ_x);
    τ_u = Vector{T}(undef, n_τ_u);
    z_encoded_views = [view(z_encoded, indices_i) for indices_i in indices];
    for (i, z_encoded_i) in enumerate(z_encoded_views)
        
        dtdτ = (transcription.intervals.t[i+1] - transcription.intervals.t[i]) / 2;
        t̄_i = (transcription.intervals.t[i+1] + transcription.intervals.t[i]) / 2;
        @. τ_x = transcription.x.τ * dtdτ + t̄_i;
        @. τ_u = transcription.u.τ * dtdτ + t̄_i;

        for m in 1:n_x
            z_encoded_i[m]                                      = warmstart.x[m](τ_x[1]);
            @. z_encoded_i[n_x+(m-1)*(n_τ_x-2)+1 : n_x+m*(n_τ_x-2)] = warmstart.x[m](τ_x[2:end-1]);
            z_encoded_i[x_offset+m]                             = warmstart.x[m](τ_x[end]);
        end
        for m in 1:n_u
            @. z_encoded_i[u_offset+(m-1)*(n_τ_u)+1 : u_offset+m*n_τ_u] = warmstart.u[m](τ_u);
        end
    end

    for m in 1:n_p
        z_encoded[p_offset+m] = warmstart.p[m];
    end

    return z_encoded
end

function encode_bounds(dop::DOProblem, metadata::Metadata, ::Transcription{T}) where {T<:Real}

    n_x = metadata.n_x;
    n_u = metadata.n_u;
    n_p = metadata.n_p;
    n_τ_x = metadata.n_τ_x;
    n_τ_u = metadata.n_τ_u;
    x_offset = metadata.x_offset;
    u_offset = metadata.u_offset;
    p_offset = metadata.p_offset;
    indices = metadata.indices;
    
    x_lower = isempty(dop.x) ? nothing : [x_m.bounds.lower for x_m in dop.x];
    x_upper = isempty(dop.x) ? nothing : [x_m.bounds.upper for x_m in dop.x];
    u_lower = isempty(dop.u) ? nothing : [u_m.bounds.lower for u_m in dop.u];
    u_upper = isempty(dop.u) ? nothing : [u_m.bounds.upper for u_m in dop.u];
    p_lower = isempty(dop.p) ? nothing : [p_m.bounds.lower for p_m in dop.p];
    p_upper = isempty(dop.p) ? nothing : [p_m.bounds.upper for p_m in dop.p];
    
    z_lower_i = Vector{T}(undef, indices[1][end]);
    z_upper_i = Vector{T}(undef, indices[1][end]);
    
    @. z_lower_i[1:n_x] = x_lower;
    @. z_upper_i[1:n_x] = x_upper;
    for m in 1:n_x
        @. z_lower_i[n_x+(m-1)*(n_τ_x-2)+1 : n_x+m*(n_τ_x-2)] = x_lower[m];
        @. z_upper_i[n_x+(m-1)*(n_τ_x-2)+1 : n_x+m*(n_τ_x-2)] = x_upper[m];
    end
    for m in 1:n_u
        @. z_lower_i[u_offset+(m-1)*(n_τ_u)+1 : u_offset+m*n_τ_u] = u_lower[m];
        @. z_upper_i[u_offset+(m-1)*(n_τ_u)+1 : u_offset+m*n_τ_u] = u_upper[m];
    end
    @. z_lower_i[x_offset+1:end] = x_lower;
    @. z_upper_i[x_offset+1:end] = x_upper;

    # Let's tessalate
    z_lower = Vector{T}(undef, indices[end][end] + n_p);
    z_upper = Vector{T}(undef, indices[end][end] + n_p);
    for indices_i in indices
        @. z_lower[indices_i] = z_lower_i;
        @. z_upper[indices_i] = z_upper_i;
    end
    for m in 1:n_p
        z_lower[p_offset+m] = p_lower[m];
        z_upper[p_offset+m] = p_upper[m];
    end
    
    @. z_lower[1:n_x] = [x_m.initial.lower for x_m in dop.x]; 
    @. z_upper[1:n_x] = [x_m.initial.upper for x_m in dop.x];
    @. z_lower[x_offset+1:end] = [x_m.final.lower for x_m in dop.x];
    @. z_upper[x_offset+1:end] = [x_m.final.upper for x_m in dop.x];

    return Progradio.Box(z_lower, z_upper)
end

function encode(warmstart::Warmstart, dop::DOProblem, transcription::Transcription)

    metadata = Metadata(dop, transcription);

    return encode_warmstart(warmstart, metadata, transcription),
        encode_bounds(dop, metadata, transcription), metadata
end

#=
function decode(z::Vector, transcription::Transcription)

    x_Warmstart = [PiecewiseInterpolant([LagrangeInterpolant()]) for m in ]


    return Warmstart(x_Warmstart, u_Warmstart)
end

function decode(problem::DOProblem{F}, M::RigidMesh{F}, z::Vector{F}) where {F<:AbstractFloat}
    n_x = length(problem.x);
    n_u = length(problem.u);
    n_τ_x = M.P_x.degree + 1;
    n_τ_u = M.P_u.degree + 1;
    x_offset = n_x * (n_τ_x - 1) + n_u * n_τ_u;
    u_offset = n_x * (n_τ_x - 1);
    
    # Slice z into views
    indices = indexing(M, n_x, n_u);
    z_views = [view(z, indices_i) for indices_i in indices];

    return [IntervalData([vcat(z_i[n], z_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)], z_i[x_offset+n]) for n in 1:n_x], 
    [z_i[u_offset+(n-1)*(n_τ_u)+1 : u_offset+n*n_τ_u] for n in 1:n_u]) for z_i in z_views]
end

function codec_indices()



    return indices
end

function X!(xz, I, n_x, n_τ_x, x_offset, z_i, j)

    for m in 1:n_x
        xz[m] =      I[j, 1]   * z_i[m];
        xz[m] += sum(I[j, l]   * z_i[n_x+(m-1)*(n_τ_x-2)+l] for l in 2:n_τ_x-1);
        xz[m] +=     I[j, end] * z_i[x_offset+m];
    end
    return nothing
end

function U!(uz, I, n_u, n_τ_u, u_offset, z_i, j)
    
    for m in 1:n_u
        uz[m] = sum(I[j, l] * z_i[u_offset+(m-1)*(n_τ_u)+l] for l in 1:n_τ_u);
    end
    return nothing
end
=#

#=
mutable struct IntervalData{F<:AbstractFloat}
    x::Vector{Vector{F}}
    u::Vector{Vector{F}}
end

function IntervalData(::UndefInitializer, n_x::I, n_u::I, M::AbstractMesh{F}) where {F<:AbstractFloat, I<:Integer}
    n_τ_x = M.P_x.degree + 1;
    n_τ_u = M.P_u.degree + 1;
    return IntervalData([Vector{F}(undef, n_τ_x) for _ in 1:n_x], [Vector{F}(undef, n_τ_u) for _ in 1:n_u])
end
=#
#function coldstart(problem::DOProblem{F}, M::RigidMesh{F}) where {F<:AbstractFloat}
#    n_x = length(problem.x);
#    n_u = length(problem.u);
#    n_τ_x = M.P_x.degree + 1;
#    n_τ_u = M.P_u.degree + 1;
#    return [IntervalData([zeros(F, n_τ_x) for _ in 1:n_x], [zeros(F, n_τ_u) for _ in 1:n_u]) for _ in 1:length(M.H)]
#end
#=
function encode(problem::DOProblem{F}, M::RigidMesh{F}, data::Vector{IntervalData{F}}) where {F<:AbstractFloat}
    n_x = length(problem.x);
    n_u = length(problem.u);
    indices = indexing(M, n_x, n_u);
    n_τ_x = M.P_x.degree + 1;
    n_τ_u = M.P_u.degree + 1;
    x_offset = n_x * (n_τ_x - 1) + n_u * n_τ_u;
    u_offset = n_x * (n_τ_x - 1);

    # Data
    z_0 = Vector{F}(undef, indices[end][end]);
    z_0_views = [view(z_0, indices_i) for indices_i in indices];
    for (i, z_0_i) in enumerate(z_0_views)
        for n in 1:n_x
            z_0_i[n]                                    = data[i].x[n][1];
            z_0_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)] .= data[i].x[n][2:end-1];
            z_0_i[x_offset+n]                           = data[i].x[n][end];
        end
        for n in 1:n_u
            z_0_i[u_offset+(n-1)*(n_τ_u)+1 : u_offset+n*n_τ_u] .= data[i].u[n];
        end
    end

    # Bounds per interval
    x_ℓ = [x_n.lower for x_n in problem.x];
    x_u = [x_n.upper for x_n in problem.x];
    u_ℓ = [u_n.lower for u_n in problem.u];
    u_u = [u_n.upper for u_n in problem.u];
    z_ℓ_i = Vector{F}(undef, indices[1][end]);
    z_u_i = similar(z_ℓ_i);
    z_ℓ_i[1:n_x] .= x_ℓ[:];
    z_u_i[1:n_x] .= x_u[:];
    for n in 1:n_x
        z_ℓ_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)] .= x_ℓ[n];
        z_u_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)] .= x_u[n];
    end
    for n in 1:n_u
        z_ℓ_i[u_offset+(n-1)*(n_τ_u)+1 : u_offset+n*n_τ_u] .= u_ℓ[n];
        z_u_i[u_offset+(n-1)*(n_τ_u)+1 : u_offset+n*n_τ_u] .= u_u[n];
    end
    z_ℓ_i[end-n_x+1:end] .= x_ℓ[:];
    z_u_i[end-n_x+1:end] .= x_u[:];

    # Let's tessellate
    z_ℓ = similar(z_0);
    z_u = similar(z_0);
    for indices_i in indices
        z_ℓ[indices_i] .= z_ℓ_i[:];
        z_u[indices_i] .= z_u_i[:];
    end

    # Boundary conditions
    z_ℓ[1:n_x] .= [x_n.initial_lower for x_n in problem.x];
    z_u[1:n_x] .= [x_n.initian_τ_upper for x_n in problem.x];
    z_ℓ[end-n_x+1:end] .= [x_n.final_lower for x_n in problem.x];
    z_u[end-n_x+1:end] .= [x_n.finan_τ_upper for x_n in problem.x];

    return z_0, z_ℓ, z_u
end
=#
