function bounds(problem::JuDOProblem{F}, transcription::InteressoTranscription{F}, indices::Vector{UnitRange{I}}) where {F<:AbstractFloat, I<:Integer}
    x_ℓ = [x_n.lower for x_n in problem.x];
    x_u = [x_n.upper for x_n in problem.x];
    u_ℓ = [u_n.lower for u_n in problem.u];
    u_u = [u_n.upper for u_n in problem.u];
    n_x = length(x_ℓ);
    n_u = length(u_ℓ);
    l_x = transcription.x.degree + 1;
    l_u = transcription.u.degree + 1;
    
    # Bounds per interval
    z_ℓ_i = Vector{F}(undef, indices[1][end]);
    z_u_i = similar(z_ℓ_i);
    z_ℓ_i[1:n_x] .= x_ℓ[:];
    z_u_i[1:n_x] .= x_u[:];
    for n in 1:n_x
        z_ℓ_i[n_x+(n-1)*(l_x-2)+1 : n_x+n*(l_x-2)] .= x_ℓ[n];
        z_u_i[n_x+(n-1)*(l_x-2)+1 : n_x+n*(l_x-2)] .= x_u[n];
    end
    offset = n_x * (l_x - 1);
    for n in 1:n_u
        z_ℓ_i[offset+(n-1)*(l_u)+1 : offset+n*l_u] .= u_ℓ[n];
        z_u_i[offset+(n-1)*(l_u)+1 : offset+n*l_u] .= u_u[n];
    end
    z_ℓ_i[end-n_x+1:end] .= x_ℓ[:];
    z_u_i[end-n_x+1:end] .= x_u[:];

    # Let's tessellate
    z_ℓ = Vector{F}(undef, indices[end][end]);
    z_u = similar(z_ℓ);
    for i in 1:length(indices)
        z_ℓ[indices[i]] .= z_ℓ_i[:];
        z_u[indices[i]] .= z_u_i[:];
    end

    # Boundary conditions
    z_ℓ[1:n_x] .= [x_n.initial[1] for x_n in problem.x];
    z_u[1:n_x] .= [x_n.initial[2] for x_n in problem.x];
    z_ℓ[end-n_x+1:end] .= [x_n.terminal[1] for x_n in problem.x];
    z_u[end-n_x+1:end] .= [x_n.terminal[2] for x_n in problem.x];

    return (z_ℓ, z_u)
end