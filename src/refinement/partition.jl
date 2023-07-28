function partition(M::AbstractMesh{F}, data::Vector{IntervalData{F}}, partitions::Vector{I}) where {F<:AbstractFloat, I<:Integer}
    n_x = length(data[1].x);
    n_u = length(data[1].u);
    
    i_partitioned = 1;
    H_partitioned = Vector{F}(undef, sum(partitions));
    data_partitioned = [IntervalData(undef, n_x, n_u, M) for _ in 1:length(H_partitioned)];
    for i in 1:length(M.H)
        if partitions[i] == 1
            H_partitioned[i_partitioned] = M.H[i];

            data_partitioned[i_partitioned] = deepcopy(data[i]);

            i_partitioned += 1;
        else
            for j in 1:partitions[i]
                # Partition Mesh Intervals
                H_partitioned[i_partitioned] = M.H[i] / partitions[i];
            
                # Interpolate Interval Data
                τartitioned = (M.P_x.nodes .+ (2 * j - 1 - partitions[i])) .* (1 / partitions[i]);
                for n in 1:n_x
                    data_partitioned[i_partitioned].x[n] .= [interpolate(M.P_x, data[i].x[n], τ) for τ in τartitioned];
                end
                τartitioned = (M.P_u.nodes .+ (2 * j - 1 - partitions[i])) .* (1 / partitions[i]);
                for n in 1:n_u
                    data_partitioned[i_partitioned].u[n] .= [interpolate(M.P_u, data[i].u[n], τ) for τ in τartitioned];
                end
                
                i_partitioned += 1;
            end
        end
    end

    return H_partitioned, data_partitioned
end