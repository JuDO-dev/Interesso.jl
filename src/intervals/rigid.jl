"""
    RigidIntervals{T}


"""
struct RigidIntervals{T} <: Intervals{T}
    n::Int
    t::Vector{T}

    """
        RigidIntervals(t::Vector{T})
    
    Constructor for RigidIntervals();
    """
    function RigidIntervals(t::Vector{T}) where {T}

        n = length(t) - 1;
        for i in 1:n
            t[i] ≤ t[i+1] ? nothing : throw(DomainError((t[i], t[i+1]), "Ensure "));
        end

        return new{T}(n, t)
    end
end

"""

"""
RigidIntervals(t_0::T, t_f::T, n::Integer) where {T} = RigidIntervals(collect(range(start=t_0, stop=t_f, length=n+1)));