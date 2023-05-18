struct DirectCollocation{T, I, PD<:PointDistribution{T}} <: Transcription{T, I}
    intervals::I
    x::PD
    u::PD
    q::PD
    
    DirectCollocation(point_distribution::PD, intervals::I) where {T,
        I<:Intervals{T}, PD<:PointDistribution{T}} = new{T, I, PD}(
        intervals,
        point_distribution,
        point_distribution,
        point_distribution
    );
end

