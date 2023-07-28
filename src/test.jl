function residuals(ẋ, x)#, u, p)
    r[1] = ẋ[1] + 1 * x[1];
    return nothing
end

z = zeros(3 + 3);

function f(z::Vector{F})::F where F<:AbstractFloat
    #for interval
    for q in 1:5 #quadrature
        rz = zeros(5);
        residuals!(rz, z[1:5], z[6:10]);
    end

#function integrate_residuals(z) =

function rigid_residuals_integration(z_start, z_mid, z_end)
    return rz
end

function flexible_residuals_integration(z1, z2, z3, t, z)









function integrate_residuals(z_n)
    w_q, t_q = FGQ.gausslegendre(5);
    
    r_q = zeros()
    
    
    #Evaluate r! at each quadrature point
    for w_q, t_q in Q
        r_q[]

        residuals!(r1, z[], z[])

    #Sum squares

    # 