Base.@kwdef struct VehicleParams
    ρ::Float64    = 1.225            # kg/m^3
    C_dA::Float64 = 1.20             # drag area (Cd*A) [m^2] -> used as CdA lump
    C_lA::Float64 = 4.50             # downforce area (|Cl|*A) [m^2]
    A_b::Float64  = 0.40             # front downforce bias [–]

    m::Float64    = 795.0            # kg
    g::Float64    = 9.81             # m/s^2
    l::Float64    = 3.135            # m
    W_b::Float64  = 0.45
    l_f::Float64  = l * (1.0 - W_b)  # m
    l_r::Float64  = l * W_b          # m
    J_zz::Float64 = 1000.0           # kg·m^2 yaw inertia

    R_e::Float64  = 0.35             # m, effective tyre radius
    J_wf::Float64 = 1.0              # kg·m^2 front wheel inertia
    J_wr::Float64 = 1.0              # kg·m^2 rear wheel inertia

    B_fx::Float64 = 20.0             # long slip "B" front
    B_rx::Float64 = 20.0             # long slip "B" rear
    B_fy::Float64 = 15.0             # lat slip "B" front
    B_ry::Float64 = 15.0             # lat slip "B" rear

    D_fx::Float64 = 1.00             # peak μ long front
    D_rx::Float64 = 1.00             # peak μ long rear
    D_fy::Float64 = 1.00             # peak μ lat  front
    D_ry::Float64 = 1.00             # peak μ lat  rear

    C_fx::Float64 = 1.50
    C_rx::Float64 = 1.50
    C_fy::Float64 = 1.20
    C_ry::Float64 = 1.20

    T_e::Float64  = 570.0            # Nm (engine torque at crank scaled by throttle)
    τ_g::Float64  = 3.00             # overall gear ratio (crank -> wheel)
    B_b::Float64  = 0.60             # brake bias to front
    B_kf::Float64 = 5000.0           # brake torque gain (front)
    B_kr::Float64 = 5000.0           # brake torque gain (rear)

    kFd::Float64  = 0.5 * ρ * C_dA
    kFlf::Float64 = 0.5 * ρ * C_lA * A_b
    kFlr::Float64 = 0.5 * ρ * C_lA * (1.0 - A_b)
    kWf::Float64  = m * g * (l_r / l)
    kWr::Float64  = m * g * (l_f / l)

    κ_c::Float64  = 0.0              # track curvature [1/m]; constant here
    n::Float64    = 5.0              # track width
end