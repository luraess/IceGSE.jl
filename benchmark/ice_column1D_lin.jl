# Solve: 0 = eta * dVx/dz - rhog(H-z)*sin(alpha)
# See https://doi.org/10.5194/tc-2020-295
using Plots, Printf, LinearAlgebra

const s2y = 3600*24*365.25

@views function ice_column1D()
    # physics
    H     = 1e3       # domain height [m]
    α     = 3.0*π/180 # slope angle   [rad]
    ρ     = 930       # ice density   [kg/m^3]
    g     = -9.81     # grav. accel.  [m/s^2]
    η0    = 1e14      # ice viscosity [Pa.s]
    T     = 253       # ice temp.     [K]
    ϵ_bg  = 1e-12     # background strain-rate [1/s]
    # numerics
    nz    = 100
    nout  = 100
    tol   = 1e-8
    iterMax = 1e4
    damp  = 1.0-2.0/nz
    # derived numerics
    dz    = H/nz
    zc    = LinRange(-H+dz/2, 0.0-dz/2, nz  )
    zv    = LinRange(-H     , 0.0,      nz+1)
    ηv    = η0.*ones(nz+1)
    Vxe   = zeros(nz+2)
    ϵxz   = zeros(nz+1)
    τxz   = zeros(nz+1)
    dVxdτ = zeros(nz)
    Fx    = zeros(nz)
    # time loop
    iter = 1; err = 2*tol
    while err>tol && iter<iterMax
        ϵxz     .= 0.5.*diff(Vxe)/dz .+ ϵ_bg
        τxz     .= 2.0*ηv.*ϵxz
        Fx      .= diff(τxz)/dz .- ρ.*g.*sin(α)
        dVxdτ   .= Fx .+ dVxdτ*damp
        dτVx     = dz*dz/η0/2.1
        Vxe[2:end-1] .= Vxe[2:end-1] .+ dτVx.*dVxdτ
        Vxe[1]   = -Vxe[2]
        Vxe[end] =  Vxe[end-1]
        if mod(iter,nout)==0  err = norm(Fx)/length(Fx); @printf("iter %d, %2.2e\n", iter, err)  end
        iter += 1
    end
    # visualisation
    p1 = plot(Vxe[2:end-1]*s2y, zc, legend=:false, linewidth=:3.0, framestyle=:box, xlabel="Vx [m/yr]", ylabel="depth [m]")
    p2 = plot(τxz*1e-6, zv, legend=:false, linewidth=:3.0, framestyle=:box, xlabel="τxz [MPa]")
    display(plot(p1, p2))
    return
end

ice_column1D()
