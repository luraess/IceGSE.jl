# Solve: 0 = eta * dVx/dz - rhog(H-z)*sin(alpha)
# See https://doi.org/10.5194/tc-2020-295
using Plots, Printf, LinearAlgebra, JLD, Interpolations

const s2yr = 3600*24*365.25

@views function compute_rheol_Tau!(η, d, Tii2, T, P, ph, R, fact)
    T_h        = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
    Tii        = sqrt.(Tii2)
    # GBS Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    m_gbs      = 1.4
    n_gbs      = 1.8
    F_gbs      = 1.0#(2.0^((n_gbs-1.0)/n_gbs)*3.0^((n_gbs+1.0)/(2.0*n_gbs)))^-1.0
    A1_gbs     = 3.9e-3*10.0^(-6.0*n_gbs)              # A_gbs: preexponential constant[m^1.4/s*Pa^1.8] in Pa^(-n_gbs).s^(-1).m^(m_gbs)
    Q1_gbs     = 49e3                                  # Q1_gbs: lower activation energy[J/mol]
    Q2_gbs     = 192e3                                 # Q2_gbs: higher activation energy[J/mol]
    Tstar_gbs  = 259.0                                 # Tstar_gbs: activation threshold[K] (-18C)
    # Dislocation Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    n_disl     = 4.0
    F_disl     = 1.0#(2.0^((n_disl-1.0)/n_disl)*3.0^((n_disl+1.0)/(2.0*n_disl)))^-1.0 # pure shear epxeriments
    A1_disl    = 4e4*10.0^(-6.0*n_disl)                # A_disl: preexponential constant[1/s*Pa^4]  in  Pa^(-n_disl).s^(-1)
    Q1_disl    = 60e3                                  # Q1_disl: lower activation energy[J/mol]
    Q2_disl    = 180e3                                 # Q2_disl: higher activation energy[J/mol]
    Tstar_disl = 259.0                                 # Tstar_gbs: activation threshold[K] (-15C)
    # Temperature dependence of activation energy ONLY for GBS and DISL (melting point)
    Tstar_gbs_h  = Tstar_gbs  .+ ph.*P                 # Tstar_h:  adjusted activation threshold[K]
    Tstar_disl_h = Tstar_disl .+ ph.*P                 # Tstar_h:  adjusted activation threshold[K]
    Q_gbs_h    = zeros(size(T))                        # Q_gbs_h:  melting point adjusted activation energy[J/mol]
    Q_disl_h   = zeros(size(T))                        # Q_disl_h: melting point adjusted activation energy[J/mol]
    Q_gbs_h[ (T_h.-Tstar_gbs_h ).<=0.0] .= Q1_gbs
    Q_gbs_h[ (T_h.-Tstar_gbs_h ).> 0.0] .= Q2_gbs
    Q_disl_h[(T_h.-Tstar_disl_h).<=0.0] .= Q1_disl
    Q_disl_h[(T_h.-Tstar_disl_h).> 0.0] .= Q2_disl
    As_gbs     =  A1_gbs .*exp.(-Q1_gbs ./(R.*Tstar_gbs_h ))
    As_disl    =  A1_disl.*exp.(-Q1_disl./(R.*Tstar_disl_h))
    # Pre-computations
    # η_p .= 0.5.*a0.^(.-1.0).*exp(Q0./(R*T)).*Tii2.^((1.0-npow)/2.0)
    B_gbs      = F_gbs  .* As_gbs .^(-1.0) .* exp.((Q_gbs_h ./R).*((1.0./T_h).-(1.0./Tstar_gbs_h ))) # B_gbs      = A1_gbs .^(-1.0/n_gbs)   .* exp.((Q1_gbs ./(n_gbs  .*R .* T_h )))
    B_disl     = F_disl .* As_disl.^(-1.0) .* exp.((Q_disl_h./R).*((1.0./T_h).-(1.0./Tstar_disl_h))) # B_disl     = A1_disl.^(-1.0/n_disl)  .* exp.((Q1_disl./(n_disl .*R .* T_h )))
    # Isolated viscosities
    η_gbs      = fact^(n_gbs-1)  .* B_gbs  .* Tii.^(1.0-n_gbs ).*d.^m_gbs
    η_disl     = fact^(n_disl-1) .* B_disl .* Tii.^(1.0-n_disl)
    # initial viscosity guess
    η         .= 1.0./(1.0./η_gbs .+ 1.0./η_disl)
    return m_gbs, n_gbs, n_disl
end

@views function load_data!(T2, dat_file, skip)
     data  = load(joinpath(@__DIR__, dat_file))
     Depth = data["depth"]
     Temp  = data["temp"]
     Depth = [0.0;    Depth[1:end-skip]]
     Temp  = [Temp[1]; Temp[1:end-skip]] .+ 273.15
     zv2   = LinRange(Depth[1], Depth[end], length(T2))
     itp1  = interpolate((Depth[:],), Temp[:], Gridded(Linear()))
     T2   .= reverse(itp1(zv2))
     return zv2[1], zv2[end]
 end

@views function ice_column1D()
    # physics
    fact   = 1.0
    # H      = 3e3         # domain height [m]
    α      = 3.0*π/180   # slope angle   [rad]
    ρ      = 920         # ice density   [kg/m^3]
    g      = -9.81       # grav. accel.  [m/s^2]
    η0     = 1e15        # ice viscosity [Pa.s]
    Ti     = 253         # ice temp.     [K]
    R      = 8.314       # gas constant
    Q0     = 6*1e4       # activation energy
    a0     = 8.75*1e-13  # prefactor
    npow   = 3.0         # power-law exponent
    mpow   = -(1.0-1.0/npow)/2.0
    P0     = 1e7*0       # [Pa]
    d_i    = 3e-3        # d: grain size[m]
    ph     = 7e-8        # ph: pressure heating coef[K/Pa]
    ϵ_bg   = 1e-14
    # numerics
    ϵ_bg   = 1e-14       # background strain-rate [1/s]
    nz     = 100
    rel    = 1e-2
    nout   = 1e3
    tol    = 1e-5
    iterMax = 2e6
    # numerics gse
    litMax = 1e2
    tolg   = 1e-15
    # load and interpolate data
    T      =  zeros(nz+1)
    ztop, zbot = load_data!(T, "./data_gisp2.jld", 2)
    H      = zbot - ztop
    # derived numerics
    dz     = H/nz
    damp   = 1.0-0.1/nz
    # array initialisation
    Vxe    = zeros(nz+2)
    ϵxz    = zeros(nz+1)
    τxz    = zeros(nz+1)
    Eii2   = zeros(nz+1)
    Tii2   =  ones(nz+1)
    dVxdτ  = zeros(nz)
    dτVx   = zeros(nz)
    Fx     = zeros(nz)
    # initial conditions
    zc     = LinRange(-zbot+dz/2, ztop-dz/2, nz  )
    zv     = LinRange(-zbot     , ztop     , nz+1)
#     zc     = LinRange(-H+dz/2, 0.0-dz/2, nz  )
#     zv     = LinRange(-H     , 0.0,      nz+1)
#     T      =  Ti.*ones(nz+1)
    ηv     =  η0.*ones(nz+1)
    η_p    =  η0.*ones(nz+1)
    d      = d_i.*ones(nz+1)
    P      =  P0.*ones(nz+1)
    global m_gbs, n_gbs, n_disl, C_gbs, C_disl
    # time loop
    iter = 1; err = 2*tol
    while err>tol && iter<iterMax
        ϵxz     .= 0.5.*diff(Vxe)/dz .+ ϵ_bg
        # Eii2    .= ϵxz.^2
        m_gbs, n_gbs, n_disl = compute_rheol_Tau!(η_p, d, Tii2, T, P, ph, R, fact)
        # η_p     .= 0.5.*a0.^(.-1.0./npow).*Eii2.^mpow.*exp.(Q0./(npow.*R.*T))
        ηv      .= rel.*η_p .+ (1.0-rel).*ηv
        τxz     .= fact*ηv.*ϵxz
        Tii2    .= τxz.^2
        Fx      .= diff(τxz)/dz .- ρ.*g.*sin(α)
        dVxdτ   .= Fx .+ dVxdτ*damp
        dτVx    .= dz.*dz./(0.5.*(ηv[1:end-1].+ηv[2:end]))./2.1
        Vxe[2:end-1] .= Vxe[2:end-1] .+ dτVx.*dVxdτ
        Vxe[1]   = -Vxe[2]
        Vxe[end] =  Vxe[end-1]
        if mod(iter, nout)==0  err = norm(Fx)/length(Fx); @printf("iter %d, %2.2e\n", iter, err) end
        iter += 1
    end
    # Behn et al.
    Txz_Behn = ρ.*g.*(zv).*sin(α); #Txz_Behn[end] = 1e-9
    # Eii_gbs  =  C_gbs  .*Txz_Behn.^n_gbs  .* d_i.^(-m_gbs)
    # Eii_disl =  C_disl .*Txz_Behn.^n_disl
    # Eii_Behn = Eii_gbs .+ Eii_disl
    n1  = 1.8
    n2  = 4
    A1o = 3.9e-3
    A2o = 4e4  # Value wrong in K&G01
    Q1  = 49e3
    Q2  = 60e3
    m1  = 1.4
    Eii_gbs  = A1o.*(d_i.^-1.4).*((Txz_Behn./1e6./fact).^n1).*exp.(-Q1./(R.*T)) #gbs
    Eii_disl = A2o.*(             (Txz_Behn./1e6./fact).^n2).*exp.(-Q2./(R.*T)) #disl
    Eii_Behn = Eii_gbs + Eii_disl
    η_Behn   = Txz_Behn./fact./Eii_Behn
    Vx_Behn  = cumsum(Eii_Behn.*dz.*2.0)
    # visualisation
    p1 = plot(Vxe[2:end-1]*s2yr, zc, linewidth=:3.0, framestyle=:box, xlabel="Vx [m/yr]", ylabel="depth [m]", label="This study")
         plot!(Vx_Behn[1:10:end]*s2yr, zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21", legend=:topleft)
    p2 = plot(log10.(ϵxz), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(ϵxz) [1/s]", label="This study")
         plot!(log10.(Eii_Behn[1:10:end]), zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21", legend=:bottomleft)
    p3 = plot(1e-6*τxz, zv, linewidth=:3.0, framestyle=:box, xlabel="τxz [MPa]", ylabel="depth [m]", label="This study")
         plot!(1e-6*Txz_Behn[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
    p4 = plot(log10.(ηv), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(η) [Pa.s]", label="This study")
         plot!(log10.(η_Behn)[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21", legend=:bottomright)
    # p5 = plot(d*1e3, zv, legend=:false, linewidth=:3.0, framestyle=:box, xlabel="d [mm]")
    display(plot(p1, p2, p3, p4, background_color=:transparent, foreground_color=:gray))
    return
end

@time ice_column1D()
