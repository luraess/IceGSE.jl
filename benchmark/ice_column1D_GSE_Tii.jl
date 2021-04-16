# Solve: 0 = eta * dVx/dz - rhog(H-z)*sin(alpha)
# See https://doi.org/10.5194/tc-2020-295
using Plots, Printf, LinearAlgebra, Statistics

const s2yr = 3600*24*365.25

@views function compute_rheol_Tau!(η, η_gbs, η_disl, d, Tii2, T, P, ph, R, fact)
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
    η_gbs     .= fact^(n_gbs-1)  .* B_gbs  .* Tii.^(1.0-n_gbs ).*d.^m_gbs
    η_disl    .= fact^(n_disl-1) .* B_disl .* Tii.^(1.0-n_disl)
    # initial viscosity guess
    η         .= 1.0./(1.0./η_gbs .+ 1.0./η_disl)
    return m_gbs, n_gbs, n_disl
end
   
@views function compute_GSE_Tau!(dd, d, r_d, d_o, η_gbs, η_disl, Tii2, τxz, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, gam, p, c, dt, R, fact, steady, rel_d, dtau_d)
    Tii        = sqrt.(Tii2)
    T_h        = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
    Eii_gbs    = Tii./(fact*η_gbs)
    Eii_disl   = Tii./(fact*η_disl)
    if steady==true
        Eii        = Eii_gbs + Eii_disl
        β          = Eii_disl./Eii; β[end] = 0.0
        d_p   = ((Kgg.*exp.(.-Qgg./(R.*T_h)).*(1.0/p).*c.*gam) ./ ((λ_gbs .- β*λ_gbs .+ β.*λ_disl).*Tii.*Eii)).^(1.0./(1.0+p))
        d    .= rel_d*d_p .+ (1-rel_d)*d
        ddn   = d_p .- d
        errd  = mean(abs.(ddn[2:end-1]))
    else
        # work
        Wxz_disl = τxz .* Eii_disl .* τxz ./ Tii # Exz = Eii* Txz/Tii that's what we do
        Wxz_gbs  = τxz .* Eii_gbs  .* τxz ./ Tii
        Q_disl   = 1.0*Wxz_disl # ACHTUNG not sure about the factor 2 here
        Q_gbs    = 1.0*Wxz_gbs
        # Reduction rate
        dred     = .-d.^2 ./c./gam .* (λ_disl.*Q_disl .+ λ_gbs.*Q_gbs)   # d grain reduction eqn. (12)
        # Growth rate
        dgg      =  p.^(-1).*d.^(1-p).*Kgg.*exp.(.-Qgg./R./T_h)           # d grain growth eq. (4)
        # Grain size evolution residual
        # d       .= d .+ dt.*(dgg .+ dred)
        # rel_d    = 0.0001
        # d_p      = d_o .+ dt.*(dgg .+ dred)
        # d       .= rel_d*d_p .+ (1-rel_d)*d
        r_d2     = d_o./dt .+ (dgg .+ dred)
        dd      .= 0.5*dd .+ 0.5*r_d2   # Crank-Nicolson here on "numerics"
        d       .= (d .+ dtau_d.*dd)./(1.0 .+ dtau_d/dt)
        r_d     .= .-(d.-d_o)./dt .+ (dgg .+ dred)
        errd  = mean(abs.(r_d*dt))
    end
    return errd
end

@views function check_Behn(T, zv, d_i, λ_disl, λ_gbs, Kgg, Qgg, R, gam, p, c, ρ, g, α, dz, fact)
    # Goldsby-Kohlstedt Constants for T <= 259K
    n1     = 1.8
    n2     = 4
    A1o    = 3.9e-3
    A2o    = 4e4  # Value wrong in K&G01
    Q1     = 49e3
    Q2     = 60e3
    m1     = 1.4
    Txz_Behn = ρ.*g.*(zv).*sin(α); #Txz_Behn[end] = 1e-9
    dB       = d_i*ones(size(T))
    η_Behn   = zeros(size(T))
    Eii_Behn = zeros(size(T))
    errd   = 2.0
    while errd>1e-10
        Eii_gbs   = A1o.*(dB.^-1.4).*((Txz_Behn./1e6./fact).^n1).*exp.(-Q1./(R.*T)) #gbs
        Eii_disl  = A2o.*(            (Txz_Behn./1e6./fact).^n2).*exp.(-Q2./(R.*T)) #disl
        Eii_Behn .= Eii_gbs + Eii_disl
        β         = Eii_disl./Eii_Behn; β[end] = 0.0
        η_Behn   .= Txz_Behn./fact./Eii_Behn
        dBi       = ((Kgg.*exp.(.-Qgg./(R.*T)).*(1.0/p).*c.*gam) ./ ((λ_gbs .- β*λ_gbs .+ β.*λ_disl).*Txz_Behn.*Eii_Behn)).^(1.0./(1.0+p))
        dd        = dBi .- dB
        ddn       = dd./dB       # normalized grain size difference
        dB       .= dB .+ dd.*0.2
        errd      = mean(abs.(ddn[2:end-1]))
    end
    Vx_Behn = cumsum(Eii_Behn.*dz.*2.0)
    return Txz_Behn, Eii_Behn, η_Behn, dB, Vx_Behn
end

@views function ice_column1D()
    steady = false
    # physics
    fact   = 1.0
    H      = 1e3         # domain height [m]
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
    ph     = 7e-8        # pressure heating coef [K/Pa]
    d_i    = 3e-3        # grain size [m]
    dt     = 1e3*s2yr    # physical time step value in seconds
    nt     = 10          # number of time steps
    # GSE - From Behn et al. 2021
    λ_disl = 0.01        # fitting param
    λ_gbs  = 0.01        # fitting param
    Qgg    = 42e3        # grain growth activation energy
    Kgg    = 1.02e-17    # grain growth prefactor !! units
    gam    = 0.065       # param
    p      = 6.03        # grain growth exponant
    c      = 3.0         # geometric factor
    # numerics
    ϵ_bg   = 1e-14       # background strain-rate [1/s]
    nz     = 100
    rel    = 1e-3
    nout   = 1e3
    tol    = 1e-5
    iterMax = 5e5
    # numerics gse
    dtau_d = 0.08*s2yr   # when steady = false
    rel_d  = 5e-1        # when steady = true
    # derived numerics
    dz     = H/nz
    damp   = 1.0-0.1/nz
    # array initialisation
    d_o    = zeros(nz+1)
    Vxe    = zeros(nz+2)
    ϵxz    = zeros(nz+1)
    τxz    = zeros(nz+1)
    Eii2   = zeros(nz+1)
    Tii2   =  ones(nz+1)
    dVxdτ  = zeros(nz)
    dτVx   = zeros(nz)
    Fx     = zeros(nz)
    η_gbs  = zeros(nz+1)
    η_disl = zeros(nz+1)
    r_d    = zeros(nz+1)
    dd     = zeros(nz+1)
    # initial conditions
    zc     = LinRange(-H+dz/2, 0.0-dz/2, nz  )
    zv     = LinRange(-H     , 0.0,      nz+1)
    ηv     =  η0.*ones(nz+1)
    η_p    =  η0.*ones(nz+1)
    T      =  Ti.*ones(nz+1)
    d      = d_i.*ones(nz+1)
    P      =  P0.*ones(nz+1)
    d_o   .= d
    global m_gbs, n_gbs, n_disl
    # time loop
    for it = 1:nt
        d_o     .= d # uncomment if steady Grain-size is wanted
        iter = 1; err = 2*tol; errd = err
        while (err>tol || errd>tol) && iter<iterMax
            ϵxz     .= 0.5.*diff(Vxe)/dz .+ ϵ_bg
            m_gbs, n_gbs, n_disl = compute_rheol_Tau!(η_p, η_gbs, η_disl, d, Tii2, T, P, ph, R, fact)
            # η_p     .= 0.5.*a0.^(.-1.0./npow).*Eii2.^mpow.*exp.(Q0./(npow.*R.*T))
            ηv      .= rel.*η_p .+ (1.0-rel).*ηv
            τxz     .= fact*ηv.*ϵxz
            Tii2    .= τxz.^2
            errd = compute_GSE_Tau!(dd, d, r_d, d_o, η_gbs, η_disl, Tii2, τxz, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, gam, p, c, dt, R, fact, steady, rel_d, dtau_d)
            Fx      .= diff(τxz)/dz .- ρ.*g.*sin(α)
            dVxdτ   .= 0.5*dVxdτ .+ 0.5.*(Fx .+ dVxdτ*damp) # Crank-Nicolson here
            dτVx    .= dz.*dz./(0.5.*(ηv[1:end-1].+ηv[2:end]))./2.1
            Vxe[2:end-1] .= Vxe[2:end-1] .+ dτVx.*dVxdτ
            Vxe[1]   = -Vxe[2]
            Vxe[end] =  Vxe[end-1]
            if mod(iter, nout)==0
                err  = norm(Fx)/length(Fx)
                @printf("it=%d, iter=%d, err=%2.2e, errd=%2.2e, gse=[%1.3f, %1.3f]\n", it, iter, err, errd, minimum(d)*1e3, maximum(d)*1e3)
            end
            iter += 1
        end
        # Behn et al.
        Txz_Behn, Eii_Behn, η_Behn, dB, Vx_Behn = check_Behn(T, zv, d_i, λ_disl, λ_gbs, Kgg, Qgg, R, gam, p, c, ρ, g, α, dz, fact)
        # visualisation
        p1 = plot(Vxe[2:end-1]*s2yr, zc, linewidth=:3.0, framestyle=:box, xlabel="Vx [m/yr]", ylabel="depth [m]", label="This study")
             plot!(Vx_Behn[1:10:end]*s2yr, zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21", legend=:topleft)
        p2 = plot(log10.(ϵxz), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(ϵxz) [1/s]", label="This study")
             plot!(log10.(Eii_Behn)[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        p3 = plot(1e-6*τxz, zv, linewidth=:3.0, framestyle=:box, xlabel="τxz [MPa]", label="This study")
             plot!(1e-6*Txz_Behn[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        p4 = plot(log10.(ηv), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(η) [Pa.s]", ylabel="depth [m]", label="This study")
             plot!(log10.(η_Behn)[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21", legend=:bottomright)
        p5 = plot(1e3*d, zv, linewidth=:3.0, framestyle=:box, xlabel="d [mm]", label="This study")
             plot!(1e3*dB[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21", legend=:bottomright)
        display(plot(p1, p2, p4, p5, background_color=:transparent, foreground_color=:gray))
    end
    # savefig("ice_column1D_GSE_Tii.png")
    return
end

@time ice_column1D()
