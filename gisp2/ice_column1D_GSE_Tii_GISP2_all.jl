# Solve: 0 = eta * dVx/dz - rhog(H-z)*sin(alpha)
# See https://doi.org/10.5194/tc-2020-295
using Plots, Printf, LinearAlgebra, Statistics, JLD, Interpolations

const s2yr = 3600*24*365.25

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

@views function load_data!(T, dat_file, A, dat_file2, skip)
    data   = load(joinpath(@__DIR__, dat_file))
    data2  = load(joinpath(@__DIR__, dat_file2))
    Depth  = data["depth"]
    Temp   = data["temp"]
    Depth2 = data2["depth"]
    Age    = data2["age"]
    Depth  = [0.0;    Depth[1:end-skip]]
    Temp   = [Temp[1]; Temp[1:end-skip]] .+ 273.15
    zv     = LinRange(Depth[1], Depth[end], length(T))
    itp1   = interpolate((Depth[:],), Temp[:], Gridded(Linear()))
    T     .= reverse(itp1(zv))
    Depth2 = [Depth[1]; Depth2; Depth[end]]
    Age    = [1       ; Age   ;   Age[end]]
    # @show(Depth2)
    itp2   = interpolate((Depth2[:],), Age[:], Gridded(Linear()))
    A     .= reverse(itp2(zv))
    return zv[1], zv[end]
end

@views function compute_rheol_Tau!(η, Eii, η_pl, η_el, η_diff, η_basal, η_gbs, η_disl, d, Tii, T, P, ph, R, Coh, ϕ, fact)
    T_h        = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
    # Diffusional Rheology (Goldsby and Kohlstedt 2001, eq.4 + Table 6)
    m_diff     = 2
    n_diff     = 1
    F_diff     = 1.0#(2^((2*n_diff-1)/n_diff))^-1          # simple shear epxeriments (not sure 100%)
    Dv         = 9.1e-4                                # Dv: exponential prefactor[m^2/s]
    Vm         = 1.97e-5                               # Vm: molar volume
    A_diff     = 42.0*Vm./(R.*T_h).*Dv
    Q_diff     = 59.4e3                                # Q1: diffusional activation energy[J/mol]
    # Basal Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    n_basal    = 2.4;
    F_basal    = 1.0#(2^((n_basal-1)/n_basal)*3^((n_basal+1)/(2*n_basal)))^-1; % pure shear epxeriments
    A_basal    = 5.5e7*10^(-6*n_basal)                 # MPa to Pa conversion
    Q_basal    = 60e3
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
    B_diff     = F_diff .* A_diff .^(-1.0) .* exp.( Q_diff  ./(R.*T_h))
    B_basal    = F_basal.* A_basal.^(-1.0) .* exp.( Q_basal ./(R.*T_h))
    B_gbs      = F_gbs  .* As_gbs .^(-1.0) .* exp.((Q_gbs_h ./R).*((1.0./T_h).-(1.0./Tstar_gbs_h ))) # B_gbs      = A1_gbs .^(-1.0/n_gbs)   .* exp.((Q1_gbs ./(n_gbs  .*R .* T_h )))
    B_disl     = F_disl .* As_disl.^(-1.0) .* exp.((Q_disl_h./R).*((1.0./T_h).-(1.0./Tstar_disl_h))) # B_disl     = A1_disl.^(-1.0/n_disl)  .* exp.((Q1_disl./(n_disl .*R .* T_h )))
    # Isolated viscosities
    η_diff    .= fact^(n_diff-1) .* B_diff .* Tii.^(1.0-n_diff ).*d.^m_diff
    η_basal   .= fact^(n_basal-1).* B_basal.* Tii.^(1.0-n_basal)
    η_gbs     .= fact^(n_gbs-1)  .* B_gbs  .* Tii.^(1.0-n_gbs  ).*d.^m_gbs
    η_disl    .= fact^(n_disl-1) .* B_disl .* Tii.^(1.0-n_disl )
    η_pl      .= (Coh*cos(ϕ) .+ P.*sin(ϕ))./fact./Eii
    # viscosity guess
    η         .= 1.0./(1.0./η_diff .+ 1.0./(η_basal .+ η_gbs) .+ 1.0./η_disl .+ 1.0./η_el .+ 0.0./η_pl)
    return
end
   
@views function compute_GSE_Tau!(dd, d, r_d, d_o, η_gbs, η_disl, Tii, τxz, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, gam, p, c, dt, R, fact, steady, rel_d, dtau_d)
    T_h        = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
    Eii_gbs    = Tii./(fact*η_gbs)
    Eii_disl   = Tii./(fact*η_disl)
    if steady==true
        Eii   = Eii_gbs + Eii_disl
        β     = Eii_disl./Eii; β[end] = 0.0
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
        d       .= (d .+ dtau_d.*dd)./(1.0 .+ dtau_d./dt)
        r_d     .= .-(d.-d_o)./dt .+ (dgg .+ dred)
        errd  = mean(abs.(r_d.*dt))
    end
    return errd
end

@views function ice_column1D()
    steady    = false
    # physics
    fact      = 2.0
    # H         = 3e3         # domain height [m]
    α         = 0.11*π/180  # slope angle   [rad]
    ρ         = 930         # ice density   [kg/m^3]
    g         = -9.81       # grav. accel.  [m/s^2]
    η0        = 1e15        # ice viscosity [Pa.s]
    Ti        = 253         # ice temp.     [K]
    R         = 8.314       # gas constant
    Q0        = 6*1e4       # activation energy
    a0        = 8.75*1e-13  # prefactor
    npow      = 3.0         # power-law exponent
    mpow      = -(1.0-1.0/npow)/2.0
    G         = 3.5e9       # Elastic modulus (E/1+nu/2 = 1e9/(1+0.3)/2)
    Coh       = 6e6         # Von Mises yield
    ϕ         = 0*π/180    # friction angel
    P0        = 1e7*0       # [Pa]
    ph        = 7e-8        # pressure heating coef [K/Pa]
    d_i       = 2e-3        # grain size [m]
    # dt        = 2e5*s2yr   # physical time step value in seconds
    nt        = 1          # number of time steps
    # GSE - From Behn et al. 2021
    λ_disl    = 0.03        # fitting param
    λ_gbs     = 0.03        # fitting param
    Qgg       = 42e3        # grain growth activation energy
    Kgg       = 1.02e-17    # grain growth prefactor !! units
    gam       = 0.065       # param
    p         = 6.03        # grain growth exponant
    c         = 3.0         # geometric factor
    # numerics
    ϵ_bg      = 1e-16       # background strain-rate [1/s]
    nz        = 300
    rel       = 1e-3
    nout      = 1e3
    tol       = 1e-5
    iterMax   = 5e5
    # numerics gse
    dtau_d    = 10.0*s2yr   # when steady = false
    rel_d     = 1e-1        # when steady = true
    # load and interpolate data
    T         =  zeros(nz+1)
    dt        =  zeros(nz+1)
    ztop, zbot = load_data!(T, "./data_gisp2.jld", dt, "./age_gisp2.jld", 0)
    dt       .= dt*s2yr
    H         = zbot - ztop
    # AL4 Azuma with T fit
    # p         = 6.9
    # k_A       = (5.37e-8.^(1/p)*1e-3).^p
    # Kgg       = k_A./exp.(.-Qgg./R./T) # m^n s^-1
    # derived numerics
    dz        = H/nz
    damp      = 1.0-0.1/nz
    # array initialisation
    d_o       = zeros(nz+1)
    Vxe       = zeros(nz+2)
    ϵxz       = zeros(nz+1)
    ϵxz1      = zeros(nz+1)
    τxz       = zeros(nz+1)
    τxz_o     = zeros(nz+1)
    Tii       =  ones(nz+1)
    Eii       = zeros(nz+1)
    dVxdτ     = zeros(nz)
    dτVx      = zeros(nz)
    Fx        = zeros(nz)
    η_diff    = zeros(nz+1)
    η_basal   = zeros(nz+1)
    η_gbs     = zeros(nz+1)
    η_disl    = zeros(nz+1)
    r_d       = zeros(nz+1)
    dd        = zeros(nz+1)
    # initial conditions
    zc        = LinRange(-zbot+dz/2, ztop-dz/2, nz  )
    zv        = LinRange(-zbot     , ztop     , nz+1)
    gisp2     = load(joinpath(@__DIR__, "gisp_data.jld"))["gisp2_gs"]
    # p1 = plot(T, zv, legend=false, linewidth=:3.0, framestyle=:box, xlabel="T [K]", ylabel="depth [m]")
    # p2 = plot(dt/s2yr/1e3, zv, legend=false, linewidth=:3.0, framestyle=:box, xlabel="Age [kyrs]")
    # p3 = plot(gisp2[:,2], .-gisp2[:,1], legend=false, linewidth=:0.0, markershape=:circle, markersize=5, framestyle=:box, xlabel="grain size [mm]")
    # display(plot(p1, p2, p3, layout=(1,3), background_color=:transparent, foreground_color=:gray))
    # savefig("ice_column1D_GSE_Tii_gisp2_initial.png")
    # error("Display initial conditions")
    ηv        =  η0.*ones(nz+1)
    η_p       =  η0.*ones(nz+1)
    η_el      =  G.*dt
    η_pl      =      ones(nz+1)
    d         = d_i.*ones(nz+1)
    P         =  P0.*ones(nz+1)
    d_o      .= d
    # time loop
    for it = 1:nt
        d_o   .= d
        τxz_o .= τxz
        iter = 1; err = 2*tol
        while err>tol && iter<iterMax
            ϵxz     .= 0.5.*diff(Vxe)/dz .+ ϵ_bg
            ϵxz1    .= ϵxz .+ τxz_o./(fact.*η_el)
            Eii     .= sqrt.(ϵxz.^2) # viscous only invariant for η_pl calculation
            compute_rheol_Tau!(η_p, Eii, η_pl, η_el, η_diff, η_basal, η_gbs, η_disl, d, Tii, T, P, ph, R, Coh, ϕ, fact)
            ηv      .= rel.*η_p .+ (1.0-rel).*ηv
            τxz     .= fact*ηv.*ϵxz1
            Tii     .= sqrt.(τxz.^2)
            errd = compute_GSE_Tau!(dd, d, r_d, d_o, η_gbs, η_disl, Tii, τxz, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, gam, p, c, dt, R, fact, steady, rel_d, dtau_d)
            Fx      .= diff(τxz)/dz .- ρ.*g.*sin(α)
            dVxdτ   .= 0.5*dVxdτ .+ 0.5.*(Fx .+ dVxdτ*damp) # Crank-Nicolson here
            dτVx    .= dz.*dz./(0.5.*(ηv[1:end-1].+ηv[2:end]))./2.1
            Vxe[2:end-1] .= Vxe[2:end-1] .+ dτVx.*dVxdτ
            Vxe[1]   = -Vxe[2]
            Vxe[end] =  Vxe[end-1]
            if mod(iter, nout)==0
                errF = norm(Fx)/length(Fx)
                errη = norm(abs.(η_p.-ηv)./η_p)./length(ηv)
                err  = maximum([errF, errη, errd])
                @printf("it=%d, iter=%d, errF=%2.2e, errη=%2.2e, errd=%2.2e, gse=[%1.3f, %1.3f]\n", it, iter, errF, errη, errd, minimum(d)*1e3, maximum(d)*1e3)
            end
            iter += 1
        end
        # Behn et al.
        # Txz_Behn, Eii_Behn, η_Behn, dB, Vx_Behn = check_Behn(T, zv, d_i, λ_disl, λ_gbs, Kgg, Qgg, R, gam, p, c, ρ, g, α, dz, fact)
        data_Behn  = load(joinpath(@__DIR__, "behn_gisp_fit.jld"))["behn_gisp2_fit"]
        # Strain rates
        ϵxz_diff    = τxz./(fact.*η_diff)
        ϵxz_bas_gbs = τxz./(fact.*(η_basal .+ η_gbs))
        ϵxz_disl    = τxz./(fact.*η_disl)
        ϵxz_el      = τxz./(fact.*η_el)
        ϵxz_pl      = τxz./(fact.*η_pl)
        τxz_basal   = fact.*η_basal.*ϵxz_bas_gbs
        τxz_gbs     = fact.*η_gbs.*ϵxz_bas_gbs
        println(maximum( (τxz .- (τxz_basal .+ τxz_gbs))./mean(τxz) ))
        println(maximum( (ϵxz1 .- (ϵxz_diff .+ ϵxz_bas_gbs .+ ϵxz_disl .+ ϵxz_el .+ 0*ϵxz_pl)) ))
        # p1 = plot(τxz[1:10:end],zv[1:10:end],label="τxz", linewidth=:0.0, markershape=:circle, markersize=4); plot!(τxz_basal,zv,label="τxz_basal"); plot!(τxz_gbs,zv,label="τxz_gbs"); 
        # display(plot(p1)); error("stop")
        Eii_diff    = sqrt.(ϵxz_diff.^2)
        Eii_bas_gbs = sqrt.(ϵxz_bas_gbs.^2)
        Eii_disl    = sqrt.(ϵxz_disl.^2)
        Eii_el      = sqrt.(ϵxz_el.^2)
        Eii_pl      = sqrt.(ϵxz_pl.^2)
        # visualisation
        # p1 = plot(Vxe[2:end-1]*s2yr, zc, linewidth=:3.0, framestyle=:box, xlabel="Vx [m/yr]", ylabel="depth [m]", label="This study", legend=:topleft, foreground_color_legend = nothing)
        #     #  plot!(Vx_Behn[1:10:end]*s2yr, zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        # p2 = plot(log10.(ϵxz), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(ϵxz) [1/s]", label="This study", legend=:bottomleft, foreground_color_legend = nothing)
        #     #  plot!(log10.(Eii_Behn)[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        # p3 = plot(1e-6*τxz, zv, linewidth=:3.0, framestyle=:box, xlabel="τxz [MPa]", label="This study")
        #     #  plot!(1e-6*Txz_Behn[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        # p4 = plot(log10.(ηv), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(η) [Pa.s]", ylabel="depth [m]", label="This study", legend=:bottomright, foreground_color_legend = nothing)
        #     #  plot!(log10.(η_Behn)[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        # p5 = plot(1e3*d, zv, linewidth=:3.0, framestyle=:box, xlabel="d [mm]", label="This study", legend=:topright, foreground_color_legend = nothing)
        #      plot!(data_Behn[:,2], .-data_Behn[:,1], linewidth=:3.0, label="Behn et al. 2021")
        #      plot!(gisp2[:,2], .-gisp2[:,1], linewidth=:0.0, markershape=:circle, markersize=3, framestyle=:box,label="Gow et al. 1997")
        #     #  plot!(1e3*dB[1:10:end], zv[1:10:end], linewidth=:0.0, markershape=:circle, markersize=4, label="Behn21")
        # display(plot(p1, p2, p4, p5, background_color=:transparent, foreground_color=:gray))
        p1 = plot(log10.(ϵxz_diff), zv, linewidth=:3.0, framestyle=:box, xlabel="log10(ϵxz) [1/s]", ylabel="depth [m]", label="ϵxz_diff", legend=:bottomleft, foreground_color_legend = nothing)
            plot!(log10.(ϵxz_bas_gbs), zv, linewidth=:3.0, thickness_scaling = 1, label="ϵxz_bas_gbs")
            plot!(log10.(ϵxz_disl), zv, linewidth=:3.0, thickness_scaling = 1, label="ϵxz_disl")
            plot!(log10.(ϵxz_el), zv, linewidth=:3.0, label="ϵxz_el")
            plot!(log10.(ϵxz_pl), zv, linewidth=:3.0, label="ϵxz_pl")
        # p1 = plot(log10.(abs.(η_p.-ηv)./η_p), zv, linewidth=:3.0, label="ϵrr_eta")
        display(plot(p1, background_color=:transparent, foreground_color=:gray))
        # savefig("ice_column1D_GSE_Tii_gisp2_exz.png")
    end
    # savefig("ice_column1D_GSE_Tii_gisp2_all.png")
    return
end

@time ice_column1D()
