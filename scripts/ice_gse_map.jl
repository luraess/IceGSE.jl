using Plots, Printf, LinearAlgebra, Statistics, Random, MAT

const s2yr = 3600*24*365.25

LogRange(astart, astop, alength) = 10.0.^(range(astart, stop=astop, length=alength))

@views function compute_rheol_Tau!(Eta, Eta_diff, Eta_basal, Eta_gbs, Eta_disl, η_el, d, Tii, T, P, ph, R, do_elast)
    fact       = 2.0
    T_h        = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
    # Diffusional Rheology (Goldsby and Kohlstedt 2001, eq.4 + Table 6)
    m_diff     = 2
    n_diff     = 1
    F_diff     = (2^((2*n_diff-1)/n_diff))^-1          # simple shear epxeriments (not sure 100%)
    Dv         = 9.1e-4                                # Dv: exponential prefactor[m^2/s]
    Vm         = 1.97e-5                               # Vm: molar volume
    A_diff     = 42.0*Vm./(R.*T_h).*Dv
    Q_diff     = 59.4e3                                # Q1: diffusional activation energy[J/mol]
    # Basal Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    n_basal    = 2.4;
    F_basal    = (2^((n_basal-1)/n_basal)*3^((n_basal+1)/(2*n_basal)))^-1; # pure shear epxeriments
    A_basal    = 5.5e7*10^(-6*n_basal)                 # MPa to Pa conversion
    Q_basal    = 60e3
    # GBS Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    m_gbs      = 1.4
    n_gbs      = 1.8
    F_gbs      = (2.0^((n_gbs-1.0)/n_gbs)*3.0^((n_gbs+1.0)/(2.0*n_gbs)))^-1.0
    A1_gbs     = 3.9e-3*10.0^(-6.0*n_gbs)              # A_gbs: preexponential constant[m^1.4/s*Pa^1.8] in Pa^(-n_gbs).s^(-1).m^(m_gbs)
    Q1_gbs     = 49e3                                  # Q1_gbs: lower activation energy[J/mol]
    Q2_gbs     = 192e3                                 # Q2_gbs: higher activation energy[J/mol]
    Tstar_gbs  = 259.0                                 # Tstar_gbs: activation threshold[K] (-18C)
    # Dislocation Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    n_disl     = 4.0
    F_disl     = (2.0^((n_disl-1.0)/n_disl)*3.0^((n_disl+1.0)/(2.0*n_disl)))^-1.0 # pure shear epxeriments
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
    Eta_diff  .= fact^(n_diff-1) .* B_diff .* Tii.^(1.0-n_diff ).*d.^m_diff
    Eta_basal .= fact^(n_basal-1).* B_basal.* Tii.^(1.0-n_basal)
    Eta_gbs   .= fact^(n_gbs-1)  .* B_gbs  .* Tii.^(1.0-n_gbs  ).*d.^m_gbs
    Eta_disl  .= fact^(n_disl-1) .* B_disl .* Tii.^(1.0-n_disl )
    # viscosity guess
    Eta       .= 1.0./(1.0./Eta_diff .+ 1.0./(Eta_basal .+ Eta_gbs) .+ 1.0./Eta_disl .+ do_elast./η_el)
    return
end

# @views function compute_GSE_Tau!(dd, d, r_d, d_o, Eta_gbs, Eta_disl, Tii, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, γ, p, c, dt, R, steady, rel_d, dtau_d)
#     fact      = 2.0
#     T_h       = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
#     Eii_gbs   = Tii./(fact*Eta_gbs)
#     Eii_disl  = Tii./(fact*Eta_disl)
#     if steady==true
#         Eii   = Eii_gbs .+ Eii_disl
#         β     = Eii_disl ./ Eii; #β[:,end] .= 0.0
#         d_p   = ((Kgg.*exp.(.-Qgg./(R.*T_h)).*(1.0/p).*c.*γ) ./ ((λ_gbs .- β*λ_gbs .+ β.*λ_disl).*Tii.*Eii)).^(1.0./(1.0+p))
#         d    .= rel_d*d_p .+ (1-rel_d)*d
#         ddn   = d_p .- d
#         errd  = mean(abs.(ddn[:,2:end-1]))
#     else
#         # work
#         W_disl = Tii .* Eii_disl # Exx = Eii* Txx/Tii that's what we do
#         W_gbs  = Tii .* Eii_gbs
#         Q_disl   = W_disl
#         Q_gbs    = W_gbs
#         # Reduction rate
#         dred     = .-d.^2 ./c./γ .* (λ_disl.*Q_disl .+ λ_gbs.*Q_gbs)   # d grain reduction eqn. (12)
#         # Growth rate
#         dgg      =  p.^(-1).*d.^(1-p).*Kgg.*exp.(.-Qgg./R./T_h)        # d grain growth eq. (4)
#         # Grain size evolution residual
#         # d       .= d .+ dt.*(dgg .+ dred)
#         # rel_d    = 0.0001
#         # d_p      = d_o .+ dt.*(dgg .+ dred)
#         # d       .= rel_d*d_p .+ (1-rel_d)*d
#         r_d2     = d_o./dt .+ (dgg .+ dred)
#         dd      .= 0.5*dd .+ 0.5*r_d2   # Crank-Nicolson here on "numerics"
#         d       .= (d .+ dtau_d.*dd)./(1.0 .+ dtau_d./dt)
#         r_d     .= .-(d.-d_o)./dt .+ (dgg .+ dred)
#         errd  = mean(abs.(r_d.*dt))
#     end
#     return errd
# end

@views function produce_map()
    do_elast = true
    do_save  = true
    # physics
    R      = 8.314
    p0     = 1e7
    ph     = 7e-8
    T      = 273
    d      = 1e-3
    dt     = 1000*0.117*s2yr
    G      = 2*3.5e7
    η_el   = G*dt
    η0     = 1e15
    # numerics
    niter  = 500
    rel    = 0.5
    nout   = 50
    nT, nE = 200, 210
    Tvec   = LinRange(170, 273, nT)
    EIIvec = LogRange(-16, -3, nE)
    (T, EII) = ([x for x=Tvec, z=EIIvec], [z for x=Tvec, z=EIIvec])
    # initial
    Eta         = η0.*ones(nT, nE)
    Eta_diff    = η0.*ones(nT, nE)
    Eta_basal   = η0.*ones(nT, nE)
    Eta_gbs     = η0.*ones(nT, nE)
    Eta_disl    = η0.*ones(nT, nE)
    TII         =     ones(nT, nE)
    P           = p0.*ones(nT, nE)
    EII_diff    = zeros(nT, nE)
    EII_bas_gbs = zeros(nT, nE)
    EII_disl    = zeros(nT, nE)
    EII_el      = zeros(nT, nE)
    TII_basal   = zeros(nT, nE)
    TII_gbs     = zeros(nT, nE)
    Mech        = zeros(nT, nE)

    for it=1:niter

        TII .= (1-rel) .* TII .+ rel .* (2.0 .* Eta .* EII)

        compute_rheol_Tau!(Eta, Eta_diff, Eta_basal, Eta_gbs, Eta_disl, η_el, d, TII, T, P, ph, R, do_elast)
        
        # compute_GSE_Tau!(dd, d, r_d, d_o, Eta_gbs, Eta_disl, TII, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, γ, p, c, dt, R, steady, rel_d, dtau_d)

        if it % nout == 0
            #Exx  = Txx ./(2.0.*Eta)
            EII_diff    .= TII ./ (2.0 .* Eta_diff)
            EII_bas_gbs .= TII ./ (2.0 .* (Eta_basal .+ Eta_gbs))
            EII_disl    .= TII ./ (2.0 .* Eta_disl)
            EII_el      .= TII ./ (2.0 .* η_el)
            TII_basal   .= 2.0 .* EII .* Eta_basal
            TII_gbs     .= 2.0 .* EII .* Eta_gbs
            # Mechanisms
            Mech .= 0
            if do_elast
                Mech[(EII_diff .> EII_bas_gbs) .& (EII_diff .> EII_disl)    .& (EII_diff .> EII_el)    ] .= 1 #.& (EII_diff .> ϵxz_pl)
                Mech[(EII_bas_gbs .> EII_diff) .& (EII_bas_gbs .> EII_disl) .& (EII_bas_gbs .> EII_el) ] .= 2 #.& (EII_bas_gbs .> ϵxz_pl)
                Mech[(EII_disl .> EII_bas_gbs) .& (EII_disl .> EII_diff)    .& (EII_disl .> EII_el)    ] .= 4 #.& (EII_disl .> ϵxz_pl)
                # Mech[(ϵxz_pl .> EII_diff)      .& (ϵxz_pl .> EII_disl)      .& (ϵxz_pl .> EII_el)      .& (ϵxz_pl .> EII_bas_gbs)] .= 5
                Mech[(Mech.==2) .& (TII_gbs .> TII_basal)] .= 3
            else
                Mech[(EII_diff .> EII_bas_gbs) .& (EII_diff .> EII_disl)    ] .= 1 #.& (EII_diff .> ϵxz_pl)
                Mech[(EII_bas_gbs .> EII_diff) .& (EII_bas_gbs .> EII_disl) ] .= 2 #.& (EII_bas_gbs .> ϵxz_pl)
                Mech[(EII_disl .> EII_bas_gbs) .& (EII_disl .> EII_diff)    ] .= 4 #.& (EII_disl .> ϵxz_pl)
                Mech[(Mech.==2) .& (TII_gbs .> TII_basal)] .= 3
            end
            # display(heatmap(Tvec, log10.(EIIvec), log10.(Eta)'))
            p1 = heatmap(Tvec, log10.(EIIvec), Mech', color=palette(:inferno, 5))
                # contour!(Tvec, log10.(EIIvec), log10.(2.0.*TII./1e6)', color = :black, levels = -5:1:3)
            display(plot(p1))
        end
    end
    if do_save
        !ispath(joinpath(@__DIR__, "../output")) && mkdir(joinpath(@__DIR__, "../output"))
        matwrite(joinpath(@__DIR__, "../output/ice_gse_map.mat"), Dict("T"=> T, "EII"=> EII, "Mech"=> Mech, "Eta"=> Eta, "TII"=> TII, "Tvec"=> Array(Tvec), "EIIvec"=> Array(EIIvec)); compress = true) end
    return
end

produce_map()
