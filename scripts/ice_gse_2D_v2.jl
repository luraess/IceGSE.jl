using Plots, Printf, LinearAlgebra, Statistics, Random

const s2yr = 3600*24*365.25

@views av_xy(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xi(A) =  0.5*(A[1:end-1,2:end-1].+A[2:end,2:end-1])
@views av_yi(A) =  0.5*(A[2:end-1,1:end-1].+A[2:end-1,2:end])

@views function compute_EII(Txx, Tzz, Txzn, Eta)
    Exx  = Txx ./(2.0.*Eta)
    Ezz  = Tzz ./(2.0.*Eta)
    Exzn = Txzn./(2.0.*Eta)
    return EII = sqrt.(0.5.*(Exx.^2 .+ Ezz.^2) .+ Exzn.^2)
end

@views function compute_TII(Exx, Ezz, Exzn, Eta)
    Txx  = 2.0.*Eta.*Exx
    Tzz  = 2.0.*Eta.*Ezz
    Txzn = 2.0.*Eta.*Exzn
    return TII = sqrt.(0.5.*(Txx.^2 .+ Tzz.^2) .+ Txzn.^2)
end

@views function initBedTopo(Lx, LzTopo, nx, nz, dz, If)
    Random.seed!(12345)
    Nh  = 10000 # inner parameter, number of harmonics
    sf  = 1     # standard deviation
    dx  = Lx/nx # numerical grid step size in x
    # preprocessing
    c   = sf/sqrt(Nh)
    Yf  = zeros(nx)
    Tmp = zeros(nx)
    for ih = 1:Nh
        ϕ = 2.0*π*rand()
        flag = true
        k = 0.0
        while flag
            k = tan(π*0.5*rand())
            d = (k*k)/(1.0+(k*k))
            if (rand()<d) flag = false  end
        end   
        θ  = acos(1.0-2.0*rand())
        v1 = k*sin(ϕ)*sin(θ)/If
        a  = randn()
        b  = randn()
        for ix=1:nx
            Tmp[ix] = dx*(ix-0.5)*v1
        end
        Yf .= Yf .+ a.*sin.(Tmp) .+ b.*cos.(Tmp)
    end
    Yf .= c.*Yf

    A  = Yf
    A .= A .- minimum(A); A .= A./maximum(A)
    A .= A.*LzTopo

    Bed = zeros(nx, nz)
    for ix = 1:size(Bed,1)
        for iz = 1:size(Bed,2)
            if (iz*dz -2*dz < A[ix]) Bed[ix,iz] = 1 end
        end
    end
    Bed[2:end-1,2:end-1] = Bed[2:end-1,2:end-1] + 1.0/8.0*( diff(diff(Bed[:,2:end-1], dims=1), dims=1) + diff(diff(Bed[2:end-1,:], dims=2), dims=2) )
    Bed[Bed.<0.99] .= 0
    return Bed
end


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

@views function compute_GSE_Tau!(dd, d, r_d, d_o, Eta_gbs, Eta_disl, Tii, Txx, Tzz, Txzn, T, P, ph, λ_disl, λ_gbs, Qgg, Kgg, γ, p, c, dt, R, steady, rel_d, dtau_d)
    fact      = 2.0
    T_h       = T .+ ph.*P  # T_h: melting point adjusted temperature[K]
    Eii_gbs   = Tii./(fact*Eta_gbs)
    Eii_disl  = Tii./(fact*Eta_disl)
    if steady==true
        Eii   = Eii_gbs .+ Eii_disl
        β     = Eii_disl ./ Eii; #β[:,end] .= 0.0
        d_p   = ((Kgg.*exp.(.-Qgg./(R.*T_h)).*(1.0/p).*c.*γ) ./ ((λ_gbs .- β*λ_gbs .+ β.*λ_disl).*Tii.*Eii)).^(1.0./(1.0+p))
        d    .= rel_d*d_p .+ (1-rel_d)*d
        ddn   = d_p .- d
        errd  = mean(abs.(ddn[:,2:end-1]))
    else
        # work
        Wxx_disl = Txx  .* Eii_disl .* Txx  ./ Tii # Exx = Eii* Txx/Tii that's what we do
        Wzz_disl = Tzz  .* Eii_disl .* Tzz  ./ Tii
        Wxz_disl = Txzn .* Eii_disl .* Txzn ./ Tii
        Wxx_gbs  = Txx  .* Eii_gbs  .* Txx  ./ Tii
        Wzz_gbs  = Tzz  .* Eii_gbs  .* Tzz  ./ Tii
        Wxz_gbs  = Txzn .* Eii_gbs  .* Txzn ./ Tii
        Q_disl   = sqrt.(0.5.*(Wxx_disl.^2 .+ Wzz_disl.^2) .+ Wxz_disl.^2)
        Q_gbs    = sqrt.(0.5.*(Wxx_gbs.^2  .+ Wzz_gbs.^2)  .+ Wxz_gbs.^2)
        # Reduction rate
        dred     = .-d.^2 ./c./γ .* (λ_disl.*Q_disl .+ λ_gbs.*Q_gbs)   # d grain reduction eqn. (12)
        # Growth rate
        dgg      =  p.^(-1).*d.^(1-p).*Kgg.*exp.(.-Qgg./R./T_h)        # d grain growth eq. (4)
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
    d[:,end] = d[:,end-1]
    return errd
end

@views function fastice2D()
    do_visu  = true
    BUG1     = 0.0
    OPEN     = 1
    MELT     = 1  # Melting if MELT = 1.0, none if =0
    do_elast = true
    do_glen  = false
    do_gse   = true
    steady   = false
    # Physics
    Lz       =  680.0
    α        =  2.0
    T0       = -11.0
    T_base   = -11.0
    npow     =  3
    Tm       = -0.1 # defined metling temperature °C
    dPdx     = 300.0
    # --
    T_K      =  273.0
    R        =  8.314
    g        =  9.8
    ρ        =  900.0
    κ        =  2.51
    cp       =  2096.9
    Q0       =  6*1e4
    a0       =  8.75*1e-13
    η0       =  1e15 #1.5e7; # cut-off viscosity
    G        =  3.5e7       # Elastic modulus (E/1+nu/2 = 1e9/(1+0.3)/2)
    # GSE
    P0       =  1e7*0       # [Pa]
    ph       =  7e-8        # pressure heating coef [K/Pa]
    d_i      =  7e-4        # grain size [m]
    # GSE - From Behn et al. 2021
    λ_disl   =  0.03        # fitting param
    λ_gbs    =  0.03        # fitting param
    Qgg      =  42e3        # grain growth activation energy
    Kgg      =  1.02e-17    # grain growth prefactor !! units
    γ        =  0.065       # param
    p        =  6.03        # grain growth exponant
    c        =  3.0         # geometric factor
    # scales
    T0       = T0 + T_K
    T_base   = T_base + T_K
    # T_sc     = npow*R*T0^2/Q0
    # τ_sc     = ρ*cp*T_sc
    # t_sc     = 2^(-npow)*a0^(-1)*τ_sc^(-npow)*exp(Q0/R/T0)
    # L_sc     = sqrt(κ/ρ/cp)*sqrt(t_sc)
    # Non-dimensional values
    # Lz       = Lz/L_sc
    # F        = ρ*g*sind(α)*L_sc/τ_sc
    F        = ρ*g*sind(α)
    # theta    = T0/T_sc;
    # q        = q0/κ*L_sc/T_sc
    # q        = q0/κ
    t_t      =  70.2*s2yr
    dt       =  0.117*s2yr #3.6933e6
    aspect   =  3.0
    # Numerics
    nx       =  255
    nz       =  125
    epsi     =  1e-4
    CN       =  0.5
    damp     =  2.0
    relV     =  1.0
    relP     =  0.5
    η_b      =  0.5
    rele     =  2e-2
    niter    =  1e5
    nout     =  1000
    # Numerics gse
    dtau_d   =  0.001*s2yr   # when steady = false
    rel_d    =  1e-2        # when steady = true
    # Pre-processing
    Lx       = aspect*Lz
    mpow     = -(1-1/npow)/2
    # q        = (T_base-T0)/Lz/T_sc
    q        = (T_base-T0)/Lz*κ
    dx, dz   = Lx/(nx-1), Lz/(nz-1)
    xn, zn   = LinRange(-dx/2, Lx+dx/2, nx+1), LinRange(-dz/2, Lz+dz/2, nz+1)
    xc, zc   = LinRange(0, Lx, nx), LinRange(0, Lz, nz)
    # Initial conditions
    nt       = cld(t_t, dt)
    η_el     = G.*dt

    amp      = Lz/7.0
    Bedc     = initBedTopo(Lx, amp, nx  , nz  , dz, 200)
    Bedv     = initBedTopo(Lx, amp, nx+1, nz+1, dz, 200)
    Bedvx    = initBedTopo(Lx, amp, nx+1, nz  , dz, 200)
    Bedvz    = initBedTopo(Lx, amp, nx  , nz+1, dz, 200)
    
    T        = repeat(LinRange((T_base-T0), (T0-T0), nz)', nx, 1) # T      = repeat(LinRange((T_base-T0)/T_sc, (T0-T0)/T_sc, nz)', nx, 1)
    T_o      = copy(T)
    Etan     = 1e10*ones(nx,nz)
    TII2     =      ones(nx,nz)
    d        = d_i.*ones(nx,nz)
    Exzn     = 0e-18*ones(nx,nz)
    P        = cotd(α).*repeat(LinRange(ρ*g*sind.(α)*Lz, 0, nz)', nx, 1)
    # display(heatmap(P')); error("stop")
    # P        = zeros(nx  ,nz  )
    Etas     = zeros(nx-1,nz-1)
    ∇V       = zeros(nx  ,nz  )
    Exx      = zeros(nx  ,nz  )
    Ezz      = zeros(nx  ,nz  )
    Exz      = zeros(nx-1,nz-1)
    Exx_v    = zeros(nx  ,nz  )
    Ezz_v    = zeros(nx  ,nz  )
    Exz_v    = zeros(nx-1,nz-1)
    Txx      = zeros(nx  ,nz  )
    Tzz      = zeros(nx  ,nz  )
    Txz      = zeros(nx-1,nz-1)
    Txzn     = zeros(nx  ,nz  )
    Txx_o    = zeros(nx  ,nz  )
    Tzz_o    = zeros(nx  ,nz  )
    Txz_o    = zeros(nx-1,nz-1)
    Txzn_o   = zeros(nx  ,nz  )
    EII2     = zeros(nx  ,nz  )
    dVxdt    = zeros(nx-1,nz-2)
    dVzdt    = zeros(nx-2,nz-1)
    ResVx    = zeros(nx-1,nz-2)
    ResVz    = zeros(nx-2,nz-1)
    dtVx     = zeros(nx-1,nz-2)
    dtVz     = zeros(nx-2,nz-1)
    Vx       = zeros(nx+1,nz  )
    Vz       = zeros(nx  ,nz+1)
    dtP      = zeros(nx  ,nz  )
    dTdt     = zeros(nx  ,nz  )
    dTdt_o   = zeros(nx  ,nz  )
    Hs       = zeros(nx  ,nz  )
    qxT      = zeros(nx+1,nz  )
    qzT      = zeros(nx  ,nz+1)
    T_res    = zeros(nx  ,nz  )
    T2       = zeros(nx  ,nz  )
    H        = zeros(nx  ,nz  )
    dd       = zeros(nx  ,nz  )
    r_d      = zeros(nx  ,nz  )
    d_o      = zeros(nx  ,nz  )
    Eta_p    = zeros(nx  ,nz  )
    Eta_it   = zeros(nx  ,nz  )
    Eta_diff = zeros(nx  ,nz  )
    Eta_basal= zeros(nx  ,nz  )
    Eta_gbs  = zeros(nx  ,nz  )
    Eta_disl = zeros(nx  ,nz  )
    Tmpx     = zeros(nx-1,nz-2)
    Tmpz     = zeros(nx-2,nz-1)
    time = 0.0; iterrorR = 2*epsi; err0 = 2*epsi; errd = 2*epsi
    # Time loop
    for it = 1:nt
        Txx_o  .= Txx
        Tzz_o  .= Tzz
        Txz_o  .= Txz
        Txzn_o .= Txzn
        d_o    .= d
        T_o    .= T
        dTdt_o .= dTdt
        iterror = 2*epsi
        iter    = 1
        while iterror>epsi && iter<niter # Pseudo-Transient cycles
            # Timesteps
            dtP     .= relP.*4.1 ./ max(nx,nz) .* Etan .* (1.0 .+ η_b)
            dtVx    .= relV.*min(dx,dz)^2 ./ av_xi(Etan) ./ (1.0 + η_b) ./ 4.1
            dtVz    .= relV.*min(dx,dz)^2 ./ av_yi(Etan) ./ (1.0 + η_b) ./ 4.1
            dtT      = 1.0 ./ (1.0 ./ (min(dx,dz)^2 ./κ.*ρ.*cp./4.1) .+ 1.0 ./ dt)
            # Pressure update and strain rates
            ∇V      .= diff(Vx, dims=1)./dx .+ diff(Vz, dims=2)./dz
            P       .= P .- ∇V.*dtP
            Exx_v   .= diff(Vx, dims=1)./dx .- 1.0/3.0 .* ∇V
            Ezz_v   .= diff(Vz, dims=2)./dz .- 1.0/3.0 .* ∇V
            Exz_v   .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dz .+ diff(Vz[:,2:end-1], dims=1)./dx)
            Exx     .= Exx_v .+ do_elast*Txx_o./(2.0.*η_el)
            Ezz     .= Ezz_v .+ do_elast*Tzz_o./(2.0.*η_el)
            Exz     .= Exz_v .+ do_elast*Txz_o./(2.0.*η_el)            
            Exzn[2:end-1,2:end-1] .= av_xy(Exz)
            if OPEN==1 Exzn[[1, end], :] = Exzn[[2, end-1], :] end # /!\ uncomment if open box on sides
            Exzn[:,1] .= Exzn[:,2]
            EII2    .= 0.5.*(Exx.^2 .+ Ezz.^2) .+ Exzn.^2 .+ 1e-18
            # Rheology
            if do_glen
                # Eta_p   .= 0.5.*a0.^(-1.0 ./ npow) .* EII2.^mpow .* exp.(Q0./(npow.*R.*(T0 .+ T))) # nondim: eta_phys = EII2.^mpow.*exp( -T./(1+T/theta) )
                Eta_p   .= 0.5.*a0.^(-1.0) .* TII2.^((1.0-npow)/2.0) .* exp.(Q0./(R.*(T0 .+ T)))
                # Eta_it  .= 1.0./(1.0./Eta_p .+ 1.0./η0)
                Eta_it  .= 1.0./(1.0./Eta_p .+ do_elast./η_el)
            else
                compute_rheol_Tau!(Eta_it, Eta_diff, Eta_basal, Eta_gbs, Eta_disl, η_el, d, sqrt.(TII2), T.+T0, P, ph, R, do_elast)
                # Eta_it  .= 1.0./(1.0./Eta_p .+ 1.0./η0)
            end
            Etan    .= (rele.*Eta_it .+ (1 .- rele).*Etan)
            Etan[:,[1,end]] .= Etan[:,[2,end-1]]
            Etan[[1,end],:] .= Etan[[2,end-1],:]
            Etas    .= av_xy(Etan)
            # Stresses
            Txx     .= 2.0 .* Etan .* (Exx .+ η_b.*∇V)
            Tzz     .= 2.0 .* Etan .* (Ezz .+ η_b.*∇V)
            Txz     .= 2.0 .* Etas .* Exz
            Txzn    .= 2.0 .* Etan .* Exzn
            TII2    .= 0.5.*(Txx.^2 .+ Tzz.^2) .+ Txzn.^2 .+ 1.0
            # if sum(isnan.(TII2))!=0 error("NaN") end
            # if do_gse && (it>1 || iter>4000)
            if do_glen==false && do_gse && it>1
                errd = compute_GSE_Tau!(dd, d, r_d, d_o, Eta_gbs, Eta_disl, sqrt.(TII2), Txx, Tzz, Txzn, T.+T0, P, ph, λ_disl, λ_gbs, Qgg, Kgg, γ, p, c, dt, R, steady, rel_d, dtau_d)
            end
            # Residuals
            ResVx   .= diff(Txz, dims=2)./dz .- diff(P[:,2:end-1], dims=1)./dx .+ diff(Txx[:,2:end-1], dims=1)./dx .+ F .+ dPdx
            ResVz   .= diff(Txz, dims=1)./dx .- diff(P[2:end-1,:], dims=2)./dz .+ diff(Tzz[2:end-1,:], dims=2)./dz .- F.*cotd(α)
            dVxdt   .= dVxdt.*(1-damp/nx) .+ ResVx
            dVzdt   .= dVzdt.*(1-damp/nz) .+ ResVz
            qxT[2:end-1,:] .= .-κ .* diff(T, dims=1)./dx
            qzT[:,2:end-1] .= .-κ .* diff(T, dims=2)./dz
            Hs      .= 4.0 .* Etan .* EII2
            # Melting
            T2      .= T .+ (T0 .- T_K)
            H       .= 1.0 .- tanh.( .-(T2 .- Tm)./0.5 )
            dTdt    .= 1.0./ρ./cp .* (.-(diff(qxT, dims=1)./dx .+ diff(qzT, dims=2)./dz) .+ Hs) .* (1.0 .- H.*MELT)
            T_res   .= .-(T .- T_o)./dt .+ (1.0 .- CN).*dTdt .+ CN.*dTdt_o # Crank Nicolson
            # Vel update
            Vx[2:end-1,2:end-1] .= Vx[2:end-1,2:end-1] .+ dVxdt.*dtVx
            Vz[2:end-1,2:end-1] .= Vz[2:end-1,2:end-1] .+ dVzdt.*dtVz
            T       .= T .+ dtT.*T_res
            # Boundary conditions
            # BC left
            Vx[1,2:end-1]   .= .-Vx[2,2:end-1] .* OPEN*(-1) # /!\ change sign if open box on sides
            Vz[1,2:end-1]   .=   Vz[2,2:end-1]
            # BC right
            Vx[end,2:end-1] .= .-Vx[end-1,2:end-1] .* OPEN*(-1) # /!\ change sign if open box on sides
            Vz[end,2:end-1] .=   Vz[end-1,2:end-1]
            # BC bottom
            Vx[:,1]  .= 0.0
            Vz[:,1]  .= .-Vz[:,2]
            qzT[:,1] .= q*0.0
            # BC top
            Vz[:,end]       .= Vz[:,end-1] .+ dz.*(P[:,end] ./ Etan[:,end]./2.0 .+ (1.0/3.0 .- η_b).*∇V[:,end] .-do_elast.*Tzz_o[:,end]./η_el./2.0)
            Vx[2:end-1,end] .= Vx[2:end-1,end-1] .+ dz.*( .-BUG1/dx.*(Vz[2:end,end-1] .- Vz[1:end-1,end-1]) .+ 1.0/3.0.*Txz[:,end-1]./Etas[:,end] .-do_elast.*Txz_o[:,end]./η_el)
            Vx[[1,end],end] .= .-Vx[[2,end-1],end] .* OPEN*(-1) # /!\ change sign if open box on sides
            T[:,end]        .= 0.0
            # Variable bed
            Tmpx .= Bedvx[2:end-1,2:end-1]
            ResVx[Tmpx.==1] .= 0.0
            Vx[Bedvx.==1]   .= 0.0
            Tmpz .= Bedvz[2:end-1,2:end-1]
            ResVz[Tmpz.==1] .= 0.0
            Vz[Bedvz.==1]   .= 0.0
            qxT[Bedvx.==1]  .= 0.0
            qzT[Bedvz.==1]  .= 0.0
            # errors
            if iter==1 || iter % nout == 0
                norm_X  = norm(ResVx)/length(ResVx)
                norm_Z  = norm(ResVz)/length(ResVz)
                norm_P  = norm(∇V)/length(∇V)
                errη    = norm(abs.(Eta_it.-Etan)./η0)./length(Etan)
                # norm_T  = norm(T_res)/length(T_res)
                if (iter==1)  err0 = max(norm_X, norm_Z, norm_P) end
                if (iter>100) iterror = max(norm_X, norm_Z, norm_P) end
                iterrorR = iterror/err0
                @printf("it=%d, iter=%d, err=%1.3e (err/err0=%1.3e), errη=%1.3e errd=%1.3e, gs=[%1.3f, %1.3f] \n", it, iter, iterror, iterrorR, errη, errd, minimum(d)*1e3, maximum(d)*1e3)
                if do_visu && iter>1
                    fontsize = 6
                    opts = (yaxis=font(fontsize, "Courier"), xaxis=font(fontsize, "Courier"),
                        ticks=nothing, framestyle=:box, titlefontsize=fontsize, titlefont="Courier", colorbar_title="", c=:davos)
                    p1 = heatmap(xn, zc, Vx' .* s2yr; title="Vx", ylabel="Ly", opts...)
                    p2 = heatmap(xc, zc, T' .+ (T0-T_K); title="T", ylabel="Ly", opts...)
                    p3 = heatmap(xc, zc, log10.(Etan)'; title="log10(η)", xlabel="Lx", ylabel="Ly", opts...)
                    p4 = heatmap(xc, zn, qzT'; title="qzT", opts...)
                    p5 = heatmap(xc, zc, Hs'; title="Hs", opts...)
                    p6 = heatmap(xc, zc, 1e3*d'; title="grain size [mm]", xlabel="Lx", opts...)
                    # p6 = heatmap(xc, zc, H'; title="H", xlabel="Lx", opts...)
                    # display(plot(p1, p4, p2, p5, p3, p6, layout=(3,2)))
                    # p1 = heatmap(Eta_p')
                    # p2 = heatmap(Eta_diff')
                    # p3 = heatmap(Eta_basal')
                    # p4 = heatmap(Eta_gbs')
                    # p5 = heatmap(Eta_disl')
                    # p6 = heatmap(TII2')
                    # display(plot(p1, p4, p2, p5, p3, p6, layout=(3,2)))
                    # Strain rates
                    EII_diff    = compute_EII(Txx, Tzz, Txzn, Eta_diff)
                    EII_bas_gbs = compute_EII(Txx, Tzz, Txzn, Eta_basal .+ Eta_gbs)
                    EII_disl    = compute_EII(Txx, Tzz, Txzn, Eta_disl)
                    EII_el      = compute_EII(Txx, Tzz, Txzn, η_el)
                    TII_basal   = compute_TII(Exx, Ezz, Exzn, Eta_basal)
                    TII_gbs     = compute_TII(Exx, Ezz, Exzn, Eta_gbs)
                    # err_τ = (sqrt.(TII2) .- (TII_basal .+ TII_gbs))./mean(sqrt.(TII2))
                    # err_ϵ = abs.((sqrt.(EII2) .- (EII_diff .+ EII_bas_gbs .+ EII_disl .+ do_elast*EII_el))./mean(sqrt.(EII2)))
                    # println("check Tii = $(round(maximum(err_τ), sigdigits=2))")
                    # println("check ϵxz = $(round(maximum(err_ϵ), sigdigits=2))")
                    # Mechanisms
                    Mech = zeros(nx, nz)
                    Mech[(EII_diff .> EII_bas_gbs) .& (EII_diff .> EII_disl)    .& (EII_diff .> EII_el)    ] .= 1 #.& (EII_diff .> ϵxz_pl)
                    Mech[(EII_bas_gbs .> EII_diff) .& (EII_bas_gbs .> EII_disl) .& (EII_bas_gbs .> EII_el) ] .= 2 #.& (EII_bas_gbs .> ϵxz_pl)
                    Mech[(EII_disl .> EII_bas_gbs) .& (EII_disl .> EII_diff)    .& (EII_disl .> EII_el)    ] .= 4 #.& (EII_disl .> ϵxz_pl)
                    # Mech[(ϵxz_pl .> EII_diff)      .& (ϵxz_pl .> EII_disl)      .& (ϵxz_pl .> EII_el)      .& (ϵxz_pl .> EII_bas_gbs)] .= 5
                    Mech[(Mech.==2) .& (TII_gbs .> TII_basal)] .= 3
                    p4 = heatmap(xc, zc, Mech'; title="mechanisms", opts...)
                    display(plot(p1, p4, p2, p5, p3, p6, layout=(3,2)))
                end
            end
            iter+=1
        end
        # if mod(it,20)==1,eval(['save run_TM2D_' int2str(it)]); end
        @printf("it=%d, iter=%d, err=%1.3e \n",it, iter, iterror)
        time = time + dt
    end
    return
end

fastice2D()
