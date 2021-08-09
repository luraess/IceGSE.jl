using Plots, Printf, LinearAlgebra, Statistics

const s2yr = 3600*24*365.25

@views av_xy(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xi(A) =  0.5*(A[1:end-1,2:end-1].+A[2:end,2:end-1])
@views av_yi(A) =  0.5*(A[2:end-1,1:end-1].+A[2:end-1,2:end])

@views function fastice2D()
    do_visu = true
    PERIOD  = 0
    # Physics
    T_K    =  273.0
    R      =  8.314
    g      =  9.8
    ρ      =  900.0
    κ      =  2.51
    cp     =  2096.9
    Q0     =  6*1e4
    a0     =  8.75*1e-13
    η0     =  1e15 #1.5e7; # cut-off viscosity
    # Input values
    Lz     = 500.0
    α      = 5.0
    T0     = -20.0
    T_base = -2.0
    npow   = 3
    MELT   = 1  # Melting if MELT = 1.0, none if =0
    Tm     = -0.2 # defined metling temperature °C
    # scales
    T0     = T0 + T_K
    T_base = T_base + T_K
    # T_sc   = npow*R*T0^2/Q0
    # τ_sc   = ρ*cp*T_sc
    # t_sc   = 2^(-npow)*a0^(-1)*τ_sc^(-npow)*exp(Q0/R/T0)
    # L_sc   = sqrt(κ/ρ/cp)*sqrt(t_sc)
    # non-dimensional values
    # Lz     = Lz/L_sc
    # F      = ρ*g*sind(α)*L_sc/τ_sc
    F      = ρ*g*sind(α)
    # theta  = T0/T_sc;
    # q      = q0/κ*L_sc/T_sc
    # q      = q0/κ
    t_t    = 2.2160e9
    dt     = 3.6933e6
    aspect = 10.0
    # Numerics
    nx     = 125 #255
    nz     = 32 #63
    epsi   = 1e-3
    CN     = 0.5
    damp   = 2.0
    relV   = 1.0
    relP   = 0.5
    η_b    = 0.5
    rele   = 1e-1
    niter  = 1e5
    nout   = 2000
    # Pre-processing
    Lx     = aspect*Lz
    mpow   = -(1-1/npow)/2
    # q      = (T_base-T0)/Lz/T_sc
    q      = (T_base-T0)/Lz*κ
    dx     = Lx/(nx-1)
    dz     = Lz/(nz-1)
    xn     = LinRange(-dx/2, Lx+dx/2, nx-1)
    zn     = LinRange(-dz/2, Lz+dz/2, nz-1)
    xc     = LinRange(0, Lx, nx)
    zc     = LinRange(0, Lz, nz)
    # T      = repeat(LinRange((T_base-T0)/T_sc, (T0-T0)/T_sc, nz)', nx, 1)
    T      = repeat(LinRange((T_base-T0), (T0-T0), nz)', nx, 1)
    T_o    = copy(T)
    Etan   = 1e10*ones(nx,nz)
    Etas   = zeros(nx-1,nz-1)
    P      = zeros(nx  ,nz  )
    ∇V     = zeros(nx  ,nz  )
    Exx    = zeros(nx  ,nz  )
    Ezz    = zeros(nx  ,nz  )
    Exz    = zeros(nx-1,nz-1)
    Exzn   = zeros(nx  ,nz  )
    Txx    = zeros(nx  ,nz  )
    Tzz    = zeros(nx  ,nz  )
    Txz    = zeros(nx-1,nz-1)
    EII2   = zeros(nx  ,nz  )
    dVxdt  = zeros(nx-1,nz-2)
    dVzdt  = zeros(nx-2,nz-1)
    ResVx  = zeros(nx-1,nz-2)
    ResVz  = zeros(nx-2,nz-1)
    dtVx   = zeros(nx-1,nz-2)
    dtVz   = zeros(nx-2,nz-1)
    Vx     = zeros(nx+1,nz  )
    Vz     = zeros(nx  ,nz+1)
    dtP    = zeros(nx  ,nz  )
    dTdt   = zeros(nx  ,nz  )
    dTdt_o = zeros(nx  ,nz  )
    Hs     = zeros(nx  ,nz  )
    qxT    = zeros(nx+1,nz  )
    qzT    = zeros(nx  ,nz+1)
    T_res  = zeros(nx  ,nz  )
    T2     = zeros(nx  ,nz  )
    H      = zeros(nx  ,nz  )
    Eta_p  = zeros(nx  ,nz  )
    Eta_it = zeros(nx  ,nz  )
    nt     = cld(t_t, dt)
    time   = 0.0
    iterrorR = 2*epsi
    err0 = 2*epsi
    # Time loop
    for it = 1:nt
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
            Exx     .= diff(Vx, dims=1)./dx .- 1.0/3.0 .* ∇V
            Ezz     .= diff(Vz, dims=2)./dz .- 1.0/3.0 .* ∇V
            Exz     .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dz .+ diff(Vz[:,2:end-1], dims=1)./dx)
            # Stresses
            Txx     .=  2.0 .* Etan .* (Exx .+ η_b.*∇V)
            Tzz     .=  2.0 .* Etan .* (Ezz .+ η_b.*∇V)
            Txz     .=  2.0 .* Etas .* Exz
            Exzn[2:end-1,2:end-1] .= av_xy(Exz)
            if PERIOD==1 Exzn[[1, end], :] = Exzn[[2, end-1], :] end # /!\ uncomment if open box on sides
            Exzn[:,1] .= Exzn[:,2]
            # Residuals
            ResVx   .= diff(Txz, dims=2)./dz .- diff(P[:,2:end-1], dims=1)./dx .+ diff(Txx[:,2:end-1], dims=1)./dx .+ F
            ResVz   .= diff(Txz, dims=1)./dx .- diff(P[2:end-1,:], dims=2)./dz .+ diff(Tzz[2:end-1,:], dims=2)./dz .- F.*cotd(α)
            dVxdt   .= dVxdt.*(1-damp/nx) .+ ResVx
            dVzdt   .= dVzdt.*(1-damp/nz) .+ ResVz
            EII2    .= 0.5.*(Exx.^2 .+ Ezz.^2) .+ Exzn.^2 .+ 1e-18
            qxT[2:end-1,:] .= .-κ .* diff(T, dims=1)./dx
            qzT[:,2:end-1] .= .-κ .* diff(T, dims=2)./dz
            Hs      .= 4.0 .* Etan .* EII2
            # Melting
            T2      .= T .+ (T0 .- T_K)
            H       .= 1.0 .- tanh.( .-(T2 .- Tm)./0.5 )
            dTdt    .= 1.0./ρ./cp .* (.-(diff(qxT, dims=1)./dx .+ diff(qzT, dims=2)./dz) .+ Hs) .* (1.0 .- H.*MELT)
            T_res   .= .-(T .- T_o)./dt .+ (1.0 .- CN).*dTdt .+ CN.*dTdt_o # Crank Nicolson
            # Rheology
            Eta_p   .= 0.5.*a0.^(-1.0 ./ npow) .* EII2.^mpow .* exp.(Q0./(npow.*R.*(T0 .+ T))) # nondim: eta_phys = EII2.^mpow.*exp( -T./(1+T/theta) );
            Eta_it  .= 1.0./(1.0./Eta_p .+ 1.0./η0)
            # Vel update
            Vx[2:end-1,2:end-1] .= Vx[2:end-1,2:end-1] .+ dVxdt.*dtVx
            Vz[2:end-1,2:end-1] .= Vz[2:end-1,2:end-1] .+ dVzdt.*dtVz
            T       .= T .+ dtT.*T_res
            Etan    .= exp.(rele.*log.(Eta_it) .+ (1 .- rele).*log.(Etan))
            Etan[:,[1,end]] .= Etan[:,[2,end-1]]
            Etan[[1,end],:] .= Etan[[2,end-1],:]
            Etas    .= av_xy(Etan)
            # Boundary conditions
            # BC left
            Vx[1,2:end-1]   .= .-Vx[2,2:end-1] .* PERIOD*(-1) # /!\ change sign if open box on sides
            Vz[1,2:end-1]   .=   Vz[2,2:end-1]
            # BC right
            Vx[end,2:end-1] .= .-Vx[end-1,2:end-1] .* PERIOD*(-1) # /!\ change sign if open box on sides
            Vz[end,2:end-1] .=   Vz[end-1,2:end-1]
            # BC bottom
            Vx[:,1]  .= 0.0
            Vz[:,1]  .= .-Vz[:,2]
            qzT[:,1] .= q
            # BC top
            Vz[:,end]       .= Vz[:,end-1] .+ dz.*P[:,end] ./ Etan[:,end]./2.0 .+ dz.*(1.0/3.0 .- η_b).*∇V[:,end]
            Vx[2:end-1,end] .=   Vx[2:end-1,end-1] .- dz./dx.*(Vz[2:end,end-1] .- Vz[1:end-1,end-1]) .+ dz./3.0 .* Txz[:,end-1] ./ Etas[:,end]
            Vx[[1,end],end] .= .-Vx[[2,end-1],end] .* PERIOD*(-1) # /!\ change sign if open box on sides
            T[:,end]        .= 0.0
            # errors
            if iter==1 || iter % nout == 0
                norm_X  = norm(ResVx)/length(ResVx)
                norm_Z  = norm(ResVz)/length(ResVz)
                norm_P  = norm(∇V)/length(∇V)
                # norm_T  = norm(T_res)/length(T_res)
                if (iter==1)  err0 = max(norm_X, norm_Z, norm_P) end
                if (iter>100) iterror = max(norm_X, norm_Z, norm_P) end
                iterrorR = iterror/err0
                @printf("it=%d, iter=%d, err=%1.3e, err/err0=%1.3e \n", it, iter, iterror, iterrorR)
                if do_visu
                    fontsize = 6
                    opts = (yaxis=font(fontsize, "Courier"), xaxis=font(fontsize, "Courier"),
                        ticks=nothing, framestyle=:box, titlefontsize=fontsize, titlefont="Courier", colorbar_title="", c=:davos)
                    p1 = heatmap(xn, zc, Vx[2:end-1,:]' .* s2yr; title="Vx", ylabel="Ly", opts...)
                    p2 = heatmap(xc, zc, T' .+ (T0-T_K); title="T", ylabel="Ly", opts...)
                    p3 = heatmap(xc, zc, log10.(Etan)'; title="log10(η)", xlabel="Lx", ylabel="Ly", opts...)
                    p4 = heatmap(xc, zn, qzT[:,2:end-1]'; title="qzT", opts...)
                    p5 = heatmap(xc, zc, Hs'; title="Hs", opts...)
                    p6 = heatmap(xc, zc, H'; title="H", xlabel="Lx", opts...)
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
