
@views function compare()
    R   = 8.314
    Qgg = 42e3
    T   = 253
    n   = 6

    Kgg_Behn_Tab1 = 1.02e-8   # m^n s^-1
    Kgg_Behn_code = 1e-17     # m^n s^-1
    # Kgg_Behn_Tab1 = (1.02e-8.^(1/n)*1e-3).^n   # m^n s^-1  # not fixing it

    K1  = Kgg_Behn_Tab1*exp(-Qgg/R/T) # m^n s^-1
    K2  = Kgg_Behn_code*exp(-Qgg/R/T) # m^n s^-1
    
    K_Azuma_AL4 = 5.37e-8                   # mm^n s^-1
    Kgg_Azuma   = K_Azuma_AL4/exp(-Qgg/R/T) # mm^n s^-1

    K_Azuma_AL4_m = (5.37e-8.^(1/n)*1e-3).^n    # m^n s^-1
    Kgg_Azuma_m   = K_Azuma_AL4_m/exp(-Qgg/R/T) # m^n s^-1

    println("K   Behn Tab1: $(round(K1, sigdigits=3)) m^n s^-1")
    println("Kgg Behn Tab1: $(round(Kgg_Behn_Tab1, sigdigits=3)) m^n s^-1")
    println("--")
    println("K   Behn code: $(round(K2, sigdigits=3)) m^n s^-1")
    println("Kgg Behn code: $(round(Kgg_Behn_code, sigdigits=3)) m^n s^-1")
    println("--")
    println("K   Azuma    : $(round(K_Azuma_AL4_m, sigdigits=3)) m^n s^-1")
    println("Kgg Azuma    : $(round(Kgg_Azuma_m,   sigdigits=3)) m^n s^-1")
    return
end

compare()

# K   Behn Tab1: 2.17e-17 m^n s^-1
# Kgg Behn Tab1: 1.02e-8 m^n s^-1
# --
# K   Behn code: 2.13e-26 m^n s^-1
# Kgg Behn code: 1.0e-17 m^n s^-1
# --
# K   Azuma    : 5.37e-26 m^n s^-1
# Kgg Azuma    : 2.52e-17 m^n s^-1