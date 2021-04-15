# Ice grain size evolution
Grain size evolution model for composite ice rheology.

This repository contains [Julia] routines implementing Goldsby composite ice rheology to assess discrepancies versus the widely used Glen's flow law. A consistent implementation of the composite rheological Goldsby model requires to resolve grain size evolution.

We implement the wattmeter (e.g. [Austin and Evans (2007)](https://doi.org/10.1130/G23244A.1)) to predict grain size evolution as a function of grain growth and grain reduction parameters for 1D and future 2D and 3D ice flow applications. In our model we reformulate the experimentally parameterized flow law relation using an effective viscosity approach as a function of the second stress or strain-rate invariants.

# Content
* [Benchmark](#benchmark)
* [References](#references)

# Benchmark
We benchmark our wattmeter-based grain-size evolution implementation against results published by [Behn et al. (2001)](https://doi.org/10.5194/tc-2020-295).

## Codes
The benchmark folder contains following 1D codes for the 1D ice column configuration (Section 3.3):
- [`ice_column1D_lin.jl`](benchmark/ice_column1D_lin.jl) a linear viscous version
- [`ice_column1D_powlawEii.jl`](benchmark/ice_column1D_powlawEii.jl) a power-law Glen version with strain-rate dependent viscosity η(ε̇II)
- [`ice_column1D_powlawEii_Tii.jl`](benchmark/ice_column1D_powlawEii_Tii.jl) a power-law Glen version with stress dependent viscosity η(τII)
- [`ice_column1D_fixGS_Eii.jl`](benchmark/ice_column1D_fixGS_Eii.jl) a composite rheology (GBS and dislocation creep) version with strain-rate dependent viscosity η(ε̇II) and fixed grain size
- [`ice_column1D_fixGS_Tii.jl`](benchmark/ice_column1D_fixGS_Tii.jl) a composite rheology (GBS and dislocation creep) version with stress dependent viscosity η(τII) and fixed grain size
- [`ice_column1D_GSE_Tii.jl`](benchmark/ice_column1D_GSE_Tii.jl) a composite rheology (GBS and dislocation creep) version with stress dependent viscosity η(τII) and grain size evolution

The symbolic math derivation of "experimental to model" formulation is done in the [`FlowLawConversionFromLab2.ipynb`](benchmark/FlowLawConversionFromLab2.ipynb) notebook.

## Results
We used the parameters and configuration from the 1D ice column case (Section 3.3, Fig. 6) to reproduce Figure 6 using the [`ice_column1D_GSE_Tii.jl`](benchmark/ice_column1D_GSE_Tii.jl) code:

![](docs/ice_column1D_GSE_Tii.png)


# References
[Austin, N. J., Evans, B. (2007). Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35 (4): 343–346, https://doi.org/10.1130/G23244A.1.](https://doi.org/10.1130/G23244A.1)

[Behn, M. D., Goldsby, D. L., and Hirth, G. (2020). The role of grain-size evolution on the rheology of ice: Implications for reconciling laboratory creep data and the Glen flow law, The Cryosphere Discuss. \[preprint\], https://doi.org/10.5194/tc-2020-295, in review.](https://doi.org/10.5194/tc-2020-295)


[Julia]: https://julialang.org