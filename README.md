# Ice Grain Size Evolution
Grain size evolution model for composite ice rheology.

This repository contains [Julia] routines implementing Goldsby composite ice rheology to assess discrepancies versus the widely used Glen's flow law. A consistent implementation of the composite rheological Goldsby model requires to resolve grain size evolution.

We implement the wattmeter (e.g. [Austin and Evans (2007)](https://doi.org/10.1130/G23244A.1)) to predict grain size evolution as a function of grain growth and grain reduction parameters for 1D and future 2D and 3D ice flow applications. In our model we reformulate the experimentally parameterized flow law relation using an effective viscosity approach as a function of the second stress or strain-rate invariants.

# Content
* [Benchmark](#benchmark)

# Benchmark
We benchmark our wattmeter-based grain-size evolution implementation against results published by [Behn et al. (2001)](https://doi.org/10.5194/tc-2020-295).

We used the parameters and configuration from the 1D ice column case (Section 3.3, Fig. 6) to reproduce Figure 6:

![](docs/ice_column1D_GSE_Tii.png)


[Julia]: https://julialang.org