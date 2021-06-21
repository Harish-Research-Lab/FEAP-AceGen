# Linear elastic (Small strain)

## How to cite

Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

## Folder information

This element in this folder is uses the small-strain linear elastic material formulation. The element has a Q4 topology. More detailed information about the formulation is available in Harish et. al. (2021) [1]. Validation of the element is based on comparison with problems presented in FEAP example manual [2] and analytical solutions from Toptan et. al. [3]. Validation with FEAP elements can be found in the *benchmarks* folder of the repository.

The folder contains:

1. AceGen Mathematica notebooks (**thermalSteady.nb** and **thermalTransient.nb**) for the steady-state and transient problems, respectively. The notebook is used to generate the FEAP user element

2. Generated FEAP user elements (**elmt13.f** and **elmt14.f**) for the steady-state and transient problems, respectively.

## References

1. Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

2. <a href=http://projects.ce.berkeley.edu/feap/example_86.pdf>FEAP example manual</a>

3. Toptan, A., Porter, N.W., Hales, J.D., Spencer, B.W., Pilch, M., Williamson, R.L., "Construction of a code verification matrix for heat conduction with finite element code applications," Journal of Verification, Validation and Uncertainty Quantification, vol. 5, pp. 041002 (2020)