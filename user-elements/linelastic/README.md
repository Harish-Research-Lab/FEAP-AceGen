# Linear elastic (Small strain)

## How to cite

Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

## Folder information

The elements in this folder uses the small-strain linear elastic material formulation. The elements have a Q4 topology. Two elements are provided, namely for the steady and the dynamic cases. More detailed information about the formulation is available in Harish et. al. (2021) [1]. This folder contains:

- AceGen Mathematica notebooks (*linelasticStatic.nb** and *linelasticDynamic.nb*) for the steady-state and  dynamic problems, respectively. These notebooks are used to generate the FEAP user element.

- Generated FEAP user elements (*elmt15.f* and *elmt16.f*) for the steady-state and dynamic problems, respectively.

## Input files for validation example

The validation of the developed elements are based on comparison with problems presented in FEAP example manual [2] and the results from the standard FEAP elements. The input files, for FEAP, for the validation examples using the above AceGen-generated elements can be found in the *examples* folder at

- **Steady state:** examples/elasticStatic/IlinElastStatic
- **Dynamic:** examples/elasticDynamic/IlinElastDyn

The input files, for FEAP, for the validation examples using standard FEAP elements can be found in the *benchmarks* folder at

- **Steady state:** benchmarks/IlinElastStaticFeap
- **Dynamic:** benchmarks/IlinElastDynFeap

## References

1. Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

2. <a href=http://projects.ce.berkeley.edu/feap/example_86.pdf>FEAP example manual</a>