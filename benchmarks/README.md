# Benchmark problems

This folder contains the benchmark problems for each of the elements available in this repository. The input files provided in this folder use standard FEAP elements [3]. The input files include:

1. **IthermalSteadyFeap**: Input file for the steady-state thermal problem. More details about the problem is provided in chapter 6 of the FEAP example manual [2] and in Sec. 5.1.3 of Harish et. al. [1]. This problem also has an analytical solution, as provided in Toptan et. al. [5].

2. **IthermalSteadyFeap3D**: Input file for the 3D steady-state thermal problem. More details about the problem is provided in chapter 6 of the FEAP example manual [2] and in Sec. 5.1.4 of Harish et. al. [1]. This problem also has an analytical solution, as provided in Toptan et. al. [5].

3. **IthermalTransientFeap**: Input file for the transient thermal problem. More details about the problem is provided in in chapter 6 of the FEAP example manual [2] and in Sec. 5.1.3 of Harish et. al. [1].

4. **IlinElastStaticFeap**: Input file for the quasi-static uniaxial test. The material model considered is linear elastic in the small deformation range. The problem description can be found in Sec 5.2.3 of Harish et. al. [1].

5. **IlinElastDynFeap**: Input file for the dynamic uniaxial test. The material model considered is linear elastic in the small deformation range. The problem description can be found in Sec 5.2.3 of Harish et. al. [1].

6. **IneoHookeFeap**: Input file for the dynamic uniaxial test. The material model considered is neo-Hookean in the large deformation range. The problem description can be found in Sec 5.4.3 of Harish et. al. [1].

## References

1. Harish, A. B., Taylor, R. L., and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

2. <a href=http://projects.ce.berkeley.edu/feap/example_86.pdf>FEAP example manual</a>

3. <a href=http://projects.ce.berkeley.edu/feap/manual_86.pdf>FEAP user manual</a>

4. Harish, A. B., Taylor, R. L., and Govindjee, S., "A Poroelastic Element for FEAP Using AceGen," Current Trends and Open Problems in Computational Mechanics, Springer (2020)

5. Toptan, A., Porter, N.W., Hales, J.D., Spencer, B.W., Pilch, M., Williamson, R.L., "Construction of a code verification matrix for heat conduction with finite element code applications," Journal of Verification, Validation and Uncertainty Quantification, vol. 5, pp. 041002 (2020)
