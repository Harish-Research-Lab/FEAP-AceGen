# Poroelastic element

## How to cite

Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

## Folder information

This element in this folder is based on small-strain poroelastic formulation. This element was developed as a part of the work presented in the Festschrift on the occasion of the 70th birthday of <a href="https://www.ikm.uni-hannover.de/de/wriggers/" target="_blank">Prof. Peter Wriggers</a>. The element has a Q9/Q4 topology as shown below and utilizes a Taylor-Hood type interpolation for the u/p dependent variables, respectively. More detailed information about the formulation is available in Harish et. al. (2020,2021) [1,2]. Detailed discussion on the validation of the element is available in Harish et. al. (2020) [1]

![Q9/Q4 Taylor-Hood element topology](./../../common/images/Q9Q4-TH_small.png "Q9/Q4 Taylor-Hood element")

The folder contains:

1. AceGen Mathematica notebook (**Poroelastic.nb**) which is required to generate the FEAP user element

2. Generated FEAP user element (**elmt16.f**).

## References

1. Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-21/01 (2021)

2. Harish, A. B. and Taylor, R. L. and Govindjee, S., "A Poroelastic Element for FEAP Using AceGen," Current Trends and Open Problems in Computational Mechanics, Springer (2020)