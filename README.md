# Elements for coupled problems

This repository is a collection of elements for <a href="http://projects.ce.berkeley.edu/feap/" target="_blank">FEAP</a> to help address coupled problems. The elements, provided here are generated using <a href = "http://symech.fgg.uni-lj.si/" target="_blank">AceGen</a>. AceGen is available as an add-on for Mathematica and facilitates the automatic generation of user elements using an automatic differentiation approach. A detailed discussion on AceGen can be found in the work of Korelc and Wriggers [1]. The FEAP programmer manual [2] outlines the usage of user element in FEAP and user issues can be addressed through the FEAP Forum [3]. 

## Folder structure

1. **Multiphysics phenomena (Poroelastic, Thermoelastic etc.):** The folder contains elements developed to address the particular phenomenon. Each element type uses a sub-folder of its own where the relates AceGen file, generated FEAP element and FEAP input file for examples are provided.

2. **dependencies:** Include the files required for generating the FEAP element with AceGen and also running the generated elements. More information on its usage is outlined below.

## Poro-elasticity

### Q9/Q4 Taylor-Hood element
This 2-D nine-noded element is developed as a part of the work [4] presented as a part of the Festschrift on the occasion of the 70th birthday of <a href="https://www.ikm.uni-hannover.de/de/wriggers/" target="_blank">Prof. Peter Wriggers</a>.

![Q9/Q4 Taylor-Hood element topology](common/images/Q9Q4-TH_small.png "Q9/Q4 Taylor-Hood element")

The element topology is as shown below. It utilizes a Q9/Q4 Taylor-Hood type interpolation for the u/p dependent variables, respectively.

# Generating the elements with AceGen

1. The mathematica files are available in the folder **Phenomena/Element-Type**. <code>*</code>

2. This is \*\

# Running the elements with FEAP

1. Check out the Youtube videos for help on compiling FEAP on <a href="https://www.youtube.com/watch?v=_ohQ__rqq3Y"> Mac or Linux</a> or <a href="http://www.youtube.com/watch?v=7QAh6QvOT6s">Windows</a>

2. 

# References
1. J. Korelc and P. Wriggers, "Automation of Finite Element Methods," Springer International Publishing, Switzerland (2016)

2. FEAP programmer manual, http://projects.ce.berkeley.edu/feap/pmanual_86.pdf

3. FEAP user forum, http://feap.berkeley.edu/forum/index.php

4. A. B. Harish, R. L. Taylor and S. Govindjee, "A Poroelastic Element for FEAP Using AceGen," Festschrift on the occasion of 70th Birthday of Peter Wriggers (Coming soon)
