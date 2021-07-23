# Automated generation of user elements for FEAP

1. [How to cite](#how-to-cite)
2. [Questions and user requests](#questions-and-user-requests)
3. [Folder information](#folder-information)
4. [Requirements](#requirements)
5. [How to generate user elements](#how-to-generate-user-elements)
6. [How to use the generated elements](#how-to-use-the-generated-elements)
7. [References](#references)

## How to cite

Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-2021/01 (2021)

[Top](#Automated-generation-of-user-elements-for-FEAP)

## Questions and user requests

If you have any questions related to the repository or would like help to explore new elements and contribute to the repository, this can be raised on the [FEAP user forum](http://feap.berkeley.edu/forum/index.php).

[Top](#Automated-generation-of-user-elements-for-FEAP)

## Folder information

This repository has the files related to the AceGen template that has been developed to automatically generate user elements for [FEAP](http://projects.ce.berkeley.edu/feap/). More information about AceGen is available in Korelc and Wriggers [6]. The folders available in this repository are as below:

- **user-elements**: This folder contains the Mathematica notebooks and the generated FEAP user elements (*elmtXX.f*)

    - *user-elements/thermal*: Steady-state and transient thermal elements. The formulation and AceGen implementation are described in Sec 5.1 of Harish et. al. [1].

    - *user-elements/linelastic*: Quasi-static and dynamic small-strain linear elastic material models. The formulation and AceGen implementation are described in Sec 5.2 of Harish et. al. [1].

    - *user-elements/poroelastic*: Dynamic small-strain poroelastic material models. The formulation and AceGen implementation are described in Sec 5.3 of Harish et. al. [1].

    - *user-elements/hyperelastic*: Large-strain quasi-static hyperelastic material models. The formulation and AceGen implementation are described in Sec 5.4 of Harish et. al. [1].

- **dependencies**: This contains the dependencies for using the Mathematica notebooks and the generated FEAP user elements. 

    - *dependencies/AceGen*: The files in this folder need to be copied to the same folder as the Mathematica notebook to generate the user element. More information can be found in the section [How to generate user elements](#how-to-generate-user-elements) 
    
    - *dependencies/FEAP*: The files in this folder need to be copied into the FEAP folder and compiled prior to usage. More information can be found in the section [How to use the generated elements](#how-to-use-the-generated-elements)

- **benchmarks**: This contains the FEAP input files used in the validation examples. These input files use the standard FEAP elements. These examples are used to compare with the results obtained from AceGen-generated FEAP user elements. The formulation and other information can be found in [2,3] and FEAP user [4] and example [5] manual.

- **examples**: This folder has the FEAP input files that use the AceGen-generated elements. The input files are the same problems as in the *benchmarks* folder.

    - *thermalSteady*: Steady-state thermal conduction in 2-D. The comparison of numerical results with the benchmark problem is described in Sec 5.1.3 and particularly in Fig. 5.7 of Harish et. al. [1].

    - *thermalTransient*: Transient thermal conduction problem in 2-D. The comparison of numerical results with the benchmark problem is described in Sec 5.1.3 and particularly in Fig. 5.8 of Harish et. al. [1].

    - *elasticStatic*: 2-D quasi-static uniaxial tensile test, in small deformation range, on an elastic material. The comparison of numerical results with the benchmark problem is described in Sec 5.2.3 and particularly in Fig. 5.15 of Harish et. al. [1].

    - *elasticDynamic*: 2-D dynamic uniaxial tensile test, in small deformation range, on an elastic material. The comparison of numerical results with the benchmark problem is described in Sec 5.2.3 and particularly in Figs. 5.16 and 5.17 of Harish et. al. [1].

    - *poroelastic*: 2-D Mandel and consolidation problems. The comparison of numerical results with the benchmark problem is described in Sec 5.3.5 and particularly in Figs. 5.24 and 5.26 of Harish et. al. [1].

    - *neoHookean*: 2-D quasi-static uniaxial tensile test, in large deformation range, on a neo-Hookean material. The comparison of numerical results with the benchmark problem is described in Sec 5.4.3 and particularly in Fig. 5.32 of Harish et. al. [1].

- **common**: This folder contains images used in the README files related to this repository.

[Top](#Automated-generation-of-user-elements-for-FEAP)

## Requirements

The following softwares are required to generate and use the user elements.

- [FEAP](http://projects.ce.berkeley.edu/feap) v8.1 or later. A free version is also available for download as [Feappv](http://projects.ce.berkeley.edu/feap/feappv).

- [Mathematica](https://www.wolfram.com/mathematica) v11.0 or later. A 15-day trial is available for download through the [Mathematica](https://www.wolfram.com/mathematica/trial) website

- [AceGen](http://symech.fgg.uni-lj.si) v.7.0 or later. An evaluation versionof AceGen is available for students and educators through the [AceGen](http://symech.fgg.uni-lj.si/Download.htm) website.

[Top](#Automated-generation-of-user-elements-for-FEAP)

## How to generate user elements

- Copy the *SMTFEAP.mf* file from the folder *dependencies/AceGen* into the same folder as the Mathematica notebook. The Mathematica noteboks are available in the folder *user-elements*.

![Copy the splice file](common/images/00_ToCopy.png "Copy the splice file")

- Open the Mathematica notebook and use *Evaluate Notebook* from the *Evaluate* menu to compile the Mathematica notebook as shown below

![Compile the Mathematica Notebook](common/images/01_Compile.png "Compiling the Mathematica Notebook")

- Upon starting to compile, often, Mathematica can ask one to evaluate the initialization cells. It is recommended to select *No* option, as shown below

![Evaluation initialization cells?](common/images/02_CompileQ.png "Evaluation initialization cells?")

- Upon successful compilation, the Mathematica notebook shows a confirmation that includes the element number. 

![Compilation successful](common/images/03_CompileSuccA.png "Compilation successful")

- In addition, upon successful compilation, the user element (*elmtXX.f*) and an *AceShare* folder are created in the same folder as the Mathematica notebook. The number *XX* is controlled by the name provided in the Mathematica notebook.

![User element generated](common/images/03_CompileSuccB.png "User element generated")

- Often, one of the common errors is that the *SMTFEAP.mf* is not in the same folder as the Mathematica notebook. This results in an error that says "*Splice file doesn't exist*" as shown below.

![Compilation fail](common/images/04_CompileFail.png "Compilation fail")

[Top](#Automated-generation-of-user-elements-for-FEAP)

## How to use the generated elements

- Detailed instructions on compiling FEAP

    - On [Windows](https://www.youtube.com/watch?v=7QAh6QvOT6s)

    - On [Linux/Mac](https://www.youtube.com/watch?v=_ohQ__rqq3Y)

- Create a new folder in the FEAP folder. Here, we have created a folder with the name *acegen* as shown below

![Copy a new folder](common/images/00_CreateFolder.png "Copy a new folder")

- Edit the *makefile*, as shown below, to add the name of the new folder (here *acegen*) into the makefile.

![Edit makefile](common/images/01_EditMake.png "Edit makefile")

- Copy all the dependency files from the *dependencies/FEAP* into the newly created folder as shown below.

![Copy dependencies](common/images/02_CopyFEAPDep.png "Copy dependencies")

- Finally copy the *makefile* into the new folder. Here, we have copied the *makefile* from the *FEAP86/modules* to the *acegen* folder as shown below .

![Copy makefile into the new folder](common/images/03_CopyMake.png "Copy makefile into the new folder")

- Compile the FEAP and use the new element generated from AceGen. Check the FEAP input files in the examples folder to learn more about calling the generated elements.

[Top](#Automated-generation-of-user-elements-for-FEAP)

## References

1. Harish, A. B. and Taylor, R. L. and Govindjee, S., "Automated Generation of User Elements (UEL) for FEAP," SEMM Report UCB/SEMM-2021/01 (2021).

2. Zienkiewicz, O., Taylor, R., and Zhu, J., The Finite Element Method: Its Basis and Fundamentals, 7th Edition, Elsevier,Oxford, 2013.

3. Zienkiewicz, O., Taylor, R., and Fox, D., The Finite Element Method for Solid and Structural Mechanics, 7th Edition,Elsevier, Oxford, 2013.

4. Taylor, R., and Govindjee, S., FEAP - A Finite Element Analysis Program, User Manual, University of California,Berkeley.,http://projects.ce.berkeley.edu/feap.

5. Taylor, R., and Govindjee, S., FEAP - A Finite Element Analysis Program, Example Manual, University of California,Berkeley.,http://projects.ce.berkeley.edu/feap.

6. Korelc, J. and Wriggers, P., Automation of Finite Element Methods, Springer International Publishing, Switzerland, 2016

[Top](#Automated-generation-of-user-elements-for-FEAP)
