# Dependencies - FEAP

This folder contains the dependency files required for using the AceGen-generated user elements with FEAP. The files in this folder need to be copied to the FEAP folder and compiled before it can be used in FEAP.

Check out the Youtube tutorial video for compiling FEAP on [Windows](https://www.youtube.com/watch?v=7QAh6QvOT6s) and [Linux/Mac](https://www.youtube.com/watch?v=_ohQ__rqq3Y).

The files in this folder include:

- **SB_Int.f**: This routine is a cover for the quadrature scheme. It initiates different quadrature schemes depending on the element topology.

- **AceFEAPQuad.f**: The routine defines the quadrature schemes (including quadrature points and weights) for different element topologies.

- **SB_ace2nst.f**: The routine transforms the residual and tangent from the AceGen storage style to FEAP storage style.

- **SB_activedofs.f**: The routine sets the active degree of freedom for the element.

- **SB_ua_set.f**: Maps the FEAP storage of elemental unknowns (including unknows, rate of change and the acceleration) to the AceGen style.

- **sms.h**: The header file with declarations required for using AceGen functions and variables. At present, a dummy header file is being used.