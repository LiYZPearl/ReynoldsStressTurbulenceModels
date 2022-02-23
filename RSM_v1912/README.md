# RSM_v1912

For any usage of RSM models please refer to  #to be updated 

	@article{li2022,
	title = {Reynolds stress turbulence modelling of surf zone breaking waves},
  	language = {eng},
  	author = {Li, Y., Larsen, B. E., and Fuhrman, D. R.},
	journal = {J. Fluid Mech.},
	VOLUME={937},
  	pages = {A7},
  	year = {2022},
	doi =  {10.1017/jfm.2022.92},
	} 

# Description

stressOmega model: 

This is the Wilcox (2006) stress-omega turbulence model, with additional buoyancy production terms added (as described in Li et al. 2021), for two-phase flow modeling of e.g. free surface waves.

stressOmegaSinglePhase model: 

This is the Wilcox (2006) stress-omega turbulence model implemented for single-phase simulations involving e.g. steady/wave boundary layers.

kOmegaWilcox2006Stab: 

This is the stabilized (as detailed in Larsen and Fuhrman, 2018) version of the Wilcox (2006) k-omega turbulence model, with the additional buoyancy production term for two-phase flow modeling e.g. free surface waves.  This model is utilized for comparison with the breaking wave cases considered in Li et al. (2021).

kOmegaWilcox2006SinglePhase: 

This is the Wilcox (2006) k-omega turbulence model for single-phase simulations involving e.g. steady/wave boundary layers.  This model is also utilized for comparison with the wave boundary layer case considered in Li et al. (2021).



## Installation

Download the repository

    git clone https://github.com/LiYZPearl/ReynoldsStressTurbulenceModels 

Create folder for turbulence model (if the folders already exist skip this part)

	mkdir -p $WM_PROJECT_USER_DIR/src/

Move the folder to the user source code

	mv ReynoldsStressTurbulenceModels $WM_PROJECT_USER_DIR/src/
	
Go to the directory and compile the turbulence models

	cd $WM_PROJECT_USER_DIR/src/ReynoldsStressTurbulenceModels/RSM_v1912
	
	wmake libso	
	
	
## Usage
Include the libary of the stabilized turbulence models in the system/controlDict folder

	libs
	(
    	"libRSM.so"
	);

Change the constant/turbulence


	simulationType  RAS;

	RAS
	{
	RASModel        stressOmega;
	//RASModel        stressOmegaSinglePhase;
	

	turbulence      on;

	printCoeffs     on;

	}


## References

Larsen, B. E., & Fuhrman, D. R. (2018). On the over-production of turbulence beneath surface waves in Reynolds-averaged Navier–Stokes models. J. Fluid Mech, 853, 419-460.

Li, Y., & Fredberg M.B., Larsen B.E. & Fuhrman, D. R. (2021). Reynolds stress turbulence modelling of breaking waves. submitted to J. Fluid Mech.

Wilcox, D. C. 2006 Turbulence Modeling for CFD, 3rd edn. DCW Industries, Inc.


	
