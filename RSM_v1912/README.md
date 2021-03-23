# RSM_v1912

For any usage of RSM models please refer to  #to be updated 

	@article{li2021,
	title = {Reynolds stress turbulence modelling of breaking waves},
  	language = {eng},
  	author = {Li, Y., Fredberg M. B., Larsen, B. E. and Fuhrman, D. R.},
	journal = {submitted J. Fluid Mech.},
	VOLUME={},
  	pages = {},
  	year = {2021},
	doi =  {},
	} 

# Description
stressOmega model: 
This is a stress-omega turbulence model with buoyancy production terms for two-phase flow modeling e.g. free surface waves

stressOmegaSinglePhase model:
This is a stress-omega turbulence model implemented for single-phase simulations e.g. steady/wave boundary layer


In addtion, the stablized two-equation models used as a comparsion study in Li et al. (2021) are included as:

kOmegaWilcox2006Stab:
This is a stabilized Wilcox (2006) k-omega turbulence model with stress-limiter and buoyancy production term for two-phase flow modeling e.g. free surface waves

kOmegaWilcox2006SinglePhase:
This is a single-phase Wilcox (2006) k-omega turbulence model with stress-limiter for single-phase simulations e.g. steady/wave boundary layer


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
Wilcox, D. C. 2006 Turbulence Modeling for CFD, 3rd edn. DCW Industries, Inc.


	
