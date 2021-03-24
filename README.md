# ReynoldsStressTurbulenceModels (RSM)
Reynolds stress turbulence closure models beyond the Boussinesq approximation

The RSM in terms of the Wilcox (2006) stress-omega model is implemented with an additional buoyancy production term for simulating multiphase flows. 

The stabilized version of the Wilcox (2006) k-omega turbulence model, with the additional buoyancy production term for two-phase flow modeling e.g. free surface waves, has also been provided. This model is utilized for comparison with the breaking wave cases considered in Li et al. (2021).

Test cases including a turbulent wave boundary layer simulation (Jensen et al. 1989, their Test 13) as well as spilling and plunging breaking wave simulations (Ting & Kirby 1994, 1996) are provided.

The paper regarding the formal analysis and applications of the Wilcox (2006) stress-omega model has been submitted to the Journal of Fluid Mechanics. 

## References (to be updated)

Li, Y., & Fredberg M.B., Larsen B.E. & Fuhrman, D. R. (2021). Reynolds stress turbulence modelling of breaking waves. submitted to J. Fluid Mech.

This library is developed for OpenFOAM v1712 to v1912.

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
