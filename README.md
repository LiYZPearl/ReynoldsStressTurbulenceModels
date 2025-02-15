# ReynoldsStressTurbulenceModels (RSM)
Reynolds stress turbulence closure models beyond the Boussinesq approximation

The RSM in terms of the Wilcox (2006) stress-omega model is implemented with an additional buoyancy production term for simulating multiphase flows. 

The stabilized version of the Wilcox (2006) k-omega turbulence model, with the additional buoyancy production term for two-phase flow modeling e.g. free surface waves, has also been provided. This model is utilized for comparison with the breaking wave cases considered in Li et al. (2021).

Test cases including a turbulent wave boundary layer simulation (Jensen et al. 1989, their Test 13) as well as spilling and plunging breaking wave simulations (Ting & Kirby 1994, 1996) are provided.

The paper regarding the formal analysis and applications of the Wilcox (2006) stress-omega model has been submitted to the Journal of Fluid Mechanics. 

This library is developed for OpenFOAM v1712 to v1912.

## References 

Li, Y., Larsen, B., & Fuhrman, D. (2022). Reynolds stress turbulence modelling of surf zone breaking waves. Journal of Fluid Mechanics, 937, A7. doi:10.1017/jfm.2022.92

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
