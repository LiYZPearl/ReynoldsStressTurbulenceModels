/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::stressOmega

Group
    grpRASTurbulence

Description
    stress-omega model.

    References:
    \verbatim
        Wilcox, D. C. (2006).
        Turbulence Modeling for CFD (Third Edition)
        (Chapter 6, pp. 333-334). La Canada, CA: DCW industries.


        Li, Y., Fredberg M. B., Larsen, B. E. and Fuhrman, D. R.
        Reynolds stress turbulence modelling of breaking waves
	submitted to J. Fluid Mech.
        2021
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        stressOmegaCoeffs
        {
            betaStar  		0.09;  // Equivalent to Cmu for epsilon-based models
            gamma		0.52;  //alpha in Wilcox Eqn
            beta0		0.0708;
            C1			1.8;
            alphaHat		0.7751;
            betaHat		0.201;
            gammaHat		0.5014;
            sigmad0		0.125;
	    alphaK		0.6;  //sigmaStar in Wilcox Eqn
	    alphaOmega		0.5;  //sigma in Wilcox Eqn
	    alphaB		1.176;  //buoyancy production coeff for REqn
            couplingFactor	0.0;    // only for incompressible
        }
    \endverbatim

Author
    Yuzhu (Pearl) Li, Technical University of Denmark

SourceFiles
    stressOmega.C
\*---------------------------------------------------------------------------*/

#ifndef stressOmega_H
#define stressOmega_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "ReynoldsStress.H"
#include "uniformDimensionedFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class stressOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class stressOmega
:
    public ReynoldsStress<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar betaStar_;
            dimensionedScalar beta0_;
            dimensionedScalar gamma_;
            dimensionedScalar alphaK_; 
            dimensionedScalar alphaOmega_; 

            dimensionedScalar C1_;
            dimensionedScalar alphaHat_;
            dimensionedScalar betaHat_;
            dimensionedScalar gammaHat_;

            dimensionedScalar sigmad0_;
            dimensionedScalar alphaB_; 

        // Fields

            volScalarField k_;
            volScalarField omega_;
            uniformDimensionedVectorField gField_; 

    // Protected Member Functions

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("stressOmega");


    // Constructors

        //- Construct from components
        stressOmega
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~stressOmega()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    alphaK_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    alphaOmega_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "stressOmega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
