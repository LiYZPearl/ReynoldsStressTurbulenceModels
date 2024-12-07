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

\*---------------------------------------------------------------------------*/

#include "kOmegaWilcox2006SinglePhase.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaWilcox2006SinglePhase<BasicTurbulenceModel>::correctNut()
{
    
   //nut as suggested by Larsen and (Fuhrman 2018)
  this->nut_ = k_/max(omega_, max( Clim_*sqrt(2.0*magSqr(symm(fvc::grad(this->U_)))/Cmu_), lambda2_*beta_/(Cmu_*gamma_)*2.0*magSqr(symm(fvc::grad(this->U_)))/(2.0*magSqr(skew(fvc::grad(this->U_)))+pOmegaSmall_)*omega_)), 

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaWilcox2006SinglePhase<BasicTurbulenceModel>::kOmegaWilcox2006SinglePhase
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta0",
            this->coeffDict_,
            0.0708
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.6
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            0.875
        )
    ),
    sigmad0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmad0",
            this->coeffDict_,
            0.125
        )
    ),
    alphaBS_ //Coefficient for the buoyancy production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaBS",
            this->coeffDict_,
            1.36
        )
    ),
    lambda2_ //New limiter suggested by Larsen and Fuhrman (2018)
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lambda2",
            this->coeffDict_,
            0.05
        )
    ),
    pOmegaSmall_("pOmegaSmall", dimless/(dimTime*dimTime), SMALL),
    pOmega_
    (
            IOobject
            (
                    "pOmega",
                   this-> runTime_.timeName(),
                   this-> mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
            ),
            this->mesh_,
	    dimensionedScalar("dime",dimensionSet(0, 0, 0, 0, 0), 1.0) 
    ),

   fBeta_
   (
	IOobject
	    (
    		"fBeta",
                this-> runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
               IOobject::AUTO_WRITE
	),
        this->mesh_,
	dimless
    ),
   beta_
   (
      IOobject
	(
	   "beta",
            this->runTime_.timeName(),
            this->mesh_,
	    IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
        this->mesh_,
	dimless
    ),
    sigmad_
    (
	IOobject
	  (
		"sigmad",
               this->runTime_.timeName(),
               this->mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
              this->mesh_, dimensionedScalar("zero", dimless, 0.0)
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
     ),
    gField_ //Needed for the buoyancy production term
    (
    IOobject
	(
	"g",
	this->db().time().constant(),
	this->db(),
	IOobject::MUST_READ,
	IOobject::NO_WRITE
	)
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaWilcox2006SinglePhase<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        Clim_.readIfPresent(this->coeffDict());
	alphaBS_.readIfPresent(this->coeffDict());
	lambda2_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaWilcox2006SinglePhase<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }


    // Local references
    const alphaField& alpha = this->alpha_;
    const volVectorField& U = this->U_;
     const rhoField& rho = this->rho_;
     const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;


    fv::options& fvOptions(fv::options::New(this->mesh_));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //fBeta according to Wilcox 2006  - YLi
    volTensorField GradU(fvc::grad(U));
    volTensorField Omij(-skew(GradU));
    volSymmTensorField Sij(symm(GradU));
    volScalarField Chi_ = (Omij & Omij) && Sij /pow((Cmu_*omega_),3);
    volScalarField absChi_ = mag(Chi_);
    fBeta_ = (1.0+85.0*absChi_)/(1.0+100.0*absChi_); //This term should be 1 for 2-D
    beta_ = 0.0708*fBeta_;


    // Calculate eddy viscosity using the stabilizing approach described in Larsen and Fuhrman (2018)
    volScalarField p0 = 2.0*magSqr(symm(GradU)); // 2S_ij S_ij
    volScalarField pOmega= 2.0*magSqr(skew(GradU)); // 2 Omega_ij Omega_ij

    volScalarField omegaTilde = max( omega_, Clim_*sqrt(p0/Cmu_) ); 
    volScalarField omegaTilde2 = max (omegaTilde + this->omegaMin_, lambda2_*(beta_/(Cmu_*gamma_))*p0*omega_/(pOmega+pOmegaSmall_));
    volScalarField nut = k_/omegaTilde2;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(),
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();


   //Crossflow diffusion term to match Wilcox 2006  - YLi
   volVectorField Gradk(fvc::grad(k_));
   volVectorField Gradomega(fvc::grad(omega_));
   volScalarField sigmadCheck_(Gradk & Gradomega); //condition to change sigmad_
  /* forAll(sigmad_,celli)
   {
	if (sigmadCheck_[celli] <= 0.00001)
	{
	sigmad_[celli]=scalar(0.0);
	}else
	{
	sigmad_[celli]=sigmad0_.value(); //scalar(0.125);
	}
    }*/
   sigmad_= pos(sigmadCheck_)*sigmad0_;
   volScalarField CDkOmega = sigmad_/omega_*(Gradk & Gradomega);
   

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*p0
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha*rho*divU, omega_)
      - fvm::Sp(beta_*alpha*rho*omega_, omega_)
      + CDkOmega*alpha*rho //Crossflow diffusion term to match Wilcox 2006 -YLi
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_) 
      - fvm::Sp(Cmu_*alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


     correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
