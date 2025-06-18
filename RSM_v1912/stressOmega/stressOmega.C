/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "stressOmega.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void stressOmega<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
stressOmega<BasicTurbulenceModel>::stressOmega
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
    ReynoldsStress<RASModel<BasicTurbulenceModel>>
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

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
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
    beta0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta0",
            this->coeffDict_,
            0.0708
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.8
        )
    ),
    alphaHat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaHat",
            this->coeffDict_,
            0.7751
        )
    ),
    betaHat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaHat",
            this->coeffDict_,
            0.201
        )
    ),
    gammaHat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gammaHat",
            this->coeffDict_,
            0.5014
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
    alphaB_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaB",
            this->coeffDict_,
            1.356
        )
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
    ),
   k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
     0.5*tr(this->R_)
     //this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
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
bool stressOmega<BasicTurbulenceModel>::read()
{
    if (ReynoldsStress<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        alphaHat_.readIfPresent(this->coeffDict());
        betaHat_.readIfPresent(this->coeffDict());
        gammaHat_.readIfPresent(this->coeffDict());
        sigmad0_.readIfPresent(this->coeffDict());
        alphaB_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void stressOmega<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho1 = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

 
    ReynoldsStress<RASModel<BasicTurbulenceModel>>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField P(-twoSymm(R & gradU));
    volScalarField G(this->GName(), 0.5*mag(tr(P)));

    // Needed coefficients
    // Coefficients for cross diffusion term
    volScalarField CD0 = (fvc::grad(k_) & fvc::grad(omega_)); 
    volScalarField sigmad = pos(CD0)*sigmad0_; // Note: pos(CD)=1 if CD>=0, else pos(CD)=0

    // beta coefficient for dissipation term
    volSymmTensorField S = symm(gradU); // Mean strain rate tensor
    volSymmTensorField Shat = S - 1.0/2.0*tr(gradU)*I; // Strain rate tensor for chiRot
    volTensorField rot = skew(gradU); // Mean rotation rate tensor
    volScalarField chiRot = mag(((rot & rot) && Shat)/pow(betaStar_ * omega_,3));
    volScalarField fBeta = (1.0+85.0*chiRot)/(1.0+100.0*chiRot);
    volScalarField beta_ = beta0_ * fBeta;

    // Needed tensor field
    volSymmTensorField D(-twoSymm(R & gradU.T())); // Anisotropic shear production tensor

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    //create a density field (rho is constant 1 for incompressible)
    const volScalarField& rhoF_ = this->db().objectRegistry::lookupObject<volScalarField>("rho");
    const surfaceScalarField& rhoPhi = this->db().objectRegistry::lookupObject<surfaceScalarField>("rhoPhi");

  
    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rhoF_, omega_)
      + fvm::div(rhoPhi, omega_)
      - fvm::laplacian(alpha*rhoF_*DomegaEff(), omega_)
     ==
        gamma_*alpha*rhoF_*G*omega_/k_
     // - fvm::SuSp(((2.0/3.0)*gamma_)*alpha*rhoF_*divU, omega_) //zero for incompressible flow due to continuity equation
      - fvm::Sp(beta_*alpha*rhoF_*omega_, omega_) //Dissipation
      - fvm::SuSp(-sigmad*rhoF_/omega_*CD0/omega_,omega_) // Cross diffusion
      + fvOptions(alpha, rhoF_, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Correct the trace of the tensorial production to be consistent
    // with the near-wall generation from the wall-functions
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label celli = curPatch.faceCells()[facei];
                P[celli] *= min
                (
                    G[celli]/(0.5*mag(tr(P[celli])) + SMALL),
                    1.0
                );
            }
        }
    }


    //Buoyancy production in Reqn
    volSymmTensorField B = alphaB_*k_/omega_* twoSymm(gField_ * fvc::grad(rhoF_)/rhoF_);

    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rhoF_, R)
      + fvm::div(rhoPhi, R)
      - fvm::laplacian(alpha*rhoF_*DkEff(), R)
      + fvm::Sp(C1_*alpha*rhoF_*betaStar_*omega_, R)
      ==
        alpha*rhoF_*P
      - (2.0/3.0*(1 - C1_)*I)*alpha*rhoF_*betaStar_*omega_*k_
      - alphaHat_*alpha*rhoF_*dev(P)
      - betaHat_*alpha*rhoF_*(D-1.0/3.0*tr(P)*I) // Mean strain effect on pressure-strain correlation term
      - gammaHat_*alpha*rhoF_*k_*dev(S) // Mean strain effect on pressure-strain correlation term
      + fvOptions(alpha, rhoF_, R)
      - alpha*rhoF_*B //Explicit buoyancy production term   
    );  


    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R);

    this->boundNormalStress(R);

    k_ = 0.5*tr(R);

    correctNut();

    // Correct wall shear-stresses when applying wall-functions
    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
