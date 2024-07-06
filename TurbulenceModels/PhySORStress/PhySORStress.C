/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "PhySORStress.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PhySORStress, 0);
addToRunTimeSelectionTable(RASModel, PhySORStress, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void PhySORStress::correctNut()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PhySORStress::PhySORStress
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<incompressible::RASModel>
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

    invScale_
    (
        Switch::lookupOrAddToDict
        (
            "invScale",
            this->coeffDict_,
            true
        )
    ),

    nut_
    (
        IOobject
        (
            "nut",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    bBot_
    (
        IOobject
        (
            "bBot",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("bBot", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor::zero)
    ),

    aBot_
    (
        IOobject
        (
            "aBot",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("aBot", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor::zero)
    ),

    a_
    (
        IOobject
        (
            "a",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("a", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor::zero)
    ),

    Rnew
    (
        IOobject
        (
            "Rnew",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("Rnew", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor::zero)
    )

{
    // bound(k_, this->kMin_);
    // bound(epsilon_, this->epsilonMin_);
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool PhySORStress::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void PhySORStress::correct()
{
    if (!turbulence_)
    {
        return;
    }

    eddyViscosity<incompressible::RASModel>::correct();
}

Foam::tmp<Foam::fvVectorMatrix>
PhySORStress::divDevRhoReff
(
    volVectorField& U
) const
{
    label timeIndex = this->mesh_.time().timeIndex();
    label startTime = this->runTime_.startTimeIndex();

    if(timeIndex - 1 == startTime){
        Info << "Entered Data-driven turbulence predictor" << endl;
        Info << "Calculating Flow Field Variables" << endl;

        volTensorField UGrad(fvc::grad(U));
        volTensorField R(-skew(UGrad));

        volTensorField S(R);
        volSymmTensorField sym = symm(UGrad);
        S.replace(0, sym.component(symmTensor::XX));
        S.replace(1, sym.component(symmTensor::XY));
        S.replace(2, sym.component(symmTensor::XZ));
        S.replace(3, sym.component(symmTensor::XY));
        S.replace(4, sym.component(symmTensor::YY));
        S.replace(5, sym.component(symmTensor::YZ));
        S.replace(6, sym.component(symmTensor::XZ));
        S.replace(7, sym.component(symmTensor::YZ));
        S.replace(8, sym.component(symmTensor::ZZ));

        volScalarField p = this->mesh_.lookupObject<volScalarField>("p");
        volScalarField d = wallDist(this->mesh_).y();

        volScalarField q1 = 0.5*(sqr(tr(UGrad)) - tr(((UGrad) & (UGrad))));
        volScalarField q2 = k_;
        volScalarField q3 = min(((sqrt(k_)*d)/(50*this->nu())), scalar(2));
        volScalarField q4 = U & fvc::grad(p);
        volScalarField q5 = k_/epsilon_;
        volScalarField q6 = sqrt(fvc::grad(p) & fvc::grad(p));
        volScalarField q7 = mag((U*U) && UGrad);
        volScalarField q8 = U & fvc::grad(k_);

        //- Now get the data-driven prediction
        Info << "Executing forward pass of the symbolic expressions" << endl;
        forAll(this->mesh_.C(), cell){       

            this->bBot_[cell].xx() = q2[cell]*sqr(q3[cell])/(sqr(1+2*q3[cell]));
            this->bBot_[cell].xy() = (q5[cell]*q8[cell])/(3*q3[cell]+1);
            this->bBot_[cell].yy() = -q2[cell]/(2+q3[cell]);
            this->bBot_[cell].zz() = q2[cell]*sqr(q3[cell])*sqr(q3[cell]-2);

            this->bBot_[cell].yx() = this->bBot_[cell].xy();
        }

        this->aBot_ = this->alpha_*this->rho_*(this->bBot_);
        Info << "Job Done" << '\n' << endl;
    }
    //- The anisotropic Reynolds stress tensor
    this->a_ = this->aBot_ - this->alpha_*this->rho_*this->nut_*(I & dev(twoSymm(fvc::grad(U))));
    
    this->Rnew = this->R() + this->bBot_;
    return
    (
        - fvm::laplacian(this->alpha_*this->rho_*(this->nu()+this->nut_), U)
        - fvc::div((this->alpha_*this->rho_*(this->nu()+this->nut_))*dev2(T(fvc::grad(U))))
        // - fvc::div((this->alpha_*this->rho_*(this->nu()+this->nut_))*T(fvc::grad(U))))
        + fvc::div(this->aBot_)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
