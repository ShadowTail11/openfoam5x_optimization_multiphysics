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

#include "tempDepPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(tempDepPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        tempDepPowerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::tempDepPowerLaw::calcNu() const
{
	const volScalarField& T= U_.mesh().lookupObject<volScalarField>("T");
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            (k_-kslope_*(T-Tbase_))*pow
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("VSMALL", dimless, VSMALL)
                ),
                n_.value() - scalar(1)
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::tempDepPowerLaw::tempDepPowerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    tempDepPowerLawCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    k_("k", dimViscosity, tempDepPowerLawCoeffs_),
    n_("n", dimless, tempDepPowerLawCoeffs_),
    kslope_(tempDepPowerLawCoeffs_.lookup("kslope")),
	Tbase_(tempDepPowerLawCoeffs_.lookup("Tbase")),
    nuMin_("nuMin", dimViscosity, tempDepPowerLawCoeffs_),
    nuMax_("nuMax", dimViscosity, tempDepPowerLawCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::tempDepPowerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    tempDepPowerLawCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    tempDepPowerLawCoeffs_.lookup("k") >> k_;
    tempDepPowerLawCoeffs_.lookup("n") >> n_;
    tempDepPowerLawCoeffs_.lookup("kslope") >> kslope_;
	tempDepPowerLawCoeffs_.lookup("Tbase") >> Tbase_;
    tempDepPowerLawCoeffs_.lookup("nuMin") >> nuMin_;
    tempDepPowerLawCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
