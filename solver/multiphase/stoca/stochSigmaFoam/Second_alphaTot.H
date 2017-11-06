#include "alphaControls.H"


if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alphaAver, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
{
    word alphaScheme("div(phi,alpha)");

    surfaceScalarField phiAlpha
    (
        fvc::flux
        (
            phi,
            alphaAver,
            alphaScheme
        )
    );

    MULES::explicitSolve(alphaAver, phi, phiAlpha, 1, 0);

    rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alphaAver.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alphaAver.name() <<") = " << min(alphaAver).value()
        << "  Max(" << alphaAver.name() <<") = " << max(alphaAver).value()
        << endl;
}
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
{
    word alphaScheme("div(phi,alpha)");

    surfaceScalarField phiAlpha
    (
        fvc::flux
        (
            phi,
            alphaAver,
            alphaScheme
        )
    );

    MULES::explicitSolve(alphaAver, phi, phiAlpha, 1, 0);

    rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alphaAver.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alphaAver.name() <<") = " << min(alphaAver).value()
        << "  Max(" << alphaAver.name() <<") = " << max(alphaAver).value()
        << endl;
}
}

rho == alphaAver*rho1 + (1-alphaAver)*rho2;

{
    fvScalarMatrix alphaAverEqn
    (
        fvm::ddt(alphaAver)
      - fvc::ddt(alphaAver)
      - fvm::laplacian
        (
            volScalarField("Dab", Dab + alphatab*turbulence->nut()),
            alphaAver
        )
    );

    alphaAverEqn.solve();

    alpha2 = 1.0 - alphaAver;
    rhoPhi += alphaAverEqn.flux()*(rho1 - rho2);
}

rho = alphaAver*rho1 + (1-alphaAver)*rho2;
