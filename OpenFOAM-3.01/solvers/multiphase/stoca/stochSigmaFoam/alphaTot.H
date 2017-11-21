#include "alphaControls.H"
//Tmp variable to store rhoPhi
surfaceScalarField tmpRhoPhiSumm
(
	IOobject
	(
            "tmpRhoPhiSumm",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
);


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
        subCycle<volScalarField> alphaSubCycle(ALPHA[counter], nAlphaSubCycles);
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
			    ALPHA[counter],
			    alphaScheme
			)
		    );

		    MULES::explicitSolve(ALPHA[counter], phi, phiAlpha, 1, 0);

		    //tmpRhoPhiSumm = phiAlpha*(rho1 - rho2) + phi*rho2;
		    rhoPhiSum = phiAlpha*(rho1 - rho2) + phi*rho2;
		}
        tmpRhoPhiSumm += (runTime.deltaT()/totalDeltaT)*rhoPhiSum;
    }

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
		    ALPHA[counter],
		    alphaScheme
		)
	    );

	    MULES::explicitSolve(ALPHA[counter], phi, phiAlpha, 1, 0);

	    tmpRhoPhiSumm += phiAlpha*(rho1 - rho2) + phi*rho2;
	}
}

rho == ALPHA[counter]*rho1 + (1.0 - ALPHA[counter])*rho2;
{

	volScalarField DADA(alphatab*turbulence->nut());//NON CI VUOLE DAB IN STOCHA

	fvScalarMatrix alpha1Eqn
	(
		fvm::ddt(ALPHA[counter])
		- fvc::ddt(ALPHA[counter])
		- fvm::laplacian
		(
			volScalarField("Dab", Dab + alphatab*turbulence->nut()),
			ALPHA[counter]
		)
		-(fvc::grad(ALPHA[counter]) & dWdt[counter])*sqrt(2*DADA)
	);

    alpha1Eqn.solve();

    tmpRhoPhiSumm += alpha1Eqn.flux()*(rho1 - rho2);

}
	ALPHA[counter].max(0.0);
	ALPHA[counter].min(1.0);

rho = ALPHA[counter]*rho1 + (1.0 - ALPHA[counter])*rho2;

if(counter == 0 )
{
alpha1 = ALPHA[counter];
rhoPhi = tmpRhoPhiSumm;
}
else
{
alpha1 += ALPHA[counter];
rhoPhi += tmpRhoPhiSumm;
}


