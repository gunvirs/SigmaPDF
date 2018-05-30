///////////////////////////////////////////////////
	// Modifica del 30 01 : Il problema e' che diverge. Diverge solo quella della superficie. 
	// questo problema e' sorto quando ho messo i termini sorgenti funzione di SIGMA[counter] stesso .
	// 
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

    volScalarField tmpElsaSigmaMin
    (
        IOobject
        (
            "tmpElsaSigmaMin",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0),
	// inserted il 30-01
       "zeroGradient"
    );

    volScalarField tmpElsaSigma
    (
        IOobject
        (
            "tmpElsaSigma",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0),
	// inserted il 30-01
       "zeroGradient"
    );

    volScalarField tmpElsaSigmaEquil
    (
        IOobject
        (
            "tmpElsaSigmaEquil",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0),
	// inserted il 30-01
       "zeroGradient"
    );

    volScalarField tmpElsaSigmaRatio
    (
        IOobject
        (
            "tmpElsaSigmaRatio",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", ElsaSigmaRatio.dimensions(), 0),
	// inserted il 30-01
       "zeroGradient"
    );


	////////////////////////////////// end tmp

