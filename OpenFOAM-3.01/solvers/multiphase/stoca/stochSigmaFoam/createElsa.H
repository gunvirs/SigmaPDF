    volScalarField ElsaRho
    (
        IOobject
        (
            "ElsaRho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho
    );

    Info<< "Reading field ElsaSigmaPrime\n" << endl;
    volScalarField ElsaSigmaPrime
    (
        IOobject
        (
            "ElsaSigmaPrime",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField ElsaSigmaMin
    (
        IOobject
        (
            "ElsaSigmaMin",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ElsaSigmaPrime
    );

    volScalarField ElsaSigma
    (
        IOobject
        (
            "ElsaSigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ElsaSigmaPrime
    );

    volScalarField ElsaSigmaEquil
    (
        IOobject
        (
            "ElsaSigmaEquil",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ElsaSigmaPrime
    );

    volScalarField ElsaD32
    (
        IOobject
        (
            "ElsaSD32",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,1,0,0,0,0,0),
       "zeroGradient"
    );



    volScalarField ElsaSigmaRatio
    (
        IOobject
        (
            "ElsaSigmaRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,0,0,0,0,0,0),
       "zeroGradient"
    );


    IOdictionary ElsaProperties
    (
        IOobject
        (
            "ElsaProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    /*
    volScalarField cellVol
    (
        IOobject
        (
            "cellVol",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE//AUTO_WRITE
        ),
	mesh,
       dimensionSet(0,3,0,0,0,0,0),
       "zeroGradient"
    );
	*/


    dimensionedScalar CstElsaSigmaMin(ElsaProperties.lookup("CstElsaSigmaMin"));

    dimensionedScalar ElsaSurfaceTensionCoef(ElsaProperties.lookup("ElsaSurfaceTensionCoef"));
    // model for equilibriume surface density choice between 
    //  Lebas IJMF 2009 ==> webEquilMod=0
    //  and Duret IJMF 2013  ==> webEquilMod=1
    dimensionedScalar webEquilMod(ElsaProperties.lookup("webEquilMod"));

    dimensionedScalar ElsaWeberEquil(ElsaProperties.lookup("ElsaWeberEquil"));

    dimensionedScalar CTurbElsaSigmaPrime(ElsaProperties.lookup("CTurbElsaSigmaPrime"));

    dimensionedScalar ElsaTiny("ElsaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);
    dimensionedScalar ElsaTinyJap("ElsaTinyJap",dimensionSet(0,2,-1,0,0,0,0),1.e-8);
    dimensionedScalar kTiny("kTiny",dimensionSet(0,2,-2,0,0,0,0),1.e-1);
    dimensionedScalar sigmaTiny("sigmaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);
    //////////	Number of Droplets   //////////

    volScalarField ndrops
    (
        IOobject
        (
            "ndrops",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,0,0,0,0,0,0),
       "zeroGradient"
    );

