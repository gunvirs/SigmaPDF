{


	for
	(
	int counter = 0 ;
	counter < nFields;
	counter++
	)
	{
       Info <<  " counter  " << counter<< endl;
///////////////////////////////////////////////////
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
        dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0)
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
        dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0)
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
        dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0)
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
        dimensionedScalar("0", ElsaSigmaRatio.dimensions(), 0)
    );

////////////////////////////////////////// end tmp
     Info << "compute Elsa SigmaMin ..." << endl;

     tmpElsaSigmaMin.internalField()=CstElsaSigmaMin*sqrt(ALPHA[counter]*(1.-ALPHA[counter]))*pow(mesh.V(),-1.0/3.0);

     Info << "compute Elsa SigmaEquil ..." << endl;

     if ( webEquilMod.value() < 0.5)
     {
       // WeberEquil Lebas 2009 IJMF
       //ElsaSigmaEquil.internalField()=rho*turbulence->k()/ElsaWeberEquil/
       //                             ElsaSurfaceTensionCoef*alpha1*(1.-alpha1);
       tmpElsaSigmaEquil.internalField()=rho*turbulence->k()/ElsaWeberEquil/
                                    ElsaSurfaceTensionCoef*ALPHA[counter]*(1.-ALPHA[counter]);
       Info << "mod Lebas" << endl;
     }
     else
     {
       // Modified WeberEquil Benjamin Duret et al IJMF 2013
       //ElsaSigmaEquil.internalField()=2.*(rho1+rho2)*turbulence->k()/ElsaWeberEquil/
       //                            ElsaSurfaceTensionCoef*alpha1*(1.-alpha1);
       tmpElsaSigmaEquil.internalField()=2.*(rho1+rho2)*turbulence->k()/ElsaWeberEquil/
                                   ElsaSurfaceTensionCoef*ALPHA[counter]*(1.-ALPHA[counter]);
       Info << "mod Duret "  << endl;
     }

     tmpElsaSigmaEquil.max(0.0);

//    Not completely satisfactory, but it is to go over the undertermination of 
//    SigmaMin that goes to infinity as DX tends to zerO  pb in average and in 2D case ...

       tmpElsaSigmaMin= min(tmpElsaSigmaMin,tmpElsaSigmaEquil*0.001);

     Info << "compute Elsa Sigma ..." << endl; 

     tmpElsaSigma.internalField()=ElsaSigmaPrime.internalField()+tmpElsaSigmaMin.internalField();

     tmpElsaSigma.internalField()=max(tmpElsaSigma.internalField(),0.0);

//   Info << ElsaSigma.boundaryField() << endl; 

   Info << "compute Elsa Ratio ..." << endl;

     tmpElsaSigmaRatio.internalField()=(max((tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny)),1.));
     tmpElsaSigmaRatio.boundaryField()=ElsaSigmaPrime.boundaryField(); 

    Info << "... done" << endl;

		word sigmaScheme("div(phi,ElsaSigmaPrime)");

		fvScalarMatrix ElsaSigmaPrimeEqn
		(
			fvm::ddt(SIGMA[counter])
			+ fvm::div(phi,SIGMA[counter],sigmaScheme)
			- fvm::laplacian
			(
			Dab + alphatab*turbulence->nut(), SIGMA[counter],
			"laplacian(Dab,ElsaSigmaPrime)"
			)
			-sqrt(2*alphatab*turbulence->nut())*(fvc::grad(SIGMA[counter]) & dWdt[counter])
			== (CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(tmpElsaSigmaRatio))*(tmpElsaSigmaEquil-tmpElsaSigmaMin )
			+fvm::Sp(-CTurbElsaSigmaPrime* (turbulence->k()/(ElsaTinyJap+turbulence->nut())) *(tmpElsaSigmaRatio),SIGMA[counter])
		);

		ElsaSigmaPrimeEqn.relax();
		ElsaSigmaPrimeEqn.solve();

		SIGMA[counter].max(0.0);
		tmpElsaSigma=SIGMA[counter]+tmpElsaSigmaMin;

		Info << "...done" <<endl;
		volScalarField tmpAlpha = ALPHA[counter];

		volScalarField tmpLimitElsa
		(
		IOobject
		(
		"tmpLimitElsa",
		runTime.timeName(),
		mesh
		),
		mesh,
		dimensionedScalar("0", ElsaSigmaPrime.dimensions(), 0)
		);

		tmpLimitElsa.internalField() =2.3249*0.020801/pow(mesh.V(),0.3333333);
		//tmpLimitElsa.internalField() =pow(4*3.14,1/3)*pow(3*1E-4,2/3)*pow(mesh.V(),-1.0/3.0);

		//tmpLimitElsa.boundaryField()=ElsaSigmaPrime.boundaryField(); 
		forAll(mesh.cells(),celli)
		{
		//tmpAlpha[celli] = ( tmpAlpha[celli] < 0.0001 ) ? 0 : tmpAlpha[celli] ;
		tmpAlpha[celli] = ((SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  tmpLimitElsa[celli]) ? 0 : tmpAlpha[celli] ;
		}

		if(counter == 0 )
		{
		ElsaSigmaPrime.internalField() = SIGMA[counter];
		ElsaD32.internalField()=6.0*(tmpAlpha)*(1.0- (tmpAlpha) )/(ElsaSigmaMin+SIGMA[counter]+ElsaTiny);
		ElsaSigma.internalField() = tmpElsaSigma;
		ElsaSigmaRatio.internalField() = tmpElsaSigmaRatio;
		ElsaSigmaEquil.internalField() = tmpElsaSigmaEquil;
		ElsaSigmaMin.internalField() = tmpElsaSigmaMin;
		}
		else
		{
		ElsaSigmaPrime.internalField() += SIGMA[counter];
		ElsaD32.internalField() += 6.0*(tmpAlpha)*(1.0- (tmpAlpha) )/(tmpElsaSigmaMin+SIGMA[counter]+ElsaTiny);
		ElsaSigma.internalField() += tmpElsaSigma;
		ElsaSigmaRatio.internalField() += tmpElsaSigmaRatio;
		ElsaSigmaEquil.internalField() += tmpElsaSigmaEquil;
		ElsaSigmaMin.internalField() += tmpElsaSigmaMin;
		}
		if(counter == 15 )
		{
		Info << "STAMPO A COUNTER\t" << counter <<endl;
		ALPHA[counter] = tmpAlpha;
		SIGMA[counter] = tmpLimitElsa;
		}

	}
		ElsaSigmaPrime /= nFields;
		ElsaSigmaPrime.max(0.0);
		ElsaD32 /= nFields;
		ElsaSigma 	/= nFields;
		ElsaSigmaRatio 	/= nFields;
		ElsaSigmaEquil 	/= nFields;
		ElsaSigmaMin 	/= nFields;
	ElsaD32.correctBoundaryConditions();
    Info << "...done" <<endl;

}
