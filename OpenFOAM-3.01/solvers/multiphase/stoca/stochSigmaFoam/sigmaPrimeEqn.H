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
		tmpElsaSigmaMin *= scalar(0);tmpElsaSigmaMin.correctBoundaryConditions();
		tmpElsaSigmaEquil *= scalar(0);tmpElsaSigmaEquil.correctBoundaryConditions();
		tmpElsaSigma *= scalar(0);tmpElsaSigma.correctBoundaryConditions();
		tmpElsaSigmaRatio *= scalar(0);tmpElsaSigmaRatio.correctBoundaryConditions();

	Info << "compute Elsa SigmaMin ..." << endl;

	tmpElsaSigmaMin.internalField()=CstElsaSigmaMin*sqrt(ALPHA[counter]*(1.-ALPHA[counter]))*pow(mesh.V(),-1.0/3.0);
	tmpElsaSigmaMin.correctBoundaryConditions();

	Info << "compute Elsa SigmaEquil ..." << endl;

	if ( webEquilMod.value() < 0.5)
	{
	// WeberEquil Lebas 2009 IJMF
		tmpElsaSigmaEquil.internalField()=rho*turbulence->k()/ElsaWeberEquil/
			ElsaSurfaceTensionCoef*ALPHA[counter]*(1.-ALPHA[counter]);
		Info << "mod Lebas" << endl;
	}
	else
	{
	// Modified WeberEquil Benjamin Duret et al IJMF 2013
		tmpElsaSigmaEquil.internalField()=2.*(rho1+rho2)*turbulence->k()/ElsaWeberEquil/
			ElsaSurfaceTensionCoef*ALPHA[counter]*(1.-ALPHA[counter]);
		Info << "mod Duret "  << endl;
	}
	tmpElsaSigmaEquil.correctBoundaryConditions();

	tmpElsaSigmaEquil.max(0.0);
	//    Not completely satisfactory, but it is to go over the undertermination of 
	//    SigmaMin that goes to infinity as DX tends to zerO  pb in average and in 2D case ...
	tmpElsaSigmaMin= min(tmpElsaSigmaMin,tmpElsaSigmaEquil*0.001);
	tmpElsaSigmaMin.correctBoundaryConditions();

	Info << "compute Elsa Sigma for the stochastic field numer\t"<< counter << "  ..." << endl; 

	tmpElsaSigma.internalField()=SIGMA[counter].internalField()+tmpElsaSigmaMin.internalField();
	tmpElsaSigma.internalField()=max(tmpElsaSigma.internalField(),0.0);
	tmpElsaSigma.correctBoundaryConditions();


	Info << "compute Elsa Ratio ..." << endl;

	tmpElsaSigmaRatio.internalField()=(max((SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny)),1.));
	tmpElsaSigmaRatio.correctBoundaryConditions();

	volScalarField  SourceElsa(CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*tmpElsaSigma*(1-(tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny)))/(1 - runTime.deltaT()*CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(1-2*((tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny))))));

	SourceElsa.correctBoundaryConditions();

	Info << "... done preparation" << endl;

		word sigmaScheme("div(phi,ElsaSigmaPrime)");
		Info << "Preparin to solve SIGMA "<< counter << " ..." ;
		fvScalarMatrix ElsaSigmaPrimeEqn
		(
			fvm::ddt(SIGMA[counter])
			+ fvm::div(phi,SIGMA[counter],sigmaScheme)
			- fvm::laplacian
			(
			Dab + alphatab*turbulence->nut(), SIGMA[counter],
			"laplacian(Dab,ElsaSigmaPrime)"
			)
			-Foam::sqrt((2/qWiener)*alphatab*turbulence->nut())*(fvc::grad(SIGMA[counter]) & dWdt[counter])
			==  SourceElsa
		);

		ElsaSigmaPrimeEqn.relax();
		Info << "... Solving SIGMA "<< counter << " ...\n";
		ElsaSigmaPrimeEqn.solve();

		Info << " SIGMA "<< counter << " Solved !!!!!\n";
		SIGMA[counter].max(0.0);

		scalar clipAlpha23=pow(3*clipAlpha,2/3);
		forAll(mesh.cells(),celli)
		{
		SIGMA[counter][celli] = ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  2.3249*clipAlpha23/pow(mesh.V()[celli],1/3)) ? 0 : SIGMA[counter][celli] ;
		scalar tmpSD32 = 0;
		tmpSD32 = ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) )  ) ? 0 : 6.0*((ALPHA[counter][celli])*(1.0-(ALPHA[counter][celli])))/(tmpElsaSigmaMin[celli]+SIGMA[counter][celli]);
		if(counter == 0)
		{
		ElsaD32[celli] = ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) ||   (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;
		}
		else
		{
		ElsaD32[celli] += ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) ||   (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;
		}
		}
		SIGMA[counter].correctBoundaryConditions();
		ElsaD32.correctBoundaryConditions();
		tmpElsaSigma=SIGMA[counter]+tmpElsaSigmaMin;
		tmpElsaSigma.correctBoundaryConditions();

		if(counter == 0 )
		{
		ElsaSigmaPrime.internalField() = SIGMA[counter];
		ElsaSigma.internalField() = SIGMA[counter].internalField()+tmpElsaSigmaMin.internalField();
		ElsaSigmaRatio.internalField() = tmpElsaSigmaRatio;
		ElsaSigmaEquil.internalField() = tmpElsaSigmaEquil;
		ElsaSigmaMin.internalField() = tmpElsaSigmaMin;
		}
		else
		{
		ElsaSigmaPrime.internalField() += SIGMA[counter];
		ElsaSigma.internalField() += SIGMA[counter].internalField()+tmpElsaSigmaMin.internalField();
		ElsaSigmaRatio.internalField() += tmpElsaSigmaRatio;
		ElsaSigmaEquil.internalField() += tmpElsaSigmaEquil;
		ElsaSigmaMin.internalField() += tmpElsaSigmaMin;
		}
	}

	ElsaSigmaPrime /= nFields;

	ElsaSigmaPrime.max(0.0);
	ElsaD32 /= nFields;
	ElsaSigma 	/= nFields;
	ElsaSigmaRatio 	/= nFields;
	ElsaSigmaEquil 	/= nFields;
	ElsaSigmaMin 	/= nFields;
//////////////////////////////
	ndrops 	/= nFields;

	ElsaD32.correctBoundaryConditions();
	ElsaSigma.correctBoundaryConditions();
	ElsaSigmaRatio.correctBoundaryConditions();
	ElsaSigmaEquil.correctBoundaryConditions();
	ElsaSigmaMin.correctBoundaryConditions();
	ElsaSigmaPrime.correctBoundaryConditions();
	ElsaSigma.correctBoundaryConditions();

	forAll(mesh.cells(),celli)
	{
	ndrops[celli] =  (ElsaD32[celli] > clipSD32) ? alpha1[celli]*6*mesh.V()[celli]/(constant::mathematical::pi*pow(ElsaD32[celli],3)) : 0 ;
	ElsaD32[celli] =  (ElsaD32[celli] > clipSD32) ? ElsaD32[celli] : 0 ;
	}
	ndrops.correctBoundaryConditions();
	ElsaD32.correctBoundaryConditions();
}
