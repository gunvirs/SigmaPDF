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
	//Info << "MAX e MIN vari \n "  <<max(tmpElsaSigmaMin).value()  << max(tmpElsaSigma).value()  << max(tmpElsaSigmaEquil).value()  << max(tmpElsaSigmaRatio).value()  << endl;

	Info << "clipSD32 \t" << clipSD32 << "\tclipAlpha \t" << clipAlpha<< endl ;

////////////////////////////////////////// end tmp
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
	//QUANDO sono in parallelo, e lavoro su buondaryField(), lavora sulle boundary del mesh totale o solo sulla 
	//tmpElsaSigma.boundaryField()=666;//SIGMA[counter].boundaryField();//Cioe' ho paura che praticamente, non inizializzando a dovere i campi anche sulle boundary, poi si sputtani il tutto

	Info << "compute Elsa Ratio ..." << endl;
	// 12 Marzo::: Tolgo il MAX
	//tmpElsaSigmaRatio.internalField()=(max((tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny)),1.));
	tmpElsaSigmaRatio.internalField()=(max((SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny)),1.));
	tmpElsaSigmaRatio.correctBoundaryConditions();
	//	28 Feb ::: Qst Un Problema?
	//tmpElsaSigmaRatio.boundaryField()=SIGMA[counter].boundaryField();
	//9 Marzo : creo un Termine sorgente del tutto esplicito e sti grandissimi cazzi
	// PS: la k usata dovrebbe essere quella risolta secondo il modello, non quella SGS. Controllare Please questa turbulence->()
	//volScalarField  SourceElsa(CTurbElsaSigmaPrime*SIGMA[counter]*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(1));//-(SIGMA[counter]/tmpElsaSigmaEquil+sigmaTiny)));
	//Questa formulazione sotto diverge....credo che vada a zero in qualche modo;
	//Pout << "Hello from processor " << Pstream::myProcNo() << "for SIgmaEquil, \n the MAX " << max(tmpElsaSigmaEquil.value()) << " the min " << min (tmpElsaSigmaEquil.value()) << endl;

	//volScalarField  SourceElsa(CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*SIGMA[counter]*(1-(SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny))));

	// Nel Caso diverga ancora, faccio come sotto
	//Lo pongo inizialmente pari a zero, e poi lo calcolo per bene solo se SigmaEquil != 0
	// Il problema credo non e' dovuto al fatto che mi diverge il calcolo di source elsa di per se, tanto che mi diverge il sistema matriciale.
	// Se SigmaEquil --> 0, sto termine tende all'infinito.
	// 
	//volScalarField  SourceElsa(0);
	//volScalarField  SourceElsa(CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*SIGMA[counter]*(1-(SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny))));
	//forAll(mesh.cells(),celli)
	//{
	//SourceElsa[celli] = (tmpElsaSigmaEquil[celli] == 0) ? 0 : SourceElsa[celli]   ;
	//}

	//SourceElsa = CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*SourceElsa;

	// Linearizzazione temporale secondo SNM \dot{S}/(1-\Delta t * \partial \dot{S} / \partial Sigma)

	//volScalarField  SourceElsa(CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*SIGMA[counter]*(1-(SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny)))/(1 - runTime.deltaT()*CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(1-2*((SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny))))));
	volScalarField  SourceElsa(CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*tmpElsaSigma*(1-(tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny)))/(1 - runTime.deltaT()*CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(1-2*((tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny))))));
	SourceElsa.correctBoundaryConditions();


	Info << "... done" << endl;

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
			-Foam::sqrt((2/qWiener)*alphatab*turbulence->nut())*(fvc::grad(SIGMA[counter]) & dWdt[counter])//Stochastic term
			//-Foam::sqrt(2*alphatab*turbulence->nut())*(fvc::grad(SIGMA[counter]) & dWdt[counter])//Stochastic term

			// IN QUESTO MODO DOVREI ESSERE CAPACE ANCHE DI FARE QUALCHE MODELLO PIU' GENERICO
			==  SourceElsa
		);
			/*
			//CLASSICI SOURCE TERMSSSSSS
			== (CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(tmpElsaSigmaRatio))*(tmpElsaSigmaEquil-tmpElsaSigmaMin )
			+fvm::Sp(-CTurbElsaSigmaPrime* (turbulence->k()/(ElsaTinyJap+turbulence->nut())) *(tmpElsaSigmaRatio),SIGMA[counter])
			//Il termine sorgente nella trasposizione classica di elsa e' stato linearizzato e diviso in parte implicita ed esplicita. Ls parte implicita 
			// Un altra idea potrebbe essere l'utilizzo del cosidetto SuSp::::: vedi link
			//    http://openfoamwiki.net/index.php/HowTo_Adding_a_new_transport_equation
			// Solo che per come l'hanno formulato I franzosi, sto SigmaRatio al sotto a uno non ci va, quindi ho sempre il termine con SP negativo (SigmaRatio non sarà minore di uno) e quindi sempre implicito. Devo capire se il termine sorgente e' corretto o meno.
			*/



		ElsaSigmaPrimeEqn.relax();
		Info << "... Solving SIGMA "<< counter << " ...\n";
		ElsaSigmaPrimeEqn.solve();
		//	HERE:::::::::: vedi note
		Info << " SIGMA "<< counter << " Solved !!!!!\n";
		SIGMA[counter].max(0.0);
		/////////	Riscontro problemi in alcuni campi. 
/*

		forAll(mesh.cells(),celli)
		{
		SIGMA[counter][celli] = (SIGMA[counter][celli] <  clipElsa) ? 0 : SIGMA[counter][celli] ;
		}
		SIGMA[counter].correctBoundaryConditions();
		//////////
		Info << "\t...done" <<endl;
		volScalarField tmpAlpha = ALPHA[counter];
*/

		/*volScalarField tmpLimitElsa
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
		//tmpLimitElsa.internalField() =2.3249*0.0044814/pow(mesh.V(),0.3333333);//Clip per alpha 1E-4
		tmpLimitElsa.internalField() =pow(4*3.14,1/3)*pow(3*clipAlpha,2/3)*pow(mesh.V(),-1.0/3.0);
		tmpLimitElsa.boundaryField()=ElsaSigmaPrime.boundaryField(); */
		//volScalarField tmpAlpha = ALPHA[counter];
		scalar clipAlpha23=pow(3*clipAlpha,2/3);
		forAll(mesh.cells(),celli)
		{
		//SIGMA[counter][celli] = (SIGMA[counter][celli] <  clipElsa) ? 0 : SIGMA[counter][celli] ;
		SIGMA[counter][celli] = ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  2.3249*clipAlpha23/pow(mesh.V()[celli],1/3)) ? 0 : SIGMA[counter][celli] ;
		scalar tmpSD32 = 0;
		//tmpSD32 = ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) || (SIGMA[counter][celli] <  clipElsa) ) ? 0 : 6.0*((ALPHA[counter][celli])*(1.0-(ALPHA[counter][celli])))/(tmpElsaSigmaMin[celli]+SIGMA[counter][celli]);
		tmpSD32 = ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) )  ) ? 0 : 6.0*((ALPHA[counter][celli])*(1.0-(ALPHA[counter][celli])))/(tmpElsaSigmaMin[celli]+SIGMA[counter][celli]);
		//tmpSD32 = ( ( ALPHA[counter][celli] <  clipAlpha ) || (SIGMA[counter][celli] <  clipElsa) ) ? 0 : 6.0*((ALPHA[counter][celli])*(1.0-(ALPHA[counter][celli])))/(tmpElsaSigmaMin[celli]+SIGMA[counter][celli]);
		if(counter == 0)
		{
		ElsaD32[celli] = ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) ||   (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;

		//ElsaD32[celli] = ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) || (SIGMA[counter][celli] <  clipElsa) ||  (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;

		//ElsaD32[celli] = ( ( ALPHA[counter][celli] <  clipAlpha  ) || (SIGMA[counter][celli] <  clipElsa) ||  (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;
		}
		else
		{
		ElsaD32[celli] += ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) ||   (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;
		//ElsaD32[celli] += ( ( (SIGMA[counter][celli]+tmpElsaSigmaMin[celli]) <  (2.3249*clipAlpha23/pow(mesh.V()[celli],1/3) ) ) || (SIGMA[counter][celli] <  clipElsa) ||  (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;


		//ElsaD32[celli] += ( ( ALPHA[counter][celli] <  clipAlpha  ) || (SIGMA[counter][celli] <  clipElsa) ||  (tmpSD32 < clipSD32)  ) ? 0 : tmpSD32;
		}
		//PercentualeCampi +=  ( ElsaD32[celli] > 0 ) ? 1 : 0;
		}
		SIGMA[counter].correctBoundaryConditions();
		ElsaD32.correctBoundaryConditions();
		tmpElsaSigma=SIGMA[counter]+tmpElsaSigmaMin;
		tmpElsaSigma.correctBoundaryConditions();
		/////	Quando lavoro con il parallelo, e faccio cicli, o meglio lavoro cella per cella, devo stare attento, che devo aggiornare le boundary del mesh relativo al songolo processore
		//SIGMA[counter].correctBoundaryConditions();
		//////////
		//tmpElsaSigma=SIGMA[counter]+tmpElsaSigmaMin;
		//tmpElsaSigma.correctBoundaryConditions();


		if(counter == 0 )
		{
		ElsaSigmaPrime.internalField() = SIGMA[counter];
		//ElsaD32.internalField()=6.0*(tmpAlpha)*(1.0- (tmpAlpha) )/(tmpElsaSigmaMin+SIGMA[counter]+ElsaTiny);
		//ElsaSigma.internalField() = tmpElsaSigma;
		ElsaSigma.internalField() = SIGMA[counter].internalField()+tmpElsaSigmaMin.internalField();
		ElsaSigmaRatio.internalField() = tmpElsaSigmaRatio;
		ElsaSigmaEquil.internalField() = tmpElsaSigmaEquil;
		ElsaSigmaMin.internalField() = tmpElsaSigmaMin;
		}
		else
		{
		ElsaSigmaPrime.internalField() += SIGMA[counter];
		//ElsaD32.internalField() += 6.0*(tmpAlpha)*(1.0- (tmpAlpha) )/(tmpElsaSigmaMin+SIGMA[counter]+ElsaTiny);
		//ElsaSigma.internalField() += tmpElsaSigma;
		ElsaSigma.internalField() += SIGMA[counter].internalField()+tmpElsaSigmaMin.internalField();
		ElsaSigmaRatio.internalField() += tmpElsaSigmaRatio;
		ElsaSigmaEquil.internalField() += tmpElsaSigmaEquil;
		ElsaSigmaMin.internalField() += tmpElsaSigmaMin;
		}
	}

	ElsaSigmaPrime /= nFields;
	//ElsaSigmaPrime.internalField() /= nFields;
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
	////////////////////////////////// FINO QUA
	//Apparentemente boudaryField() non cambia il valore alla boundary cioè al confine del dominio, ma solo al confine del dominio decomposto.
	//Cioe se faccio in parallelo, con 2 cpu, e vedo il risultato, troverò che il comando cambia il valore del campo non a tutte le boudary del dominio, ma a quello tra il 
	//il mesh0 e mesh1.
	ElsaSigma.correctBoundaryConditions();
	//IMP:::: BOUDARYFIELD corregge il valore della boundary ma NON IL VALORE DI CELLA, IL VALORE DI CENTRO FACCIA!!!!!!!!!!!!!!!!!
	// Se giro in seriale (no parallel) pare che non abbia effetto
		//ElsaSigma.boundaryField()=777;//.internalField() = tmpElsaSigma;
	// Se Faccio correctBoundaryConditions() dopo boudaryField, quello corregge il valore alla BOUDARY tra i 2 mesh decomposti
	// e la modifica fatta sotto scompare.

	//nDropps
	//ndrops 	/= nFields;
	//PercentualeCampi/= nFields;
	//ndrops.correctBoundaryConditions();
	Info << "...done" <<endl;
/////////////////		N Droplets
	Info << "Inizio Numb" <<endl;

		///////	Printsssssss
/*

	if (printSourceTerms)
	{
	//SnonLin = SIGMA[counter]*(1-(SIGMA[counter]/(tmpElsaSigmaEquil+sigmaTiny)));
	SnonLin= tmpElsaSigma*(1-(tmpElsaSigma/(tmpElsaSigmaEquil+sigmaTiny)));
	SLin	=  SourceElsa;
	//SLinCappato	=  CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*SIGMA[counter]*(1-tmpElsaSigmaRatio)/(1 - runTime.deltaT()*CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(1-2*tmpElsaSigmaRatio));
	//SLinCappato	=  CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*tmpElsaSigma*(1-tmpElsaSigmaRatio)/(1 - runTime.deltaT()*CTurbElsaSigmaPrime*(turbulence->k()/(ElsaTinyJap+turbulence->nut()))*(1-2*tmpElsaSigmaRatio));
	}
*/




	//PerctSigm/=nFields;
	/*forAll(mesh.cells(),celli)
	{
	ElsaD32[celli] =  (  PerctSigm[celli] < 0.51  ) ? 0 : ElsaD32[celli] ;
	}*/
	forAll(mesh.cells(),celli)
	{
	ndrops[celli] =  (ElsaD32[celli] > clipSD32) ? alpha1[celli]*6*mesh.V()[celli]/(constant::mathematical::pi*pow(ElsaD32[celli],3)) : 0 ;
	ElsaD32[celli] =  (ElsaD32[celli] > clipSD32) ? ElsaD32[celli] : 0 ;
	//ElsaD32[celli] =  ( (ElsaD32[celli] > 1E-6) && (PercentualeCampi > 0.51) ) ? ElsaD32[celli] : 0 ;
	}
	Info << "Finisco Numb" <<endl;
	ndrops.correctBoundaryConditions();
	ElsaD32.correctBoundaryConditions();
	//ElsaD32.correctBoundaryConditions();
}
