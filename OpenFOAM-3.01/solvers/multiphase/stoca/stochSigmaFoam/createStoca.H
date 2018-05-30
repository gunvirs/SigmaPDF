	//Dizionario dove specifico caratteristiche per risolvere il problema stocastico
	IOdictionary PDFProperties
	(
		IOobject
		(
		    "PDFProperties",
		    runTime.constant(),
		    mesh,
		    IOobject::MUST_READ_IF_MODIFIED,
		    IOobject::NO_WRITE
		)
	);

	scalar clipElsa(PDFProperties.lookupOrDefault<scalar>("clipElsa", 100));

	scalar clipAlpha(PDFProperties.lookupOrDefault<scalar>("clipAlpha", 1e-5));

	scalar clipSD32(PDFProperties.lookupOrDefault<scalar>("clipSD32", 1e-6));
	// Premultiplimative factor for tuebolent diffusion
	// Default value is 1
	// Effect of density (or better momentum flux) ratio should be checked.
	// This parameter has this goal.
	scalar qWiener(PDFProperties.lookupOrDefault<scalar>("qWiener", 1));
	//scalar qWiener(PDFProperties.lookupOrDefault<scalar>("qWiener", rho1.value()/rho2.value()));
	//BOOLEAN operator to decide if calculate or not dissipation turbulence rate EPS (CHECK Ruiz et al. ,2015,Physics of Fluids, 27 (4) , EQ 16 )
	// This is done to evaluate Integral Length Scale l_t = c_0 * TurbulentKineticEnergy^(3/2)/EPS ;
	bool	calcTurb(PDFProperties.lookupOrDefault<bool>("calcTurb", false));


	//BOOLEAN operator to decide if evaluet second statistical moment from stochastic fields
	bool	calcTurb(PDFProperties.lookupOrDefault<bool>("calcSM", false));

	bool	printSourceTerms(PDFProperties.lookupOrDefault<bool>("printSourceSigma", false));



	Info << "clipSD32 \t" << clipSD32 << "\tclipAlpha \t" << clipAlpha<< endl;
	//Number of stocastic Fields

        scalar nFields(readScalar(PDFProperties.lookup("nFields")));


	PtrList <volScalarField> ALPHA(nFields);
	PtrList <volScalarField> SIGMA(nFields);

	int counter;
	for
	(
	counter = 0 ;
	counter < nFields;
	counter++
	)
	{
		word fieldName(name(counter));
		//Alfa 
		ALPHA.set
		(
			counter,
			new volScalarField
			(
				IOobject
				(
					IOobject::groupName("alpha1_",fieldName),
					runTime.timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
				alpha1
			)
		);
		//Sigma 
		SIGMA.set
		(
			counter,
			new volScalarField
			(
				IOobject
				(
					IOobject::groupName("ElsaSigmaPrime",fieldName),
					runTime.timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
				ElsaSigmaPrime
			)
		);
	}

	//The Pointer List For the Wiener Increment is created now, and initialisated as random variable at each time step
	PtrList <uniformDimensionedVectorField> dWdt(nFields);
	for
	(
	counter = 0 ;
	counter < nFields;
	counter++
	)
	{
		word fieldName(name(counter));
		dWdt.set
		(
			counter,
			new uniformDimensionedVectorField
			(
				IOobject
				(
					IOobject::groupName("Wiener_",fieldName),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE//AUTO_WRITE
				),
				dimensionedVector
				(
					"0",
					dimensionSet(0, 0, -0.5, 0, 0),
					vector::zero
				)
			)
		);
	}
	////////////////////	2ND MOMENTS
	volScalarField alphaVariance
	(
	IOobject
	(
		"alphaVariance",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,//NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0, 0, 0, 0, 0,0,0),
	"zeroGradient"
	);

	volScalarField sigmaPrimeVariance
	(
	IOobject
	(
		"sigmaPrimeVariance",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,//NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0, -2, 0, 0, 0,0,0),
	"zeroGradient"
	);

	/////////////////////	Calculation Turbulence Properties
	volScalarField dissip
	(
	IOobject
	(
		"dissip",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,//NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0, 2, -3, 0, 0,0,0),
	"zeroGradient"
	);




