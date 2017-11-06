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

	//Number of stocastic Fields

        scalar nFields(readScalar(PDFProperties.lookup("nFields")));

	//Preso spunto da https://www.cfd-online.com/Forums/openfoam-programming-development/119642-define-array-volvectrofield-type.html

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
					IOobject::AUTO_WRITE
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

