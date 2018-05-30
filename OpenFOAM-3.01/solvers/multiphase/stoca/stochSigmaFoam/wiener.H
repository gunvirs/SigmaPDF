	// Takes a random number between 0 and 1, after takes the nearest integer to this number. All this with the function intRand
	// Then, if the result is 0, xrand is equal to -1, otherwise xrand equal is to 1. This for the first half of the number fields.
	// The other half is the opposite(e.g. for nFields = 4, if xrand(1) = -1 and  xrand(2) = 1, then xrand(3) = -xrand(1) = 1 and xrand(4) = -xrand(3) = -1 )
	// Finally the Weiner increment term is equal the xrand/sqrt(deltaT). In such a way the Wiener increment has mean 0 and variance equal to dt. 
	// Formarly, the qiener term should be a Gaussian random variable N(0,dT), but to obtain a proper gaussian distribution the number of stochastic fields 
	// must be high. This procedure is cheaper because allows to obtain the same N(0,dT) with few nFields.

	// Initialise the Wiener Term	 
	//- New object of class Random

	// - time step

	dimensionedScalar overSqrtDeltaT = 1/sqrt(runTime.deltaT());

	//Random obj(1);
	Random obj(mesh.time().value()*10000000000);
	
	int counter;
	if (nFields==1)
	{
		dWdt[0] = Foam::vector(0,0,0)*overSqrtDeltaT;
	}
	else
	{
		for
		(
		counter = 0 ;
		counter < nFields/2;
		counter++
		)
		// vedi link sotto per esempio su generatore random number in OF
		// https://bitbucket.org/shor-ty/grayscottmodel/src/6dc3187186c2aaec2927d78b9d82e77506482081/applications/solvers/grayScottFoam/perputation.H?at=grayScottFoam-2.3.x&fileviewer=file-view-default
		{
			dWdt[counter] =	Foam::vector(1,0,0)*(((obj.scalar01())>= 0.5 ? (1) : (0))*2 -1)*overSqrtDeltaT+Foam::vector(0,1,0)*(((obj.scalar01())>= 0.5 ? (1) : (0))*2 -1)*overSqrtDeltaT+Foam::vector(0,0,1)*(((obj.scalar01())>= 0.5 ? (1) : (0))*2 -1)*overSqrtDeltaT;

			dWdt[(nFields/2)+counter] = -dWdt[counter];

		}
	}
