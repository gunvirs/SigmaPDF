if (calcSM==true)
{
	for
	(
	int counter = 0 ;
	counter < nFields;
	counter++
	)
	{
	if(counter == 0 )
	{
		alphaVariance = pow((ALPHA[counter]-alpha1),2);
		sigmaPrimeVariance = pow((SIGMA[counter]-ElsaSigmaPrime),2);
	}

	else
	{
		alphaVariance += pow((ALPHA[counter]-alpha1),2);
		sigmaPrimeVariance += pow((SIGMA[counter]-ElsaSigmaPrime),2);
	}
	}
	alphaVariance/= nFields;
	sigmaPrimeVariance/= nFields;
}
