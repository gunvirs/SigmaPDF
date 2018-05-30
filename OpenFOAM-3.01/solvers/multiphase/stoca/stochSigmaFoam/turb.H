if (calcTurb==true)
{
	//scalarField modulostrain( magSqr(symm(fvc::grad(U))) );//S:S = S_{ij}*S_{ij} , dove S = 0.5*((\partial u_i /\partial x_ j)+(\partial u_j /\partial x_ i)) = symm( grad(U) ) --> Strain Rate Tensor
	//totVisc = turbulence->nut()+mixture.mu();
	// gently reminder:::: mu e' la viscosita' cinematica, non dinamica. Controllando dove e' definita, ha le dimensioni di  quella cinematica: dimViscosity (dimArea/dimTime) --> lui sottointende viscosita' cinematica
	//dissip	= ( turbulence->nut()+mixture.mu() )*modulostrain;
	dissip	= mixture.nu()*magSqr(symm(fvc::grad(U)))+turbulence->nut()*magSqr(symm(fvc::grad(U)));
	//integralL = c_0*pow(TKE,2/3)/dissip;
}
