// Parametrized case by Giovanni Tretola

//Run using:
//m4 -P system/inletL.m4 > system/blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926)

//Geometry
m4_define(dinj,0.1)	// [mm] -->210micron

m4_define(dTot,calc(dinj*20))
m4_define(rTot,calc(dTot/2))
m4_define(rRel,0.9)//PROVOOOOOOOOO
//m4_define(rL,3)
m4_define(rL,4)

//Geometry
m4_define(rinj,calc(dinj/2))

m4_define(r1,calc((dinj/5)/2))
m4_define(r2,calc(r1*3))
m4_define(r3,calc(rinj))
m4_define(rTot,calc(dTot/2))

m4_define(devLength, calc(dinj*5))	//Length inlet
m4_define(length, calc(dinj*25))//50))


//Grid points (integers!):

//Vedi schema 
m4_define(nC, 10)	//numero di celle tra 0 e R1	-- Liquid Inlet	//al massimo 6
m4_define(nR2, 1)	//numero di celle tra R1 e R2
m4_define(nR3, 8)	//numero di celle tra rinj e rL	-- Gas Inlet	//in questo modo qua è sotto i 120 micron, meglio 2 forse
//m4_define(nR3, calc((rL-1)*nC))	//numero di celle tra rinj e rL -- Questo secondo cerchio per non avere asimmetria subito dopo jet. Siccome questo cerchio ha un raggio che vale rL*raggioInjector, con questa espressione avro la stessa dimensione di cella che ho per l iniettore centrale //--> QUESTO SE rGrading2 = 1
m4_define(nS, 25)	//numero di celle tra fine iniettore gas e fine dominio laterale

// Number of cells along cylinder axis
m4_define(zAENC, 50)//256)
m4_define(zEFNC, 128)//256)

//Grading
m4_define(rGrading1, 1)
//m4_define(rGrading2, 1)//--> Ho un secondo cerchio nella sezione di ingresso per non forzare una asimmetria nel getto. 
m4_define(rGrading2, 0.24)
m4_define(nsGrading, 1.8)
m4_define(reverNSG, calc(1/nsGrading))	//0.2)
m4_define(zGrading1, 5)
m4_define(zGrading2, 0.2)
m4_define(zGradingAE, 1)
m4_define(zGradingEF, 7.2)

//Plane E:
m4_define(zE, devLength)//calc(0))
m4_define(r1E, r1)
m4_define(r2E, r2)
m4_define(r3E, r3)
m4_define(rLE, calc(rL*rinj))	//rL = radius Large
m4_define(rTotHalfE, rTotHalf)
m4_define(rTotE, rTot)
m4_define(rRelE, rRel)
m4_define(rRelEc, 0.9)

//Plane F:
m4_define(zF, calc(length + devLength))
m4_define(r1F, r1)
m4_define(r2F, r2)
m4_define(r3F, r3)
m4_define(rLF, calc(rL*rinj))	//rL = radius Large
m4_define(rTotHalfF, rTotHalf)
m4_define(rTotF, rTot)
m4_define(rRelF, rRel)
m4_define(rRelFc, 0.9)

//Plane A:
m4_define(zA, calc(0))//-devLength)
m4_define(r1A, r1)
m4_define(r2A, r2)
m4_define(r3A, r3)
m4_define(rRelA, rRel)
m4_define(rRelAc, 0.9)

/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.6                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(

//Plane A:

//Central O-grid
(calc(rRelA*r3A*cos(pi/4)) -calc(rRelA*r3A*sin(pi/4)) zA) vlabel(A0)
(calc(rRelA*r3A*cos(pi/4)) calc(rRelA*r3A*sin(pi/4)) zA) vlabel(A1)
(calc(-rRelA*r3A*cos(pi/4)) calc(rRelA*r3A*sin(pi/4)) zA) vlabel(A2)
(calc(-rRelA*r3A*cos(pi/4)) -calc(rRelA*r3A*sin(pi/4)) zA) vlabel(A3)
//4 block rounding
(calc(r3A*cos(pi/4)) -calc(r3A*sin(pi/4)) zA) vlabel(A4)
(calc(r3A*cos(pi/4)) calc(r3A*sin(pi/4)) zA) vlabel(A5)
(calc(-r3A*cos(pi/4)) calc(r3A*sin(pi/4)) zA) vlabel(A6)
(calc(-r3A*cos(pi/4)) -calc(r3A*sin(pi/4)) zA) vlabel(A7)


//Plane E:
//Central O-grid
(calc(rRelE*r3E*cos(pi/4)) -calc(rRelE*r3E*sin(pi/4)) zE) vlabel(E0)
(calc(rRelE*r3E*cos(pi/4)) calc(rRelE*r3E*sin(pi/4)) zE) vlabel(E1)
(calc(-rRelE*r3E*cos(pi/4)) calc(rRelE*r3E*sin(pi/4)) zE) vlabel(E2)
(calc(-rRelE*r3E*cos(pi/4)) -calc(rRelE*r3E*sin(pi/4)) zE) vlabel(E3)
//4 block rounding
(calc(r3E*cos(pi/4)) -calc(r3E*sin(pi/4)) zE) vlabel(E4)
(calc(r3E*cos(pi/4)) calc(r3E*sin(pi/4)) zE) vlabel(E5)
(calc(-r3E*cos(pi/4)) calc(r3E*sin(pi/4)) zE) vlabel(E6)
(calc(-r3E*cos(pi/4)) -calc(r3E*sin(pi/4)) zE) vlabel(E7)
//EXTRA block rounding
(calc(rLE*cos(pi/4)) -calc(rLE*sin(pi/4)) zE) vlabel(Eext4)
(calc(rLE*cos(pi/4)) calc(rLE*sin(pi/4)) zE) vlabel(Eext5)
(calc(-rLE*cos(pi/4)) calc(rLE*sin(pi/4)) zE) vlabel(Eext6)
(calc(-rLE*cos(pi/4)) -calc(rLE*sin(pi/4)) zE) vlabel(Eext7)

// Vertices to define block North (6 5 13 14):
(calc(rLE*cos(pi/4))	calc(rTotE)		zE) vlabel(E13)
(-calc(rLE*cos(pi/4))		calc(rTotE)	zE) vlabel(E14)

// Vertices to define block North West: (6 14 15 16)
(-calc(rTotE)		calc(rTotE)		zE) vlabel(E15)
(-calc(rTotE)	calc(rLE*sin(pi/4))		zE) vlabel(E16)
// Vertices to define block West : (7 6 16 17)
(-calc(rTotE)		-calc(rLE*sin(pi/4))	zE) vlabel(E17)
// Vertices to define block South West : (19 7 17 18)
(-calc(rLE*cos(pi/4))	-calc(rTotE)		zE) vlabel(E19)
(-calc(rTotE)		-calc(rTotE)		zE) vlabel(E18)
// Vertices to define block South 4: (8 4 7 19)
(calc(rLE*cos(pi/4))	-calc(rTotE)		zE) vlabel(E8)
// Vertices to define block South Est: (9 10 4 8)
(calc(rTotE)		-calc(rLE*sin(pi/4))	zE) vlabel(E10)
(calc(rTotE)		-calc(rTotE)		zE) vlabel(E9)
// Vertices to define block Est: (10 11 5 4)
(calc(rTotE)		calc(rLE*sin(pi/4))	zE) vlabel(E11)
// Vertices to define block North Est: (11 12 13 5)
(calc(rTotE)		calc(rTotE)		zE) vlabel(E12)


//	plane F
//Central O-grid
(calc(rRelF*r3F*cos(pi/4)) -calc(rRelF*r3F*sin(pi/4)) zF) vlabel(F0)
(calc(rRelF*r3F*cos(pi/4)) calc(rRelF*r3F*sin(pi/4)) zF) vlabel(F1)
(calc(-rRelF*r3F*cos(pi/4)) calc(rRelF*r3F*sin(pi/4)) zF) vlabel(F2)
(calc(-rRelF*r3F*cos(pi/4)) -calc(rRelF*r3F*sin(pi/4)) zF) vlabel(F3)
//4 block rounding
(calc(r3F*cos(pi/4)) -calc(r3F*sin(pi/4)) zF) vlabel(F4)
(calc(r3F*cos(pi/4)) calc(r3F*sin(pi/4)) zF) vlabel(F5)
(calc(-r3F*cos(pi/4)) calc(r3F*sin(pi/4)) zF) vlabel(F6)
(calc(-r3F*cos(pi/4)) -calc(r3F*sin(pi/4)) zF) vlabel(F7)
//EXTRA block rounding
(calc(rLF*cos(pi/4)) -calc(rLF*sin(pi/4)) zF) vlabel(Fext4)
(calc(rLF*cos(pi/4)) calc(rLF*sin(pi/4)) zF) vlabel(Fext5)
(calc(-rLF*cos(pi/4)) calc(rLF*sin(pi/4)) zF) vlabel(Fext6)
(calc(-rLF*cos(pi/4)) -calc(rLF*sin(pi/4)) zF) vlabel(Fext7)

// Vertices to define block North (6 5 13 14):
(calc(rLF*cos(pi/4))	calc(rTotF)		zF) vlabel(F13)
(-calc(rLF*cos(pi/4))		calc(rTotF)	zF) vlabel(F14)

// Vertices to define block North West: (6 14 15 16)
(-calc(rTotF)		calc(rTotF)		zF) vlabel(F15)
(-calc(rTotF)	calc(rLF*sin(pi/4))		zF) vlabel(F16)
// Vertices to define block West : (7 6 16 17)
(-calc(rTotF)		-calc(rLF*sin(pi/4))	zF) vlabel(F17)
// Vertices to define block South West : (19 7 17 18)
(-calc(rLF*cos(pi/4))	-calc(rTotF)		zF) vlabel(F19)
(-calc(rTotF)		-calc(rTotF)		zF) vlabel(F18)
// Vertices to define block South 4: (8 4 7 19)
(calc(rLF*cos(pi/4))	-calc(rTotF)		zF) vlabel(F8)
// Vertices to define block South Est: (9 10 4 8)
(calc(rTotF)		-calc(rLF*sin(pi/4))	zF) vlabel(F10)
(calc(rTotF)		-calc(rTotF)		zF) vlabel(F9)
// Vertices to define block Est: (10 11 5 4)
(calc(rTotF)		calc(rLF*sin(pi/4))	zF) vlabel(F11)
// Vertices to define block North Est: (11 12 13 5)
(calc(rTotF)		calc(rTotF)		zF) vlabel(F12)

);

// Defining blocks:
blocks
(

    // ------------------  Blocks between plane A and plane E  --------------------------

    //Blocks between plane A and plane E:
    // Central Block
    // block0 - positive x O-grid block
    hex (A5 A1 A0 A4 E5 E1 E0 E4 ) AE (nR2 nC zAENC) simpleGrading (rGrading1 1 zGradingAE)
    // block1 - positive y O-grid block
    hex (A6 A2 A1 A5 E6 E2 E1 E5 ) AE (nR2 nC zAENC) simpleGrading (rGrading1 1 zGradingAE)
    // block2 - negative x O-grid block
    hex (A7 A3 A2 A6 E7 E3 E2 E6 ) AE (nR2 nC zAENC) simpleGrading (rGrading1 1 zGradingAE)
    // block3 - negative y O-grid block
    hex (A4 A0 A3 A7 E4 E0 E3 E7 ) AE (nR2 nC zAENC) simpleGrading (rGrading1 1 zGradingAE)
    // block4 - central O-grid block
    hex (A0 A1 A2 A3 E0 E1 E2 E3 ) AE (nC nC zAENC) simpleGrading (1 1 zGradingAE)

    // ------------------  Blocks between plane E and plane F  --------------------------

    //Blocks between plane E and plane F:
    // Central Block
    // block0 - positive x O-grid block
    hex (E5 E1 E0 E4 F5 F1 F0 F4 ) EF (nR2 nC zEFNC) simpleGrading (rGrading1 1 zGradingEF)
    // block1 - positive y O-grid block
    hex (E6 E2 E1 E5 F6 F2 F1 F5 ) EF (nR2 nC zEFNC) simpleGrading (rGrading1 1 zGradingEF)
    // block2 - negative x O-grid block
    hex (E7 E3 E2 E6 F7 F3 F2 F6 ) EF (nR2 nC zEFNC) simpleGrading (rGrading1 1 zGradingEF)
    // block3 - negative y O-grid block
    hex (E4 E0 E3 E7 F4 F0 F3 F7 ) EF (nR2 nC zEFNC) simpleGrading (rGrading1 1 zGradingEF)
    // block4 - central O-grid block
    hex (E0 E1 E2 E3 F0 F1 F2 F3 ) EF (nC nC zEFNC) simpleGrading (1 1 zGradingEF)
		// Second Ring
    hex (Eext5 E5 E4 Eext4 Fext5 F5 F4 Fext4 ) EF (nR3 nC zEFNC) simpleGrading (rGrading2 1 zGradingEF)
    hex (Eext6 E6 E5 Eext5 Fext6 F6 F5 Fext5 ) EF (nR3 nC zEFNC) simpleGrading (rGrading2 1 zGradingEF)
    hex (Eext7 E7 E6 Eext6 Fext7 F7 F6 Fext6 ) EF (nR3 nC zEFNC) simpleGrading (rGrading2 1 zGradingEF)
    hex (Eext4 E4 E7 Eext7 Fext4 F4 F7 Fext7 ) EF (nR3 nC zEFNC) simpleGrading (rGrading2 1 zGradingEF)
 		//--------- ------- BLOCCHI TUTTI ATTORNI  
    // block North:
    hex (Eext6 Eext5 E13 E14 Fext6 Fext5 F13 F14 ) EF (nC nS zEFNC) simpleGrading (1 nsGrading zGradingEF)
    // block North EST:
    hex (Eext5 E11 E12 E13 Fext5 F11 F12 F13 ) EF (nS nS zEFNC) simpleGrading (nsGrading nsGrading zGradingEF)
    // block EST:
    hex (Eext4 E10 E11 Eext5 Fext4 F10 F11 Fext5 ) EF (nS nC zEFNC) simpleGrading (nsGrading 1 zGradingEF)
    // block SOUTH EST:
    //hex (E8 E9 E10 Eext4 F8 F9 F10 Fext4 ) EF (nS nS zEFNC) simpleGrading (nsGrading nsGrading zGradingEF)
    hex (E8 E9 E10 Eext4 F8 F9 F10 Fext4 ) EF (nS nS zEFNC) simpleGrading (nsGrading reverNSG zGradingEF)
    // block SOUTH:
    //hex (E19 E8 Eext4 Eext7 F19 F8 Fext4 Fext7 ) EF (nC nS zEFNC) simpleGrading (1 nsGrading zGradingEF)
    hex (E19 E8 Eext4 Eext7 F19 F8 Fext4 Fext7 ) EF (nC nS zEFNC) simpleGrading (1 reverNSG zGradingEF)
    // block SOUTH WEST:
    //hex (E18 E19 Eext7 E17 F18 F19 Fext7 F17 ) EF (nS nS zEFNC) simpleGrading (nsGrading nsGrading zGradingEF)
    hex (E18 E19 Eext7 E17 F18 F19 Fext7 F17 ) EF (nS nS zEFNC) simpleGrading (reverNSG reverNSG zGradingEF)
    // block WEST:
    //hex (E17 Eext7 Eext6 E16 F17 Fext7 Fext6 F16 ) EF (nS nC zEFNC) simpleGrading (nsGrading 1 zGradingEF)
    hex (E17 Eext7 Eext6 E16 F17 Fext7 Fext6 F16 ) EF (nS nC zEFNC) simpleGrading (reverNSG 1 zGradingEF)
    // block North WEST:
    //hex (E16 Eext6 E14 E15 F16 Fext6 F14 F15 ) EF (nS nS zEFNC) simpleGrading (nsGrading nsGrading zGradingEF)
    hex (E16 Eext6 E14 E15 F16 Fext6 F14 F15 ) EF (nS nS zEFNC) simpleGrading (reverNSG nsGrading zGradingEF)
);

edges
(

    //  - - - -  Plane A  - - - -
    arc A0 A1 (calc(rRelEc*rRelE*r3A) 0 zA)
    arc A1 A2 (0 calc(rRelEc*rRelE*r3A) zA)
    arc A2 A3 (-calc(rRelEc*rRelE*r3A) 0 zA)
    arc A3 A0 (0 -calc(rRelEc*rRelE*r3A) zA)
    arc A4 A5 (r3A 0 zA)
    arc A5 A6 (0 r3A zA)
    arc A6 A7 (-r3A 0 zA)
    arc A7 A4 (0 -r3A zA)

    //  - - - -  Plane E  - - - -
    arc E0 E1 (calc(rRelEc*rRelE*r3E) 0 zE)
    arc E1 E2 (0 calc(rRelEc*rRelE*r3E) zE)
    arc E2 E3 (-calc(rRelEc*rRelE*r3E) 0 zE)
    arc E3 E0 (0 -calc(rRelEc*rRelE*r3E) zE)
    arc E4 E5 (r3E 0 zE)
    arc E5 E6 (0 r3E zE)
    arc E6 E7 (-r3E 0 zE)
    arc E7 E4 (0 -r3E zE)
    arc Eext4 Eext5 (rLE 0 zE)
    arc Eext5 Eext6 (0 rLE zE)
    arc Eext6 Eext7 (-rLE 0 zE)
    arc Eext7 Eext4 (0 -rLE zE)

    //  - - - -  Plane F  - - - -
    arc F0 F1 (calc(rRelFc*rRelF*r3F) 0 zF)
    arc F1 F2 (0 calc(rRelFc*rRelF*r3F) zF)
    arc F2 F3 (-calc(rRelFc*rRelF*r3F) 0 zF)
    arc F3 F0 (0 -calc(rRelFc*rRelF*r3F) zF)
    arc F4 F5 (r3F 0 zF)
    arc F5 F6 (0 r3F zF)
    arc F6 F7 (-r3F 0 zF)
    arc F7 F4 (0 -r3F zF)
    arc Fext4 Fext5 (rLF 0 zF)
    arc Fext5 Fext6 (0 rLF zF)
    arc Fext6 Fext7 (-rLF 0 zF)
    arc Fext7 Fext4 (0 -rLF zF)

);

// Defining patches:
patches
(
    wall inletL //INLET nDec --> fissate come wall 
	(
       (A1 A5 A4 A0)
       (A2 A6 A5 A1)
       (A3 A7 A6 A2)
       (A0 A4 A7 A3)
       (A3 A2 A1 A0)
   )
    patch wallInlet //seconAirInlet  Blocchi attorno
    (
	(Eext5 E5 E4 Eext4)
	(Eext6 E6 E5 Eext5)
	(Eext7 E7 E6 Eext6)
	(Eext4 E4 E7 Eext7)
	//////
	(Eext6 Eext5 E13 E14)
	(Eext5 E11 E12 E13)
	(Eext4 E10 E11 Eext5)
	(E8 E9 E10 Eext4)
	(E19 E8 Eext4 Eext7)
	(E18 E19 Eext7 E17)
	(E17 Eext7 Eext6 E16)
	(E16 Eext6 E14 E15)
	///////
    )



    patch outlet 
    (
       (F1 F5 F4 F0)
       (F2 F6 F5 F1)
       (F3 F7 F6 F2)
       (F0 F4 F7 F3)
       (F3 F2 F1 F0)
	//////
	(Fext5 F5 F4 Fext4)
	(Fext6 F6 F5 Fext5)
	(Fext7 F7 F6 Fext6)
	(Fext4 F4 F7 Fext7)
	//////
	(Fext6 Fext5 F13 F14)
	(Fext5 F11 F12 F13)
	(Fext4 F10 F11 Fext5)
	(F8 F9 F10 Fext4)
	(F19 F8 Fext4 Fext7)
	(F18 F19 Fext7 F17)
	(F17 Fext7 Fext6 F16)
	(F16 Fext6 F14 F15)
    )
    wall lateralInj
    (
       (A4 A5 E5 E4)
       (A5 A6 E6 E5)
       (A6 A7 E7 E6)
       (A7 A4 E4 E7)
    )

    wall lateral
    (	
       (E10 E11 F11 F10)
       (E11 E12 F12 F11)
       (E12 E13 F13 F12)
       (E13 E14 F14 F13)
       (E14 E15 F15 F14)
       (E15 E16 F16 F15)
       (E16 E17 F17 F16)
       (E17 E18 F18 F17)
       (E18 E19 F19 F18)
       (E19 E8 F8 F19)
       (E8 E9 F9 F8)
       (E9 E10 F10 F9)
    )
);

// ************************************************************************* //
