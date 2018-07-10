// Parametrized case by Giovanni Tretola

//Run using:
//m4 -P system/blockMeshDict.m4 > system/blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926)

//Geometry
m4_define(dinj,0.1)	// [mm] -->100micron

//m4_define(dTot,calc(dinj*20))
m4_define(dTot,calc(dinj*20))
m4_define(rTot,calc(dTot/2))
m4_define(rRel,0.8)
m4_define(rL,3)
//m4_define(rL,4)

//Geometry
m4_define(rinj,calc(dinj/2))

//m4_define(r1,calc((dinj/5)/2))

m4_define(r2,calc(r1*3))
m4_define(r3,calc(rinj))
m4_define(rTot,calc(dTot/2))

m4_define(length, calc(dinj*50))
//m4_define(length, calc(dinj*35))

//Grid points (integers!):

//Vedi schema 
//m4_define(nC, 10)//5)	//numero di Liquid Inlet	//al massimo 6
m4_define(nC, 10)//5)	//numero di Liquid Inlet -- Dx 0.01 mm
m4_define(nR2, 30)	//
m4_define(zEFNC, 250)

//Grading
//m4_define(rGrading1, 1)
m4_define(rGrading1, 0.1394)
//m4_define(rGrading2, 1)//--> Ho un secondo cerchio nella sezione di ingresso per non forzare una asimmetria nel getto. 
m4_define(rGrading2, 7.1753)//1/7.1753	//	1/0.1394
m4_define(zGradingEF, 1)

//Plane E:
m4_define(zE, calc(0))
m4_define(r1E, rinj)
m4_define(r2E, rTot)
m4_define(rRelEc, 1)
//m4_define(rRelEc, 0.8)

//Plane F:
m4_define(zF, calc(length))
m4_define(r1F, rinj)
m4_define(r2F, rTot)
m4_define(rRelFc, 1)
//m4_define(rRelFc, 0.8)

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


//Plane E:
//Central O-grid
(calc(r1E*cos(pi/4)) -calc(r1E*sin(pi/4)) zE) vlabel(E0)
(calc(r1E*cos(pi/4)) calc(r1E*sin(pi/4)) zE) vlabel(E1)
(calc(-r1E*cos(pi/4)) calc(r1E*sin(pi/4)) zE) vlabel(E2)
(calc(-r1E*cos(pi/4)) -calc(r1E*sin(pi/4)) zE) vlabel(E3)
//4 block rounding
(calc(r2E*cos(pi/4)) -calc(r2E*sin(pi/4)) zE) vlabel(E4)
(calc(r2E*cos(pi/4)) calc(r2E*sin(pi/4)) zE) vlabel(E5)
(calc(-r2E*cos(pi/4)) calc(r2E*sin(pi/4)) zE) vlabel(E6)
(calc(-r2E*cos(pi/4)) -calc(r2E*sin(pi/4)) zE) vlabel(E7)

//	plane F
//Central O-grid
(calc(r1F*cos(pi/4)) -calc(r1F*sin(pi/4)) zF) vlabel(F0)
(calc(r1F*cos(pi/4)) calc(r1F*sin(pi/4)) zF) vlabel(F1)
(calc(-r1F*cos(pi/4)) calc(r1F*sin(pi/4)) zF) vlabel(F2)
(calc(-r1F*cos(pi/4)) -calc(r1F*sin(pi/4)) zF) vlabel(F3)
//4 block rounding
(calc(r2F*cos(pi/4)) -calc(r2F*sin(pi/4)) zF) vlabel(F4)
(calc(r2F*cos(pi/4)) calc(r2F*sin(pi/4)) zF) vlabel(F5)
(calc(-r2F*cos(pi/4)) calc(r2F*sin(pi/4)) zF) vlabel(F6)
(calc(-r2F*cos(pi/4)) -calc(r2F*sin(pi/4)) zF) vlabel(F7)

);

// Defining blocks:
blocks
(

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

);

edges
(

    //  - - - -  Plane E  - - - -
    arc E0 E1 (calc(rRelEc*r1E) 0 zE)
    arc E1 E2 (0 calc(rRelEc*r1E) zE)
    arc E2 E3 (-calc(rRelEc*r1E) 0 zE)
    arc E3 E0 (0 -calc(rRelEc*r1E) zE)
    arc E4 E5 (r2E 0 zE)
    arc E5 E6 (0 r2E zE)
    arc E6 E7 (-r2E 0 zE)
    arc E7 E4 (0 -r2E zE)

    //  - - - -  Plane F  - - - -
    arc F0 F1 (calc(rRelFc*r1F) 0 zF)
    arc F1 F2 (0 calc(rRelFc*r1F) zF)
    arc F2 F3 (-calc(rRelFc*r1F) 0 zF)
    arc F3 F0 (0 -calc(rRelFc*r1F) zF)
    arc F4 F5 (r2F 0 zF)
    arc F5 F6 (0 r2F zF)
    arc F6 F7 (-r2F 0 zF)
    arc F7 F4 (0 -r2F zF)

);

// Defining patches:
patches
(
    wall inletL //INLET nDec --> fissate come wall 
	(
      /* (E1 E5 E4 E0)
       (E2 E6 E5 E1)
       (E3 E7 E6 E2)
       (E0 E4 E7 E3)*/
       (E3 E2 E1 E0)
   )
    patch wallInlet //seconAirInlet  Blocchi attorno
    (
       (E1 E5 E4 E0)
       (E2 E6 E5 E1)
       (E3 E7 E6 E2)
       (E0 E4 E7 E3)
    )



    patch outlet 
    (
       (F1 F5 F4 F0)
       (F2 F6 F5 F1)
       (F3 F7 F6 F2)
       (F0 F4 F7 F3)
       (F3 F2 F1 F0)
    )

    wall lateral
    (	
       (E4 E5 F5 F4)
       (E5 E6 F6 F5)
       (E6 E7 F7 F6)
       (E7 E4 F4 F7)
    )
);

// ************************************************************************* //
