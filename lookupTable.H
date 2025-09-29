volScalarField _Dt(turbulence->nut()/0.4);
volScalarField DEff(turbulence->nuEff());
volVectorField GradZ(fvc::grad(Z));
//const volScalarField &delta = mesh.objectRegistry::lookupObject<volScalarField>("delta");

scalar C_phi = 0.2;

Kai = 2.0 * _Dt * magSqr(fvc::grad(Z));
Zvar1 = C_phi * sqr(delta) * magSqr(fvc::grad(Z)) ;
forAll(ZCells,celli)
{
  Zvar1[celli] = Zvar1[celli]/(Z[celli]*(1-Z[celli])+1e-4);
}
volScalarField Yt(0.0 * Y[0]);

forAll(ZCells, celli)
{
       Kai[celli] = (Kai[celli]>1e-10)? Kai[celli]:1e-10;
       Kai[celli] = std::log10(Kai[celli]);
	   
    if (ignition)
    { 

        pa = {ZCells[celli], Zvar1[celli], C[celli]};
        // pa = {ZCells[celli], 0.0, C[celli]};
        flameletTable.interp(pa);
        Yc[celli] = (flameletTable.lookupC()>1e-5)?C[celli]/flameletTable.lookupC():0.0;
	    pa = {ZCells[celli], Zvar1[celli], Yc[celli]};
	    // pa = {ZCells[celli], 0.0, Yc[celli]};
  
	    flameletTable.interp(pa);
    }
    else
    {
        pa = {ZCells[celli], 0.0, 0.0};
        flameletTable.interp(pa);
    }

    heCells[celli] = 0.0;
    TCells[celli] = flameletTable.lookupT();

    rho_lookCells[celli] = flameletTable.lookupRho();
    Source_CCells[celli] = flameletTable.lookupSourceC();
    forAll(Y, i)
    {
        scalarField &YCells = Y[i].primitiveFieldRef();
        YCells[celli] = flameletTable.lookupY(i);
        Yt[celli] += Y[i][celli];
    }

    forAll(Y, i)
    {
        scalarField &YCells = Y[i].primitiveFieldRef();
        if(Yt[celli] == 0.0)
        {
          Info << " Yt[celli]==0.0" << endl;
        }
        YCells[celli] /= Yt[celli];
        heCells[celli] += YCells[celli] * thermo.composition().Hs(i, pCells[celli], TCells[celli]);
    }
}
Zvar1.correctBoundaryConditions();
Zvar.correctBoundaryConditions();
T.correctBoundaryConditions(); 
/*
forAll(Y, i)
{
    Y[i].correctBoundaryConditions();
}
*/
forAll(mesh.boundary(), patchI)
{
    const fvPatchScalarField &pp = p.boundaryField()[patchI];

    fvPatchScalarField &phe = he.boundaryFieldRef()[patchI];
    fvPatchScalarField &pT = T.boundaryFieldRef()[patchI];
    fvPatchScalarField &pYt = Yt.boundaryFieldRef()[patchI];
    fvPatchScalarField &pRho = rho_look.boundaryFieldRef()[patchI];
    fvPatchScalarField &pSourceC = Source_C.boundaryFieldRef()[patchI];
    fvPatchScalarField &pZ = Z.boundaryFieldRef()[patchI];
    fvPatchScalarField &pC = C.boundaryFieldRef()[patchI];
    fvPatchScalarField &pYc = Yc.boundaryFieldRef()[patchI];
    fvPatchScalarField &pZvar = Zvar1.boundaryFieldRef()[patchI];
    

    forAll(mesh.boundary()[patchI], facei)
    {
        phe[facei] = 0.0;
          pa = {pZ[facei], 0.0, 0.0};
          flameletTable.interp(pa);
        pRho[facei] = flameletTable.lookupRho();
        pSourceC[facei] = flameletTable.lookupSourceC();      
        forAll(Y, i)
        {
           
            fvPatchScalarField &pY = Y[i].boundaryFieldRef()[patchI];
            pYt[facei] += pY[facei];
            phe[facei] += pY[facei] * thermo.composition().Hs(i, pp[facei], pT[facei]);
        }

    }
}

thermo.correct();
T_inter=T;
/*
forAll(ZCells, celli)
{
TCells[celli] = flameletTable.lookupT();
}
*/
Info << "==================================================================" << endl;
Info << "U max: " << gMax(U) << endl;
Info << "Z Vol: " << Z.weightedAverage(mesh.V()).value() << endl;
Info << "Zvar min/max :" << gMin(Zvar) << "," << gMax(Zvar) << endl;
Info << "C min/max :" << gMin(C) << "," << gMax(C) << endl;
Info << "rho min/max : " << gMin(rho) << ", " << gMax(rho) << endl;
Info << "T min/max: " << gMin(T) << ", " << gMax(T) << endl;
forAll (Y,i)
{
Info << "Y[" << i << "] min/max: " << gMin(Y[i]) << ", " << gMax(Y[i]) << endl;
}
Info << "p min/max: " << gMin(p) << ", " << gMax(p) << endl;
// Info << "Dt min/max: " << std::log10(gMin(_Dt)) << ", " << std::log10(gMax(_Dt)) << endl;
Info << "==================================================================" << endl;
