// Dt = thermophysicalTransport->alphaEff();

volScalarField _Dt1(turbulence->nut());
volVectorField GradZ1(fvc::grad(Z));
volScalarField NGradZ1(magSqr(fvc::grad(Z)));


forAll (zVarTmp,celli)
{
eta[celli] = Z[celli]*(1.0-Z[celli]);
zVarTmp[celli] = (Zvar[celli] > eta[celli])? eta[celli] : Zvar[celli];
zVarTmp[celli] = (zVarTmp[celli] > 0.0)? zVarTmp[celli] : 0.0;
}

scalar cPhi = 2.0;
scalar cEps_cU = 2.0;
scalar ScTurb = 0.4 ; 



forAll (zVarTmp,celli)
{
  Source1_Z[celli] = 2.0*rho[celli]*_Dt1[celli]*NGradZ1[celli];

  Source2_Z[celli] = rho[celli]*cPhi * cEps_cU / ScTurb * _Dt1[celli] * zVarTmp[celli]  \
                    / std::pow(mesh.V()[celli],2.0/3.0);
}

Dt = turbulence->muEff();
/*
fvScalarMatrix ZvarEqn
(
    fvm::ddt(rho, Zvar)
  + fvm::div(phi, Zvar)
  - fvm::laplacian(Dt, Zvar)
  == Source1_Z - Source2_Z
);
*/
        //const volScalarField& delta = turbulence->delta();
        //const volScalarField &delta = mesh.objectRegistry::lookupObject<volScalarField>("delta");
        volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("delta", dimLength, 0.0)
);

// Calculate delta (typically as cell size measure)
delta.ref() = pow(mesh.V(), 1.0/3.0);  // Cube root of cell volume

Zvar =0;
/*
        fvScalarMatrix ZvarEqn(
            fvm::ddt(rho, Zvar) + fvm::div(phi, Zvar) - fvm::laplacian(turbulence->alphaEff(), Zvar)
            + 4.0*turbulence->mut()*Zvar/(delta*delta)
            - 2.0*turbulence->alphat()*magSqr(fvc::grad(Z))
            );

Info<< "OK4"<< nl << endl;
ZvarEqn.relax();

ZvarEqn.solve();

Zvar.max(0.0);


forAll (zVarTmp,celli)
{
Zvar[celli] = (Zvar[celli] > eta[celli])? eta[celli] : Zvar[celli];
}

Zvar.correctBoundaryConditions();

Zvar =0;
*/

 
