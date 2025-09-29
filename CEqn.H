// Dt = thermophysicalTransport->alphaEff();
Dt = thermo.alpha() + turbulence->mut()/0.7;

// Dt = turbulence->muEff();
fvScalarMatrix CEqn
(
    fvm::ddt(rho, C)
  + fvm::div(phi, C)
  - fvm::laplacian(Dt, C)
  == Source_C
);


CEqn.relax();
fvOptions.constrain(CEqn);
CEqn.solve();

C.max(0.0);
C.min(1.0);

C.correctBoundaryConditions();


