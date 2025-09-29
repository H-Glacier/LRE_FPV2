// Dt = thermophysicalTransport->alphaEff();
Dt = thermo.alpha() + turbulence->mut()/0.7;

// Dt = turbulence->muEff();
fvScalarMatrix ZEqn
(
    fvm::ddt(rho, Z)
  + fvm::div(phi, Z)
  - fvm::laplacian(Dt, Z)
);


ZEqn.relax();
fvOptions.constrain(ZEqn);
ZEqn.solve();

Z.max(0.0);
Z.min(1.0);

Z.correctBoundaryConditions();

rhoUZ = rho * U * Z;

