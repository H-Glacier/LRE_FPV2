volScalarField &he = thermo.he();
const scalarField &ZCells = Z.internalField();
const scalarField &pCells = p.internalField();
scalarField &heCells = he.primitiveFieldRef();
scalarField &TCells = T.primitiveFieldRef();
//scalarField &T_lookCells = T_look.primitiveFieldRef();
scalarField &rhoCells = rho.primitiveFieldRef();
scalarField &rho_lookCells = rho_look.primitiveFieldRef();
scalarField &Source_CCells = Source_C.primitiveFieldRef();
// scalarField &DtCells = Dt.primitiveFieldRef();

bool ignition = true;

size_t Table_D = 3;
std::vector<double> pa(Table_D, 0.0);

Info<< "Reading flamelet tables" << endl;
// flameletLib flameletTable(turbulenceModelType);
std::string path = "tables";
std::ifstream fin(path+"/chemTable",std::ios::in);
if(!fin)
{
    throw std::runtime_error("cannot open chemTable");
}

flameletLib flameletTable(path+"/chemTable");
