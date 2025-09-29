#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <cstdlib>

class flameletLib
{
public:

    flameletLib(std::string loc_table);

    ~flameletLib();

    void getFileContent(std::string fileName);

    void interp(std::vector<double> pa);

    void interp_C(std::vector<double> pa_C);

    double lookupT();

    double lookupRho();

    double lookupAlpha();
    
    double lookupC();
    
    double lookupSourceC();

    double lookupSourceE();
    
    double lookupY(size_t Yi);

    void print_index(size_t n);
    
private:

    size_t nD, n, nS,n_C;
    std::vector<std::vector<size_t>> _D, iD, nW,_D_C;
    std::vector<std::vector<double>> pD, wD;
    std::vector<size_t> sD, siD, D0, s_prod;
    std::vector<double> pS, Table, weight, Max_C,weight_C;

};
