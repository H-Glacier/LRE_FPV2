#include "flameletLib.H"

/* 
 * Constructor of the class
 */
flameletLib::flameletLib(std::string loc_table)
{
   // Get the contents of file in a vector
   std::cout << "Reading NTF table.........." << std::endl;
   getFileContent(loc_table);
   
   // Initialize the vectors
   _D.resize(nD, std::vector<size_t>(n, 0));
   nW.resize(nD, std::vector<size_t>(n, 0));
   wD.resize(nD, std::vector<double>(2, 0.0));
   D0.resize(nD, 0);
   weight.resize(n, 1.0);
 //hq  
   _D_C.resize(nD-1, std::vector<size_t>(n_C, 0));
   weight_C.resize(n_C, 1.0);
 //hq  

   // prepare a binary list for the interpolating
   for (size_t i = 0; i < n; ++i)
   {
      size_t index = 0;
      size_t tmp = i;
      while (index < nD)
      {
         tmp % 2 ? nW[index][i] = 1 : nW[index][i]= 0;
         tmp /= 2;
         ++index;
      }
   }
}

/* 
 * Destructor of the class
 */
flameletLib::~flameletLib() {}

/*
 * It will iterate through all the lines and
 * put the properties and values of the table into given vectors
 */
void flameletLib::getFileContent(std::string fileName)
{
   // Open the File
   std::ifstream in(fileName.c_str(),std::ios::in|std::ios::binary);
   // Check if object is valid
   if (!in)
   {
      std::cerr << "Cannot open the File : " << fileName << std::endl;
      exit (0);
   }

   std::string str;

   // Read the dimensions of the table
   getline(in, str);
   nD = int(stod(str)) - 1;
   n = std::pow(2.0, nD);
   n_C = std::pow(2.0, nD-1);

   // Read the numbers of each variable
   sD.resize(nD, 0);
   s_prod.resize(nD + 1, 1);
   for (size_t i = 0; i < nD; ++i)
   {
      getline(in, str);
      sD[i] = int(stod(str));
      if (i == 0)
      {
         s_prod[i] *= sD[i];
      }
      else
      {
         s_prod[i] *= sD[i] * s_prod[i - 1];
      }   
   }
   getline(in, str);
   nS = int(stod(str));
   s_prod[nD] *= nS * s_prod[nD - 1];

   // Read the distribution of each variable
   for (size_t i = 0; i < nD; ++i)
   {
      iD.push_back(std::vector<size_t>());
      pD.push_back(std::vector<double>());
      iD[i].push_back(0);

      for (size_t j = 0; j < sD[i]; ++j)
      {
         getline(in, str);
         pD[i].push_back(stod(str));
         
         if (j > 1 && fabs(pD[i][j] - 2.0 * pD[i][j - 1] + pD[i][j - 2]) > 1e-10)
         {
            iD[i].push_back(j - 1);
         }
      }
      iD[i].push_back(sD[i] - 1);
      siD.push_back(iD[i].size());
   }

   for (size_t i = 0; i < nS; ++i)
   {
      getline(in, str);
      pS.push_back(stod(str));
   }
//hq   
   Max_C.reserve(s_prod[nD-2]);
//   std::cout << "s_prod[nD-2].........." << s_prod[nD-2]<<std::endl;
   for (size_t i = 0; i < s_prod[nD-2]; ++i)
   {
      getline(in, str);
      Max_C.push_back(stod(str));      
   }
//hq
   // Read the values of the scalars corresponds to the variables
   Table.reserve(s_prod[nD]);
   while (getline(in, str))
   {
      Table.push_back(stod(str));
   }

   //Close The File
   in.close();
}

/*
 * Find the nearby points and their interpolating weights
 * in the table corresponds to the given variables
 */
void flameletLib::interp(std::vector<double> pa)
{
   if (pa.size() != nD)
   {
      std::cerr << "Input dimension is " << pa.size() 
                << ", which disagree with " << nD << " of the table." 
                << std::endl;
      exit (0);
   }

   std::fill(weight.begin(), weight.end(), 1.0);
   std::fill(weight_C.begin(), weight_C.end(), 1.0);

   size_t tmpI;
   double tmpx, tmpd;
   for (size_t i = 0; i < nD; ++i)
   {
      pa[i] = (pa[i] > pD[i].back()) ? pD[i].back() : pa[i];
      pa[i] = (pa[i] < pD[i][0]) ? pD[i][0] : pa[i];

      for (size_t j = 1; j < siD[i]; ++j)
      {
         if (pa[i] <= pD[i][iD[i][j]])
         {
            tmpI = j - 1;
            break;
         }
      }
      tmpd = (pD[i][iD[i][tmpI + 1]] - pD[i][iD[i][tmpI]]) 
           / (iD[i][tmpI + 1] - iD[i][tmpI]) + 1e-10;
      tmpx = (pa[i] - pD[i][iD[i][tmpI]]) / tmpd;
      D0[i] = floor(tmpx) + iD[i][tmpI];
      wD[i][1] = tmpx - floor(tmpx);
      wD[i][0] = 1.0 - wD[i][1];
   }

   for (size_t i = 0; i < nD; ++i)
   {
    //	  std::cout << "D0 = " << D0[i] <<std::endl;
    // std::cout << "D0 = " << D0[i] << std::endl;
      for (size_t j = 0; j < n; ++j)
      {
         _D[i][j] = D0[i] + nW[i][j];
         weight[j] *= wD[i][nW[i][j]];
      }
   }

//hq   
   for (size_t i = 0; i < nD-1; ++i)
   {
      for (size_t j = 0; j < n_C; ++j)
      {
         _D_C[i][j] = D0[i] + nW[i][j];
         weight_C[j] *= wD[i][nW[i][j]];
      }
   }
//hq   
}

double flameletLib::lookupC()
{
   double tmpC = 0.0;
   for (size_t i = 0; i < n_C; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD-1; ++j)
      {      
         loc += _D_C[j][i] * s_prod[nD-2] / s_prod[j];
      }
      tmpC += weight_C[i] * Max_C[loc];
    //  std::cout << "Max_C = " << Max_C[loc] <<std::endl;
   }
    //  std::cout << "tmpC = " << tmpC <<std::endl;
   return tmpC;
}

/* 
 * Interpolating the temperature corresponds to the given variables
 */
double flameletLib::lookupT()
{
   double tmpT = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {      
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
 //        std::cout << "i =  " << i<<"j="<< j<< "  "<<_D[j][i] << std::endl;
      }
      tmpT += weight[i] * Table[loc];
 //     std::cout << "TT  " << i << "  weight[i] = "<< weight[i] << "  TT = "<< Table[loc] << std::endl;
   }

   return tmpT;
}

/* 
 * Interpolating the temperature corresponds to the given variables
 */

double flameletLib::lookupRho()
{
   double tmpRho = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {      
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
      }
      loc += pS[nS-3];
      tmpRho += weight[i] * Table[loc];
   }

   return tmpRho;
}

double flameletLib::lookupAlpha()
{
   double tmpAlpha = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {      
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
      }
      loc += pS[nS-2];
      tmpAlpha += weight[i] * Table[loc];
   }

   return tmpAlpha;
}
//Interpolating the source term of C corresponds to the given variables
 
double flameletLib::lookupSourceC()
{
   double tmpSourceC = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {      
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
      }
      loc += pS[nS-4];
      tmpSourceC += weight[i] * Table[loc];
   }

   return tmpSourceC;
}

// Interpolating the evaporation source term corresponds to the given variables
 
double flameletLib::lookupSourceE()
{
   double tmpSourceE = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {      
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
      }
      loc += pS[nS-4];
      tmpSourceE += weight[i] * Table[loc];
   }

   return tmpSourceE;
}

/* 
 * Interpolating the specie concentration corresponds to the given variables
 */
double flameletLib::lookupY(size_t Yi)
{
   double tmpY = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {      
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
      }
      loc += pS[Yi];
      tmpY += weight[i] * Table[loc];
   }

   if (tmpY < 0.0)
   {
      return 0.0;
   }
   else
   {
      return tmpY;
   }
}

/* 
 * print the location, weight and value of each neibouring point, for testing purpose.
 */
void flameletLib::print_index(size_t S)
{
   double tmpS = 0.0;
   for (size_t i = 0; i < n; ++i)
   {
      size_t loc = 0;
      for (size_t j = 0; j < nD; ++j)
      {
         loc += _D[j][i] * s_prod[nD] / s_prod[j];
    //     std::cout << _D[j][i] << " ";
      }
      loc += pS[S] - 1;
      tmpS += weight[i] * Table[loc];
    //  std::cout << weight[i] << " " << Table[loc] << std::endl;
   }
    //  std::cout << tmpS << std::endl;
}

// int main(int argc, char *argv[])
// {
   
//    if (argc == 2)
//    {
//       flameletLib NTF("../../Table_flamelet_4D.txt");
//       std::vector<double> p{0.31, 1.88, -3.55, -4.12};
//       NTF.interp(p);
//       // double T = NTF.lookupT();
//       size_t S = std::stoi(argv[1]);
//       NTF.print_index(S);
//    }

//    return 0;
// }
