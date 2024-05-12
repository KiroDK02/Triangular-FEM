#include <iostream>
#include "structs.h"

int main()
{
   FEM fem{ };
   TimeMesh time{ };
   SLAE slae{ }, LU{ };
   vector<Conditions> conds{ };
   LOS v{ };

   readMesh(fem);
   readBoundaryCondition(conds);
   
   readTimeMesh(time);
   readSplitTimeMesh(time);

   portraitSparseMatrix(fem, slae);

   getSolution(fem, time, slae, conds);

   vector<double> xSet = { 1.5, 2.5 };
   vector<double> ySet = { 1.5, 2.5 };

   outputResult(fem, time, xSet, ySet, 4);

   //calcGlobalMatrixAndVector(fem, time, slae, cond, 0);
   //
   //calcLU(slae, LU);
   //localOptimalSchemeLU(slae, LU, v, 10000, 1e-14);
   
//   double uNum = uNumerical(fem, 1.5, 1.5, time.q[0]);

   return 0;
}