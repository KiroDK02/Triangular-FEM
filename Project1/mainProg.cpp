#include <iostream>
#include "structs.h"

int main()
{
   FEM fem { };
   TimeMesh time { };
   SLAE slae { }, LU { };
   vector<Conditions> conds { };
   LOS v { };

   readMesh(fem);
   readBoundaryCondition(conds);

   readTimeMesh(time);
   readSplitTimeMesh(time);

   portraitSparseMatrix(fem, slae);

   //getSolutionParabolicProblem(fem, time, slae, conds);
   getSolutionHyperbolicProblem(fem, time, slae, conds);

   vector<double> xSet { };
   vector<double> ySet { };
   vector<double> tSet { };

   double x0 = 1.1;
   double x1 = 2.9;
   double y0 = 1.1;
   double y1 = 2.9;
   double t0 = 1.1;
   double t1 = 3.9;

   int xSize = 100;
   int ySize = 100;
   int tSize = 100;

   double hx = (x1 - x0) / xSize;
   double hy = (y1 - y0) / ySize;
   double ht = (t1 - t0) / tSize;

   xSet.resize(xSize + 1);
   ySet.resize(ySize + 1);
   tSet.resize(tSize + 1);

   xSet[0] = x0;
   for (int i = 0; i < xSize; i++)
      xSet[i] = x0 + i * hx;
   xSet[xSize] = x1;

   ySet[0] = y0;
   for (int i = 0; i < ySize; i++)
      ySet[i] = y0 + i * hy;
   ySet[ySize] = y1;

   tSet[0] = t0;
   for (int i = 0; i < tSize; i++)
      tSet[i] = t0 + i * ht;
   tSet[tSize] = t1;

   //outputResult(fem, time, xSet, ySet, 3);
   //std::cout << uNumericalAnyTime(fem, time, 1.5, 2.5, 1.5);

   std::cout << "Error: " << calcErrorNumSolution(fem, time, xSet, ySet, tSet);

   //calcGlobalMatrixAndVector(fem, time, slae, cond, 0);
   //
   //calcLU(slae, LU);
   //localOptimalSchemeLU(slae, LU, v, 10000, 1e-14);

//   double uNum = uNumerical(fem, 1.5, 1.5, time.q[0]);

   return 0;
}