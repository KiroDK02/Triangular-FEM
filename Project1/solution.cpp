#include "structs.h"

void getSolution(FEM &fem, TimeMesh &time, SLAE &slae, vector<Conditions> &conds)
{
   const int tSize = time.tSplitting.size();

   getWeightsInitU(time, fem);
   time.q[0] = time.qti_2;
   time.q[1] = time.qti_1;

   for (int ti = 2; ti < tSize; ti++)
   {
      SLAE LU{ };
      LOS LOSvectors{ };

      double tValue = time.tSplitting[ti];

      calcGlobalMatrixAndVector(fem, time, slae, conds, ti);
      calcLU(slae, LU);
      localOptimalSchemeLU(slae, LU, LOSvectors, 10000, 1e-14);

      time.qti = slae.q;
      time.saveWeightsTi(ti);
      time.swapVectors();

      clearSLAE(slae);
   }
}