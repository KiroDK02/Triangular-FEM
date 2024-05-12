#include "structs.h"
#include "matrices.h"

// Начальное условие
// их два: на 0ом и на 1ом временных слоях
// t[0] - нулевой временной слой
// t[1] - первый временной слой
// В эту функцию передавайте прям настоящую
// По ней начальное условие заберется
std::function<double(double, double, double)> uInit
{
   [](double x, double y, double t) { return x * x + y * y + t * t; } 
};

// Функция правой части
std::function<double(double, double, double)> F
{
   [](double x, double y, double t) { return -4 + 2 * (x + y) * t; }
};

// Лямбда. Только от времени.
std::function<double(double)> lambda
{
   [](double t) { return 1; }
};

// гамма из прошлого курсача превратилась в сигму
// я оставил, чтоб она раскладывалась по лин-му базису
// в пространстве, время просто передается с временного слоя
std::function<double(double, double, double)> sigma
{
   [](double x, double y, double t) { return x + y; }
};

// Функции первых краевых
vector<std::function<double(double, double, double)>> ug 
{
   [](double x, double y, double t) { return 1 + y * y + t * t; }, // 1
   [](double x, double y, double t) { return x * x + 9 + t * t; }, // 2
   [](double x, double y, double t) { return 9 + y * y + t * t; }, // 3
   [](double x, double y, double t) { return x * x + 1 + t * t; } // 4
};

// Функции вторых краевых
vector<std::function<double(double, double, double)>> tetta
{
   [](double x, double y, double t) { return -2; }, // 1
   [](double x, double y, double t) { return x * x + 4; }, // 2
   [](double x, double y, double t) { return 4 + y * y; }, // 3
   [](double x, double y, double t) { return x * x; } // 4
};

void outputResult(FEM &fem, TimeMesh &time, vector<double> &xSet, vector<double> &ySet, int ti)
{
   const int tSize = time.tSplitting.size();
   const int xSize = xSet.size();
   const int ySize = ySet.size();

   vector<double> q = time.q[ti];
   double tValue = time.tSplitting[ti];

   FILE *result = NULL;
   fopen_s(&result, "result.txt", "w");

   if (result == NULL)
      exit(EXIT_FAILURE);

   fprintf_s(result, "\t\tx\t\t\t\t  y\t\t\t\t\tt\t\t\t\t   u  \t\t\t\t  u*\t\t\t|u - u*|\n");

   for (int s = 0; s < ySize; s++)
      for (int p = 0; p < xSize; p++)
      {
         double x = xSet[p];
         double y = ySet[s];

         double uNum = uNumerical(fem, x, y, q);
         double uReal = uInit(x, y, tValue);

         double difference = abs(uNum - uReal);

         fprintf_s(result, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", x, y, tValue, uNum, uReal, difference);
      }

   fclose(result);
}

void calcGlobalMatrixAndVector(FEM &fem, TimeMesh &time, SLAE &slae, vector<Conditions> &cond, int ti)
{
   auto &elements = fem.elements;
   auto &vertexCoord = fem.vertexCoord;
   
   auto &A = slae.A;
   auto &b = slae.b;

   double tValue = time.tSplitting[ti];
   double deltaT = time.tSplitting[ti] - time.tSplitting[ti - 2];
   double deltaT0 = time.tSplitting[ti] - time.tSplitting[ti - 1];
   double deltaT1 = time.tSplitting[ti - 1] - time.tSplitting[ti - 2];

   const int sizeSlae = vertexCoord.size();
   
   vector<vector<double>> stiffness{ };
   vector<vector<double>> mass{ };
   vector<vector<double>> tempMatrix{ };
   vector<double> bLocal{ };
   vector<double> tempVector{ };

   initSizeLocalMatrix(tempMatrix);
   initSizeLocalMatrix(stiffness);
   initSizeLocalMatrix(mass);
   bLocal.resize(6);
   tempVector.resize(6);
   
   for (auto &elem : fem.elements)
   {
      calcSigmaWeights(fem.vertexCoord, elem, tValue);

      calcLocalStiffnessMatrix(elem, stiffness, tValue);
      calcLocalMassMatrix(elem, mass);
      calcLocalb(vertexCoord, elem, bLocal, tValue);

      addLocalMatrixInGlobal(A, elem.localVertex, stiffness);
      addLocalVectorInGlobal(b, elem.localVertex, bLocal);

      multiplyMatrixToCoef(mass, (deltaT + deltaT0) / (deltaT * deltaT0), tempMatrix);
      addLocalMatrixInGlobal(A, elem.localVertex, tempMatrix);

      multiplyMatrixToVector(mass, time.qti_2, tempVector, elem.localVertex);
      multiplyVectorCoef(tempVector, -deltaT0 / (deltaT * deltaT1));
      addLocalVectorInGlobal(b, elem.localVertex, tempVector);

      multiplyMatrixToVector(mass, time.qti_1, tempVector, elem.localVertex);
      multiplyVectorCoef(tempVector, deltaT / (deltaT1 * deltaT0));
      addLocalVectorInGlobal(b, elem.localVertex, tempVector);
   }

   addCondSecond(slae, cond, vertexCoord, tValue);
   addCondFirst(slae, cond, vertexCoord, tValue);
}

void addCondFirst(SLAE &slae, vector<Conditions> &cond, vector<Point> &coord, double tValue)
{
   auto &A = slae.A;
   auto &b = slae.b;

   for (auto &condition : cond)
      if (condition.typeOfCond == 1)
      {
         auto &globalNumVertex = condition.globalNumVertex;
         int numOfFunction = condition.numOfFunction;

         const int numVertex = globalNumVertex.size();

         for (int i = 0; i < numVertex; i++)
         {
            int globalNum = globalNumVertex[i];
            stringMatrixInNull(A, globalNum);
            b[globalNum] = ug[numOfFunction](coord[globalNum].x, coord[globalNum].y, tValue);
         }
   }
}

void addCondSecond(SLAE &slae, vector<Conditions> &conds, vector<Point> &coords, double tValue)
{
   auto &A = slae.A;
   auto &b = slae.b;

   vector<double> bLocal(3);

   for (auto &condition : conds)
   {
      if (condition.typeOfCond == 2)
      {
         auto &globalNumVertex = condition.globalNumVertex;
         int numOfFunction = condition.numOfFunction;

         vector<int> localVertex{ globalNumVertex[0], globalNumVertex[1], globalNumVertex[2] };

         Point point1 = coords[globalNumVertex[0]];
         Point point2 = coords[globalNumVertex[1]];
         Point point3 = coords[globalNumVertex[2]];

         double hm = getLengthOfEdge(point1.x, point3.x, point1.y, point3.y);

         bLocal[0] = AS3[0][0] * tetta[numOfFunction](point1.x, point1.y, tValue) +
                     AS3[0][1] * tetta[numOfFunction](point2.x, point2.y, tValue) +
                     AS3[0][2] * tetta[numOfFunction](point3.x, point3.y, tValue);
         
         bLocal[1] = AS3[1][0] * tetta[numOfFunction](point1.x, point1.y, tValue) +
                     AS3[1][1] * tetta[numOfFunction](point2.x, point2.y, tValue) +
                     AS3[1][2] * tetta[numOfFunction](point3.x, point3.y, tValue);
         
         bLocal[2] = AS3[2][0] * tetta[numOfFunction](point1.x, point1.y, tValue) +
                     AS3[2][1] * tetta[numOfFunction](point2.x, point2.y, tValue) +
                     AS3[2][2] * tetta[numOfFunction](point3.x, point3.y, tValue);

         multiplyVectorCoef(bLocal, hm / 30.);

         addLocalVectorInGlobal(b, localVertex, bLocal);
      }
   }
}

void calcDetD(vector<Point> &coords, FiniteElement &elem)
{
   auto &localVertex = elem.localVertex;
   auto &detD = elem.detD;

   auto point1 = coords[localVertex[0]];
   auto point2 = coords[localVertex[1]];
   auto point3 = coords[localVertex[2]];

   detD = (point2.x - point1.x) * (point3.y - point1.y) - (point3.x - point1.x) * (point2.y - point1.y);
}

void calcAlpha(vector<Point> &coords, FiniteElement &elem)
{
   auto &alpha = elem.alpha;
   auto &localVertex = elem.localVertex;

   auto point1 = coords[localVertex[0]];
   auto point2 = coords[localVertex[1]];
   auto point3 = coords[localVertex[2]];

   alpha[0] = point2.x * point3.y - point3.x * point2.y;
   alpha[1] = point2.y - point3.y;
   alpha[2] = point3.x - point2.x;
   alpha[3] = point3.x * point1.y - point1.x * point3.y;
   alpha[4] = point3.y - point1.y;
   alpha[5] = point1.x - point3.x;
   alpha[6] = point1.x * point2.y - point2.x * point1.y;
   alpha[7] = point1.y - point2.y;
   alpha[8] = point2.x - point1.x;

   multiplyVectorCoef(alpha, 1. / elem.detD);
}

void calcSigmaWeights(vector<Point> &coords, FiniteElement &elem, double tValue)
{
   auto &localVertex = elem.localVertex;
   auto &sigmaWeights = elem.sigmaWeights;

   vector<Point> points = { coords[localVertex[0]], coords[localVertex[1]], coords[localVertex[2]] };

   for (int i = 0; i < 3; i++)
      sigmaWeights[i] = sigma(points[i].x, points[i].y, tValue);
}

void portraitSparseMatrix(FEM &fem, SLAE &slae)
{
   auto &vertexCoord = fem.vertexCoord;
   auto &elements = fem.elements;
   auto &ig = slae.A.ig, &jg = slae.A.jg;
   auto &b = slae.b;

   vector<std::set<int>> rowCount{ };

   const int sizeSlae = vertexCoord.size();
   
   slae.A.di.resize(sizeSlae);
   slae.b.resize(sizeSlae);
   slae.q.resize(sizeSlae);
   ig.resize(sizeSlae + 1);
   rowCount.resize(sizeSlae);

   for (auto &elem : elements)
      for (auto i : elem.localVertex)
         for (auto j : elem.localVertex)
            if (j < i) rowCount[i].insert(j);
   
   ig[0] = 0;
   for (int i = 0; i < sizeSlae; i++)
   {
      ig[i + 1] = ig[i] + rowCount[i].size();
      for (auto j : rowCount[i])
         jg.push_back(j);
   }

   slae.A.ggl.resize(ig[sizeSlae]);
   slae.A.ggu.resize(ig[sizeSlae]);
}

void addLocalMatrixInGlobal(SparseMatrix &A, vector<int> &localVertex, vector<vector<double>> &localMatrix)
{
   const int sizeLocal = 6;

   for (int i = 0; i < sizeLocal; i++)
      for (int j = 0; j < sizeLocal; j++)
      {
         double elem = localMatrix[i][j];
         addElemInGlobalMatrix(A, localVertex[i], localVertex[j], elem);
      }
}

void addLocalVectorInGlobal(vector<double> &b, vector<int> &localVertex, vector<double> &bLocal)
{
   const int sizeLocal = bLocal.size();

   for (int i = 0; i < sizeLocal; i++)
      b[localVertex[i]] += bLocal[i];
}

void addElemInGlobalMatrix(SparseMatrix &A, int i, int j, double elem)
{
   auto &ig = A.ig, &jg = A.jg;
   auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di;

   if (i == j) di[i] += elem;
   else if (i > j)
   {
      int beg = ig[i], end = ig[i + 1] - 1;
      while (jg[beg] != j)
      {
         int ind = (beg + end) / 2;
         if (jg[ind] < j)
            beg = ind + 1;
         else
            end = ind;
      }
      ggl[beg] += elem;
   }
   else
   {
      int beg = ig[j], end = ig[j + 1] - 1;
      while (jg[beg] != i)
      {
         int ind = (beg + end) / 2;
         if (jg[ind] < i)
            beg = ind + 1;
         else
            end = ind;
      }
      ggu[beg] += elem;
   }
}

void calcLocalStiffnessMatrix(FiniteElement &elem, vector<vector<double>> &G, double tValue)
{
   auto &alpha = elem.alpha;

   double a1 = alpha[0 * 3 + 1] * alpha[0 * 3 + 1] + alpha[0 * 3 + 2] * alpha[0 * 3 + 2];
   double a2 = alpha[0 * 3 + 1] * alpha[1 * 3 + 1] + alpha[0 * 3 + 2] * alpha[1 * 3 + 2];
   double a3 = alpha[0 * 3 + 1] * alpha[2 * 3 + 1] + alpha[0 * 3 + 2] * alpha[2 * 3 + 2];
   double a4 = alpha[1 * 3 + 1] * alpha[1 * 3 + 1] + alpha[1 * 3 + 2] * alpha[1 * 3 + 2];
   double a5 = alpha[1 * 3 + 1] * alpha[2 * 3 + 1] + alpha[1 * 3 + 2] * alpha[2 * 3 + 2];
   double a6 = alpha[2 * 3 + 1] * alpha[2 * 3 + 1] + alpha[2 * 3 + 2] * alpha[2 * 3 + 2];

   for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
         G[i][j] = lambda(tValue) * abs(elem.detD) * (
            a1 * G1[i][j] +
            a2 * G2[i][j] +
            a3 * G3[i][j] +
            a4 * G4[i][j] +
            a5 * G5[i][j] +
            a6 * G6[i][j]);
}

void calcLocalMassMatrix(FiniteElement &elem, vector<vector<double>> &M)
{
   auto &localVertex = elem.localVertex;
   auto &sigmaWeights = elem.sigmaWeights;

   for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
         M[i][j] = abs(elem.detD) * 
                  (sigmaWeights[0] * M1[i][j] +
                   sigmaWeights[1] * M2[i][j] +
                   sigmaWeights[2] * M3[i][j]);
}

void calcLocalb(vector<Point> &coords, FiniteElement &elem, vector<double> &b, double tValue)
{
   auto &localVertex = elem.localVertex;

   for (int i = 0; i < 6; i++)
   {
      double sumbi = 0;

      for (int j = 0; j < 6; j++)
         sumbi += F(coords[localVertex[j]].x, coords[localVertex[j]].y, tValue) * Mass[i][j];

      b[i] = abs(elem.detD) * sumbi;
   }
}

void getWeightsInitU(TimeMesh &time, FEM &fem)
{
   auto &qti = time.qti;
   auto &qti_1 = time.qti_1;
   auto &qti_2 = time.qti_2;

   const int sizeMatrix = fem.vertexCoord.size();

   qti.resize(sizeMatrix);
   qti_1.resize(sizeMatrix);
   qti_2.resize(sizeMatrix);

   for (int i = 0; i < sizeMatrix; i++)
   {
      qti_2[i] = uInit(fem.vertexCoord[i].x, fem.vertexCoord[i].y, time.tSplitting[0]);
      qti_1[i] = uInit(fem.vertexCoord[i].x, fem.vertexCoord[i].y, time.tSplitting[1]);
   }
}