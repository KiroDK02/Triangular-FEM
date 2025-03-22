#include "structs.h"

double uNumerical(FEM &fem, double x, double y, vector<double> &q)
{
   auto &elements = fem.elements;
   auto &vertexCoord = fem.vertexCoord;
   
   bool isFound = false;

   FiniteElement element{ };

   for (auto &elem : elements)
   {
      auto &point1 = vertexCoord[elem.localVertex[0]];
      auto &point2 = vertexCoord[elem.localVertex[1]];
      auto &point3 = vertexCoord[elem.localVertex[2]];

      double realTriangleS = abs(elem.detD);

      double xyTriangleS = 2 *
         (triangleS(x, y, point2.x, point2.y, point3.x, point3.y) +
          triangleS(point1.x, point1.y, x, y, point3.x, point3.y) +
          triangleS(point1.x, point1.y, point2.x, point2.y, x, y));

      if (abs(realTriangleS - xyTriangleS) < 1e-10)
      {
         isFound = true;
         element = elem;
         break;
      }
   }

   if (!isFound)
   {
      std::cout << "Element was not found.\n";
      return -10000000;
   }

   double uNum = q[element.localVertex[0]] * element.L1(x, y) * (2. * element.L1(x, y) - 1.) +
                 q[element.localVertex[1]] * element.L2(x, y) * (2. * element.L2(x, y) - 1.) + 
                 q[element.localVertex[2]] * element.L3(x, y) * (2. * element.L3(x, y) - 1.) + 
                 q[element.localVertex[3]] * 4. * element.L1(x, y) * element.L2(x, y) +
                 q[element.localVertex[4]] * 4. * element.L2(x, y) * element.L3(x, y) +
                 q[element.localVertex[5]] * 4. * element.L1(x, y) * element.L3(x, y);

   return uNum;
}

double uNumericalAnyTime(FEM &fem, TimeMesh &time, double x, double y, double tValue)
{
   auto &timeSplitting = time.tSplitting;
   auto &q = time.q;

   int ti = 0;

   for (int i = 0; i < timeSplitting.size(); i++)
      if (tValue <= timeSplitting[i])
      {
         ti = i;
         break;
      }

   auto &qti_2 = q[ti - 2];
   auto &qti_1 = q[ti - 1];
   auto &qti = q[ti];

   double uValue = uNumerical(fem, x, y, qti_2) * time.etta2(ti, tValue) +
                   uNumerical(fem, x, y, qti_1) * time.etta1(ti, tValue) +
                   uNumerical(fem, x, y, qti)   * time.etta0(ti, tValue);

   return uValue;
}

void initSizeLocalMatrix(vector<vector<double>> &matrix)
{
   matrix.resize(6);
   for (int i = 0; i < 6; i++)
      matrix[i].resize(6);
}

void stringMatrixInNull(SparseMatrix &A, int i)
{
   auto &ig = A.ig, &jg = A.jg;
   auto &ggl = A.ggl, &ggu = A.ggu, &di = A.di;
   const int size = ig[di.size()];

   int i0 = ig[i], i1 = ig[i + 1];

   for (i0; i0 < i1; i0++)
      ggl[i0] = 0.;

   int j0 = i1;

   for (j0; j0 < size; j0++)
      if (jg[j0] == i) ggu[j0] = 0.;

   di[i] = 1.;
}

void multiplyVectorCoef(vector<double> &vector, double coef)
{
   const int size = vector.size();

   for (int i = 0; i < size; i++)
      vector[i] *= coef;
}

void findBaseNodes(FiniteElement &elem, vector<Point> &vertexCoord)
{
   auto &localVertex = elem.localVertex;
   double maxS = 0;
   vector<int> indexes(3);

   for (int a1 = 0; a1 < 6; a1++)
      for (int a2 = a1 + 1; a2 < 6; a2++)
         for (int a3 = a2 + 1; a3 < 6; a3++)
         {
            double curS = triangleS(vertexCoord[localVertex[a1]].x, vertexCoord[localVertex[a1]].y,
                                    vertexCoord[localVertex[a2]].x, vertexCoord[localVertex[a2]].y,
                                    vertexCoord[localVertex[a3]].x, vertexCoord[localVertex[a3]].y);

            if (curS > maxS)
            {
               maxS = curS;
               indexes[0] = a1;
               indexes[1] = a2;
               indexes[2] = a3;
            }
         }

   std::swap(localVertex[0], localVertex[indexes[0]]);
   std::swap(localVertex[1], localVertex[indexes[1]]);
   std::swap(localVertex[2], localVertex[indexes[2]]);

   for (int idx = 3; idx < 6; idx++)
   {
      double xm3 = (vertexCoord[localVertex[0]].x + vertexCoord[localVertex[1]].x) / 2.;
      double ym3 = (vertexCoord[localVertex[0]].y + vertexCoord[localVertex[1]].y) / 2.;
      if (abs(vertexCoord[localVertex[idx]].x - xm3) < 1e-12 &&
         abs(vertexCoord[localVertex[idx]].y - ym3) < 1e-12)
      {
         std::swap(localVertex[3], localVertex[idx]);
         break;
      }
   }

   double xm4 = (vertexCoord[localVertex[1]].x + vertexCoord[localVertex[2]].x) / 2.;
   double ym4 = (vertexCoord[localVertex[1]].y + vertexCoord[localVertex[2]].y) / 2.;
   if (abs(vertexCoord[localVertex[5]].x - xm4) < 1e-12 &&
      abs(vertexCoord[localVertex[5]].y - ym4) < 1e-12)
      std::swap(localVertex[4], localVertex[5]);
}

void multiplyMatrixToCoef(vector<vector<double>> &matrix, double coef, vector<vector<double>> &resultMatrix)
{
   const int sizeMatrix = matrix.size();

   for (int i = 0; i < sizeMatrix; i++)
      for (int j = 0; j < sizeMatrix; j++)
         resultMatrix[i][j] = coef * matrix[i][j];
}

void multiplyMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec, vector<double> &result, vector<int> &localNum)
{
   const int sizeMatrix = matrix.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      double sum = 0;
      
      for (int j = 0; j < sizeMatrix; j++)
         sum += matrix[i][j] * vec[localNum[j]];

      result[i] = sum;
   }
}

void clearSLAE(SLAE &slae)
{
   auto &A = slae.A;
   auto &b = slae.b;
   auto &q = slae.q;

   const int sizeSLAE = A.di.size();
   const int countNonZeroElems = A.ig[sizeSLAE];

   for (int i = 0; i < sizeSLAE; i++)
   {
      A.di[i] = 0;
      q[i] = 0;
      b[i] = 0;
   }

   for (int i = 0; i < countNonZeroElems; i++)
   {
      A.ggl[i] = 0;
      A.ggu[i] = 0;
   }
}

void TimeMesh::swapVectors()
{
   std::swap(qti_2, qti_1);
   std::swap(qti_1, qti);
}

void TimeMesh::saveWeightsTi(int ti)
{
   q[ti] = qti;
}

double getLengthOfEdge(double x1, double x2, double y1, double y2)
{
   return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double triangleS(double x1, double y1, double x2, double y2, double x3, double y3)
{
   return 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
}

double FiniteElement::L1(double x, double y)
{
   return alpha[0] + alpha[1] * x + alpha[2] * y;
}

double FiniteElement::L2(double x, double y)
{
   return alpha[3] + alpha[4] * x + alpha[5] * y;
}

double FiniteElement::L3(double x, double y)
{
   return alpha[6] + alpha[7] * x + alpha[8] * y;
}

double TimeMesh::etta2(int ti, double t)
{
   double deltaT = tSplitting[ti] - tSplitting[ti - 2];
   double deltaT0 = tSplitting[ti] - tSplitting[ti - 1];
   double deltaT1 = tSplitting[ti - 1] - tSplitting[ti - 2];
   
   return (t - tSplitting[ti - 1]) * (t - tSplitting[ti]) / (deltaT * deltaT1);
}

double TimeMesh::etta1(int ti, double t)
{
   double deltaT = tSplitting[ti] - tSplitting[ti - 2];
   double deltaT0 = tSplitting[ti] - tSplitting[ti - 1];
   double deltaT1 = tSplitting[ti - 1] - tSplitting[ti - 2];

   return -(t - tSplitting[ti - 2]) * (t - tSplitting[ti]) / (deltaT1 * deltaT0);
}

double TimeMesh::etta0(int ti, double t)
{
   double deltaT = tSplitting[ti] - tSplitting[ti - 2];
   double deltaT0 = tSplitting[ti] - tSplitting[ti - 1];
   double deltaT1 = tSplitting[ti - 1] - tSplitting[ti - 2];

   return (t - tSplitting[ti - 2]) * (t - tSplitting[ti - 1]) / (deltaT * deltaT0);
}