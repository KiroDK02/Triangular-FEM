#include "structs.h"

void readTimeMesh(TimeMesh &time)
{
   int countTimeLayer = 0;

   std::ifstream timeMesh("timeMesh.txt");

   timeMesh >> countTimeLayer;
   time.t.resize(countTimeLayer);

   for (int ti = 0; ti < countTimeLayer; ti++)
   {
      double t = 0;

      timeMesh >> t;
      time.t[ti] = t;
   }

   timeMesh.close();
}

void readSplitTimeMesh(TimeMesh &time)
{
   auto &tSplitting = time.tSplitting;

   const int tSize = time.t.size();

   std::ifstream splittingTimeMesh("splittingTimeMesh.txt");

   int nk = 0;
   tSplitting.resize(1, time.t[0]);

   for (int ti = 0, j = 1; ti < tSize - 1; ti++, j++)
   {
      int countIntervals = 0;

      double coef = 0;
      double step = 0;

      splittingTimeMesh >> countIntervals >> coef;
      nk += countIntervals;
      tSplitting.resize(nk + 1);

      if (coef != 1)
      {
         double sumProgression = (pow(coef, countIntervals) - 1.) / (coef - 1.);
         step = (time.t[ti + 1] - time.t[ti]) / sumProgression;

         int jk = 1;
         for (j; j < nk; j++, jk++)
            tSplitting[j] = time.t[ti] + step * (pow(coef, jk) - 1.) / (coef - 1.);
      }
      else
      {
         step = (time.t[ti + 1] - time.t[ti]) / countIntervals;

         int jk = 1;
         for (j; j < nk; j++, jk++)
            tSplitting[j] = time.t[ti] + step * jk;
      }

      tSplitting[j] = time.t[ti + 1];
   }

   time.q.resize(time.tSplitting.size());
}

void readMesh(FEM &fem)
{
   auto &vertexCoord = fem.vertexCoord;
   auto &elements = fem.elements;

   int numVertices = 0, numElements = 0;

   std::ifstream mesh("mesh.txt");
   mesh >> numVertices;
   vertexCoord.resize(numVertices);

   for (int i = 0; i < numVertices; i++)
   {
      double x = 0, y = 0;
      mesh >> x >> y;
      vertexCoord[i].x = x;
      vertexCoord[i].y = y;
   }

   mesh.close();

   std::ifstream element("elements.txt");
   element >> numElements;
   elements.resize(numElements);

   for (int i = 0; i < numElements; i++)
   {
      elements[i].alpha.resize(9);
      elements[i].localVertex.resize(6);
      elements[i].sigmaWeights.resize(3);

      int p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 = 0;
      element >> p1 >> p2 >> p3 >> p4 >> p5 >> p6;
     
      elements[i].localVertex[0] = p1;
      elements[i].localVertex[1] = p2;
      elements[i].localVertex[2] = p3;
      elements[i].localVertex[3] = p4;
      elements[i].localVertex[4] = p5;
      elements[i].localVertex[5] = p6;

      findBaseNodes(elements[i], vertexCoord);
      calcDetD(vertexCoord, elements[i]);
      calcAlpha(vertexCoord, elements[i]);
   }

   element.close();
}

void readBoundaryCondition(vector<Conditions> &cond)
{
   int numEdgeConditions = 0;
   
   std::ifstream conditions("conditions.txt");
   conditions >> numEdgeConditions;
   cond.resize(numEdgeConditions);

   for (int i = 0; i < numEdgeConditions; i++)
   {
      auto &globalNumVertex = cond[i].globalNumVertex;
      int numVertex = 0, typeOfCond = 0, numOfFunction = 0;

      conditions >> typeOfCond >> numOfFunction >> numVertex;
      
      cond[i].typeOfCond = typeOfCond;
      cond[i].numOfFunction = numOfFunction - 1;
      globalNumVertex.resize(numVertex);

      for (int k = 0; k < numVertex; k++)
      {
         int vertex = 0;

         conditions >> vertex;
         globalNumVertex[k] = vertex;
      }
   }
}