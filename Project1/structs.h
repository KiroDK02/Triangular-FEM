#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <functional>

using std::vector;

struct Point
{
   double x{ }, y{ };
};

struct FiniteElement
{
   vector<int> localVertex{ }; // Глобальные номера из массива
   double detD{ };             // vertexCoord локальных узлов элемента,
   vector<double> alpha{ };    
   vector<double> sigmaWeights{ };   

   double L1(double x, double y);
   double L2(double x, double y);
   double L3(double x, double y);
};                             

// локальные номера узлов элемента в массиве localVertex 
   /*     5     3-4-5
   *     /|     |  /
   *    3 4 или 1 2
   *   /  |     |/
   *  0-1-2     0
   */
// т.е. localVertex[0] - глобальный номер узла 0 из vertexCoord, 
// localVertex[1] - узла 1 из vertexCoord и т.д.

struct FEM
{
   vector<Point> vertexCoord{ };
   vector<FiniteElement> elements{ };
};

struct TimeMesh
{
   vector<double> t{ };
   vector<double> tSplitting{ };

   vector<vector<double>> q{ };
   vector<double> qti{ };
   vector<double> qti_1{ };
   vector<double> qti_2{ };

   void swapVectors();
   void saveWeightsTi(int ti);

   double etta2(int ti, double t);
   double etta1(int ti, double t);
   double etta0(int ti, double t);
};

// globalNumVertex:
// Для 1го краевого:
// Глобальные координаты всех узлов ребра,
// на котором задано краевое.
// Для 2го краевого:
// Задается каждое ребро, входящее в границу.
// Для одного ребра будут 3 глобальные номера
struct Conditions
{
   int typeOfCond{ }, numOfFunction{ };
   vector<int> globalNumVertex{ };
};

struct SparseMatrix
{
   vector<int> ig{ }, jg{ };
   vector<double> ggl{ }, ggu{ }, di{ };
};

struct SLAE
{
   SparseMatrix A{ };
   vector<double> q{ }, b{ };
};

struct LOS
{
   vector <double> r1{ }, rk{ }, z1{ }, p1{ }, Ar{ }, p{ }, mult{ };
};

// Для расчета на элементе
void calcDetD(vector<Point> &coords, FiniteElement &elem);
void calcAlpha(vector<Point> &coords, FiniteElement &elem);
void calcSigmaWeights(vector<Point> &coords, FiniteElement &elem, double tValue);

// Считывание сетки и ее разбиения по времени
void readTimeMesh(TimeMesh &time);
void readSplitTimeMesh(TimeMesh &time);

// Считываение сетки
void readMesh(FEM &fem);
void readBoundaryCondition(vector<Conditions> &cond);

// Для получения решения по времени
void getSolutionParabolicProblem(FEM &fem, TimeMesh &time, SLAE &slae, vector<Conditions> &conds);
void getSolutionHyperbolicProblem(FEM &fem, TimeMesh &time, SLAE &slae, vector<Conditions> &conds);

// Для глобальной матрицы
void portraitSparseMatrix(FEM &fem, SLAE &slae);
void calcGlobalMatrixAndVectorParabolicProblem(FEM &fem, TimeMesh &time, SLAE &slae, vector<Conditions> &cond, int ti);
void calcGlobalMatrixAndVectorHyperbolicProblem(FEM &fem, TimeMesh &time, SLAE &slae, vector<Conditions> &cond, int ti);
void addCondFirst(SLAE &slae, vector<Conditions> &cond, vector<Point> &coord, double tValue);
void addCondSecond(SLAE &slae, vector<Conditions> &conds, vector<Point> &coords, double tValue);
void calcLocalStiffnessMatrix(FiniteElement &elem, vector<vector<double>> &G, double tValue);
void calcLocalMassSigmaMatrix(FiniteElement &elem, vector<vector<double>> &M);
void calcLocalMassChiMatrix(FiniteElement &elem, vector<vector<double>> &M, double tValue);
void calcLocalb(vector<Point> &coords, FiniteElement &elem, vector<double> &b, double tValue);
void addLocalMatrixInGlobal(SparseMatrix &A, vector<int> &localVertex, vector<vector<double>> &localMatrix);
void addElemInGlobalMatrix(SparseMatrix &A, int i, int j, double elem);
void addLocalVectorInGlobal(vector<double> &b, vector<int> &localVertex, vector<double> &bLocal);
void getWeightsInitU(TimeMesh &time, FEM &fem);

// Рассчет численного значения u в рандомной точки области
double uNumerical(FEM &fem, double x, double y, vector<double> &q);
double uNumericalAnyTime(FEM &fem, TimeMesh &time, double x, double y, double tValue);
double calcErrorNumSolution(FEM &fem, TimeMesh &time, vector<double> &xSet, vector<double> &ySet, vector<double> &timeSet);

// вывод результата в наборе точек
void outputResult(FEM &fem, TimeMesh &time, vector<double> &xSet, vector<double> &ySet, int ti);

// Вспомогательные функции
void initSizeLocalMatrix(vector<vector<double>> &matrix);
void stringMatrixInNull(SparseMatrix &A, int i);
void multiplyVectorCoef(vector<double> &vector, double coef);
void findBaseNodes(FiniteElement &elem, vector<Point> &vertexCoord);
double triangleS(double x1, double x2, double x3, double y1, double y2, double y3);
double getLengthOfEdge(double x1, double x2, double y1, double y2);
void multiplyMatrixToCoef(vector<vector<double>> &matrix, double coef, vector<vector<double>> &resultMatrix);
void multiplyMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec, vector<double> &result, vector<int> &localNum);
void clearSLAE(SLAE &slae);

// Функции для ЛОСа(LU)
void localOptimalSchemeLU(SLAE &slae, SLAE &LU, LOS &v, int maxIter, double eps);
void calcLU(SLAE &slae, SLAE &LU);
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);
void multOfMatrix(SparseMatrix &A, vector<double> &x, vector<double> &F);
void calcDiscrepancy(SLAE &slae, LOS &v, vector<double> &x, double &normb);
void calcVectorMultCoef(vector <double> &a, double coef, vector <double> &res);
void calcSumVectors(vector <double> &a, vector <double> &b, vector <double> &res);
double scalarMult(vector <double> &a, vector <double> &b);