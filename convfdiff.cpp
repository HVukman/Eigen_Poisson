/* Solving CD equation -d*u''+c*u'=0 on 0,2*pi grid
with n points in 1d
f(0)=0 f(1)=1
exact solution (exp(Rl/x)-1)/(exp(Rl)-1)
RL= cl/d
Dirichlet Boundary left and right */


#include<Eigen/Eigen/SparseCholesky>
#include <Eigen/Eigen/IterativeLinearSolvers>
#include <iostream>
#include <vector>
#include <fstream>


#define _USE_MATH_DEFINES
#include <cmath>

// declaring Pi (maybe not needed)
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::SparseMatrix;


using std::cout;
using std::endl;

// Output format clean with precision of 4
const Eigen::IOFormat CleanFmt(9, 0, ", ", "\n", "", "");

// define sparse matrix template
typedef SparseMatrix<double> SpMat;


void set_A(SpMat& A,double h,double diffusion, double velocity) {
 // define sparse matrix template triplet
 typedef Eigen::Triplet<double> T;
 // h=1/n
 int n=A.rows();
 // list of non-zeros coefficients
 std::vector<T> aa;



 // put stencil 1/h^2 [1 -2 1]
 for(int i=0 ; i < n ; i++) {

      if(i==0)         {
        aa.push_back( T(i,i, 1 ));
                       }

      if(i>0 && i<n-1) {

        aa.push_back( T(i,i-1,  -diffusion/(h*h)));
        aa.push_back( T(i,i,   +diffusion*2/(h*h) -velocity/h) );
        aa.push_back( T(i,i+1,  -diffusion/(h*h) + velocity/h));
                       }

     if (i==n-1)       {
        aa.push_back( T(i,i, 1));
                       }

                            }
 // put triplets of a into A
 A.setFromTriplets(aa.begin(), aa.end());


 /* For checking a row
 VectorXd Af;
 Af = VectorXd(A.row(3));
 std::cout  << Af.format(CleanFmt) << std::endl;*/

 return;
}
// Fill matrix for convdiff

// Fill right hand side  b
void set_f(VectorXd& f ,int n)
 {

 // Right hand side
 for(int i=0 ; i<n-1 ; i++) {

 f(i)=  0;

 }

 f(n-1)=1;

 return;
}


int main(void) {
  int n=100;

 // Geschw und Diffusion
  double diffusion = 0.05;
  double velocity  = 1;

  VectorXd  f(n);


 // Set up grid
 // Set grid from 0 to 2*Pi

  auto beginning = 0;
  auto ending = 1;

  VectorXd grid(n),Solution(n);
  grid.setLinSpaced(n,beginning,ending);
 // Length of step
  double h=grid(1);
 // set A
  SpMat A(n,n);
  set_A(A,h,velocity,diffusion);

 // fill f
  set_f(f,n);

 /* solve Ax = b
   https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html
    Direct LDLT (doesn't really work)
  Eigen::SimplicialLDLT<SpMat> solver;
  solving iteratively with CG (also doen't work,maybe with good guess)*/
  VectorXd x0(n);
  x0(n-1)=1;
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setMaxIterations(10000);
  solver.setTolerance	( 10^-23);
  solver.compute(A);

  if(solver.info()!=Eigen::Success) {
   std::cout << "Failure decomposing\n";
  return 1;
}

  Solution  = solver.solveWithGuess(f,x0);
  if(solver.info()!=Eigen::Success) {
  std::cout << "Failure solving\n";
  return 1;
 }
  std::cout << "#iteration: " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error() << std::endl;
  std::cout  << Solution.format(CleanFmt) << std::endl;

  std::ofstream out("approx_cd.txt");
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
  std::cout  << Solution.format(CleanFmt) << std::endl;

 return 0;
}
