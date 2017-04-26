/* Solving Poisson equation -u''=b on 0,2*pi grid
with n points in 1d
h is width of grid
here u=sin */



#include<Eigen/Eigen/SparseCholesky>
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
const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "", "");

// define sparse matrix template
typedef SparseMatrix<double> SpMat;


void set_A(SpMat& A,double h) {
 // define sparse matrix template triplet
 typedef Eigen::Triplet<double> T;
 // h=1/n
 int n=A.rows();
 // list of non-zeros coefficients
 std::vector<T> aa;

 auto beginning = 0;
 auto ending = 2*M_PI;
 // Set grid from 0 to Pi
 VectorXd grid;
 grid.setLinSpaced(n,beginning,ending);

 // put stencil 1/h^2 [-1 2 -1]
 for(int i=0 ; i < n ; i++) {

      if(i==0)         {
        aa.push_back( T(i,i, 1 ));
                       }

      if(i>0 && i<n-1) {

        aa.push_back( T(i,i-1, -1/(h*h)));
        aa.push_back( T(i,i, 2/(h*h) ) );
        aa.push_back( T(i,i+1, -1/(h*h)));
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
 std::cout  << Af.format(CleanFmt) << std::endl; */

 return;
}
// Fill right hand side  b
void set_b(VectorXd& b, VectorXd& grid ,int n) {

 // Right hand side sin(x)
 for(int i=0 ; i<n ; i++) {

 double x=grid(i);
 b(i)=  std::sin(x);

 }

 return;
}


int main(void) {
 int n=100;
 VectorXd solution(n), b(n);
 SpMat A(n,n);

 // Set up grid
 // Set grid from 0 to 2*Pi

 auto beginning = 0;
 auto ending = 2*M_PI;

 VectorXd grid;
 grid.setLinSpaced(n,beginning,ending);
 // Length of step
 double h=grid(1);
 // fill A

 set_A(A,h);

 // fill b
 set_b(b,grid, n);


 /* solve Ax = b
   https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html
    Direct LDLT */
 Eigen::SimplicialLDLT<SpMat> solver;

  solver.compute(A);
  if(solver.info()!=Eigen::Success) {
   std::cout << "Failure decomposing\n";
  return 1;
}
  solution = solver.solve(b);

  if(solver.info()!=Eigen::Success) {
  std::cout << "Failure solving\n";
  return 1;
 }
 std::ofstream out("approx_sin.txt");
 std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
 std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
 std::cout  << solution.format(CleanFmt) << std::endl;

 return 0;
}
