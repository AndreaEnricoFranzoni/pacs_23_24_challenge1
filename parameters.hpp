#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
#include <string>
#include "types_def.hpp"

//struct that defines the structure of the parameters

struct parameters
{ //analytical expression of the function R^2 --> R
  std::string f="x1*x2+4*x1*x1*x1*x1+x2*x2+3*x1";
  //analytical expression of the partial derivative wrt x1
  std::string df_x1="x2+16*x1*x1*x1+3";
  //analytical expression of the partial derivative wrt x2
  std::string df_x2="x1+2*x2";
  // first component of the initial point
  double x1_0 = 0.;
  // second component of the initial point
  double x2_0 = 0.;
  // initial value of the learning rate
  double alpha_0 = 0.5;
  // tolerance for residual 
  double eps_r = 1.e-6;
  // tolerance for step length
  double eps_s = 1.e-6;
  // value of the step paramter
  double mu = 0.2;
  // value of the Armijo's rule coefficient
  double sigma = 0.3;
  // value of the coefficient fot the momentum's method
  double eta = 0.9;
  // max number of iterations
  int itermax = 1000000;
  // rule for the gradient descent method
  int rule = 0;
  // if we want to use the exact formulartion of the gradient or an approximation of it
  bool exact_grad = true;
  // spacing for finite difference method for gradient approximation
  double h = 0.001;
};
//! Prints parameters
std::ostream &operator<<(std::ostream &, const parameters &);
#endif