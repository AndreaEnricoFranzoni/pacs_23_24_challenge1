#include "parameters.hpp"
#include <iostream>
std::ostream &
operator<<(std::ostream &out, const parameters &p)
{
  out << "PARAMETER VALUES:"
      << "\n";
  out << "f= " << p.f << "\n";
  out << "df_x1= " << p.df_x1 << "\n";
  out << "df_x2= " << p.df_x2 << "\n";
  out << "x1_0= " << p.x1_0 << "\n";
  out << "x2_0= " << p.x2_0 << "\n";
  out << "alpha_0= " << p.alpha_0 << "\n";
  out << "eps_r= " << p.eps_r << "\n";
  out << "eps_s= " << p.eps_s << "\n";
  out << "mu= " << p.mu << "\n";
  out << "sigma= " << p.sigma << "\n";
  out << "eta= " << p.eta << "\n";
  out << "itermax= " << p.itermax << "\n";
  out << "rule= " << p.rule << "\n";
  out << "Exact gradient= " << p.exact_grad << "\n\n";
  out << "h= " << p.h << "\n\n";
  return out;
}