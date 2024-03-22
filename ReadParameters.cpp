#include "ReadParameters.hpp"
#include "GetPot"
#include "json.hpp"
#include <fstream>

//challenge1
parameters
readParameters(std::string const &filename, bool verbose)
{
  #include "GetPot"
  // Parameter default constructor fills it with the defaults values
  parameters defaults;
  // checks if file exixts and is readable
  std::ifstream check(filename);
  if(!check)
    {
      std::cerr << "ERROR: Parameter file " << filename << " does not exist"
                << std::endl;
      std::cerr << "Reverting to default values." << std::endl;
      if(verbose)
        std::cout << defaults;
      check.close();
      return defaults;
    }
  else
    check.close();

  GetPot     ifile(filename.c_str());
  parameters values;
  // Read parameters from getpot ddata base
  values.f = ifile("f", defaults.f);
  values.df_x1 = ifile("df_x1", defaults.df_x1);
  values.df_x2 = ifile("df_x2", defaults.df_x2);
  values.x1_0 = ifile("x1_0", defaults.x1_0);
  values.x2_0 = ifile("x2_0", defaults.x2_0);
  values.alpha_0 = ifile("alpha_0", defaults.alpha_0);
  values.eps_r = ifile("eps_r", defaults.eps_r);
  values.eps_s = ifile("eps_s", defaults.eps_r);
  values.mu = ifile("mu", defaults.mu);
  values.sigma = ifile("sigma", defaults.sigma);
  values.eta = ifile("eta", defaults.eta);
  values.itermax = ifile("itermax", defaults.itermax);
  values.rule = ifile("rule", defaults.rule);
  values.exact_grad = ifile("exact_grad", defaults.exact_grad);
  values.h = ifile("h", defaults.h);
  
  if(verbose)
    {
      std::cout << "PARAMETER VALUES IN GETPOT FILE"
                << "\n";
      ifile.print();
      std::cout << std::endl;
      std::cout << "ACTUAL VALUES"
                << "\n"
                << values;
    }
  return values;
  

}

parameters
readParameters_json(std::string const &filename, bool verbose)
{
  // Parameter default constructor fills it with the defaults values
  parameters defaults;
  // checks if file exixts and is readable
  std::ifstream check(filename);
  if(!check)
    {
      std::cerr << "ERROR: Parameter file " << filename << " does not exist"
                << std::endl;
      std::cerr << "Reverting to default values." << std::endl;
      if(verbose)
        std::cout << defaults;
      check.close();
      return defaults;
    }
  else
    check.close();

  std::ifstream jfile(filename);
  nlohmann::json ifile;
  jfile>>ifile;
  parameters values;
  // Read parameters from getpot ddata base
  values.f = ifile.value("f", defaults.f);
  values.df_x1 = ifile.value("df_x1", defaults.df_x1);
  values.df_x2 = ifile.value("df_x2", defaults.df_x2);
  values.x1_0 = ifile.value("x1_0", defaults.x1_0);
  values.x2_0 = ifile.value("x2_0", defaults.x2_0);
  values.alpha_0 = ifile.value("alpha_0", defaults.alpha_0);
  values.eps_r = ifile.value("eps_r", defaults.eps_r);
  values.eps_s = ifile.value("eps_s", defaults.eps_r);
  values.mu = ifile.value("mu", defaults.mu);
  values.sigma = ifile.value("sigma", defaults.sigma);
  values.eta = ifile.value("eta", defaults.eta);
  values.itermax = ifile.value("itermax", defaults.itermax);
  values.rule = ifile.value("rule", defaults.rule);
  values.exact_grad = ifile.value("exact_gradient", defaults.exact_grad);
  values.h = ifile.value("h", defaults.h);

  if(verbose)
    {
      std::cout << "PARAMETER VALUES IN JSON FILE"
                << "\n";
      std::cout<<std::setw(4)<<ifile;
      std::cout << std::endl;
      std::cout << "ACTUAL VALUES"
                << "\n"
                << values;
    }
  return values;
}