#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include "GetPot"
#include "readParameters.hpp"
#include "Challenge1.hpp"
#include "muparser_fun_2D.hpp"
#include "types_def.hpp"


void
printHelp()
{
  std::cout
    << "USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"
    << std::endl;
  std::cout << "-h this help" << std::endl;
  std::cout << "-v verbose output" << std::endl;
}


int 
main(int argc, char **argv)
{   
    int  status{0}; 
    using namespace std;

    bool jsonfile=false;
    GetPot cl(argc, argv);
    if(cl.search(2, "-h", "--help"))
    {
      printHelp();
      return 0;
    }
    bool verbose = cl.search(1, "-v");

    // Get file with parameter values
    string filename = cl.follow("parameters.pot", "-p");
    auto pos = filename.find(".json");
    if(pos != std::string::npos)
    {
        jsonfile=true;
        std::cout<<"Json input file\n";
    }
    else
    {
        jsonfile=false;
        std::cout<<"Getpot input file\n";
    }
    cout << "Reading parameters from " << filename << std::endl;
    parameters param;
        if(jsonfile)
            param = readParameters_json(filename, verbose);
        else
            param = readParameters(filename, verbose);

    //parameters
    const auto &[f, df_x1, df_x2, x1_0, x2_0, alpha_0, eps_r, eps_s, mu, sigma, eta, itermax, rule, exact_grad, h] = param;

    //initial point
    vector<double> x0{x1_0,x2_0};
    
    //Muparser ===> wrapper: done in the wrapper for muparser muparser_fun_2D.hpp
    //f: is a function that takes two double and returns one double
    MuparserFun2 f_pars(f);
    function<double(double,double)> fun = f_pars;   //wrapper for the function
    //grad: is a vector of two function, each one of them takes two double and returns one double
    MuparserFun2 df_x1_pars(df_x1);
    function<double(double,double)> dfun_x1 = df_x1_pars;
    MuparserFun2 df_x2_pars(df_x2);
    function<double(double,double)> dfun_x2 = df_x2_pars;
    vector<function<double(double,double)>> dfun;   //wrapper for the gradient
    dfun.reserve(2);
    dfun.push_back(dfun_x1);
    dfun.push_back(dfun_x2);

    //storage of the approximation of the min
    std::vector<double> result;

    //the idea is to have two constant expression, one that indicates if we are using the exact gradient
    //or its approximation due to finite differnces, and one indicating which algorithm we are actually using.
    //We want know this two values once and for all since we need these informations at each iterations.
    // I could not being able to create the constant expression through a single constexpr function: the motivation
    // resides in the fact that I had to pass the parameteres read from the files. Explicitly, I had this error,
    // in passing in the two functions below rule and exact_grad:
    // error: the value of '<structured bindings>' was not initialized with a constant expression
    // I do not actually really understand the motivation behind this, since rule and exact_grad are actually constant
    //So I had to relay on this if sequence, that is not really elegant, in order to pass the int and the bool as actual constants.
    // In case other algo are coded, other than in the types_def.hpp files, they have to be added in this if-sequence,
    // updating the number from 5 on
    if ( exact_grad == true){

      constexpr bool exact_grad_chosen = ex_grad_chosen(true);
      cout << "Solving using exact gradient" << endl;

      if (rule==0)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(0);
        cout << "Solving using Armijo rule" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==1)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(1);
        cout << "Solving using exponential decay rule" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);


      } else if (rule==2)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(2);
        cout << "Solving using inverse decay rule" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==3)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(3);
        cout << "Solving using heavy ball method" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==4)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(4);
        cout << "Solving using Nesterov method" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else
      {
        cerr << "Wrong algorithm requested: expected an integer between 0 and 4" << endl;
        status = 1;
        return status;

      }
      
    } else
    {
      constexpr bool exact_grad_chosen = ex_grad_chosen(false);
      cout << "Solving using finite differnce method" << endl;
      

      if (rule==0)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(0);
        cout << "Solving using Armijo rule" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==1)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(1);
        cout << "Solving using exponential decay rule" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==2)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(2);
        cout << "Solving using inverse decay rule" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==3)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(3);
        cout << "Solving using heavy ball method" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else if (rule==4)
      {
        constexpr AlgoAvail algo_chosen = rule_chosen(4);
        cout << "Solving using Nesterov method" << endl;
        result = function_min(fun,dfun,x0,alpha_0,itermax,eps_s,eps_r,sigma,mu,eta,algo_chosen,exact_grad_chosen,h);

      } else
      {
        cerr << "Wrong algorithm requested: expected an integer between 0 and 4" << endl;
        status = 1;
        return status;

      }
    } 
 
   cout << "Approximation of the minimum found: (" << result[0] << "," << result[1] << ")" << endl;
        
    return status;

}