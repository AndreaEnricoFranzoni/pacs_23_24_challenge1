#ifndef HPP_CHALLENGE_1_HPP
#define HPP_CHALLENGE_1_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <concepts>
#include "types_def.hpp"
#include "utils.hpp"


//CREATING CONSTEXPR for WHICH ALGO and IF EXACT GRAD
constexpr AlgoAvail rule_chosen(const int &rule){
    AlgoAvail algo_used;

    switch (rule)
    {
    case 0:
        algo_used = Armijo;
        return algo_used;
        break;

    case 1:
        algo_used = exp_decay;
        return algo_used;
        break;
    
    case 2:
        algo_used = inv_decay;
        return algo_used;
        break;
    
    case 3:
        algo_used = heavy_ball;
        return algo_used;
        break;

    case 4:
        algo_used = nesterov;
        return algo_used;
        break;

}
}

constexpr bool ex_grad_chosen(const bool &choice){
    if(choice){
        return true;
    } else{
        return false;
    }
}

//GRADIENT EVALUATION
//In case we'd like to evaluate the gradient using finite differences
//(taken from lecture7)
template <unsigned N>
inline auto numDeriv = [](auto const &f , double const &h){ 
    if constexpr (N==0u)
        return [f](auto x){ return f(x); }; 
    else
        {
            return [f, h](auto x) {
                auto const center = (N%2 == 0) ? x : x + h ; 
                auto prev=numDeriv<N-1u>(f, h);
                return (prev(center) - prev(center - h)) / h; 
                };
        } 
};

//to decide on the fly if evaluating with exact gradient or by finite differences
//The grad has been defined as a vector of two scalar functions with 2-dim input
template<class T>
std::vector<T> grad_evaluation(std::vector<T> &x_k,
                               const std::vector<std::function<T(T,T)>> &grad,
                               const std::function<T(T,T)> &f,
                               double const &h, 
                               const bool &exact){
                                
    std::vector<T> grad_eval;
    grad_eval.reserve(x_k.size());   //but only for 2-D input space functions
    if (exact == true)
    {   
        grad_eval.push_back(grad[0](x_k[0],x_k[1]));
        grad_eval.push_back(grad[1](x_k[0],x_k[1]));
        return grad_eval;

    } else
    {   
        
        //Evaluation of df_dx1
        auto f_givenx2 = [&f,&x_k](T x){return f(x,x_k[1]);};
        auto df_dx1 = numDeriv<1>(f_givenx2,h);
        grad_eval.push_back(df_dx1(x_k[0]));

        //Evaluation of df_dx1
        auto f_givenx1 = [&f,&x_k](T x){return f(x_k[0],x);};
        auto df_dx2 = numDeriv<1>(f_givenx1,h);
        grad_eval.push_back(df_dx2(x_k[1]));

        return grad_eval;  
    }
       
}


//CONVERGENCE CRITERIA
//check on step length: checking if the euclidean norm of two consecutive iterations is under a certain fixed tolerance
template<class T>
bool control_on_the_step_length(std::vector<T> const &x_k_1, std::vector<T>  const &x_k, double const &eps_s){

    std::vector<T> diff = vector_diff(x_k_1,x_k);

    if (euclidean_norm(diff)<eps_s)
    {
        return true;
    } else{
        return false;
    }
};

//check on residuals: checking if the euclidean norm of the gradient or the absolute value between two 
// consecutive iterations is under a certain tolerance
template<class T>
bool control_on_the_residual(T const &f_k_1, T const &f_k, std::vector<T> const &grad_k, double const &eps_r){

    const T norm_grad = euclidean_norm(grad_k);
    const T diff_update = std::abs(f_k_1-f_k);
    if ( norm_grad < eps_r || diff_update < eps_r)
    {
        return true;
    } else{
        return false;
    }
};


//RULES FOR UPDATING THE LEARNING RATE
//Armijo rule for updating the learning rate: everything is passed as parameter since is already calculated during the gradient descent algo
template<class T>
double Armijo_rule(const double &starting_alpha,
                   const std::vector<T> &x_k,
                   const double &sigma,
                   const std::function<T(T,T)> &f,
                   T const &f_ev,
                   const std::vector<T> &grad_ev,
                   T const &df_ev_norm_squares,
                   const std::integral auto &dim_input){
              
              //evaluating the improvement wrt the current point
              std::vector<double> temp_adjustment;
              temp_adjustment.reserve(dim_input);
              for (size_t j = 0; j < dim_input; ++j){   
                temp_adjustment.push_back(starting_alpha*grad_ev[j]);
              }
              
              
              //two iterations value difference
              std::vector<T> updated_x_k = vector_diff(x_k, temp_adjustment);
              //checking Armijo condition
              if (f_ev - f(updated_x_k[0],updated_x_k[1]) >= sigma*starting_alpha*df_ev_norm_squares)
              {
                return starting_alpha;
              } else
              { 
                return Armijo_rule(starting_alpha/2,x_k,sigma,f,f_ev,grad_ev,df_ev_norm_squares,dim_input);
              }
};

//Exponential decay rule for updating the learning rate
inline double exponential_decay(const double &alpha,const double &mu){ 
    return alpha*std::exp(-mu);
};

//Inverse decay rule for updating the learning rate
inline double inverse_decay(const double &alpha_0, const double &mu, const size_t &k){
    return alpha_0/(1+mu*k);
}


//LEARNING ALGO
//Gradient descent
template<class T>
std::vector<T> gradient_descent(
    const std::function<T(T,T)> &f,
    const std::vector<std::function<T(T,T)>> &df,                    
    const std::vector<T> &x0,
    const double &alpha0,
    const int &itermax,
    const double &eps_s,
    const double &eps_r,
    const double &sigma,
    const double &mu,
    const AlgoAvail &rule,
    const bool &grad_ex,
    const double &h
){

  //initialization of the algorithm
  auto n = x0.size();
  std::vector<T> x_k = x0;
  double alpha_k = alpha0;
  std::vector<std::vector<T>> iterations;
  iterations.push_back(x_k);  

  for (size_t i = 0; i < static_cast<decltype(i)>(itermax); ++i)
  { 
    //function evaluation in the current point
    T f_eval = f(x_k[0],x_k[1]);
    //gradient evaluation im the current point (exact or by finite differences)
    std::vector<T> df_eval= grad_evaluation(x_k,df,f,h,grad_ex);

    //evaluating the adjustment in the descent
    std::vector<T> adjustment;
    adjustment.reserve(n);
    for (size_t j = 0; j < n; ++j)
    {   
        adjustment.push_back(alpha_k*df_eval[j]);
    }
    
    //next canididate point and its evaluation
    std::vector<T> x_k_1 = vector_diff(x_k,adjustment);
    T f_eval_next = f(x_k_1[0],x_k_1[1]);
    iterations.push_back(x_k_1);
    
    //check if convergence reached
    if (control_on_the_step_length(x_k_1,x_k,eps_s) ||
        control_on_the_residual(f_eval_next, f_eval, df_eval, eps_r))
    {   
        //if convergence reached: returns the approximated minimum
        std::cout << "Convergence reached on the " << i+1 << "-th iteration" << std::endl;
        return x_k_1;
    } else  
    {
        //if convergence not reached: update of the learning rate
        switch (rule)
        {
        case Armijo: //Armijo
            alpha_k = Armijo_rule(alpha_k,x_k,sigma,f,f_eval,df_eval,pow_integer(euclidean_norm(df_eval),2),n);
            break;

        case exp_decay: //Exponential decay
            alpha_k = exponential_decay(alpha_k,mu);
            break;

        case inv_decay: //Inverse decay
            alpha_k = inverse_decay(alpha0,mu,i);
            break;
        }
    }

    //swapping the point for the next iteration
    x_k.swap(x_k_1);
  }
    
  //If convergence not reached: error is printed and null vector is returned
  std::cerr << "No convergence" << std::endl;
  return std::vector<T>(n,0.);  
};


template<class T>
std::vector<T> HEAVYBALL(
    const std::function<T(T,T)> &f,
    const std::vector<std::function<T(T,T)>> &df,
    const std::vector<T> &x0,
    const double &alpha,
    const int &itermax,
    const double &eps_s,
    const double &eps_r,
    const double &eta,
    const bool &grad_ex,
    const double &h
){
    //initialization of the algorithm
    auto n = x0.size();
    std::vector<T> x_k = x0;
    std::vector<T> df_0 = grad_evaluation(x_k,df,f,h,grad_ex);

    std::vector<T> descent;
    descent.reserve(n);
    for (size_t j = 0; j < n; ++j)
    {   
        descent.push_back(alpha*df_0[j]);
    }

    x_k = vector_diff(x0,descent);

    std::vector<std::vector<T>> iterations;
    iterations.push_back(x0);
    iterations.push_back(x_k);

    
    //the first iteration has been already done
    for (size_t i = 1; i < static_cast<decltype(i)>(itermax); ++i)
    { 
        //function evaluation in the current point
        T f_eval = f(x_k[0],x_k[1]);
        //gradient evaluation im the current point (exact or by finite differences)
        std::vector<T> df_eval = grad_evaluation(x_k,df,f,h,grad_ex);
    
        //evaluating the adjustment: descent and momentum factors
        //descent
        std::vector<T> descent;
        descent.reserve(n);
        for (size_t j = 0; j < n; ++j)
        {   
            descent.push_back(alpha*df_eval[j]);
        }
        //momentum
        std::vector<T> momentum;
        momentum.reserve(n);
        for (size_t j = 0; j < n; ++j)
        {   
            momentum.push_back(eta*vector_diff(iterations[i-1],iterations[i])[j]);
            //it is actually the opposite of the momentum
        }
        //adjustment
        std::vector<T> adjustment = vector_diff(descent,momentum);

        //next canididate point and its evaluation
        std::vector<T> x_k_1 = vector_diff(x_k,adjustment);
        T f_eval_next = f(x_k_1[0],x_k_1[1]);

        iterations.push_back(x_k_1);

    
        //check if convergence reached
        if (control_on_the_step_length(x_k_1,x_k,eps_s) ||
            control_on_the_residual(f_eval_next, f_eval, df_eval, eps_r))
        {   
            //if convergence reached: returns the approximated minimum
            std::cout << "Convergence reached on the " << i+1 << "-th iteration" << std::endl;
            return x_k_1;
        } 
        //if convergence not reached: next iteration
        //swapping the point for the next iteration
        x_k.swap(x_k_1);
    }

    //If convergence not reached: error is printed and null vector is returned
    std::cerr << "No convergence" << std::endl;
    return std::vector<T>(n,0.);
};


template<class T>
std::vector<T> NESTEROV(
    const std::function<T(T,T)> &f,
    const std::vector<std::function<T(T,T)>> &df,
    const std::vector<T> &x0,
    const double &alpha,
    const int &itermax,
    const double &eps_s,
    const double &eps_r,
    const double &eta,
    const bool &grad_ex,
    const double &h
){
    //initialization of the algorithm
    auto n = x0.size();
    std::vector<T> x_k = x0;
    std::vector<T> df_0 = grad_evaluation(x_k,df,f,h,grad_ex);

    std::vector<T> descent;
    descent.reserve(n);
    for (size_t j = 0; j < n; ++j)
    {   
        descent.push_back(alpha*df_0[j]);
    }

    x_k = vector_diff(x0,descent);

    std::vector<std::vector<T>> iterations;
    iterations.push_back(x0);
    iterations.push_back(x_k);

    
    //the first iteration has been already done
    for (size_t i = 1; i < static_cast<decltype(i)>(itermax); ++i)
    { 
        //function evaluation in the current point
        T f_eval = f(x_k[0],x_k[1]);
        //gradient evaluation im the current point (exact or by finite differences)
        std::vector<T> df_eval = grad_evaluation(x_k,df,f,h,grad_ex);
    

        //evaluating the adjustment: descent and momentum factors 
        //momentum
        std::vector<T> momentum;
        momentum.reserve(n);
        for (size_t j = 0; j < n; ++j)
        {   
            momentum.push_back(eta*vector_diff(iterations[i-1],iterations[i])[j]);
            //it is actually the opposite of the momentum
        }
        std::vector<T> y = vector_diff(x_k,momentum);
        //descent
        std::vector<T> descent;
        descent.reserve(n);
        std::vector<T> grad_f_y = grad_evaluation(y,df,f,h,grad_ex);
        for (size_t j = 0; j < n; ++j)
        {   
            descent.push_back(alpha*grad_f_y[j]);
        }

        //next canididate point and its evaluation
        std::vector<T> x_k_1 = vector_diff(y,descent);
        T f_eval_next = f(x_k_1[0],x_k_1[1]);

        iterations.push_back(x_k_1);

        //check if convergence reached
        if (control_on_the_step_length(x_k_1,x_k,eps_s) ||
            control_on_the_residual(f_eval_next, f_eval, df_eval, eps_r))
        {   
            //if convergence reached: returns the approximated minimum
            std::cout << "Convergence reached on the " << i+1 << "-th iteration" << std::endl;
            return x_k_1;
        } 
        //if convergence not reached: next iteration
        //swapping the point for the next iteration
        x_k.swap(x_k_1);
    }

    //If convergence not reached: error is printed and null vector is returned
    std::cerr << "No convergence" << std::endl;
    return std::vector<T>(n,0.);
};


//This is the function called actually, that simply calls the correct 
//method to be used
template<class T>
std::vector<T> function_min(
        const std::function<T(T,T)> &f,
        const std::vector<std::function<T(T,T)>> &df,                    
        const std::vector<T> &x0,
        const double &alpha0,
        const int &itermax,
        const double &eps_s,
        const double &eps_r,
        const double &sigma,
        const double &mu,
        const double &eta,
        const AlgoAvail &rule,
        const bool &grad_ex,
        const double &h
){

    
    switch (rule)
    {
    case Armijo:
        return gradient_descent(f,df,x0,alpha0,itermax,eps_s,eps_r,sigma,mu,rule,grad_ex,h);
        break;
    
    case exp_decay:
        return gradient_descent(f,df,x0,alpha0,itermax,eps_s,eps_r,sigma,mu,rule,grad_ex,h);
        break;

    case inv_decay:
        return gradient_descent(f,df,x0,alpha0,itermax,eps_s,eps_r,sigma,mu,rule,grad_ex,h);
        break;

    case heavy_ball:
        return HEAVYBALL(f,df,x0,alpha0,itermax,eps_s,eps_r,eta,grad_ex,h);
        break;
    
    case nesterov:
        return NESTEROV(f,df,x0,alpha0,itermax,eps_s,eps_r,eta,grad_ex,h);
        break;


    }

}
    

#endif // HPP_CHALLENGE_1_HPP