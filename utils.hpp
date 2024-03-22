#ifndef HPP_UTILS_CHALLENGE_1_HPP
#define HPP_UTILS_CHALLENGE_1_HPP

#include<cmath>
#include<vector>

//some functions that are useful

//function to do the power with integer exponent
inline double pow_integer(double base, unsigned int exp){
  double res = 1.0;
  while (exp > 0) {
    if (exp & 1)
      res *= base;
    base *= base;
    exp >>= 1;
  }
  return res;
};

//function to evaluate the Euclidean norm of a vector
template<typename T>
double euclidean_norm(std::vector<T> const &vec){
    
    double norm = 0.;
    for (const auto i: vec){
        norm += pow_integer(i,2);
    }    
    return std::sqrt(norm);
};

//function to handle vector differences
template<typename T>
std::vector<T> vector_diff(std::vector<T> const &v1, std::vector<T> const &v2){

    const size_t n = v1.size();

    std::vector<T> differences;
    differences.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        differences.emplace_back(v1[i] - v2[i]);
    }

    return differences;
};


#endif