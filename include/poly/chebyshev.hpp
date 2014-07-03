
#ifndef POLY_CHEBYSHEV_HPP
#define POLY_CHEBYSHEV_HPP

#include <vector>

namespace poly
{
void chebyshev(double x, int n, std::vector<double>& results);

void chebyshev_to_poly(const std::vector<double>& chebychev, std::vector<double>& poly);
void correct_poly(const std::vector<double>& poly, std::vector<double>& correct, double min, double max);

};

#endif // POLY_CHEBYSHEV_HPP
