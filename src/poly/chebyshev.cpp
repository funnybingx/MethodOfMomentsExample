
#include "poly/chebyshev.hpp"
#include <assert.h>
#include <TMath.h>
using std::vector;

//calculates all chebyshev polynomials up to including degree n
void poly::chebyshev(double x, int n, vector<double>& results)
{
  assert(n >= 0);
  results.clear();
  results.reserve(n+1);
  results.push_back(1.0);
  results.push_back(x);
  for (int i=2; i<n+1; i++)
    results.push_back(2.0*x*results.at(i-1) - results.at(i-2));
}

//get polynomial coefficients from chebyshev coefficients
void poly::chebyshev_to_poly(const vector<double>& chebyshev, vector<double>& poly)
{
  poly.clear();
  poly.resize(chebyshev.size(), 0.0);
  if (chebyshev.size() > 0)
    poly.at(0) = chebyshev.at(0);
  if (chebyshev.size() > 1)
    poly.at(1) = chebyshev.at(1);
  for (unsigned int n =2; n<chebyshev.size(); n++)
    {
      for (unsigned int r=0; r<=n/2; r++)
	{
	  poly.at(n-2*r) += chebyshev.at(n)
	    *n/2.0*pow(-1.0,r)/double(n-r)*TMath::Binomial(n-r,r)*pow(2.0,n-2*r);
	}
    }
}


//correct the polynomial coefficients from -1..1 to min..max
void poly::correct_poly(const vector<double>& poly, vector<double>& correct, double min, double max)
{
  correct.clear();
  
  double c = 2.0/(max-min);
  double d = -2.0*min/(max-min)-1.0;
  double e = d/c;
  //unsigned int the_bin, new_bin;

  unsigned int order = poly.size();
  correct.resize(order, 0.0);
  for (unsigned int i = 0; i<order; i++)
    {
      double c_dash = poly.at(i);
      for (unsigned int m=0; m<=i; m++)
	{
	  //new_bin = m;
	  correct.at(m) += 
	    c_dash
	    *pow(c,int(i))
	    *TMath::Binomial(i,m)
	    *pow(e,int(i-m));
	}
    }
}

