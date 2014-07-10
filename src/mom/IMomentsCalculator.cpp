#include "mom/IMomentsCalculator.hpp"

IMomentsCalculator::IMomentsCalculator(unsigned int order, RooRealVar* x, 
                                       RooRealVar* w) :
    _order(order),
    _xvar(x),
    _wvar(w),
    _debug(false)
{
}
