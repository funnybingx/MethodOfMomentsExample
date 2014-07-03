#ifndef UTILS_HELP_HPP
#define UTILS_HELP_HPP

#include <TString.h>

namespace utils
{
  namespace help
  {
    bool mkdir(TString dir);
    void reportOnLoop(unsigned int i, unsigned int t);
  }
}

#endif // UTILS_HELP_HPP
