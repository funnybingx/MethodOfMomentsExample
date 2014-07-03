
#include "utils/help.hpp"

#include <iostream>
#include <TSystem.h>

using std::cout;

bool utils::help::mkdir(TString dir)
{
  if (dir.EndsWith("/")) dir.Chop();
  int returncode = gSystem->Exec(TString("mkdir -p ") + dir);
  if (returncode ==0) return true;
  else return false;
}

void utils::help::reportOnLoop(unsigned int i, unsigned int t)
{
  cout << " entry : " << i << " of : " << t
      << " ( " << 100.*(double(i)/double(t)) << " %) " << "\r" ;
  cout.flush() ;
  return;
}

