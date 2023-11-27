//---------------------------------------------------------------------------
#include "MainUnit.h"
//---------------------------------------------------------------------------
void usage(std::string fname)
{
  std::cout << "usage : thermo-hydro-mechanical analysis program" << std::endl;
  std::cout << "      : " << fname << " \"data file\"" << std::endl;
}

int main(int argc, char *argv[])
{
  TBEMInput FBEMInput;
  // char fname[256];
  std::string fname;
  if(argc==1) {
    // printf("Data file: \n");
    std::cout << "Put your \"data file here\": " << std::endl;
    std::cin >> fname;
    // gets(fname);
    if(fname.length()==0)
      fname="THM_TeST2.dat";
  }
  else if(argc==2) {
    fname = argv[1];
  }
  else{
    usage(argv[0]);
    exit(0);
  }

  if(FBEMInput.ReadFile(fname)){
    TTHMBEM FBEM(FBEMInput);
    FBEM.Execute();
  }
  else{
    usage(argv[0]);
  }
  return 0;
}

