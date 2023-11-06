//---------------------------------------------------------------------------
#include "MainUnit.h"
//---------------------------------------------------------------------------
void usage(char *fname)
{
  printf("usage : thermo-hydro-mechanical analysis program\n",fname);
  printf("      : %s data_file.dat\n",fname);
}

int main(int argc, char *argv[])
{
    char fname[256];
    if(argc==1) {
      printf("Data file: \n");
      gets(fname);
      if(strlen(fname)==0)
        strcpy(fname,"THM_TeST2.dat");
    }
    else if(argc==2) {
      strcpy(fname, argv[1]);
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
#ifdef __BCPLUSPLUS__
    getc(stdin);
#endif
    return 0;
}

