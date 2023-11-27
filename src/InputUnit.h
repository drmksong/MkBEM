//---------------------------------------------------------------------------
#ifndef InputUnitH
#define InputUnitH
//---------------------------------------------------------------------------
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "MkArray.hpp"
#include "MkMatrix.hpp"
#include "BEMData.h"

class TBEMInput
{
private:
   std::string FileName; // input file name
   TBEMData FData;
   bool readStatus;
   int _i0 = 0, _i1 = 1, _i2 = 2, _i3 = 3;

public:
   TBEMInput()
   {
      FileName = "";
      FData.Clear(); //TODO: need to be implemented, remove todo when implemented
   }

   TBEMInput(std::string fname) : FileName(fname)
   {
      readStatus = ReadFile();
   };
   // void SetFileName(char *fname) { strcpy(FileName, fname); } // obsolete and need to be removed
   void SetFileName(std::string fname) { FileName = fname; }
   bool ReadFile();
   bool ReadFile(std::string fname)
   {
      SetFileName(fname);
      return readStatus = ReadFile();
   }

   TBEMData &GetData() { return FData; }
   TBEMInput &operator=(TBEMInput &inp);
   // char *GetName() { return FileName; } // obsolete and need to be removed
   std::string GetName() const { return FileName; }
};

//---------------------------------------------------------------------------
#endif
