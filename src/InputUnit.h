//---------------------------------------------------------------------------
#ifndef InputUnitH
#define InputUnitH

#ifndef __BCPLUSPLUS__
 #ifndef __fastcall
  #define __fastcall
 #endif
#endif

//---------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "MkInt.h"
#include "MkDouble.h"
#include "MkMatrix.h"

//#ifndef expl 
//#define expl exp
//#endif

class TBEMInput;

struct TBEMData {
// for input
   int N,L,M;
   //N : Number of boundary node(input),
   //L : Number of internal point for evaluation(calculated by the domain integral type and NE)
   //M : Number of different boundaries(input)
   int NI,NT,NE,NS,NSH,NSF;
   //NI : Number of node for internal cell(input)
   //NT : Number of time step for Time(input)
   //NE : Number of element for internal cell(input)
   //NS : Number of Source
   //NSH : Number of Heat Source
   //NSF : Number of Fluid Srouce
   double Mu,Nu,Lambda,K_T,K_P,Alpha_T,Alpha_P,T_0,Eta,Tau,Rho,Cp,K,B,Beta_TS;
   //Mu : Coeff. of elasticity, Shear modulus(input)
   //Nu : Poisson's ratio(input)
   //Lambda : Coeff. of elasticity(calculated from Mu, Nu)
   //K_T : thermal conductivity of rockmass(input)
   //K_P : hydraulic conductivity of rockmass(input)
   //Alpha_T : Thermal expansion coefficient(input)
   //Alpha_P : Hydraulic expansion coefficient(input)
   //T_0 : absolute initial temperature(input)
   //Eta : Thermal expansion coefficient for pore fluid(input)
   //Tau : Time for step calucalation(allocated from Time array)
   //Rho : Density of rockmass(input)
   //Cp : Specific heat(input)
   //K : Bulk modulus(Calculated from Mu and Nu)
   //B : Skempton's coeff.(input)
   MkInt Code,Last;
   //Code : Index to identify the node is primary or secondary variable(input)
   //Last : Last node number of each boundary(input)
   MkDouble X,Y,Bc,BcBk,Xi,Yi;//Xi,Yi -> Node[NI,2]
   //X,Y : Coordinate of boundaries(input)
   //Bc : Boundary condition(input)
   //BcBk : Back up of Bc
   MkDouble Time,Node,Area;
   MkInt Elem;
   MkDouble neXi,neYi;
   //Time : Times for caculating(input)
   //Node : Coordinate of internal points(input)
   //Elem : Element of internal cell(input)
   //Area : Array of areas of internal cell(input)
   //neXi : x coord. of center of internal cell(computed)
   //neYi : y coord. of center of internal cell(computed)
   MkDouble HeatSrc, HeatX, HeatY;
   MkDouble FluidSrc, FluidX, FluidY;
// for computation
   double D;

   MkDouble Xm,Ym,F; // F : Counter part of Bc, unknowns

   MkDouble Displ, Temp, Press;            // internal primary
   MkDouble OriVolStrain, OriTemp, OriPress;
   MkDouble CurVolStrain, CurTemp, CurPress;
   MkDouble Stress, ThermFlux, FlowRate;   // internal secondary
   MkMatrix G,H;
   MkMatrix g,h;
// member function
   TBEMData();
   TBEMData & operator =(TBEMData &data);
   TBEMData & operator =(TBEMInput &inp);
};

class TBEMInput {
private:
   char FileName[256];
   TBEMData FData;
public:
   TBEMInput();
   void SetFileName(char *fname){strcpy(FileName,fname);}
   bool ReadFile();
   bool ReadFile(char *fname){SetFileName(fname);return ReadFile();}
   TBEMData &GetData(){return FData;}
   TBEMInput & operator =(TBEMInput &inp);
   char *GetName(){return FileName;}
};

//---------------------------------------------------------------------------
#endif
