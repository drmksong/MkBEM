//---------------------------------------------------------------------------
#ifndef BEMUnitH
#define BEMUnitH

#ifndef __BCPLUSPLUS__
 #ifndef __fastcall
  #define __fastcall
 #endif
#else
 #ifndef PROPOUT
  #define PROPOUT
 #endif
#endif

//---------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "InputUnit.h"
#include "BEMData.h"
#include "InputUnit.h"
#include "MkMisc.hpp"

enum ExecType {etSingle,etStep};
//---------------------------------------------------------------------------
class TBEM 
{
private:

  TBEMData FData;
  double B1,B2; // B1 = lambda_1*s, B2 = lambda_2*s, s is Laplace transform coeff.
  double B3,B4;
  bool isInputed;
  MkDouble Z,W;
  double a1_ij,a2_ij,a3_ij,a4_ij,a5_ij;
  double a1_it,a2_it,a3_it;
  double a1_ip,a2_ip,a3_ip;
  double a1_tj,a2_tj,a3_tj;
  double a1_pj,a2_pj,a3_pj;
  double a1_tp,a1_pt;
  double a1_tt,a2_tt;
  double a1_pp,a2_pp;
  double Tau,DTau;
  double Q;
  ExecType fExecType;

  double G11, G12, G13, G14, G21, G22, G23, G24, G31, G32, G33, G34, G41, G42, G43, G44;
  double H11, H12, H13, H14, H21, H22, H23, H24, H31, H32, H33, H34, H41, H42, H43, H44;
  double g111,g112,g121,g122,g131,g132,g141,g142,
        g211,g212,g221,g222,g231,g232,g241,g242,
        g311,g312,g321,g322,g331,g332,g341,g342,
        g411,g412,g421,g422,g431,g432,g441,g442;
  double h111,h112,h121,h122,h131,h132,h141,h142,
        h211,h212,h221,h222,h231,h232,h241,h242,
        h311,h312,h321,h322,h331,h332,h341,h342,
        h411,h412,h421,h422,h431,h432,h441,h442;
  double Bi1,Bi2,Bi3,Bi4;
  double Bi1Src,Bi2Src,Bi3Src,Bi4Src;
  double Bi1Inter,Bi2Inter,Bi3Inter,Bi4Inter;

protected:
      //   char FileName[256];
      std::string FileName;

public:
        TBEM();
        __fastcall void Execute();        
        __fastcall void CalcParam();
        __fastcall void OutProp();
        __fastcall void Sys();
        __fastcall void SysStep();
        __fastcall void Inter();

        __fastcall void InterDisp();
        __fastcall void InterTemp();
        __fastcall void InterPressure();

        __fastcall void InterStress();
        __fastcall void InterThermFlux();
        __fastcall void InterFlowRate();

        __fastcall void InterPrimary();
        __fastcall void InterSecondary();

        __fastcall void Solve(MkMatrix<double> A,MkDouble B,double &D,int N);
        __fastcall void Solve();

        __fastcall void Quad(const int i,const int j,const int k);
        __fastcall void Diag(const int i,const int j,const int k);
        __fastcall void QuadStep(const int i,const int j,const int k);
        __fastcall void DiagStep(const int i,const int j,const int k);

        __fastcall void Output();
        __fastcall void BackupInternal();
        __fastcall void OutDeform();
        __fastcall void SetInput(TBEMInput &input){FData = input.GetData();}
        __fastcall void SetData(TBEMData &data){FData = data;}
        __fastcall void ClearG();
        __fastcall void ClearH();
        __fastcall void ClearGH(){ClearG(); ClearH(); }
        __fastcall void Clearg();
        __fastcall void Clearh();
        __fastcall void Cleargh(){Clearg(); Clearh(); }
        __fastcall void CalcG(double xp,double yp, double x1, double y1,double x2, double y2);
        __fastcall void CalcG(double x1, double y1,double x2, double y2);
        __fastcall void CalcH(double xp,double yp, double x1, double y1,double x2, double y2);
        __fastcall void CalcH(double x1, double y1,double x2, double y2);
        __fastcall void Calcg(double xp,double yp, double x1, double y1,double x2, double y2);
        __fastcall void Calch(double xp,double yp, double x1, double y1,double x2, double y2);
        __fastcall void CalcBi_Ori(double x,double y);
        __fastcall void CalcBi(double x,double y,int e);
        __fastcall void CalcBiSrc(double x,double y);
        __fastcall void CalcBiInter(double x,double y,int e);
        __fastcall void CalcBiStep(double x,double y,int e);
        __fastcall void CalcGi(double xp,double yp, double x1, double y1,double x2, double y2);
        __fastcall void CalcGi(double x1, double y1,double x2, double y2);
        __fastcall void CalcHi(double xp,double yp, double x1, double y1,double x2, double y2);
        __fastcall void CalcHi(double x1, double y1,double x2, double y2);
        __fastcall void TempOut();

};
//---------------------------------------------------------------------------
double E1(double);

class TTHMBEM : public TBEM {
public:
   TTHMBEM(TBEMInput &input);
   TTHMBEM(TBEMData &data);
};
#endif




