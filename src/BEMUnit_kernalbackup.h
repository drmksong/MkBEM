//---------------------------------------------------------------------------
#ifndef BEMUnit_kernalbackupH
#define BEMUnit_kernalbackupH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <ComCtrls.hpp>
#include <stdlib.h>
#include <stdio.h>
#include "InputUnit.h"
//---------------------------------------------------------------------------
class TBEM : public TThread
{
private:
  TBEMData FData;
  float B1,B2; // B1 = lambda_1*s, B2 = lambda_2*s, s is Laplace transform coeff.
  float B3,B4;
  bool isInputed;
  Float Z,W;
  TProgressBar *FProgBar;
  TLabel *FLabel;
  float a1_ij,a2_ij,a3_ij,a4_ij,a5_ij;
  float a1_it,a2_it,a3_it;
  float a1_ip,a2_ip,a3_ip;
  float a1_tj,a2_tj,a3_tj;
  float a1_pj,a2_pj,a3_pj;
  float a1_tp,a1_pt;
  float a1_tt,a2_tt;
  float a1_pp,a2_pp;
  float Tau;
  float Q;

  float G11, G12, G13, G14, G21, G22, G23, G24, G31, G32, G33, G34, G41, G42, G43, G44;
  float H11, H12, H13, H14, H21, H22, H23, H24, H31, H32, H33, H34, H41, H42, H43, H44;
  float g111,g112,g121,g122,g131,g132,g141,g142,
        g211,g212,g221,g222,g231,g232,g241,g242,
        g311,g312,g321,g322,g331,g332,g341,g342,
        g411,g412,g421,g422,g431,g432,g441,g442;
  float h111,h112,h121,h122,h131,h132,h141,h142,
        h211,h212,h221,h222,h231,h232,h241,h242,
        h311,h312,h321,h322,h331,h332,h341,h342,
        h411,h412,h421,h422,h431,h432,h441,h442;
  float Bi1,Bi2,Bi3,Bi4;

protected:
        void __fastcall Execute();
public:
        __fastcall TBEM(bool CreateSuspended);
        __fastcall void CalcParam();
        __fastcall void Sys();
        __fastcall void Inter();

        __fastcall void InterDisp();
        __fastcall void InterTemp();
        __fastcall void InterPressure();

        __fastcall void InterStress();
        __fastcall void InterThermFlux();
        __fastcall void InterFlowRate();

        __fastcall void InterPrimary();
        __fastcall void InterSecondary();

        __fastcall void Solve(TMatrix A,Float B,float &D,int N);
        __fastcall void Solve();

        __fastcall void Quad(const int i,const int j,const int k);
        __fastcall void Diag(const int i,const int j,const int k);

        __fastcall void Output();
        __fastcall void SetInput(TBEMInput &input){FData = input;}
        __fastcall void SetData(TBEMData &data){FData = data;}
        __fastcall void SetProgBar(TProgressBar *prog_bar){FProgBar = prog_bar;}
        __fastcall void SetLabel(TLabel *label){FLabel = label;}
        __fastcall void ClearG();
        __fastcall void ClearH();
        __fastcall void ClearGH(){ClearG(); ClearH(); }
        __fastcall void Clearg();
        __fastcall void Clearh();
        __fastcall void Cleargh(){Clearg(); Clearh(); }
        __fastcall void CalcG(float xp,float yp, float x1, float y1,float x2, float y2);
        __fastcall void CalcG(float x1, float y1,float x2, float y2);
        __fastcall void CalcH(float xp,float yp, float x1, float y1,float x2, float y2);
        __fastcall void CalcH(float x1, float y1,float x2, float y2);
        __fastcall void Calcg(float xp,float yp, float x1, float y1,float x2, float y2);
        __fastcall void Calch(float xp,float yp, float x1, float y1,float x2, float y2);
        __fastcall void CalcBi(float x,float y);
};
//---------------------------------------------------------------------------
float E1(float);

class TTHMBEM : public TBEM {
public:
   TTHMBEM(TProgressBar *progbar, TLabel *label, TBEMInput &input);
   TTHMBEM(TProgressBar *progbar, TLabel *label, TBEMData &data);
};
#endif
