#include "BEMData.h"

TBEMData::TBEMData()
{
  D = Mu = Nu = 0;
  N = L = M = 0;
  NI = NE = NT = NS = 0;
}

TBEMData &TBEMData::operator=(TBEMData &data)
{
  // for input
  N = data.N;
  L = data.L;
  M = data.M;
  NI = data.NI;
  NE = data.NE;
  NT = data.NT;
  NS = data.NS;
  NSH = data.NSH;
  NSF = data.NSF;

  Mu = data.Mu;
  Nu = data.Nu;
  Lambda = data.Lambda;
  K_T = data.K_T;
  K_P = data.K_P;
  Beta_TS = data.Beta_TS;
  Alpha_T = data.Alpha_T;
  Alpha_P = data.Alpha_P;
  T_0 = data.T_0;
  Eta = data.Eta;
  Rho = data.Rho;
  Cp = data.Cp;
  K = data.K;
  B = data.B; 

  Code = (data.Code);                 //.CopyFrom is replaced with =operator
  Last = (data.Last);                 //.CopyFrom is replaced with =operator
  X = (data.X);                       //.CopyFrom is replaced with =operator
  Y = (data.Y);                       //.CopyFrom is replaced with =operator
  Xi = (data.Xi);                     //.CopyFrom is replaced with =operator
  Yi = (data.Yi);                     //.CopyFrom is replaced with =operator
  neXi = (data.neXi);                 //.CopyFrom is replaced with =operator
  neYi = (data.neYi);                 //.CopyFrom is replaced with =operator
  Bc = (data.Bc);                     //.CopyFrom is replaced with =operator
  BcBk = (data.BcBk);                 //.CopyFrom is replaced with =operator
  Time = (data.Time);                 //.CopyFrom is replaced with =operator
  Node = (data.Node);                 //.CopyFrom is replaced with =operator
  Elem = (data.Elem);                 //.CopyFrom is replaced with =operator
  Area = (data.Area);                 //.CopyFrom is replaced with =operator
  OriVolStrain = (data.OriVolStrain); //.CopyFrom is replaced with =operator
  OriTemp = (data.OriTemp);           //.CopyFrom is replaced with =operator
  OriPress = (data.OriPress);         //.CopyFrom is replaced with =operator

  CurVolStrain = (data.CurVolStrain); //.CopyFrom is replaced with =operator
  CurTemp = (data.CurTemp);           //.CopyFrom is replaced with =operator
  CurPress = (data.CurPress);         //.CopyFrom is replaced with =operator

  HeatSrc = (data.HeatSrc);   //.CopyFrom is replaced with =operator
  HeatX = (data.HeatX);       //.CopyFrom is replaced with =operator
  HeatY = (data.HeatY);       //.CopyFrom is replaced with =operator
  FluidSrc = (data.FluidSrc); //.CopyFrom is replaced with =operator
  FluidX = (data.FluidX);     //.CopyFrom is replaced with =operator
  FluidY = (data.FluidY);     //.CopyFrom is replaced with =operator

  // for calculation
  D = data.D;
  Xm = (data.Xm); //.CopyFrom is replaced with =operator
  Ym = (data.Ym); //.CopyFrom is replaced with =operator
  F = (data.F);   //.CopyFrom is replaced with =operator

  Displ = (data.Displ); //.CopyFrom is replaced with =operator
  Temp = (data.Temp);   //.CopyFrom is replaced with =operator
  Press = (data.Press); //.CopyFrom is replaced with =operator

  Stress = (data.Stress);       //.CopyFrom is replaced with =operator
  ThermFlux = (data.ThermFlux); //.CopyFrom is replaced with =operator
  FlowRate = (data.FlowRate);   //.CopyFrom is replaced with =operator

  G = data.G;
  H = data.H;
  return *this;
}

