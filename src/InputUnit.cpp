//---------------------------------------------------------------------------
#include <stdio.h>
#include "InputUnit.h"

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

TBEMData &TBEMData::operator=(TBEMInput &inp)
{
  // for input
  N = inp.GetData().N;
  L = inp.GetData().L;
  M = inp.GetData().M;
  NI = inp.GetData().NI;
  NE = inp.GetData().NE;
  NT = inp.GetData().NT;
  NS = inp.GetData().NS;
  NSH = inp.GetData().NSH;
  NSF = inp.GetData().NSF;

  Mu = inp.GetData().Mu;
  Nu = inp.GetData().Nu;
  Lambda = inp.GetData().Lambda;
  K_T = inp.GetData().K_T;
  K_P = inp.GetData().K_P;
  Beta_TS = inp.GetData().Beta_TS;
  Alpha_T = inp.GetData().Alpha_T;
  Alpha_P = inp.GetData().Alpha_P;
  T_0 = inp.GetData().T_0;
  Eta = inp.GetData().Eta;
  Rho = inp.GetData().Rho;
  Cp = inp.GetData().Cp;
  K = inp.GetData().K;
  B = inp.GetData().B;

  Code = (inp.GetData().Code);                 //.CopyFrom is replaced with =operator
  Last = (inp.GetData().Last);                 //.CopyFrom is replaced with =operator
  X = (inp.GetData().X);                       //.CopyFrom is replaced with =operator
  Y = (inp.GetData().Y);                       //.CopyFrom is replaced with =operator
  Xi = (inp.GetData().Xi);                     //.CopyFrom is replaced with =operator
  Yi = (inp.GetData().Yi);                     //.CopyFrom is replaced with =operator
  neXi = (inp.GetData().neXi);                 //.CopyFrom is replaced with =operator
  neYi = (inp.GetData().neYi);                 //.CopyFrom is replaced with =operator
  Bc = (inp.GetData().Bc);                     //.CopyFrom is replaced with =operator
  BcBk = (inp.GetData().BcBk);                 //.CopyFrom is replaced with =operator
  Time = (inp.GetData().Time);                 //.CopyFrom is replaced with =operator
  Node = (inp.GetData().Node);                 //.CopyFrom is replaced with =operator
  Elem = (inp.GetData().Elem);                 //.CopyFrom is replaced with =operator
  Area = (inp.GetData().Area);                 //.CopyFrom is replaced with =operator
  OriVolStrain = (inp.GetData().OriVolStrain); //.CopyFrom is replaced with =operator
  OriTemp = (inp.GetData().OriTemp);           //.CopyFrom is replaced with =operator
  OriPress = (inp.GetData().OriPress);         //.CopyFrom is replaced with =operator
  CurVolStrain = (inp.GetData().CurVolStrain); //.CopyFrom is replaced with =operator
  CurTemp = (inp.GetData().CurTemp);           //.CopyFrom is replaced with =operator
  CurPress = (inp.GetData().CurPress);         //.CopyFrom is replaced with =operator

  HeatSrc = (inp.GetData().HeatSrc);   //.CopyFrom is replaced with =operator
  HeatX = (inp.GetData().HeatX);       //.CopyFrom is replaced with =operator
  HeatY = (inp.GetData().HeatY);       //.CopyFrom is replaced with =operator
  FluidSrc = (inp.GetData().FluidSrc); //.CopyFrom is replaced with =operator
  FluidX = (inp.GetData().FluidX);     //.CopyFrom is replaced with =operator
  FluidY = (inp.GetData().FluidY);     //.CopyFrom is replaced with =operator

  // for calculation
  D = inp.GetData().D;
  Xm = (inp.GetData().Xm); //.CopyFrom is replaced with =operator
  Ym = (inp.GetData().Ym); //.CopyFrom is replaced with =operator
  F = (inp.GetData().F);   //.CopyFrom is replaced with =operator

  Displ = (inp.GetData().Displ); //.CopyFrom is replaced with =operator
  Temp = (inp.GetData().Temp);   //.CopyFrom is replaced with =operator
  Press = (inp.GetData().Press); //.CopyFrom is replaced with =operator

  Stress = (inp.GetData().Stress);       //.CopyFrom is replaced with =operator
  ThermFlux = (inp.GetData().ThermFlux); //.CopyFrom is replaced with =operator
  FlowRate = (inp.GetData().FlowRate);   //.CopyFrom is replaced with =operator

  G = inp.GetData().G;
  H = inp.GetData().H;
  return *this;
}

TBEMInput::TBEMInput()
{
  memset(FileName, '\0', 255);
}

bool TBEMInput::ReadFile()
{
  FILE *fp;
  char line1[256], str[256];
  int i, k;
  //   int dummy,dummy2;
  //   double sdummy,sdummy2;

  fp = fopen(FileName, "r");

  if (fp == NULL)
    return false;

  fgets(line1, 255, fp);
  snprintf(str, 256, "Input::%s ", line1);
  puts(str);

  fgets(str, 255, fp);
  sscanf(str, "%d %d %d %d", &FData.N, &FData.L, &FData.M, &FData.NS);
  // # of Boundary elem, Interior elem, Different bound.
  if (!FData.N)
    exit(-1);

  snprintf(str, 256, " Number of Boundary Elements = %d ", FData.N);
  puts(str);
  snprintf(str, 256, "Number of Interior Points = %d ", FData.L);
  puts(str);
  snprintf(str, 256, "Number of Source Points = %d ", FData.NS);
  puts(str);

  if (FData.M > 0)
  {
    snprintf(str, 256, "Number of different Boundaries = %d  ", FData.M);
    puts(str);
    if (FData.M)
      FData.Last.Initialize(FData.M);
    for (i = 0; i < FData.M; i++)
      FData.Last(i) = 0;
    fgets(str, 255, fp);
    char *str2;
    str2 = str;
    for (i = 0; i < FData.M; i++)
    {
      sscanf(str2, "%d ", &FData.Last(i));
      str2 = strchr(str2, ' ');
      // TrimLeft(str2,str2); // TrimLeft may not be needed, I can't find the reason why it is used. 2023 11 07 by MK
      snprintf(line1, 256, "Last node on  boundary %2d = %2d ", i, FData.Last(i));
      puts(line1);
    };
  };

  FData.X.Initialize(FData.N + 1);
  FData.Y.Initialize(FData.N + 1);
  FData.Code.Initialize(4 * FData.N);

  FData.Bc.Initialize(4 * FData.N);
  FData.BcBk.Initialize(4 * FData.N);
  FData.F.Initialize(4 * FData.N);

  if (FData.L)
    FData.Xi.Initialize(FData.L);
  if (FData.L)
    FData.Yi.Initialize(FData.L);

  FData.Xm.Initialize(FData.N);
  FData.Ym.Initialize(FData.N);

  //     u.Initialize(L);
  //     q1.Initialize(L);
  //     q2.Initialize(L);
  FData.G.Initialize(4 * FData.N, 4 * FData.N);
  FData.H.Initialize(4 * FData.N, 4 * FData.N);
  if (FData.L)
    FData.Displ.Initialize(2 * FData.L);
  if (FData.L)
    FData.Temp.Initialize(FData.L);
  if (FData.L)
    FData.Press.Initialize(FData.L);

  if (FData.L)
    FData.Stress.Initialize(3 * FData.L);
  if (FData.L)
    FData.ThermFlux.Initialize(2 * FData.L);
  if (FData.L)
    FData.FlowRate.Initialize(2 * FData.L);

  fgets(str, 255, fp);
  sscanf(str, "%d %d %d ", &FData.NI, &FData.NE, &FData.NT);
  // # of NI(Internal Node) NT(Time) NE(Internal Cell) NS(N of Src)

  snprintf(str, 256, " Number of Internal Node = %d ", FData.NI);
  puts(str);
  snprintf(str, 256, " Number of Time Step = %d ", FData.NT);
  puts(str);
  snprintf(str, 256, " Number of Internal Elem = %d ", FData.NE);
  puts(str);

  if (FData.NI)
    FData.Node.Initialize(FData.NI, 2);
  if (FData.NE)
    FData.Elem.Initialize(FData.NE, 4);
  if (FData.NE)
    FData.Area.Initialize(FData.NE);
  if (FData.NT)
    FData.Time.Initialize(FData.NT);
  if (FData.NE)
    FData.OriVolStrain.Initialize(FData.NE);
  if (FData.NE)
    FData.OriTemp.Initialize(FData.NE);
  if (FData.NE)
    FData.OriPress.Initialize(FData.NE);
  if (FData.NE)
    FData.CurVolStrain.Initialize(FData.NE);
  if (FData.NE)
    FData.CurTemp.Initialize(FData.NE);
  if (FData.NE)
    FData.CurPress.Initialize(FData.NE);
  if (FData.NE)
    FData.neXi.Initialize(FData.NE);
  if (FData.NE)
    FData.neYi.Initialize(FData.NE);

  fgets(str, 255, fp);
  sscanf(str, "%lf %lf", &FData.Mu, &FData.Nu);
  snprintf(str, 256, "Shear modulus = %lf ", FData.Mu);
  puts(str);
  snprintf(str, 256, "Poisson ratio = %lf ", FData.Nu);
  puts(str);

  fgets(str, 255, fp);
  sscanf(str, "%lf %lf", &FData.K_T, &FData.K_P);
  snprintf(str, 256, "Thermal Conductivity = %lf ", FData.K_T);
  puts(str);
  snprintf(str, 256, "Hydraulic Conductivity = %lf ", FData.K_P);
  puts(str);

  fgets(str, 255, fp);
  sscanf(str, "%lf %lf %lf", &FData.Beta_TS, &FData.Alpha_P, &FData.T_0);
  snprintf(str, 256, "Coeff of elasticity = %lf ", FData.Alpha_T);
  puts(str);
  snprintf(str, 256, "Thermal Conductivity = %lf ", FData.Alpha_P);
  puts(str);
  snprintf(str, 256, "Hydraulic Conductivity = %lf ", FData.T_0);
  puts(str);

  fgets(str, 255, fp);
  sscanf(str, "%lf %lf %lf", &FData.Eta, &FData.Rho, &FData.Cp);
  snprintf(str, 256, "Thermal expansion coeff = %lf ", FData.Eta);
  puts(str);
  snprintf(str, 256, "Density = %lf ", FData.Rho);
  puts(str);
  snprintf(str, 256, "Specific heat = %lf ", FData.Cp);
  puts(str);

  fgets(str, 255, fp);
  sscanf(str, "%lf", &FData.B);
  snprintf(str, 256, "Skemton coefficient = %lf ", FData.B);
  puts(str);

  strcpy(str, " Cordinates of The Extreme");
  puts(str);
  strcpy(str, " Points of the Boundary Elements");
  puts(str);

  strcpy(str, "Point    X          Y ");
  puts(str);
  for (i = 0; i < FData.N; i++)
  {
    fgets(str, 255, fp);
    sscanf(str, "%lf %lf", &FData.X(i), &FData.Y(i));
    snprintf(str, 256, "%2d %10.4lf %10.4lf ", i, FData.X(i), FData.Y(i));
    puts(str);
  }
  FData.X(FData.N) = FData.X(0);
  FData.Y(FData.N) = FData.Y(0);

  puts(" ");

  /* Read boundary conditions in Bc[i] vector; if Code[i] = 0,  the Bc[i]
    value is a known potential; if Code[i] = 1, the Bc[i] value is a known
    potential derivative (flux).*/

  puts("Boundary Conditions ");
  puts("Node		Code	Prescribed Value ");

  for (i = 0; i < FData.N; i++)
  {
    fgets(str, 255, fp);
    sscanf(str, "%d %lf %d %lf %d %lf %d %lf", &FData.Code(4 * i), &FData.Bc(4 * i),
           &FData.Code(4 * i + 1), &FData.Bc(4 * i + 1),
           &FData.Code(4 * i + 2), &FData.Bc(4 * i + 2),
           &FData.Code(4 * i + 3), &FData.Bc(4 * i + 3));
  }

  for (i = 0; i < FData.N; i++)
  {
    snprintf(str, 256, "%3d : %2d  %10.4lf %2d  %10.4lf %2d  %10.4lf %2d  %10.4lf",
             i, FData.Code(4 * i), FData.Bc(4 * i), FData.Code(4 * i + 1), FData.Bc(4 * i + 1),
             FData.Code(4 * i + 2), FData.Bc(4 * i + 2), FData.Code(4 * i + 3), FData.Bc(4 * i + 3));
    puts(str);
  }

  puts(" ");

  /*  Read coordinates of the interior points   */

  if (FData.L != 0)
  {
    puts(" Interior Point Coordinates ");
    puts("Point    Xi         Yi");
    for (i = 0; i < FData.L; i++)
    {
      fgets(str, 255, fp);
      sscanf(str, "%lf %lf", &FData.Xi(i), &FData.Yi(i));
      snprintf(str, 256, "%2d %10.4lf %10.4lf", i, FData.Xi(i), FData.Yi(i));
      puts(str);
    }
  }

  if (FData.NT != 0)
  {
    puts(" Time Step ");
    puts(" Index Time");
    for (i = 0; i < FData.NT; i++)
    {
      fgets(str, 255, fp);
      sscanf(str, "%lf ", &FData.Time(i));
      snprintf(str, 256, "%2d %10.4lf", i, FData.Time(i));
      puts(str);
    }
  }

  if (FData.NS != 0)
  {
    MkDouble TempSrc(FData.NS), TempX(FData.NS), TempY(FData.NS);
    char(*SrcType)[10] = new char[FData.NS][10];
    for (i = 0; i < FData.NS; i++)
    {
      memset(SrcType[i], '\0', 9);
      fgets(str, 255, fp);
      sscanf(str, "%s %lf %lf %lf", (SrcType[i]), &TempSrc(i), &TempX(i), &TempY(i));
      snprintf(str, 256, "%s %10.4lf %10.4lf %10.4lf", SrcType[i], TempSrc(i), TempX(i), TempY(i));
      puts(str);
    }
    FData.NSH = FData.NSF = 0;
    for (i = 0; i < FData.NS; i++)
    {
      if (!strcmp(SrcType[i], "HEAT"))
        FData.NSH++;
      else if (!strcmp(SrcType[i], "FLUID"))
        FData.NSF++;
    }
    if (FData.NSH)
    {
      FData.HeatSrc.Initialize(FData.NSH);
      FData.HeatX.Initialize(FData.NSH);
      FData.HeatY.Initialize(FData.NSH);
    }
    if (FData.NSF)
    {
      FData.FluidSrc.Initialize(FData.NSF);
      FData.FluidX.Initialize(FData.NSF);
      FData.FluidY.Initialize(FData.NSF);
    }
    FData.NSH = FData.NSF = 0;
    for (i = 0; i < FData.NS; i++)
    {
      if (!strcmp(SrcType[i], "HEAT"))
      {
        double s, x, y;
        FData.HeatSrc(FData.NSH) = TempSrc(i);
        FData.HeatX(FData.NSH) = TempX(i);
        FData.HeatY(FData.NSH) = TempY(i);
        FData.NSH++;
      }
      else if (!strcmp(SrcType[i], "FLUID"))
      {
        FData.FluidSrc(FData.NSF) = TempSrc(i);
        FData.FluidX(FData.NSF) = TempX(i);
        FData.FluidY(FData.NSF) = TempY(i);
        FData.NSF++;
      }
    }
    delete[] SrcType;
  }

  puts("Cell Node    X          Y ");
  for (i = 0; i < FData.NI; i++)
  {
    fgets(str, 255, fp);
    sscanf(str, "%lf %lf", &FData.Node(i, _i0), &FData.Node(i, _i1));
    snprintf(str, 256, "%2d %10.4lf %10.4lf ", i, FData.Node(i, _i0), FData.Node(i, _i1));
    puts(str);
  }

  puts(" ");

  puts("Cell Elem    X       Y ");
  for (i = 0; i < FData.NE; i++)
  {
    fgets(str, 255, fp);
    sscanf(str, "%d %d %d %d %lf %lf %lf", &FData.Elem(i, _i0), &FData.Elem(i, _i1),
           &FData.Elem(i, _i2), &FData.Elem(i, _i3),
           &FData.OriVolStrain(i), &FData.OriTemp(i),
           &FData.OriPress(i));
    FData.Elem(i, _i0)--;
    FData.Elem(i, _i1)--;
    FData.Elem(i, _i2)--;
    FData.Elem(i, _i3)--;

    snprintf(str, 256, "%2d %5d %5d %5d %5d %10.3lf %10.3lf %10.3lf", i, FData.Elem(i, _i0),
             FData.Elem(i, _i1),
             FData.Elem(i, _i2),
             FData.Elem(i, _i3),
             FData.OriVolStrain(i),
             FData.OriTemp(i),
             FData.OriPress(i));
    puts(str);
  }

  for (i = 0; i < FData.NE; i++)
  {
    FData.neXi(i) = 0;
    FData.neYi(i) = 0;
    for (k = 0; k < 4; k++)
    {
      FData.neXi(i) = FData.neXi(i) + FData.Node(FData.Elem(i, k), _i0) / 4.0;
      FData.neYi(i) = FData.neYi(i) + FData.Node(FData.Elem(i, k), _i1) / 4.0;
    }
    snprintf(str, 256, "%2d %10.4lf %10.4lf", i, FData.neXi(i), FData.neYi(i));
    puts(str);
  }

  fclose(fp);
  return true;
}

TBEMInput &TBEMInput::operator=(TBEMInput &inp)
{
  strcpy(FileName, inp.FileName);
  FData = inp.FData;
  return *this;
}
//---------------------------------------------------------------------------
