//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "BEMUnit_kernalbackup.h"
#pragma package(smart_init)
//---------------------------------------------------------------------------
//   Important: Methods and properties of objects in VCL can only be
//   used in a method called using Synchronize, for example:
//
//      Synchronize(UpdateCaption);
//
//   where UpdateCaption could look like:
//
//      void __fastcall Unit1::UpdateCaption()
//      {
//        Form1->Caption = "Updated in a thread";
//      }
//---------------------------------------------------------------------------

float sqr(float a){return a*a;}
float cube(float a){return a*a*a;}
float E1(float x)
{
  float result;
  if ((x > 0) && (x <= 1)) {
    result = - 0.57721566 + 0.99999193*x
              - 0.24991055*x*x + 0.05519968 * x*x*x
              - 0.00976004*x*x*x*x + 0.00107857*pow(x,5) - log(x);
  }
  else if (x > 1) {
    result = (pow(x,4)+8.57332874*pow(x,3)+18.05901697*pow(x,2)+
             8.63476089*x+0.26777373) /
             ( (pow(x,4)+9.57332234*pow(x,3)+25.63295614*pow(x,2)+
                21.09965308*x+3.95849692)*x*exp(x) );
  }
//  else raise exception;
  return result;
}

__fastcall TBEM::TBEM(bool CreateSuspended)
        : TThread(CreateSuspended)
{
  Z.Initialize(8);
  W.Initialize(8);

  Z(0) = 0.183434642495650;
  Z(2) = 0.525532409916329;
  Z(4) = 0.796666477413627;
  Z(6) = 0.960289856497536;
  Z(1) = -0.183434642495650;
  Z(3) = -0.525532409916329;
  Z(5) = -0.796666477413627;
  Z(7) = -0.960289856497536;

  W(0) = 0.362683783378362;
  W(2) = 0.313706645877887;
  W(4) = 0.222381034453374;
  W(6) = 0.101228536290376;
  W(1) = 0.362683783378362;
  W(3) = 0.313706645877887;
  W(5) = 0.222381034453374;
  W(7) = 0.101228536290376;

  a1_ij=a2_ij=a3_ij=a4_ij=a5_ij=0;
  a1_it=a2_it=a3_it=0;
  a1_ip=a2_ip=a3_ip=0;
  a1_tj=a2_tj=a3_tj=0;
  a1_pj=a2_pj=a3_pj=0;
  a1_tp=a1_pt=0;
  a1_tt=a2_tt=0;
  a1_pp=a2_pp=0;

  FreeOnTerminate = True;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::CalcParam()
{
  float lm,lmm,at,ap,t0,kp,kt,eta,mu,nu,rho,cp;
  at = FData.Alpha_T;
  ap = FData.Alpha_P;
  t0 = FData.T_0;
  kp = FData.K_P;
  kt = FData.K_T;
  eta = FData.Eta;
  mu = FData.Mu;
  nu = FData.Nu;
  rho = FData.Rho;
  cp = FData.Cp;
  FData.Lambda =2*FData.Mu*FData.Nu/(1-2*FData.Nu);
  FData.K = 3*FData.Lambda+2*FData.Mu;
  lm = FData.Lambda+FData.Mu;
  lmm = FData.Lambda+2*FData.Mu;
  float a,b,a1,b1,b2;
  B1 = 2;
  B2 = 1;
  B3 =sqrt(FData.Cp/FData.K_P);
  B4 =sqrt(ap/(kp*FData.K*FData.B));
  a = rho*cp/kt+ap/(kp*FData.K*FData.B)+(at*at*t0*kp+ap*ap*kt)/(lmm*kp*kt);
  b = rho*cp*ap/(kt*kp*FData.K*FData.B)+(-lmm*eta*eta*t0
      +at*at*t0*B4*B4*kp+ap*ap*B3*B3*kt-2*at*eta*t0*ap)/(lmm*kt*kp);
  a1 = a*a - 4*b;
  if (a1 > 0) {
    b1 = (a+sqrt(a1))/2;
    b2 = (a-sqrt(a1))/2;
  }
  if (b1>0 && b2>0) {
     B1 =sqrt(b1);
     B2 =sqrt(b2);
  }

  a1_ij=-mu*(at*at*t0*kp*(B1*B1-B4*B4)+ap*ap*kt*(B1*B1-B3*B3)-2*at*eta*t0*ap)
         /(B1*B1*(B1*B1-B2*B2)*lm*lmm);
  a2_ij=-mu*(at*at*t0*kp*(B2*B2-B4*B4)+ap*ap*kt*(B2*B2-B3*B3)-2*at*eta*t0*ap)
         /(B2*B2*(B1*B1-B2*B2)*lm*lmm);
  a3_ij=mu*(at*at*t0*kp*(B1*B1*B2*B2-(B1*B1+B2*B2)*B4*B4)
           +ap*ap*kt*(B1*B1*B2*B2-(B1*B1+B2*B2)*B3*B3)
           -2*at*eta*t0*ap*(B1*B1+B2*B2))/(B1*B1*B1*B2*B2*B2*lm*lmm);
  Q = mu*(at*at*t0*kp*B4*B4+ap*ap*kt*B3*B3+2*at*eta*t0*ap)/(B1*B1*B2*B2*lm*lmm);
                                          // changed - to +  a4_ij=(3-4*nu)*kt*kp/4-Q/4;             // changed  a5_ij=kt*kp/2+Q/2;                      // changed
  a1_it=(eta*ap+at*kp*(B1*B1-B4*B4))/(B1*B1-B2*B2);
  a2_it=(eta*ap+at*kp*(B2*B2-B4*B4))/(B1*B1-B2*B2);
  a3_it=(eta*ap-at*kp*B4*B4)/(B1*B2);
  a1_ip=(eta*at*t0+ap*kt*(B1*B1-B3*B3))/(B1*B1-B2*B2);
  a2_ip=(eta*at*t0+ap*kt*(B2*B2-B3*B3))/(B1*B1-B2*B2);
  a3_ip=(eta*at*t0-ap*kt*B3*B3)/(B1*B2);
  a1_tj=(eta*t0*ap+at*t0*kp*(B1*B1-B4*B4))/(B1*B1-B2*B2);
  a2_tj=(eta*t0*ap+at*t0*kp*(B2*B2-B4*B4))/(B1*B1-B2*B2);
  a3_tj=(eta*t0*ap-at*t0*kp*B4*B4)/(B1*B2);
  a1_pj=a1_ip;
  a2_pj=a2_ip;
  a3_pj=a3_ip;
  a1_tp=-(at*t0*ap+lmm*eta*t0)/(B1*B1-B2*B2);
  a1_pt=-(at*ap+lmm*eta)/(B1*B1-B2*B2);
  a1_tt=(ap*ap+lmm*kp*(B1*B1-B4*B4))/(B1*B1-B2*B2);
  a2_tt=(ap*ap+lmm*kp*(B2*B2-B4*B4))/(B1*B1-B2*B2);
  a1_pp=(at*at+lmm*kt*(B1*B1-B3*B3))/(B1*B1-B2*B2);
  a2_pp=(at*at+lmm*kt*(B2*B2-B3*B3))/(B1*B1-B2*B2);
}
//---------------------------------------------------------------------------
void __fastcall TBEM::Execute()
{
        //---- Place thread code here ----
  int NN,i;
  char str[256];

  FProgBar->Position = 0;
  CalcParam();
  FData.BcBk.CopyFrom(FData.Bc);
  for (int nt=0;nt<FData.NT;nt++) {
      Tau = FData.Time(nt);
      NN = FData.N;
      for (i = 0 ; i<NN;i++) {
          sprintf(str," %d-th Traction X is %f ",i, FData.F(4*i));
      }
      
      if (nt!=0) FData.Bc.CopyFrom(FData.BcBk);

//  { Compute matrices G and H, and form the system A X = F }
      FLabel->Caption = "Compute matrices G and H, and form the system A X = F ";
      Sys();
      FProgBar->StepIt();

//  {Solve the system AX = F }

      FLabel->Caption = "Compute the system A X = F ";
      Solve();

      for (i = 0 ; i<NN ;i++) {
          sprintf(str," %d-th Traction X is %f ",i, FData.F(4*i));
      }
      FProgBar->StepIt();
//  { Compute the potential at interior points }
      FLabel->Caption = "Compute the potential at interior points";
      Inter();
      InterDisp();
      InterPressure();
      InterTemp();
//      InterPrimary();
//      InterSecondary();
      InterStress();
      InterThermFlux();
      InterFlowRate();
      FProgBar->StepIt();

//  { Output solution at boundary nodes and interior points }
      FLabel->Caption = "Output solution at boundary nodes and interior points";
      Output();
      FProgBar->StepIt();
  }
//  FProgBar->Position = 0;
  FLabel->Caption = "Processing Terminated";

}
//---------------------------------------------------------------------------
void __fastcall TBEM::Sys()
{

  int i,j,k,kk,NN,found;
  float qx,qy,ux,uy,temp;
  float dummy,dummy2;
  float dumarr[11];

  found = 0;
 // This function computes the matrices G and H. and forms the system A X = F }

  for (i=0 ; i < FData.N;i++) {
    FData.Xm(i) = (FData.X(i) + FData.X(i+1))/2.0;
    FData.Ym(i) = (FData.Y(i) + FData.Y(i+1))/2.0;
  }

  if((FData.M-1) > 0) {
    FData.Xm(FData.Last(0)) = (FData.X(FData.Last(0))  + FData.X(0))/2.0;
    FData.Ym(FData.Last(0)) = (FData.Y(FData.Last(0))  + FData.Y(0))/2.0;
    for (k=1 ; k < FData.M ; k++ ) {
      FData.Xm(FData.Last(k))=(FData.X(FData.Last(k))+FData.X(FData.Last(k-1)+1))/2.0;
      FData.Ym(FData.Last(k))=(FData.Y(FData.Last(k))+FData.Y(FData.Last(k-1)+1))/2.0;
    }
  }

// Compute the coefficients of G and H matrices }

  for (i=0 ; i < FData.N ; i++) {
    for (j=0 ; j < FData.N ; j++) {
      if(FData.M > 1) {
        if (j==FData.Last(0)) kk = 0;
        else {
          found = 0;
          for (k=1 ; k < FData.M ; k++) {
            if (j==FData.Last(k))  {
              kk = FData.Last(k-1) + 1;
              found = 1;
              break;
            }
          }
          if (found == 0)  kk = j + 1;
        }
      }
      else {
        kk = j + 1;
        if (i !=j )  {
          Quad(i,j,kk);
        }

        else {
          Diag(i,j,kk);
        }
      }
    }
  }
  //	{ Reorder the columns of the system of equations as in (5.28)
  //	 and form the system matrix A which is stored G}

  NN = 4*FData.N;
  for (j=0 ; j < NN ;j++) {
    if(FData.Code(j) > 0)  {
      for (i=0 ; i<NN ; i++)  {
    	  temp = FData.G(i,j);
    	  FData.G(i,j) = -FData.H(i,j);
    	  FData.H(i,j) = -temp;
      }
    }
    else if (NN%4 == 0 || NN%4 == 1)
      for (i = 0 ; i < NN ; i++)
        FData.G(i,j) = FData.G(i,j)*FData.Mu;
  }

  //	{ Form the right-side vector F which is stored in F }

  for (i=0 ; i < NN ; i++)   {
    FData.F(i) = 0.0;
    for (j=0 ; j < NN ; j++){
      FData.F(i) = FData.F(i) + FData.H(i,j) * FData.Bc(j);
    }
  }

  for (i=0; i < FData.N ; i++) {
    float x,y;
    x = FData.Xm(i);
    y = FData.Ym(i);
    CalcBi(x,y);
    FData.F(i+0) += Bi1;
    FData.F(i+1) += Bi2;
    FData.F(i+2) += Bi3;
    FData.F(i+3) += Bi4;
  }
}
//---------------------------------------------------------------------------
// This subroutine computes the G and H matrices coefficients tha
// relate a collocation pint with a different element using Gauss quadrature
//---------------------------------------------------------------------------
void __fastcall TBEM::Quad(int i,int j,int k)
{
  float Xp,Yp,X1,Y1,X2,Y2;

  Xp = FData.Xm(i);
  Yp = FData.Ym(i);
  X1 = FData.X(j);
  Y1 = FData.Y(j);
  X2 = FData.X(k);
  Y2 = FData.Y(k);

  CalcG(Xp, Yp, X1, Y1, X2, Y2);
  CalcH(Xp, Yp, X1, Y1, X2, Y2);

  FData.G(4*i,4*j)     = G11;
  FData.G(4*i,4*j+1)   = G12;
  FData.G(4*i,4*j+2)   = G13;
  FData.G(4*i,4*j+3)   = G14;
  FData.G(4*i+1,4*j+0) = G21;
  FData.G(4*i+1,4*j+1) = G22;
  FData.G(4*i+1,4*j+2) = G23;
  FData.G(4*i+1,4*j+3) = G24;
  FData.G(4*i+2,4*j+0) = G31;
  FData.G(4*i+2,4*j+1) = G32;
  FData.G(4*i+2,4*j+2) = G33;
  FData.G(4*i+2,4*j+3) = G34;
  FData.G(4*i+3,4*j+0) = G41;
  FData.G(4*i+3,4*j+1) = G42;
  FData.G(4*i+3,4*j+2) = G43;
  FData.G(4*i+3,4*j+3) = G44;

  FData.H(4*i,4*j)     = H11;
  FData.H(4*i,4*j+1)   = H12;
  FData.H(4*i,4*j+2)   = H13;
  FData.H(4*i,4*j+3)   = H14;
  FData.H(4*i+1,4*j+0) = H21;
  FData.H(4*i+1,4*j+1) = H22;
  FData.H(4*i+1,4*j+2) = H23;
  FData.H(4*i+1,4*j+3) = H24;
  FData.H(4*i+2,4*j+0) = H31;
  FData.H(4*i+2,4*j+1) = H32;
  FData.H(4*i+2,4*j+2) = H33;
  FData.H(4*i+2,4*j+3) = H34;
  FData.H(4*i+3,4*j+0) = H41;
  FData.H(4*i+3,4*j+1) = H42;
  FData.H(4*i+3,4*j+2) = H43;
  FData.H(4*i+3,4*j+3) = H44;
}
//---------------------------------------------------------------------------
// This subroutine computes the values of the matrix G coefficient
// that relate an element with itself
//---------------------------------------------------------------------------
void __fastcall TBEM::Diag(int i,int j,int k)
{
  float X1,Y1,X2,Y2;

  X1 = FData.X(j);
  Y1 = FData.Y(j);
  X2 = FData.X(k);
  Y2 = FData.Y(k);

  CalcG(X1, Y1,X2, Y2);
  CalcH(X1, Y1,X2, Y2);

  FData.G(4*i,4*j)     = G11;
  FData.G(4*i,4*j+1)   = G12;
  FData.G(4*i,4*j+2)   = G13;
  FData.G(4*i,4*j+3)   = G14;
  FData.G(4*i+1,4*j+0) = G21;
  FData.G(4*i+1,4*j+1) = G22;
  FData.G(4*i+1,4*j+2) = G23;
  FData.G(4*i+1,4*j+3) = G24;
  FData.G(4*i+2,4*j+0) = G31;
  FData.G(4*i+2,4*j+1) = G32;
  FData.G(4*i+2,4*j+2) = G33;
  FData.G(4*i+2,4*j+3) = G34;
  FData.G(4*i+3,4*j+0) = G41;
  FData.G(4*i+3,4*j+1) = G42;
  FData.G(4*i+3,4*j+2) = G43;
  FData.G(4*i+3,4*j+3) = G44;

  FData.H(4*i,4*j)     = H11;
  FData.H(4*i,4*j+1)   = H12;
  FData.H(4*i,4*j+2)   = H13;
  FData.H(4*i,4*j+3)   = H14;
  FData.H(4*i+1,4*j+0) = H21;
  FData.H(4*i+1,4*j+1) = H22;
  FData.H(4*i+1,4*j+2) = H23;
  FData.H(4*i+1,4*j+3) = H24;
  FData.H(4*i+2,4*j+0) = H31;
  FData.H(4*i+2,4*j+1) = H32;
  FData.H(4*i+2,4*j+2) = H33;
  FData.H(4*i+2,4*j+3) = H34;
  FData.H(4*i+3,4*j+0) = H41;
  FData.H(4*i+3,4*j+1) = H42;
  FData.H(4*i+3,4*j+2) = H43;
  FData.H(4*i+3,4*j+3) = H44;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::Solve()
{
  TMatrix A(FData.G);
  TVector X,B(FData.F);

  A.Solve(B,stLUD);
  FData.F.CopyFrom(B.GetFloat());
}
//---------------------------------------------------------------------------
void __fastcall TBEM::Solve(TMatrix A, Float B,float &D,int N)
{

  int N1,i,j,k,i1,l,k1,found;
  float c;

//  { found is a flag which is used to check if any non zero coeff is found }

  N1 = N - 1;
  for ( k=1 ; k<=N1 ; k++ )  {
    k1 = k + 1;
    c = A(k,k);

    if( (abs(c) - EPS) <=0 )  {
      found = 0;
      for ( j=k1 ; k<=N ; k++ )  { // {Interchange rows to get Nonzero }
        if( (abs(A(j,k)) - EPS) > 0.0 )  {
          for ( l=k ; l<=N ; l++ )  {
            c = A(k,l);
            A(k,l) = A(j,l);
            A(j,l) = c;
          }
          c = B(k);
          B(k) = B(j);
          B(j) = c;
          c = A(k,k);
          found = 1;		//{ coeff is found }
          break;
        }
      }
    }

    if(found==0)  {
//      fputs("Singularity in Row %d   ",k);
      D = 0.0;
      exit(EXIT_FAILURE);   //{ If no coeff is found the control is transferred to main }
    }

//    { Divide row by diagonal coefficient }

    c = A(k,k);
    for ( j = k1 ; j<=N ; j++ )  A(k,j) = A(k,j) /c;
    B(k) = B(k) /c;

//    { Eliminate unknown X[k] from row i }

    for ( i = k1 ; i<=N ; i++ )  {
      c = A(i,k);
      for ( j = k1 ; j<=N ; j++ )  A(i,j) = A(i,j) - c*A(k,j);
      B(i) = B(i) -c*B(k);
    }
  }

//  { Compute the last unknown }

  if((abs(A(N,N)) - EPS) > 0.0) {
    B(N) = B(N) /A(N,N);

//  {Apply back substitution to compute the remaining unknowns  }
    for ( l = 1 ; l<=N1 ; l++ )  {
      k = N - l;
      k1 = k +1;
      for ( j = k1 ; j<=N ; j++ )  B(k) = B(k) -A(k,j)*B(j);
    }

//    {Compute the value of the determinent }

    D = 1.0;
    for ( i=1 ; i<=N ; i++ )  D = D *A(i,i);
  }
  else {
//    fputs("Singularity in Row %d   ",k);
    D = 0.0;
  }

}

//---------------------------------------------------------------------------
void __fastcall TBEM::Inter()
{
  int i, NN;
  float temp;
//  { Rearrange the Bc and F arrays to store all values of the
//    potential in Bc and all potential derivatives in F }
  NN = 4*FData.N;
  for (int i=0 ; i < NN ; i++ )  {
    if (FData.Code(i) > 0) {
      temp = FData.Bc(i);
      FData.Bc(i) = FData.F(i);
      FData.F(i) = temp;
    }
    else if (NN%4==0 || NN%4==1) FData.F(i) = FData.F(i)*FData.Mu;
  }

  NN = 4*FData.N;
  for (int j=0 ; j < NN ;j++) {
    if(FData.Code(j) > 0)  {
      for (i=0 ; i<NN ; i++)  {
    	  temp = FData.G(i,j);
    	  FData.G(i,j) = -FData.H(i,j);
    	  FData.H(i,j) = -temp;
      }
    }
    else if (NN%4 == 0 || NN%4 == 1)
      for (i = 0 ; i < NN ; i++)
        FData.G(i,j) = FData.G(i,j)/FData.Mu;
  }

}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterDisp()
{
  int i,j,k, kk,found;
  float temp;

  float d1,d2; // delete it
  float f1,f2;
  float b1,b2;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.Displ(2*i) = 0.0;
      FData.Displ(2*i+1) = 0.0;
      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j == FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j == FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        CalcG(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        CalcH(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.Displ(2*i) = FData.Displ(2*i) + FData.F(4*j)*G11+FData.F(4*j+1)*G12                          +FData.F(4*j+2)*G13+FData.F(4*j+3)*G14-FData.Bc(4*j)*H11
                          -FData.Bc(4*j+1)*H12-FData.Bc(4*j+2)*H13-FData.Bc(4*j+3)*H14;
        FData.Displ(2*i+1) = FData.Displ(2*i+1) + FData.F(4*j)*G21+FData.F(4*j+1)*G22
                          +FData.F(4*j+2)*G23+FData.F(4*j+3)*G24-FData.Bc(4*j)*H21
                          -FData.Bc(4*j+1)*H22-FData.Bc(4*j+2)*H23-FData.Bc(4*j+3)*H24;
        float x,y;
        x = FData.Xi(i);        y = FData.Yi(i);        CalcBi(x,y);        FData.Displ(2*i) -= Bi1;
        FData.Displ(2*i+1) -= Bi2;
      }
    }
  }
}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterTemp()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.Temp(i) = 0.0;
      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j == FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j == FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        CalcG(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        CalcH(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.Temp(i) = FData.Temp(i) + FData.F(4*j)*G31+FData.F(4*j+1)*G32
                        +FData.F(4*j+2)*G33+FData.F(4*j+3)*G34-FData.Bc(4*j)*H31
                        -FData.Bc(4*j+1)*H32-FData.Bc(4*j+2)*H33-FData.Bc(4*j+3)*H34;

        float x,y;
        x = FData.Xi(i);        y = FData.Yi(i);        CalcBi(x,y);        FData.Temp(i) -= Bi3;
      }
    }
  }
}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterPressure()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.Press(i) = 0.0;
      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j == FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j == FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        CalcG(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        CalcH(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.Press(i) = FData.Press(i) + FData.F(4*j)*G41+FData.F(4*j+1)*G42
                         +FData.F(4*j+2)*G43+FData.F(4*j+3)*G44-FData.Bc(4*j)*H41
                         -FData.Bc(4*j+1)*H42-FData.Bc(4*j+2)*H43-FData.Bc(4*j+3)*H44;
        float x,y;
        x = FData.Xi(i);        y = FData.Yi(i);        CalcBi(x,y);        FData.Press(i) -= Bi4;
      }
    }
  }
}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterPrimary()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.Displ(2*i) = 0.0;
      FData.Displ(2*i+1) = 0.0;
      FData.Temp(i) = 0.0;
      FData.Press(i) = 0.0;
      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j == FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j == FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        CalcG(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        CalcH(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.Displ(2*i) = FData.Displ(2*i) + FData.F(4*j)*G11+FData.F(4*j+1)*G12+ FData.F(4*j+2)*G13+FData.F(4*j+3)*G14-FData.Bc(4*j)*H11-FData.Bc(4*j+1)*H12-FData.Bc(4*j+2)*H13-FData.Bc(4*j+3)*H14;
        FData.Displ(2*i+1) = FData.Displ(2*i+1) + FData.F(4*j)*G21+FData.F(4*j+1)*G22+ FData.F(4*j+2)*G23+FData.F(4*j+3)*G24-FData.Bc(4*j)*H21-FData.Bc(4*j+1)*H22-FData.Bc(4*j+2)*H23-FData.Bc(4*j+3)*H24;
        FData.Temp(i) = FData.Temp(i) + FData.F(4*j)*G31+FData.F(4*j+1)*G32+ FData.F(4*j+2)*G33+FData.F(4*j+3)*G34-FData.Bc(4*j)*H31-FData.Bc(4*j+1)*H32-FData.Bc(4*j+2)*H33-FData.Bc(4*j+3)*H34;
        FData.Press(i) = FData.Press(i) + FData.F(4*j)*G41+FData.F(4*j+1)*G42+ FData.F(4*j+2)*G43+FData.F(4*j+3)*G44-FData.Bc(4*j)*H41-FData.Bc(4*j+1)*H42-FData.Bc(4*j+2)*H43-FData.Bc(4*j+3)*H44;
      }
    }
  }

}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterStress()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.Stress(3*i) = 0.0;//
      FData.Stress(3*i+1) = 0.0;//
      FData.Stress(3*i+2) = 0.0;//

      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j==FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j==FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        Calcg(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        Calch(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.Stress(3*i) = FData.Stress(3*i) + FData.F(4*j)*g111+FData.F(4*j+1)*g112+FData.F(4*j+2)*g311+FData.F(4*j+3)*g411-FData.Bc(4*j)*h111-FData.Bc(4*j+1)*h112-FData.Bc(4*j+2)*h311-FData.Bc(4*j+3)*h411;
        FData.Stress(3*i+1) = FData.Stress(3*i+1) + FData.F(4*j)*g121+FData.F(4*j+1)*g122+FData.F(4*j+2)*g312+FData.F(4*j+3)*g412-FData.Bc(4*j)*h121-FData.Bc(4*j+1)*h122-FData.Bc(4*j+2)*h312-FData.Bc(4*j+3)*h412;
        FData.Stress(3*i+2) = FData.Stress(3*i+2) + FData.F(4*j)*g221+FData.F(4*j+1)*g222+FData.F(4*j+2)*g322+FData.F(4*j+3)*g422-FData.Bc(4*j)*h221-FData.Bc(4*j+1)*h222-FData.Bc(4*j+2)*h322-FData.Bc(4*j+3)*h422;
      }
    }
  }
}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterThermFlux()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.ThermFlux(2*i) = 0.0;
      FData.ThermFlux(2*i+1) = 0.0;

      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j==FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j==FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        Calcg(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        Calch(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.ThermFlux(2*i) = FData.ThermFlux(2*i) + FData.F(4*j)*g141+FData.F(4*j+1)*g241+FData.F(4*j+2)*g341+FData.F(4*j+3)*g441-FData.Bc(4*j)*h141-FData.Bc(4*j+1)*h241-FData.Bc(4*j+2)*h341-FData.Bc(4*j+3)*h441;
        FData.ThermFlux(2*i+1) = FData.ThermFlux(2*i+1) + FData.F(4*j)*g142+FData.F(4*j+1)*g242+FData.F(4*j+2)*g342+FData.F(4*j+3)*g442-FData.Bc(4*j)*h142-FData.Bc(4*j+1)*h242-FData.Bc(4*j+2)*h342-FData.Bc(4*j+3)*h442;
      }
    }
  }
}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterFlowRate()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.FlowRate(2*i) = 0.0;
      FData.FlowRate(2*i+1) = 0.0;

      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j==FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j==FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        Calcg(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        Calch(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.FlowRate(2*i) = FData.FlowRate(2*i) + FData.F(4*j)*g131+FData.F(4*j+1)*g231+FData.F(4*j+2)*g331+FData.F(4*j+3)*g431-FData.Bc(4*j)*h131-FData.Bc(4*j+1)*h231-FData.Bc(4*j+2)*h331-FData.Bc(4*j+3)*h431;
        FData.FlowRate(2*i+1) = FData.FlowRate(2*i+1) + FData.F(4*j)*g132+FData.F(4*j+1)*g232+FData.F(4*j+2)*g332+FData.F(4*j+3)*g432-FData.Bc(4*j)*h132-FData.Bc(4*j+1)*h232-FData.Bc(4*j+2)*h332-FData.Bc(4*j+3)*h432;
      }
    }
  }
}
//---------------------------------------------------------------------------
void __fastcall TBEM::InterSecondary()
{
  int i,j,k,kk,found;

//     {Initialization }
  found = 0;
//     { Compute the potential and fluxes at the interior points }

  if(FData.L!=0) {
    for ( i=0 ; i < FData.L ; i++ )  {
      FData.Stress(3*i) = 0.0;//
      FData.Stress(3*i+1) = 0.0;//
      FData.Stress(3*i+2) = 0.0;//
      FData.ThermFlux(2*i) = 0.0;
      FData.ThermFlux(2*i+1) = 0.0;
      FData.FlowRate(2*i) = 0.0;
      FData.FlowRate(2*i+1) = 0.0;

      for ( j=0 ; j < FData.N ; j++ )  {
        if(FData.M > 1) {
          if(j==FData.Last(0)) kk = 0;
          else {
            found = 0;
            for ( k=1 ; k < FData.M ; k++ )  {
              if(j==FData.Last(k)) {
                kk = FData.Last(k-1) + 1;
                found = 1;
                break;
              }
            }
            if(found==0) kk = j + 1;
          }
        }
        else kk = j + 1;
        Calcg(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));
        Calch(FData.Xi(i),FData.Yi(i),FData.X(j),FData.Y(j),FData.X(kk),FData.Y(kk));

        FData.Stress(3*i) = FData.Stress(3*i) + FData.F(4*j)*g111+FData.F(4*j+1)*g112+FData.F(4*j+2)*g311+FData.F(4*j+3)*g411-FData.Bc(4*j)*h111-FData.Bc(4*j+1)*h112-FData.Bc(4*j+2)*h311-FData.Bc(4*j+3)*h411;
        FData.Stress(3*i+1) = FData.Stress(3*i+1) + FData.F(4*j)*g121+FData.F(4*j+1)*g122+FData.F(4*j+2)*g312+FData.F(4*j+3)*g412-FData.Bc(4*j)*h121-FData.Bc(4*j+1)*h122-FData.Bc(4*j+2)*h312-FData.Bc(4*j+3)*h412;
        FData.Stress(3*i+2) = FData.Stress(3*i+2) + FData.F(4*j)*g221+FData.F(4*j+1)*g222+FData.F(4*j+2)*g322+FData.F(4*j+3)*g422-FData.Bc(4*j)*h221-FData.Bc(4*j+1)*h222-FData.Bc(4*j+2)*h322-FData.Bc(4*j+3)*h422;

        FData.ThermFlux(2*i) = FData.ThermFlux(2*i) + FData.F(4*j)*g141+FData.F(4*j+1)*g241+FData.F(4*j+2)*g341+FData.F(4*j+3)*g441-FData.Bc(4*j)*h141-FData.Bc(4*j+1)*h241-FData.Bc(4*j+2)*h341-FData.Bc(4*j+3)*h441;
        FData.ThermFlux(2*i+1) = FData.ThermFlux(2*i+1) + FData.F(4*j)*g142+FData.F(4*j+1)*g242+FData.F(4*j+2)*g342+FData.F(4*j+3)*g442-FData.Bc(4*j)*h142-FData.Bc(4*j+1)*h242-FData.Bc(4*j+2)*h342-FData.Bc(4*j+3)*h442;

        FData.FlowRate(2*i) = FData.FlowRate(2*i) + FData.F(4*j)*g131+FData.F(4*j+1)*g231+FData.F(4*j+2)*g331+FData.F(4*j+3)*g431-FData.Bc(4*j)*h131-FData.Bc(4*j+1)*h231-FData.Bc(4*j+2)*h331-FData.Bc(4*j+3)*h431;
        FData.FlowRate(2*i+1) = FData.FlowRate(2*i+1) + FData.F(4*j)*g132+FData.F(4*j+1)*g232+FData.F(4*j+2)*g332+FData.F(4*j+3)*g432-FData.Bc(4*j)*h132-FData.Bc(4*j+1)*h232-FData.Bc(4*j+2)*h332-FData.Bc(4*j+3)*h432;
      }
    }
  }

}
//---------------------------------------------------------------------------
void __fastcall TBEM::Output()
{
    FILE *outfile;
    char outname[256];
    int i,j,k;
    char str[256];
    static int isfirst=0;

    if (isfirst==0) {
       outfile = fopen("333.out","w");
       fputs("                       Results                \n",outfile);
       fputs("------------------------------------------------------\n",outfile);
       fputs("      (X,Y)     Disp X    Disp Y  Temp Pressure Tract. X   Tract. Y Thermflux Flowrate\n",outfile);

       isfirst = 1;
    }
    else outfile = fopen("333.out","a");

    fprintf(outfile, "Time = %f\n", Tau);

    for ( i=0 ; i<FData.N ; i++ ) {
      sprintf(str,"(%10.5f, %10.5f) %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5e %15.5e \n",
          FData.Xm(i),FData.Ym(i),FData.Bc(4*i),FData.Bc(4*i+1),FData.Bc(4*i+2),FData.Bc(4*i+3),FData.F(4*i),FData.F(4*i+1),FData.F(4*i+2),FData.F(4*i+3));
      fputs(str,outfile);
    }
    if(FData.L!=0)  {
      fputs("",outfile);
      fputs("       	Interior Points Displacement    \n",outfile);
      fputs("      (Xi,Yi)   Disp X   Disp Y  \n",outfile);
      for ( k=0 ; k<FData.L ; k++ ) {
        sprintf(str,"(%10.5f,%10.5f) %15.5f %15.5f %15.5f %15.5f\n",FData.Xi(k),FData.Yi(k),FData.Displ(2*k),FData.Displ(2*k+1),FData.Temp(k),FData.Press(k));
        fputs(str,outfile);
      }
      fputs("",outfile);
      fputs("       	Interior Points Stress    \n",outfile);
      fputs("      (Xi,Yi)   Stress X   Stress Y  \n",outfile);
      for ( k=0 ; k<FData.L ; k++ ) {
        sprintf(str,"(%10.5f,%10.5f) %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f\n",
                       FData.Xi(k),FData.Yi(k),FData.Stress(3*k),FData.Stress(3*k+1),FData.Stress(3*k+2),FData.ThermFlux(2*k),FData.ThermFlux(2*k+1),FData.FlowRate(2*k),FData.FlowRate(2*k+1));
        fputs(str,outfile);
      }

      int NN;
      NN = FData.N;      for (int j = 0 ; j < NN ; j++) {        for (int i = 0 ; i < NN ; i++) {           if (i==0) fputs("|",outfile);           sprintf(str," %9.5f ",FData.G(4*i+2,4*j+2));           fputs(str,outfile);           if (i==NN-1) fputs("|",outfile);        }        sprintf(str," %9.5f = ",FData.F(4*j+2));        fputs(str,outfile);        for (int i = 0 ; i < NN ; i++) {           if (i==0) fputs("|",outfile);           sprintf(str," %9.5f ",FData.H(4*i+2,4*j+2));           fputs(str,outfile);           if (i==NN-1) fputs("|",outfile);        }        sprintf(str," %9.5f \n",FData.Bc(4*j+2));        fputs(str,outfile);      }    }
    fclose(outfile);

}
//---------------------------------------------------------------------------
void __fastcall TBEM::ClearG()
{
  G11=G12=G13=G14=G21=G22=G23=G24=G31=G32=G33=G34=G41=G42=G43=G44=0.0;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::ClearH()
{
  H11=H12=H13=H14=H21=H22=H23=H24=H31=H32=H33=H34=H41=H42=H43=H44=0.0;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::Clearg()
{
  g111=g112=g121=g122=g131=g132=g141=g142=
  g211=g212=g221=g222=g231=g232=g241=g242=
  g311=g312=g321=g322=g331=g332=g341=g342=
  g411=g412=g421=g422=g431=g432=g441=g442=0;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::Clearh()
{
  h111=h112=h121=h122=h131=h132=h141=h142=
  h211=h212=h221=h222=h231=h232=h241=h242=
  h311=h312=h321=h322=h331=h332=h341=h342=
  h411=h412=h421=h422=h431=h432=h441=h442=0;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::CalcBi(float x, float y)
{
  float Ra,eta_1,eta_2,eta_11,eta_12,eta_22,rx,ry;
  float G31, G32, G33, G34, G41, G42, G43, G44;
  Ra = sqrt(sqr(FData.HeatX-x)+sqr(FData.HeatY-y));  rx = (FData.HeatX-x)/Ra;  ry = (FData.HeatY-y)/Ra;
  eta_1 = B1*Ra/sqrt(4*Tau);  eta_2 = B2*Ra/sqrt(4*Tau);  eta_12 = eta_1*eta_2*sqrt(Tau);  eta_11 = eta_1*eta_1*sqrt(Tau);  eta_22 = eta_2*eta_2*sqrt(Tau);  G31 = Ra*rx/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))           *(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_tj/(eta_1*eta_2));
  G32 = Ra*ry/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))           *(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_tj/(eta_1*eta_2));
  G33 = -1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))           *(a1_tt*E1(eta_1*eta_1)-a2_tt*E1(eta_2*eta_2));
  G34 = -1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))           *(a1_tp*E1(eta_1*eta_1)-a1_tp*E1(eta_2*eta_2));  Ra = sqrt(sqr(FData.FluidX-x)+sqr(FData.FluidY-y));  rx = (FData.FluidX-x)/Ra;  ry = (FData.FluidY-y)/Ra;
  eta_1 = B1*Ra/sqrt(4*Tau);  eta_2 = B2*Ra/sqrt(4*Tau);  eta_12 = eta_1*eta_2*sqrt(Tau);  eta_11 = eta_1*eta_1*sqrt(Tau);  eta_22 = eta_2*eta_2*sqrt(Tau);  G41 = Ra*rx/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))           *(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_pj/(eta_1*eta_2));
  G42 = Ra*ry/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))           *(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_pj/(eta_1*eta_2));  G43 = -1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))           *(a1_pt*E1(eta_1*eta_1)-a1_pt*E1(eta_2*eta_2));  G44 = -1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))         *(a1_pp*E1(eta_1*eta_1)-a2_pp*E1(eta_2*eta_2));

  Bi1 = FData.HeatSrc*G31+FData.FluidSrc*G41;
  Bi2 = FData.HeatSrc*G32+FData.FluidSrc*G42;
  Bi3 = FData.HeatSrc*G33+FData.FluidSrc*G43;
  Bi4 = FData.HeatSrc*G34+FData.FluidSrc*G44;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::CalcG(float Xp, float Yp, float X1, float Y1,float X2, float Y2)
{
  float Xg[8],Yg[8];
  float Ax,Bx,Ay,By,HL,nx,ny,rn,sgn,Denom,Ra,rx,ry,slope,Perp;
  int I;

  Ax = (X2 - X1)/2;
  Bx = (X2 + X1)/2;
  Ay = (Y2 - Y1)/2;
  By = (Y2 + Y1)/2;
  HL = sqrt(sqr(Ax) + sqr(Ay));
  nx =  Ay/HL;
  ny = -Ax/HL;

  if (Ax != 0)  { // Warning !!!!!!!!!!!!
    slope = Ay/Ax;
    Perp = abs((slope*Xp-Yp+Y1-slope*X1)/sqrt(sqr(slope)+1));
  }
  else Perp = abs(Xp-X1);

  sgn = (X1-Xp)*(Y2-Yp)-(X2-Xp)*(Y1-Yp);

  if (sgn < 0)  Perp = - Perp;

  ClearGH();

//  { Compute coefficients of G[x,y], H[x,y] }

  for (I=0 ;I < 8;I++ ) {
    float eta_1,eta_2;
    Xg[I] = Ax*Z(I) + Bx;
    Yg[I] = Ay*Z(I) + By;
    Ra = sqrt( sqr(Xp - Xg[I]) + sqr(Yp - Yg[I]));
    rx = (Xg[I]-Xp)/Ra;
    ry = (Yg[I]-Yp)/Ra;
    rn = rx*nx + ry*ny;
    eta_1 = B1*Ra/sqrt(4*Tau);
    eta_2 = B2*Ra/sqrt(4*Tau);

    G11 += -1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
             *((a1_ij*(1/(4*eta_1*eta_1)*expl(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
                  a2_ij*(1/(4*eta_2*eta_2)*expl(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
                  a3_ij/(4*eta_1*eta_2)+
                  a4_ij)*(2*log(Ra))
              +(rx*rx)*(a1_ij/(2*eta_1*eta_1)*expl(-eta_1*eta_1)-
                        a2_ij/(2*eta_2*eta_2)*expl(-eta_2*eta_2)-
                        a3_ij/(2*eta_1*eta_2)-
                        a5_ij))*W(I)*HL; // Gariginal

    G22 += -1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))             *((a1_ij*(1/(4*eta_1*eta_1)*expl(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
                  a2_ij*(1/(4*eta_2*eta_2)*expl(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
                  a3_ij/(4*eta_1*eta_2)+
                  a4_ij)*(2*log(Ra))
              +(ry*ry)*(a1_ij/(2*eta_1*eta_1)*expl(-eta_1*eta_1)-
                        a2_ij/(2*eta_2*eta_2)*expl(-eta_2*eta_2)-
                        a3_ij/(2*eta_1*eta_2)-
                        a5_ij))*W(I)*HL; // Gariginal

/*  // Original
    G11 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *((a1_ij*(1/(4*eta_1*eta_1)*exp(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
              a2_ij*(1/(4*eta_2*eta_2)*exp(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
              a3_ij/(4*eta_1*eta_2)-
              a4_ij*(2*log(Ra)+1))
            -(rx*rx)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                      a3_ij/(2*eta_1*eta_2)-
                      a5_ij))*W(I)*HL;
*/
    G12 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *(-(rx*ry)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                       a3_ij/(2*eta_1*eta_2)-
                       a5_ij))*W(I)*HL;

    G13 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_it*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_it*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_it/(eta_1*eta_2))*W(I)*HL;

    G14 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_ip*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_ip*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_ip/(eta_1*eta_2))*W(I)*HL;

    G21 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *(-(ry*rx)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                       a3_ij/(2*eta_1*eta_2)-
                       a5_ij))*W(I)*HL;

/*  // Original
    G22 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *((a1_ij*(1/(4*eta_1*eta_1)*exp(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
              a2_ij*(1/(4*eta_2*eta_2)*exp(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
              a3_ij/(4*eta_1*eta_2)-
              a4_ij*(2*log(Ra)+1))
            -(ry*ry)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                      a3_ij/(2*eta_1*eta_2)-
                      a5_ij))*W(I)*HL;
*/

    G23 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_it*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_it*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_it/(eta_1*eta_2))*W(I)*HL;

    G24 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_ip*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_ip*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_ip/(eta_1*eta_2))*W(I)*HL;

    G31 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_tj/(eta_1*eta_2))*W(I)*HL;

    G32 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_tj/(eta_1*eta_2))*W(I)*HL;

    G33 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))  //changed - to +
           *(a1_tt*E1(eta_1*eta_1)-a2_tt*E1(eta_2*eta_2))*W(I)*HL;

    G34 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))  //changed - to +
           *(a1_tp*E1(eta_1*eta_1)-a1_tp*E1(eta_2*eta_2))*W(I)*HL;

    G41 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_pj/(eta_1*eta_2))*W(I)*HL;

    G42 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_pj/(eta_1*eta_2))*W(I)*HL;

    G43 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))  //changed - to +
           *(a1_pt*E1(eta_1*eta_1)-a1_pt*E1(eta_2*eta_2))*W(I)*HL;

    G44 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu)) //changed - to +
           *(a1_pp*E1(eta_1*eta_1)-a2_pp*E1(eta_2*eta_2))*W(I)*HL;
  }

  G13= 0;
  G23= 0;
  G43= 0;
//  G14= 0;
//  G24= 0;
//  G34= 0;
  G31= 0;
  G41= 0;
  G32= 0;
  G42= 0;

}
//---------------------------------------------------------------------------
void __fastcall TBEM::CalcG(float X1, float Y1,float X2, float Y2)
{
  float Ax,Bx,Ay,By,HL,Ra,rx,ry;
  int I;

  float Xg[8],Yg[8];
  float Xp,Yp;

  Xp = (X1+X2)/2;
  Yp = (Y1+Y2)/2;

  Ax = (X2 - X1)/2;
  Bx = (X2 + X1)/2;
  Ay = (Y2 - Y1)/2;
  By = (Y2 + Y1)/2;
  HL = sqrt(sqr(Ax) + sqr(Ay));

  ClearG();

  for (I=0 ;I < 8;I++ ) {
    float eta_1,eta_2;
    Xg[I] = Ax*Z(I) + Bx;
    Yg[I] = Ay*Z(I) + By;
    Ra = sqrt( sqr(Xp - Xg[I]) + sqr(Yp - Yg[I]));
    rx = (Xg[I]-Xp)/Ra;
    ry = (Yg[I]-Yp)/Ra;
    eta_1 = B1*Ra/sqrt(4*Tau);
    eta_2 = B2*Ra/sqrt(4*Tau);

    G11 += -1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
             *((a1_ij*(1/(4*eta_1*eta_1)*expl(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
                  a2_ij*(1/(4*eta_2*eta_2)*expl(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
                  a3_ij/(4*eta_1*eta_2)+
                  a4_ij)*(2*log(Ra))
              +(rx*rx)*(a1_ij/(2*eta_1*eta_1)*expl(-eta_1*eta_1)-
                        a2_ij/(2*eta_2*eta_2)*expl(-eta_2*eta_2)-
                        a3_ij/(2*eta_1*eta_2)-
                        a5_ij))*W(I)*HL; // Gariginal

    G22 += -1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))             *((a1_ij*(1/(4*eta_1*eta_1)*expl(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
                  a2_ij*(1/(4*eta_2*eta_2)*expl(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
                  a3_ij/(4*eta_1*eta_2)+
                  a4_ij)*(2*log(Ra))
              +(ry*ry)*(a1_ij/(2*eta_1*eta_1)*expl(-eta_1*eta_1)-
                        a2_ij/(2*eta_2*eta_2)*expl(-eta_2*eta_2)-
                        a3_ij/(2*eta_1*eta_2)-
                        a5_ij))*W(I)*HL; // Gariginal

/*  // Original
    G11 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *((a1_ij*(1/(4*eta_1*eta_1)*exp(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
              a2_ij*(1/(4*eta_2*eta_2)*exp(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
              a3_ij/(4*eta_1*eta_2)-
              a4_ij*(2*log(Ra)+1))
            -(rx*rx)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                      a3_ij/(2*eta_1*eta_2)-
                      a5_ij))*W(I)*HL;
*/
    G12 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *(-(rx*ry)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                       a3_ij/(2*eta_1*eta_2)-
                       a5_ij))*W(I)*HL;

    G13 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_it*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_it*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_it/(eta_1*eta_2))*W(I)*HL;

    G14 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_ip*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_ip*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_ip/(eta_1*eta_2))*W(I)*HL;

    G21 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *(-(ry*rx)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                       a3_ij/(2*eta_1*eta_2)-
                       a5_ij))*W(I)*HL;

/*  // Original
    G22 += 1/(4*M_PI*FData.K_T*FData.K_P*FData.Mu*(1-FData.Nu))
           *((a1_ij*(1/(4*eta_1*eta_1)*exp(-eta_1*eta_1)-1/4*E1(eta_1*eta_1))-
              a2_ij*(1/(4*eta_2*eta_2)*exp(-eta_2*eta_2)-1/4*E1(eta_2*eta_2))-
              a3_ij/(4*eta_1*eta_2)-
              a4_ij*(2*log(Ra)+1))
            -(ry*ry)*(a1_ij/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ij/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-
                      a3_ij/(2*eta_1*eta_2)-
                      a5_ij))*W(I)*HL;
*/

    G23 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_it*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_it*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_it/(eta_1*eta_2))*W(I)*HL;

    G24 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(a1_ip*(1/(eta_1*eta_1)*exp(-eta_1*eta_1)-E1(eta_1*eta_1))-
             a2_ip*(1/(eta_2*eta_2)*exp(-eta_2*eta_2)-E1(eta_2*eta_2))+
             a3_ip/(eta_1*eta_2))*W(I)*HL;

    G31 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_tj/(eta_1*eta_2))*W(I)*HL;

    G32 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_tj/(eta_1*eta_2))*W(I)*HL;

    G33 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))   //changed - to +
           *(a1_tt*E1(eta_1*eta_1)-a2_tt*E1(eta_2*eta_2))*W(I)*HL;

    G34 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))   //changed - to +
           *(a1_tp*E1(eta_1*eta_1)-a1_tp*E1(eta_2*eta_2))*W(I)*HL;

    G41 += Ra*rx/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_pj/(eta_1*eta_2))*W(I)*HL;

    G42 += Ra*ry/(8*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
           *(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
             a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
             a3_pj/(eta_1*eta_2))*W(I)*HL;

    G43 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))  //changed - to +
           *(a1_pt*E1(eta_1*eta_1)-a1_pt*E1(eta_2*eta_2))*W(I)*HL;

    G44 += 1/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))   //changed - to +
           *(a1_pp*E1(eta_1*eta_1)-a2_pp*E1(eta_2*eta_2))*W(I)*HL;
  }

  G13= 0;
  G23= 0;
  G43= 0;
//  G14= 0;
//  G24= 0;
//  G34= 0;
  G31= 0;
  G41= 0;
  G32= 0;
  G42= 0;

}
//---------------------------------------------------------------------------
void __fastcall TBEM::CalcH(float Xp,float Yp, float X1, float Y1,float X2, float Y2)
{
  float Xg[8],Yg[8];
  float Ax,Bx,Ay,By,HL,nx,ny,rn,sgn,Denom,Ra,rx,ry,slope,Perp;
  int I;

  Ax = (X2 - X1)/2;
  Bx = (X2 + X1)/2;
  Ay = (Y2 - Y1)/2;
  By = (Y2 + Y1)/2;
  HL = sqrt(sqr(Ax) + sqr(Ay));
  nx =  Ay/HL;
  ny = -Ax/HL;

  if (Ax != 0)  { // Warning !!!!!!!!!!!!
    slope = Ay/Ax;
    Perp = abs((slope*Xp-Yp+Y1-slope*X1)/sqrt(sqr(slope)+1));
  }
  else Perp = abs(Xp-X1);

  sgn = (X1-Xp)*(Y2-Yp)-(X2-Xp)*(Y1-Yp);

  if (sgn < 0)  Perp = - Perp;

  ClearH();

//  { Compute coefficients of G[x,y], H[x,y] }

  for (I=0 ;I < 8;I++ ) {
    float eta_1,eta_2;
    Xg[I] = Ax*Z(I) + Bx;
    Yg[I] = Ay*Z(I) + By;
    Ra = sqrt( sqr(Xp - Xg[I]) + sqr(Yp - Yg[I]));
    rx = (Xg[I]-Xp)/Ra;
    ry = (Yg[I]-Yp)/Ra;
    rn = rx*nx + ry*ny;
    eta_1 = B1*Ra/sqrt(4*Tau);
    eta_2 = B2*Ra/sqrt(4*Tau);

    H11 += 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(Perp/Ra*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                 +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                 +a3_ij/(eta_1*eta_2)
                 -a4_ij*2
                 +a5_ij)+
             2*Perp/Ra*rx*rx*(a1_ij*((2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1))
                        -a2_ij*((2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2))
                        -a3_ij*2/(eta_1*eta_2)
                        -a5_ij*2)+
             rx*nx*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                    +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)
                    +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                                -a2_ij*exp(-eta_2*eta_2)))+
             rx*nx*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*2
                    +a5_ij))*W(I)*HL;

    H12 += 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(2*Perp/Ra*rx*ry*(a1_ij*((2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1))
                        -a2_ij*((2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2))
                        -a3_ij*2/(eta_1*eta_2)
                        -a5_ij*2)+
             rx*ny*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                    +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)
                    +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                                -a2_ij*exp(-eta_2*eta_2)))+
             ry*nx*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*2
                    +a5_ij))*W(I)*HL;

    H13 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(nx*2*(a1_it*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1)-1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_1*eta_1))
                -a2_it*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2)-1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_2*eta_2))
                +a3_it/(2*eta_1*eta_2))
             -2*rx*Perp/Ra*(a1_it/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                            a2_it/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                            a3_it/(eta_1*eta_2)));

    H14 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(nx*2*(a1_ip*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1)
                       -1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_1*eta_1))-
                 a2_ip*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2)
                       -1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_2*eta_2))+
                 a3_ip/(2*eta_1*eta_2))
             -2*rx*Perp/Ra*(a1_ip/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_ip/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                       a3_ip/(eta_1*eta_2)));

    H21 += 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(2*Perp/Ra*ry*rx*(a1_ij*((2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1))
                        -a2_ij*((2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2))
                        -a3_ij*2/(eta_1*eta_2)
                        -a5_ij*2)+
             ry*nx*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                    +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)
                    +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                                -a2_ij*exp(-eta_2*eta_2)))+
             rx*ny*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*2
                    +a5_ij))*W(I)*HL;

    H22 += 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(Perp/Ra*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                 +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                 +a3_ij/(eta_1*eta_2)
                 -a4_ij*2
                 +a5_ij)+
             2*Perp/Ra*ry*ry*(a1_ij*((2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1))
                        -a2_ij*((2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2))
                        -a3_ij*2/(eta_1*eta_2)
                        -a5_ij*2)+
             ry*ny*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                    +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)
                    +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                                -a2_ij*exp(-eta_2*eta_2)))+
             ry*ny*(-a1_ij/(eta_1*eta_1)*exp(-eta_1*eta_1)
                    +a2_ij/(eta_2*eta_2)*exp(-eta_2*eta_2)
                    +a3_ij/(eta_1*eta_2)
                    -a4_ij*2
                    +a5_ij))*W(I)*HL;

    H23 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(ny*2*(a1_it*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1)-1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_1*eta_1))-
                 a2_it*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2)-1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_2*eta_2))+
                 a3_it/(2*eta_1*eta_2))
             -2*ry*Perp/Ra*(a1_it/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_it/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                       a3_it/(eta_1*eta_2)));

    H24 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*(FData.Lambda+2*FData.Mu))
           *(ny*2*(a1_ip*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1)-1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_1*eta_1))-
                   a2_ip*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2)-1/2*2*FData.Nu/(1-2*FData.Nu)*E1(eta_2*eta_2))+
                   a3_ip/(2*eta_1*eta_2))
             -2*ry*Perp/Ra*(a1_ip/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                       a2_ip/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                       a3_ip/(eta_1*eta_2)));

    H31 -= 1/(8*M_PI*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *(nx*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                  a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                  a3_tj/(eta_1*eta_2))
              -2*rx*Perp/Ra*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_tj/(eta_1*eta_2)))*W(I)*HL;

    H32 -= 1/(8*M_PI*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *(ny*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                  a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                  a3_tj/(eta_1*eta_2))
              -2*ry*Perp/Ra*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_tj/(eta_1*eta_2)))*W(I)*HL;

    H33 -= 1/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(Perp/Ra*(a1_tt*exp(-eta_1*eta_1)-a2_tt*exp(-eta_2*eta_2)))*W(I)*HL;

    H34 -= 1/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(Perp/Ra*(a1_tp*exp(-eta_1*eta_1)-a1_tp*exp(-eta_2*eta_2)))*W(I)*HL;

    H41 -= 1/(8*M_PI*FData.K_T*Tau*(FData.Lambda+2*FData.Mu))
            *(nx*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                  a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                  a3_pj/(eta_1*eta_2))
              -2*rx*Perp/Ra*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_pj/(eta_1*eta_2)))*W(I)*HL;

    H42 -= 1/(8*M_PI*FData.K_T*Tau*(FData.Lambda+2*FData.Mu))
            *(ny*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                  a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                  a3_pj/(eta_1*eta_2))
              -2*ry*Perp/Ra*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_pj/(eta_1*eta_2)))*W(I)*HL;

    H43 -= 1/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(Perp/Ra*(a1_pt*exp(-eta_1*eta_1)-a1_pt*exp(-eta_2*eta_2)))*W(I)*HL;

    H44 -= 1/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(Perp/Ra*(a1_pp*exp(-eta_1*eta_1)-a2_pp*exp(-eta_2*eta_2)))*W(I)*HL;
  }

  H13= 0;
  H23= 0;
  H43= 0;
//  H14= 0;
//  H24= 0;
//  H34= 0;
  H31= 0;
  H41= 0;
  H32= 0;
  H42= 0;

}
//---------------------------------------------------------------------------
void __fastcall TBEM::CalcH(float X1,float Y1,float X2,float Y2)
{
  float Ax,Bx,Ay,By,HL,Ra,rx,ry;
  int I;

  float Xg[8],Yg[8];

  Ax = (X2 - X1)/2;
  Bx = (X2 + X1)/2;
  Ay = (Y2 - Y1)/2;
  By = (Y2 + Y1)/2;
  HL = sqrt(sqr(Ax) + sqr(Ay));

  H11 = 0.5;
  H12 = 0.0;
  H13 = 0.0;
  H14 = 0.0;
  H21 = 0.0;
  H22 = 0.5;
  H23 = 0.0;
  H24 = 0.0;
  H31 = 0.0;
  H32 = 0.0;
  H33 = 0.5;
  H34 = 0.0;
  H41 = 0.0;
  H42 = 0.0;
  H43 = 0.0;
  H44 = 0.5;
}
//---------------------------------------------------------------------------
void __fastcall TBEM::Calcg(float Xp,float Yp, float X1,float Y1,float X2,float Y2)
{
  float Xg[8],Yg[8];
  float Ax,Bx,Ay,By,HL,nx,ny,rn,sgn,Denom,Ra,rx,rx2,rx3,ry,ry2,ry3,rxy,ryx,slope,Perp;
  int I;

  Ax = (X2 - X1)/2;
  Bx = (X2 + X1)/2;
  Ay = (Y2 - Y1)/2;
  By = (Y2 + Y1)/2;
  HL = sqrt(sqr(Ax) + sqr(Ay));
  nx =  Ay/HL;
  ny = -Ax/HL;

  if (Ax != 0)  { // Warning !!!!!!!!!!!!
    slope = Ay/Ax;
    Perp = abs((slope*Xp-Yp+Y1-slope*X1)/sqrt(sqr(slope)+1));
  }
  else Perp = abs(Xp-X1);

  sgn = (X1-Xp)*(Y2-Yp)-(X2-Xp)*(Y1-Yp);

  if (sgn < 0)  Perp = - Perp;

  Clearg();

//  { Compute coefficients of G[x,y], H[x,y] }

  for (I=0 ;I < 8;I++ ) {
    float eta_1,eta_2,eta_12,eta_11,eta_22;
    Xg[I] = Ax*Z(I) + Bx;
    Yg[I] = Ay*Z(I) + By;
    Ra = sqrt( sqr(Xp - Xg[I]) + sqr(Yp - Yg[I]));
    rx = (Xg[I]-Xp)/Ra;
    rx2 = rx*rx;
    rx3 = rx2*rx;
    ry = (Yg[I]-Yp)/Ra;
    ry2 = ry*ry;
    ry3 = ry2*ry;
    rxy = rx*ry;
    ryx = ry*rx;
    rn = rx*nx + ry*ny;
    eta_1 = B1*Ra/sqrt(4*Tau);
    eta_2 = B2*Ra/sqrt(4*Tau);
    eta_11 = eta_1*eta_1;
    eta_22 = eta_2*eta_2;
    eta_12 = eta_1*eta_2;

    g111 -= 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(2*rx*(-a1_ij/(eta_11)*exp(-eta_11)
                   +a2_ij/(eta_22)*exp(-eta_22)
                   +a3_ij/(eta_12)
                   -a4_ij*2
                   +a5_ij)+
             rx*(-a1_ij/(eta_11)*exp(-eta_11)
                 +a2_ij/(eta_22)*exp(-eta_22)
                 +a3_ij/(eta_12)
                 -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                 +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)  //4를 2로 바꿔야 한다.
                 +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_11)
                                 -a2_ij*exp(-eta_22)))+
             2*rx3*( a1_ij*(2+eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*(2+eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*2/(eta_12)
                    -2*a5_ij))*W(I)*HL;            // 2를 빼먹었었음.

    g112 -= 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(ry*(-a1_ij/(eta_11)*exp(-eta_11)
                 +a2_ij/(eta_22)*exp(-eta_22)
                 +a3_ij/(eta_12)
                 -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                 +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)   // 4를 2로 바꿔야 한다.
                 +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_11)
                                             -a2_ij*exp(-eta_22)))+
             2*rx2*ry*( a1_ij*(2+eta_11)/(eta_11)*exp(-eta_11)
                       -a2_ij*(2+eta_22)/(eta_22)*exp(-eta_22)
                       -a3_ij*2/(eta_12)
                       -a5_ij*2))*W(I)*HL;

    g121 -= 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(ry*(-a1_ij/(eta_11)*exp(-eta_11)
                 +a2_ij/(eta_22)*exp(-eta_22)
                 +a3_ij/(eta_12)
                 -a4_ij*2
                 +a5_ij)+
             2*rx2*ry*( a1_ij*(2+eta_11)/(eta_11)*exp(-eta_11)
                       -a2_ij*(2+eta_22)/(eta_22)*exp(-eta_22)
                       -a3_ij*2/(eta_12)
                       -a5_ij*2))*W(I)*HL;

    g122 -= 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(rx*(-a1_ij/(eta_11)*exp(-eta_11)
                 +a2_ij/(eta_22)*exp(-eta_22)
                 +a3_ij/(eta_12)     // eta_2 로 잘못되어있음.
                 -a4_ij*2
                 +a5_ij)+
             2*ry2*rx*( a1_ij*(2+eta_11)/(eta_11)*exp(-eta_11)
                       -a2_ij*(2+eta_22)/(eta_22)*exp(-eta_22)
                       -a3_ij*2/(eta_12)
                       -a5_ij*2))*W(I)*HL;

    g221 -= 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(rx*(-a1_ij/(eta_11)*exp(-eta_11)
                 +a2_ij/(eta_22)*exp(-eta_22)
                 +a3_ij/(eta_12)
                 -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                 +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)
                 +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_11)
                                             -a2_ij*exp(-eta_22)))+
             2*ry2*rx*( a1_ij*(2+eta_11)/(eta_11)*exp(-eta_11)
                       -a2_ij*(2+eta_22)/(eta_22)*exp(-eta_22)
                       -a3_ij*2/(eta_12)
                       -a5_ij*2))*W(I)*HL;

    g222 -= 1/(4*M_PI*FData.K_T*FData.K_P*Ra*(1-FData.Nu))
           *(2*ry*(-a1_ij/(eta_11)*exp(-eta_11)
                   +a2_ij/(eta_22)*exp(-eta_22)
                   +a3_ij/(eta_12)
                   -a4_ij*2
                   +a5_ij)+
             ry*(-a1_ij/(eta_11)*exp(-eta_11)
                 +a2_ij/(eta_22)*exp(-eta_22)
                 +a3_ij/(eta_12)
                 -a4_ij*4*FData.Nu/(1-2*FData.Nu)
                 +a5_ij*2*(1-FData.Nu)/(1-2*FData.Nu)
                 +2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_11)
                                             -a2_ij*exp(-eta_22)))+
             2*ry3*( a1_ij*(2+eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*(2+eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*2/(eta_12)
                    -a5_ij*2))*W(I)*HL;


    g311 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *((a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_tj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*rx2*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2)))*W(I)*HL;

    g312 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *(-2*rxy*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2)))*W(I)*HL;

    g321 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *(-2*ryx*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2)))*W(I)*HL;

    g322 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *((a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_tj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*ry2*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2)))*W(I)*HL;


    g411 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *((a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_pj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*rx2*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2)))*W(I)*HL;

    g412 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *(-2*rxy*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2)))*W(I)*HL;

    g421 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *(-2*ryx*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2)))*W(I)*HL;

    g422 -= FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Tau*(FData.Lambda+2*FData.Mu))
            *((a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_pj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*ry2*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2)))*W(I)*HL;

    g131 += FData.Mu/(8*M_PI*FData.K_P*(FData.Lambda+2*FData.Mu))
            *(2*(a1_it*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1))-
                 a2_it*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2))+
                 a3_it/(2*eta_1*eta_2))-
               2*rx2*(a1_it/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_it/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_it/(eta_1*eta_2)))*W(I)*HL;

    g132 += FData.Mu/(8*M_PI*FData.K_P*(FData.Lambda+2*FData.Mu))
            *(-2*rxy*(a1_it/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_it/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_it/(eta_1*eta_2)))*W(I)*HL;

    g231 += FData.Mu/(8*M_PI*FData.K_P*(FData.Lambda+2*FData.Mu))
            *(-2*ryx*(a1_it/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_it/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_it/(eta_1*eta_2)))*W(I)*HL;

    g232 += FData.Mu/(8*M_PI*FData.K_P*(FData.Lambda+2*FData.Mu))
            *(2*(a1_it*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1))-
                 a2_it*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2))+
                 a3_it/(2*eta_1*eta_2))-
               2*ry2*(a1_it/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_it/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_it/(eta_1*eta_2)))*W(I)*HL;

    g141 += FData.Mu/(8*M_PI*FData.K_T*(FData.Lambda+2*FData.Mu))
            *(2*(a1_ip*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1))-
                 a2_ip*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2))+
                 a3_ip/(2*eta_1*eta_2))-
               2*rx2*(a1_ip/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ip/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_ip/(eta_1*eta_2)))*W(I)*HL;

    g142 += FData.Mu/(8*M_PI*FData.K_T*(FData.Lambda+2*FData.Mu))
            *(-2*rxy*(a1_ip/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ip/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_ip/(eta_1*eta_2)))*W(I)*HL;

    g241 += FData.Mu/(8*M_PI*FData.K_T*(FData.Lambda+2*FData.Mu))
            *(-2*ryx*(a1_ip/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ip/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_ip/(eta_1*eta_2)))*W(I)*HL;

    g242 += FData.Mu/(8*M_PI*FData.K_T*(FData.Lambda+2*FData.Mu))
            *(2*(a1_ip*(1/(2*eta_1*eta_1)*exp(-eta_1*eta_1)-1/2*E1(eta_1*eta_1))-
                 a2_ip*(1/(2*eta_2*eta_2)*exp(-eta_2*eta_2)-1/2*E1(eta_2*eta_2))+
                 a3_ip/(2*eta_1*eta_2))-
               2*ry2*(a1_ip/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_ip/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_ip/(eta_1*eta_2)))*W(I)*HL;

    g331 -= FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(rx*(a1_tt*exp(-eta_1*eta_1)-
                  a2_tt*exp(-eta_2*eta_2)))*W(I)*HL;

    g332 -= FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(ry*(a1_tt*exp(-eta_1*eta_1)-
                  a2_tt*exp(-eta_2*eta_2)))*W(I)*HL;

    g441 -= FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(rx*(a1_pp*exp(-eta_1*eta_1)-
                  a2_pp*exp(-eta_2*eta_2)))*W(I)*HL;

    g442 -= FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(ry*(a1_pp*exp(-eta_1*eta_1)-
                  a2_pp*exp(-eta_2*eta_2)))*W(I)*HL;

    g341 -= FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(rx*(a1_tp*exp(-eta_1*eta_1)-
                  a1_tp*exp(-eta_2*eta_2)))*W(I)*HL;

    g342 -= FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(ry*(a1_tp*exp(-eta_1*eta_1)-
                  a1_tp*exp(-eta_2*eta_2)))*W(I)*HL;

    g431 -= FData.Mu/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(rx*(a1_pt*exp(-eta_1*eta_1)-
                  a1_pt*exp(-eta_2*eta_2)))*W(I)*HL;

    g432 -= FData.Mu/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(ry*(a1_pt*exp(-eta_1*eta_1)-
                  a1_pt*exp(-eta_2*eta_2)))*W(I)*HL;
  }

}
//---------------------------------------------------------------------------
void __fastcall TBEM::Calch(float Xp,float Yp, float X1,float Y1,float X2,float Y2)
{
  float Xg[8],Yg[8];
  float Ax,Bx,Ay,By,HL,nx,ny,rn,sgn,Denom,Ra,Ra2,rx,rx2,rx3,ry,ry2,ry3,rxy,ryx,slope,Perp;
  int I;

  Ax = (X2 - X1)/2;
  Bx = (X2 + X1)/2;
  Ay = (Y2 - Y1)/2;
  By = (Y2 + Y1)/2;
  HL = sqrt(sqr(Ax) + sqr(Ay));
  nx =  Ay/HL;
  ny = -Ax/HL;

  if (Ax != 0)  { // Warning !!!!!!!!!!!!
    slope = Ay/Ax;
    Perp = abs((slope*Xp-Yp+Y1-slope*X1)/sqrt(sqr(slope)+1));
  }
  else Perp = abs(Xp-X1);

  sgn = (X1-Xp)*(Y2-Yp)-(X2-Xp)*(Y1-Yp);

  if (sgn < 0)  Perp = - Perp;

  Clearh();

//  { Compute coefficients of G[x,y], H[x,y] }

  for (I=0 ;I < 8;I++ ) {
    float eta_1,eta_2,eta_12,eta_11,eta_22;
    Xg[I] = Ax*Z(I) + Bx;
    Yg[I] = Ay*Z(I) + By;
    Ra = sqrt( sqr(Xp - Xg[I]) + sqr(Yp - Yg[I]));
    Ra2 = Ra*Ra;
    rx = (Xg[I]-Xp)/Ra;
    rx2 = rx*rx;
    rx3 = rx2*rx;
    ry = (Yg[I]-Yp)/Ra;
    ry2 = ry*ry;
    ry3 = ry2*ry;
    rxy = rx*ry;
    ryx = ry*rx;
    rn = rx*nx + ry*ny;
    eta_1 = B1*Ra/sqrt(4*Tau);
    eta_2 = B2*Ra/sqrt(4*Tau);
    eta_11 = eta_1*eta_1;
    eta_22 = eta_2*eta_2;
    eta_12 = eta_1*eta_2;

    h111 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Ra2*(1-FData.Nu))
         *(rx*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  -a4_ij*4
                  +a5_ij*6)+
           rx*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  +a5_ij*8
                  +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1+eta_11)*exp(-eta_11)
                                 -a2_ij*(1+eta_22)*exp(-eta_22)
                                 -a4_ij*2
                                 +a5_ij))+
           rx*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  -a4_ij*4
                  +a5_ij*6)+
           4*rx3*rn*(a1_ij*2*(6+4*eta_11+eta_11*eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*2*(6+4*eta_22+eta_22*eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*12/(eta_12)
                    -a5_ij*8)+
           rx2*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)

                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)+
           rx2*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*16*FData.Nu/(1-2*FData.Nu)
                   +a5_ij*8*(1-FData.Nu)/(1-2*FData.Nu)
                   +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1-eta_11)*exp(-eta_11)-a2_ij*(1-eta_22)*exp(-eta_22)))+
           rx2*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)-
           nx*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*4
               +a5_ij*2)-
           nx*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*4
               +a5_ij*2)-
           nx*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*8*FData.Nu/(1-2*FData.Nu)
               +a5_ij*4*(1-FData.Nu)/(1-2*FData.Nu)
               +4*FData.Nu /(1-2*FData.Nu)*(a1_ij*2*exp(-eta_11)
                               -a2_ij*2*exp(-eta_22)
                               -a4_ij*2
                               +a5_ij)
               +8*FData.Nu*FData.Nu/(1-2*FData.Nu)/(1-2*FData.Nu)*(a1_ij*eta_11*exp(-eta_11)
                                          -a2_ij*eta_22*exp(-eta_22))))*W(I)*HL;

    h112 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Ra2*(1-FData.Nu))
         *(ry*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  +a5_ij*8
                  +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1+eta_11)*exp(-eta_11)
                                 -a2_ij*(1+eta_22)*exp(-eta_22)
                                 -a4_ij*2
                                 +a5_ij))+
           4*rx2*ry*rn*(a1_ij*2*(6+4*eta_11+eta_11*eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*2*(6+4*eta_22+eta_22*eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*12/(eta_12)
                    -a5_ij*8)+
           rxy*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)+
           rx2*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*16*FData.Nu/(1-2*FData.Nu)
                   +a5_ij*8*(1-FData.Nu)/(1-2*FData.Nu)
                   +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1-eta_11)*exp(-eta_11)-a2_ij*(1-eta_22)*exp(-eta_22)))+
           rxy*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)-
           ny*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*8*FData.Nu/(1-2*FData.Nu)
               +a5_ij*4*(1-FData.Nu)/(1-2*FData.Nu)
               +4*FData.Nu /(1-2*FData.Nu)*(a1_ij*2*exp(-eta_11)
                               -a2_ij*2*exp(-eta_22)
                               -a4_ij*2
                               +a5_ij)
               +8*FData.Nu*FData.Nu/(1-2*FData.Nu)/(1-2*FData.Nu)*(a1_ij*eta_11*exp(-eta_11)
                                          -a2_ij*eta_22*exp(-eta_22))))*W(I)*HL;

    h121 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Ra2*(1-FData.Nu))
         *(ry*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  -a4_ij*4
                  +a5_ij*6)+
           4*rx2*ry*rn*(a1_ij*2*(6+4*eta_11+eta_11*eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*2*(6+4*eta_22+eta_22*eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*12/(eta_12)
                    -a5_ij*8)+
           rxy*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)+
           rxy*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*16*FData.Nu/(1-2*FData.Nu)
                   +a5_ij*8*(1-FData.Nu)/(1-2*FData.Nu)
                   +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1-eta_11)*exp(-eta_11)-a2_ij*(1-eta_22)*exp(-eta_22)))+
           rx2*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)-
           ny*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*4
               +a5_ij*2))*W(I)*HL;

    h122 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Ra2*(1-FData.Nu))
         *(rx*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  -a4_ij*4
                  +a5_ij*6)+
           4*rx*ry2*rn*(a1_ij*2*(6+4*eta_11+eta_11*eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*2*(6+4*eta_22+eta_22*eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*12/(eta_12)
                    -a5_ij*8)+
           ry2*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)+
           rxy*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*16*FData.Nu/(1-2*FData.Nu)
                   +a5_ij*8*(1-FData.Nu)/(1-2*FData.Nu)
                   +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1-eta_11)*exp(-eta_11)-a2_ij*(1-eta_22)*exp(-eta_22)))+
           rxy*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)-
           nx*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*4
               +a5_ij*2))*W(I)*HL;

    h221 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Ra2*(1-FData.Nu))
           *(rx*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  +a5_ij*8
                  +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1+eta_11)*exp(-eta_11)
                                 -a2_ij*(1+eta_22)*exp(-eta_22)
                                 -a4_ij*2
                                 +a5_ij))+
           4*rx*ry2*rn*(a1_ij*2*(6+4*eta_11+eta_11*eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*2*(6+4*eta_22+eta_22*eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*12/(eta_12)
                    -a5_ij*8)+
           rxy*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)+
           ry2*nx*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*16*FData.Nu/(1-2*FData.Nu)
                   +a5_ij*8*(1-FData.Nu)/(1-2*FData.Nu)
                   +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1-eta_11)*exp(-eta_11)-a2_ij*(1-eta_22)*exp(-eta_22)))+
           rxy*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)-
           nx*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*8*FData.Nu/(1-2*FData.Nu)
               +a5_ij*4*(1-FData.Nu)/(1-2*FData.Nu)
               +4*FData.Nu /(1-2*FData.Nu)*(a1_ij*2*exp(-eta_11)
                               -a2_ij*2*exp(-eta_22)
                               -a4_ij*2
                               +a5_ij)
               +8*FData.Nu*FData.Nu/(1-2*FData.Nu)/(1-2*FData.Nu)*(a1_ij*eta_11*exp(-eta_11)
                                          -a2_ij*eta_22*exp(-eta_22))))*W(I)*HL;

    h222 += FData.Mu/(4*M_PI*FData.K_T*FData.K_P*Ra2*(1-FData.Nu))
         *(ry*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  -a4_ij*4
                  +a5_ij*6)+
           ry*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  +a5_ij*8
                  +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1+eta_11)*exp(-eta_11)
                                 -a2_ij*(1+eta_22)*exp(-eta_22)
                                 -a4_ij*2
                                 +a5_ij))+
           ry*rn*(-a1_ij*((8+4*eta_11)/(eta_11)*exp(-eta_11))
                  +a2_ij*((8+4*eta_22)/(eta_22)*exp(-eta_22))
                  +a3_ij*8/(eta_12)
                  -a4_ij*4
                  +a5_ij*6)+
           4*ry3*rn*(a1_ij*2*(6+4*eta_11+eta_11*eta_11)/(eta_11)*exp(-eta_11)
                    -a2_ij*2*(6+4*eta_22+eta_22*eta_22)/(eta_22)*exp(-eta_22)
                    -a3_ij*12/(eta_12)
                    -a5_ij*8)+
           ry2*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)+
           ry2*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*16*FData.Nu/(1-2*FData.Nu)
                   +a5_ij*8*(1-FData.Nu)/(1-2*FData.Nu)
                   +8*FData.Nu/(1-2*FData.Nu)*(a1_ij*(1-eta_11)*exp(-eta_11)-a2_ij*(1-eta_22)*exp(-eta_22)))+
           ry2*ny*(-a1_ij*2*(4+2*eta_11)/(eta_11)*exp(-eta_11)
                   +a2_ij*2*(4+2*eta_22)/(eta_22)*exp(-eta_22)
                   +a3_ij*8/(eta_12)
                   -a4_ij*4
                   +a5_ij*6)-
           ny*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*4
               +a5_ij*2)-
           ny*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*4
               +a5_ij*2)-
           ny*(-a1_ij*2/(eta_11)*exp(-eta_11)
               +a2_ij*2/(eta_22)*exp(-eta_22)
               +a3_ij*2/(eta_12)
               -a4_ij*8*FData.Nu/(1-2*FData.Nu)
               +a5_ij*4*(1-FData.Nu)/(1-2*FData.Nu)
               +4*FData.Nu /(1-2*FData.Nu)*(a1_ij*2*exp(-eta_11)
                               -a2_ij*2*exp(-eta_22)
                               -a4_ij*2
                               +a5_ij)
               +8*FData.Nu*FData.Nu/(1-2*FData.Nu)/(1-2*FData.Nu)*(a1_ij*eta_11*exp(-eta_11)
                                          -a2_ij*eta_22*exp(-eta_22))))*W(I)*HL;

    h311 -= FData.Mu/(2*M_PI*FData.K_T*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *((a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_tj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*rx2*(a1_tj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj*2/(eta_1*eta_2))+
               2*nx*rx*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h312 -= FData.Mu/(2*M_PI*FData.K_T*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *(-2*rxy*(a1_tj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj*2/(eta_1*eta_2))+
               nx*ry*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2))+
               ny*rx*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h321 -= FData.Mu/(2*M_PI*FData.K_T*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *(-2*ryx*(a1_tj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj*2/(eta_1*eta_2))+
               ny*rx*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2))+
               nx*ry*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h322 -= FData.Mu/(2*M_PI*FData.K_T*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *((a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_tj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*ry2*(a1_tj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_tj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_tj*2/(eta_1*eta_2))+
               2*ny*ry*(a1_tj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_tj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h411 -= FData.Mu/(2*M_PI*FData.K_P*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *((a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_pj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*rx2*(a1_pj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj*2/(eta_1*eta_2))+
               2*nx*rx*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h412 -= FData.Mu/(2*M_PI*FData.K_P*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *(-2*rxy*(a1_pj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj*2/(eta_1*eta_2))+
               nx*ry*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2))+
               ny*rx*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h421 -= FData.Mu/(2*M_PI*FData.K_P*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *(-2*ryx*(a1_pj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj*2/(eta_1*eta_2))+
               ny*rx*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2))+
               nx*ry*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h422 -= FData.Mu/(2*M_PI*FData.K_P*Tau*Ra*(FData.Lambda+2*FData.Mu))
            *((a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
               a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
               a3_pj/(eta_1*eta_2)-
               2*FData.Nu /(1-2*FData.Nu)*(a1_ij*exp(-eta_1*eta_1)
                                          -a2_ij*exp(-eta_2*eta_2)))-
               2*ry2*(a1_pj*(2+2*eta_1*eta_1+eta_1*eta_1*eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                      a2_pj*(2+2*eta_2*eta_2+eta_2*eta_2*eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                      a3_pj*2/(eta_1*eta_2))+
               2*ny*ry*(a1_pj*(1+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_pj*(1+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                        a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h131 += FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(2*rx2*(a1_tj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_tj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_tj*2/(eta_1*eta_2))-
             nx*rx*(a1_tj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_tj*exp(-eta_1*eta_1)
                                               -a2_tj*exp(-eta_2*eta_2)))-
             nx*rx*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj/(eta_1*eta_2))-
             rn*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                 a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                 a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h132 += FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(2*rxy*(a1_tj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_tj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_tj*2/(eta_1*eta_2))-
             nx*ry*(a1_tj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_tj*exp(-eta_1*eta_1)
                                               -a2_tj*exp(-eta_2*eta_2)))-
             ny*rx*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h231 += FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(2*ryx*(a1_tj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_tj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_tj*2/(eta_1*eta_2))-
             ny*rx*(a1_tj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_tj*exp(-eta_1*eta_1)
                                               -a2_tj*exp(-eta_2*eta_2)))-
             nx*ry*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h232 += FData.Mu/(2*M_PI*FData.K_P*Ra*(FData.Lambda+2*FData.Mu))
            *(2*ry2*(a1_tj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_tj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_tj*2/(eta_1*eta_2))-
             ny*ry*(a1_tj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_tj*exp(-eta_1*eta_1)
                                               -a2_tj*exp(-eta_2*eta_2)))-
             ny*ry*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_tj/(eta_1*eta_2))-
             rn*(a1_tj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                 a2_tj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                 a3_tj/(eta_1*eta_2)))*W(I)*HL;

    h141 += FData.Mu/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(2*rx2*(a1_pj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_pj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_pj*2/(eta_1*eta_2))-
             nx*rx*(a1_pj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_pj*exp(-eta_1*eta_1)
                                               -a2_pj*exp(-eta_2*eta_2)))-
             nx*rx*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj/(eta_1*eta_2))-
             rn*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                 a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                 a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h142 += FData.Mu/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(2*rxy*(a1_pj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_pj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_pj*2/(eta_1*eta_2))-
             nx*ry*(a1_pj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_pj*exp(-eta_1*eta_1)
                                               -a2_pj*exp(-eta_2*eta_2)))-
             ny*rx*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h241 += FData.Mu/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(2*ryx*(a1_pj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_pj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_pj*2/(eta_1*eta_2))-
             ny*rx*(a1_pj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_pj*exp(-eta_1*eta_1)
                                               -a2_pj*exp(-eta_2*eta_2)))-
             nx*ry*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h242 += FData.Mu/(2*M_PI*FData.K_T*Ra*(FData.Lambda+2*FData.Mu))
            *(2*ry2*(a1_pj*(2+eta_1*eta_1)/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                     a2_pj*(2+eta_2*eta_2)/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                     a3_pj*2/(eta_1*eta_2))-
             ny*ry*(a1_pj*2/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj*2/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj*2/(eta_1*eta_2)+
                    2*FData.Nu /(1-2*FData.Nu)*(a1_pj*exp(-eta_1*eta_1)
                                               -a2_pj*exp(-eta_2*eta_2)))-
             ny*ry*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                    a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                    a3_pj/(eta_1*eta_2))-
             rn*(a1_pj/(eta_1*eta_1)*exp(-eta_1*eta_1)-
                 a2_pj/(eta_2*eta_2)*exp(-eta_2*eta_2)+
                 a3_pj/(eta_1*eta_2)))*W(I)*HL;

    h331 += FData.K_T/(2*M_PI*FData.K_P*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(nx*(a1_tt*exp(-eta_1*eta_1)-
                  a2_tt*exp(-eta_2*eta_2))
              -2*rx*rn*(a1_tt*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_tt*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h332 += FData.K_T/(2*M_PI*FData.K_P*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(ny*(a1_tt*exp(-eta_1*eta_1)-
                  a2_tt*exp(-eta_2*eta_2))
              -2*ry*rn*(a1_tt*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_tt*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h441 += FData.K_P/(2*M_PI*FData.K_T*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(nx*(a1_pp*exp(-eta_1*eta_1)-
                  a2_pp*exp(-eta_2*eta_2))
              -2*rx*rn*(a1_pp*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_pp*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h442 += FData.K_P/(2*M_PI*FData.K_T*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(ny*(a1_pp*exp(-eta_1*eta_1)-
                  a2_pp*exp(-eta_2*eta_2))
              -2*ry*rn*(a1_pp*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a2_pp*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h341 += 1/(2*M_PI*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(nx*(a1_tp*exp(-eta_1*eta_1)-
                  a1_tp*exp(-eta_2*eta_2))
              -2*rx*rn*(a1_tp*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a1_tp*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h342 += 1/(2*M_PI*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(ny*(a1_tp*exp(-eta_1*eta_1)-
                  a1_tp*exp(-eta_2*eta_2))
              -2*ry*rn*(a1_tp*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a1_tp*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h431 += 1/(2*M_PI*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(nx*(a1_pt*exp(-eta_1*eta_1)-
                  a1_pt*exp(-eta_2*eta_2))
              -2*rx*rn*(a1_pt*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a1_pt*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;

    h432 += 1/(2*M_PI*Ra*Ra*(FData.Lambda+2*FData.Mu))
            *(ny*(a1_pt*exp(-eta_1*eta_1)-
                  a1_pt*exp(-eta_2*eta_2))
              -2*ry*rn*(a1_pt*(1+eta_1*eta_1)*exp(-eta_1*eta_1)-
                        a1_pt*(1+eta_2*eta_2)*exp(-eta_2*eta_2)))*W(I)*HL;
  }
}
//---------------------------------------------------------------------------
TTHMBEM::TTHMBEM(TProgressBar *progbar, TLabel *label, TBEMInput &input)
:TBEM(False)
{
    SetProgBar(progbar);
    SetLabel(label);
    SetInput(input);

}
//---------------------------------------------------------------------------
TTHMBEM::TTHMBEM(TProgressBar *progbar, TLabel *label, TBEMData &data)
:TBEM(False)
{
    SetProgBar(progbar);
    SetLabel(label);
    SetData(data);
}
//---------------------------------------------------------------------------
 