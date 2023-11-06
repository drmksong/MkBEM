//---------------------------------------------------------------------------
#include <stdio.h>
#include "InputUnit.h"

TBEMData::TBEMData()
{
   D = Mu = Nu = 0;
   N = L = M = 0;
   NI = NE = NT = NS = 0;
}

TBEMData & TBEMData::operator =(TBEMData &data)
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
   Lambda   = data.Lambda;
   K_T      = data.K_T;
   K_P      = data.K_P;
   Beta_TS  = data.Beta_TS;
   Alpha_T  = data.Alpha_T;
   Alpha_P  = data.Alpha_P;
   T_0      = data.T_0;
   Eta      = data.Eta;
   Rho      = data.Rho;
   Cp       = data.Cp;
   K        = data.K;
   B        = data.B;

   Code.CopyFrom(data.Code);
   Last.CopyFrom(data.Last);
   X.CopyFrom(data.X);
   Y.CopyFrom(data.Y);
   Xi.CopyFrom(data.Xi);
   Yi.CopyFrom(data.Yi);
   neXi.CopyFrom(data.neXi);
   neYi.CopyFrom(data.neYi);
   Bc.CopyFrom(data.Bc);
   BcBk.CopyFrom(data.BcBk);
   Time.CopyFrom(data.Time);
   Node.CopyFrom(data.Node);
   Elem.CopyFrom(data.Elem);
   Area.CopyFrom(data.Area);
   OriVolStrain.CopyFrom(data.OriVolStrain);
   OriTemp.CopyFrom(data.OriTemp);
   OriPress.CopyFrom(data.OriPress);

   CurVolStrain.CopyFrom(data.CurVolStrain);
   CurTemp.CopyFrom(data.CurTemp);
   CurPress.CopyFrom(data.CurPress);

   HeatSrc.CopyFrom(data.HeatSrc);
   HeatX.CopyFrom(data.HeatX);
   HeatY.CopyFrom(data.HeatY);
   FluidSrc.CopyFrom(data.FluidSrc);
   FluidX.CopyFrom(data.FluidX);
   FluidY.CopyFrom(data.FluidY);

// for calculation
   D = data.D;
   Xm.CopyFrom(data.Xm);
   Ym.CopyFrom(data.Ym);
   F.CopyFrom(data.F);

   Displ.CopyFrom(data.Displ);
   Temp.CopyFrom(data.Temp);
   Press.CopyFrom(data.Press);

   Stress.CopyFrom(data.Stress);
   ThermFlux.CopyFrom(data.ThermFlux);
   FlowRate.CopyFrom(data.FlowRate);

   G = data.G;
   H = data.H;
   return *this;
}

TBEMData & TBEMData::operator =(TBEMInput &inp)
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
   Lambda   = inp.GetData().Lambda;
   K_T      = inp.GetData().K_T;
   K_P      = inp.GetData().K_P;
   Beta_TS  = inp.GetData().Beta_TS;
   Alpha_T  = inp.GetData().Alpha_T;
   Alpha_P  = inp.GetData().Alpha_P;
   T_0      = inp.GetData().T_0;
   Eta      = inp.GetData().Eta;
   Rho      = inp.GetData().Rho;
   Cp       = inp.GetData().Cp;
   K        = inp.GetData().K;
   B        = inp.GetData().B;

   Code.CopyFrom(inp.GetData().Code);
   Last.CopyFrom(inp.GetData().Last);
   X.CopyFrom(inp.GetData().X);
   Y.CopyFrom(inp.GetData().Y);
   Xi.CopyFrom(inp.GetData().Xi);
   Yi.CopyFrom(inp.GetData().Yi);
   neXi.CopyFrom(inp.GetData().neXi);
   neYi.CopyFrom(inp.GetData().neYi);
   Bc.CopyFrom(inp.GetData().Bc);
   BcBk.CopyFrom(inp.GetData().BcBk);
   Time.CopyFrom(inp.GetData().Time);
   Node.CopyFrom(inp.GetData().Node);
   Elem.CopyFrom(inp.GetData().Elem);
   Area.CopyFrom(inp.GetData().Area);
   OriVolStrain.CopyFrom(inp.GetData().OriVolStrain);
   OriTemp.CopyFrom(inp.GetData().OriTemp);
   OriPress.CopyFrom(inp.GetData().OriPress);
   CurVolStrain.CopyFrom(inp.GetData().CurVolStrain);
   CurTemp.CopyFrom(inp.GetData().CurTemp);
   CurPress.CopyFrom(inp.GetData().CurPress);

   HeatSrc.CopyFrom(inp.GetData().HeatSrc);
   HeatX.CopyFrom(inp.GetData().HeatX);
   HeatY.CopyFrom(inp.GetData().HeatY);
   FluidSrc.CopyFrom(inp.GetData().FluidSrc);
   FluidX.CopyFrom(inp.GetData().FluidX);
   FluidY.CopyFrom(inp.GetData().FluidY);

// for calculation
   D = inp.GetData().D;
   Xm.CopyFrom(inp.GetData().Xm);
   Ym.CopyFrom(inp.GetData().Ym);
   F.CopyFrom(inp.GetData().F);

   Displ.CopyFrom(inp.GetData().Displ);
   Temp.CopyFrom(inp.GetData().Temp);
   Press.CopyFrom(inp.GetData().Press);

   Stress.CopyFrom(inp.GetData().Stress);
   ThermFlux.CopyFrom(inp.GetData().ThermFlux);
   FlowRate.CopyFrom(inp.GetData().FlowRate);

   G = inp.GetData().G;
   H = inp.GetData().H;
   return *this;
}

TBEMInput::TBEMInput()
{
   memset(FileName,'\0',255);
}

bool TBEMInput::ReadFile()
{
     FILE *fp;
     char line1[256],str[256];
     int i,k;
//   int dummy,dummy2;
//   double sdummy,sdummy2;

     fp = fopen(FileName,"r");

     if(fp==NULL) return false;
     
     fgets(line1,255,fp);
     sprintf(str,"Input::%s ",line1);
     puts(str);

     fgets(str,255,fp);
     sscanf(str,"%d %d %d %d",&FData.N,&FData.L,&FData.M,&FData.NS);
     // # of Boundary elem, Interior elem, Different bound.
     if (!FData.N) exit(-1);

     sprintf(str," Number of Boundary Elements = %d ",FData.N);
     puts(str);
     sprintf(str,"Number of Interior Points = %d ",FData.L);
     puts(str);
     sprintf(str,"Number of Source Points = %d ",FData.NS);
     puts(str);

     if (FData.M > 0) {
          sprintf(str,"Number of different Boundaries = %d  ",FData.M);
          puts(str);
          if (FData.M) FData.Last.Initialize(FData.M);
          for(i=0;i<FData.M;i++) FData.Last(i) = 0;
          fgets(str,255,fp);
          char *str2;
          str2 = str;
          for (i=0 ; i< FData.M;i++) {
            sscanf(str2,"%d ",&FData.Last(i));
            str2 = strchr(str2,' ');
            TrimLeft(str2,str2);
            sprintf(line1,"Last node on  boundary %2d = %2d ",i,FData.Last(i));
            puts(line1);
          };
     };

     FData.X.Initialize(FData.N+1);
     FData.Y.Initialize(FData.N+1);
     FData.Code.Initialize(4*FData.N);

     FData.Bc.Initialize(4*FData.N);
     FData.BcBk.Initialize(4*FData.N);
     FData.F.Initialize(4*FData.N);

     if (FData.L) FData.Xi.Initialize(FData.L);
     if (FData.L) FData.Yi.Initialize(FData.L);

     FData.Xm.Initialize(FData.N);
     FData.Ym.Initialize(FData.N);

//     u.Initialize(L);
//     q1.Initialize(L);
//     q2.Initialize(L);
     FData.G.Initialize(4*FData.N,4*FData.N);
     FData.H.Initialize(4*FData.N,4*FData.N);
     if (FData.L) FData.Displ.Initialize(2*FData.L);
     if (FData.L) FData.Temp.Initialize(FData.L);
     if (FData.L) FData.Press.Initialize(FData.L);

     if (FData.L) FData.Stress.Initialize(3*FData.L);
     if (FData.L) FData.ThermFlux.Initialize(2*FData.L);
     if (FData.L) FData.FlowRate.Initialize(2*FData.L);

     fgets(str,255,fp);
     sscanf(str,"%d %d %d ",&FData.NI,&FData.NE,&FData.NT);
     // # of NI(Internal Node) NT(Time) NE(Internal Cell) NS(N of Src)

     sprintf(str," Number of Internal Node = %d ",FData.NI);
     puts(str);
     sprintf(str," Number of Time Step = %d ",FData.NT);
     puts(str);
     sprintf(str," Number of Internal Elem = %d ",FData.NE);
     puts(str);

     if (FData.NI) FData.Node.Initialize(FData.NI,2);
     if (FData.NE) FData.Elem.Initialize(FData.NE,4);
     if (FData.NE) FData.Area.Initialize(FData.NE);
     if (FData.NT) FData.Time.Initialize(FData.NT);
     if (FData.NE) FData.OriVolStrain.Initialize(FData.NE);
     if (FData.NE) FData.OriTemp.Initialize(FData.NE);
     if (FData.NE) FData.OriPress.Initialize(FData.NE);
     if (FData.NE) FData.CurVolStrain.Initialize(FData.NE);
     if (FData.NE) FData.CurTemp.Initialize(FData.NE);
     if (FData.NE) FData.CurPress.Initialize(FData.NE);
     if (FData.NE) FData.neXi.Initialize(FData.NE);
     if (FData.NE) FData.neYi.Initialize(FData.NE);

     fgets(str,255,fp);
     sscanf(str,"%lf %lf",&FData.Mu,&FData.Nu);
     sprintf(str,"Shear modulus = %lf ",FData.Mu);
     puts(str);
     sprintf(str,"Poisson ratio = %lf ",FData.Nu);
     puts(str);

     fgets(str,255,fp);
     sscanf(str,"%lf %lf",&FData.K_T,&FData.K_P);
     sprintf(str,"Thermal Conductivity = %lf ",FData.K_T);
     puts(str);
     sprintf(str,"Hydraulic Conductivity = %lf ",FData.K_P);
     puts(str);

     fgets(str,255,fp);
     sscanf(str,"%lf %lf %lf",&FData.Beta_TS,&FData.Alpha_P,&FData.T_0);
     sprintf(str,"Coeff of elasticity = %lf ",FData.Alpha_T);
     puts(str);
     sprintf(str,"Thermal Conductivity = %lf ",FData.Alpha_P);
     puts(str);
     sprintf(str,"Hydraulic Conductivity = %lf ",FData.T_0);
     puts(str);

     fgets(str,255,fp);
     sscanf(str,"%lf %lf %lf",&FData.Eta,&FData.Rho,&FData.Cp);
     sprintf(str,"Thermal expansion coeff = %lf ",FData.Eta);
     puts(str);
     sprintf(str,"Density = %lf ",FData.Rho);
     puts(str);
     sprintf(str,"Specific heat = %lf ",FData.Cp);
     puts(str);

     fgets(str,255,fp);
     sscanf(str,"%lf",&FData.B);
     sprintf(str,"Skemton coefficient = %lf ",FData.B);
     puts(str);

     strcpy(str," Cordinates of The Extreme");
     puts(str);
     strcpy(str," Points of the Boundary Elements");
     puts(str);

     strcpy(str,"Point    X          Y ");
     puts(str);
     for( i=0 ; i <FData.N;i++){
       fgets(str,255,fp);
       sscanf(str,"%lf %lf",&FData.X(i),&FData.Y(i));
       sprintf(str,"%2d %10.4lf %10.4lf ",i,FData.X(i),FData.Y(i));
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

     for (i=0 ; i < FData.N ; i++) {
       fgets(str,255,fp);
       sscanf(str,"%d %lf %d %lf %d %lf %d %lf", &FData.Code(4*i),&FData.Bc(4*i),
                                             &FData.Code(4*i+1),&FData.Bc(4*i+1),
                                             &FData.Code(4*i+2),&FData.Bc(4*i+2),
                                             &FData.Code(4*i+3),&FData.Bc(4*i+3));
     }
     
     for (i=0 ; i< FData.N;i++) {
       sprintf(str,"%3d : %2d  %10.4lf %2d  %10.4lf %2d  %10.4lf %2d  %10.4lf",
              i,FData.Code(4*i),FData.Bc(4*i),FData.Code(4*i+1),FData.Bc(4*i+1),
                FData.Code(4*i+2),FData.Bc(4*i+2),FData.Code(4*i+3),FData.Bc(4*i+3));
       puts(str);
     }

     puts(" ");

     /*  Read coordinates of the interior points   */

     if ( FData.L != 0 ) {
       puts(" Interior Point Coordinates ");
       puts("Point    Xi         Yi");
       for (i=0 ; i< FData.L;i++) {
         fgets(str,255,fp);
         sscanf(str,"%lf %lf",&FData.Xi(i),&FData.Yi(i));
         sprintf(str,"%2d %10.4lf %10.4lf",i,FData.Xi(i),FData.Yi(i));
         puts(str);
       }
     }

     if ( FData.NT != 0 ) {
       puts(" Time Step ");
       puts(" Index Time");
       for (i=0 ; i< FData.NT;i++) {
         fgets(str,255,fp);
         sscanf(str,"%lf ",&FData.Time(i));
         sprintf(str,"%2d %10.4lf",i,FData.Time(i));
         puts(str);
       }
     }

     if (FData.NS != 0) {
        MkDouble TempSrc(FData.NS),TempX(FData.NS),TempY(FData.NS);
        char (*SrcType)[10] = new char[FData.NS][10];
        for (i = 0 ; i < FData.NS ; i++) {
            memset(SrcType[i],'\0',9) ;
            fgets(str,255,fp);
            sscanf(str,"%s %lf %lf %lf",&SrcType[i],&TempSrc(i),&TempX(i),&TempY(i));
            sprintf(str,"%s %10.4lf %10.4lf %10.4lf",SrcType[i],TempSrc(i),TempX(i),TempY(i));
            puts(str);
        }
        FData.NSH = FData.NSF = 0;
        for (i = 0 ; i < FData.NS ; i++) {
            if (!strcmp(SrcType[i],"HEAT")) FData.NSH++;
            else if (!strcmp(SrcType[i],"FLUID")) FData.NSF++;
        }
        if (FData.NSH) {
           FData.HeatSrc.Initialize(FData.NSH);
           FData.HeatX.Initialize(FData.NSH);
           FData.HeatY.Initialize(FData.NSH);
        }
        if (FData.NSF) {
           FData.FluidSrc.Initialize(FData.NSF);
           FData.FluidX.Initialize(FData.NSF);
           FData.FluidY.Initialize(FData.NSF);
        }
        FData.NSH = FData.NSF = 0;
        for (i = 0 ; i < FData.NS ; i++) {
            if (!strcmp(SrcType[i],"HEAT")) {
               double s,x,y;
               FData.HeatSrc(FData.NSH) = TempSrc(i);
               FData.HeatX(FData.NSH) = TempX(i);
               FData.HeatY(FData.NSH) = TempY(i);
               FData.NSH++;
            }
            else if (!strcmp(SrcType[i],"FLUID")){
               FData.FluidSrc(FData.NSF) = TempSrc(i);
               FData.FluidX(FData.NSF) = TempX(i);
               FData.FluidY(FData.NSF) = TempY(i);
               FData.NSF++;
            }
        }
        delete[] SrcType;
     }

     puts("Cell Node    X          Y ");
     for (i= 0 ; i < FData.NI;i++) {
       fgets(str,255,fp);
       sscanf(str,"%lf %lf",&FData.Node(i,0),&FData.Node(i,1));
       sprintf(str,"%2d %10.4lf %10.4lf ",i,FData.Node(i,0),FData.Node(i,1));
       puts(str);
     }

     puts(" ");

     puts("Cell Elem    X       Y ");
     for ( i= 0 ; i< FData.NE;i++){
       fgets(str,255,fp);
       sscanf(str,"%d %d %d %d %lf %lf %lf",&FData.Elem(i,0),&FData.Elem(i,1),
                                         &FData.Elem(i,2),&FData.Elem(i,3),
                                         &FData.OriVolStrain(i),&FData.OriTemp(i),
                                         &FData.OriPress(i));
       FData.Elem(i,0)--;
       FData.Elem(i,1)--;
       FData.Elem(i,2)--;
       FData.Elem(i,3)--;

       sprintf(str,"%2d %5d %5d %5d %5d %10.3lf %10.3lf %10.3lf",i,FData.Elem(i,0),
                                                                FData.Elem(i,1),
                                                                FData.Elem(i,2),
                                                                FData.Elem(i,3),
                                                                FData.OriVolStrain(i),
                                                                FData.OriTemp(i),
                                                                FData.OriPress(i));
       puts(str);
     }

     for (i=0 ;i < FData.NE;i++ ) {
           FData.neXi(i) =0; FData.neYi(i) =0;
           for (k = 0;k< 4;k++) {
             FData.neXi(i) = FData.neXi(i)+FData.Node(FData.Elem(i,k),0)/4.0;
             FData.neYi(i) = FData.neYi(i)+FData.Node(FData.Elem(i,k),1)/4.0;
           }
           sprintf(str,"%2d %10.4lf %10.4lf",i,FData.neXi(i),FData.neYi(i));
           puts(str);
     }

     fclose(fp);
     return true;
}

TBEMInput &TBEMInput::operator=(TBEMInput &inp)
{
   strcpy(FileName,inp.FileName);
   FData = inp.FData;
   return *this;
}
//---------------------------------------------------------------------------
