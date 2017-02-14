#include <iostream>
//#include <math>
#include "TVector3.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"

 const Double_t R3 =-3.243410198;
 const Double_t Rn =-2.91646954;
 const Double_t DT = 0.000001;
 const Double_t TMax = 10.;
 const Double_t Y = 0.01; 	 // critical dressing parameter 
 const Double_t X = 1.189;	 // critical dressing parameter
 const Double_t B0 = 10.;	 // holding field
 const Double_t taubeta = 885.;
 const Double_t tau3 = 500.;
 const Double_t taucell = 1150.;

 const Int_t neutronnumber = 1680000;

using namespace std;


TVector3 cross(TVector3 S, TVector3 B){
  return S.Cross(B);
}

TVector3 BField(Double_t Time, Double_t Freq, Double_t B0, Double_t B1, Double_t Bm, Double_t Arg_m, Double_t Ba, TVector3 S){
 Double_t Arg;
 Arg = TMath::TwoPi()*Time*Freq;
 Double_t Bx,By,Bz;
 Bx = (B1+Bm*Arg_m)*TMath::Sin(Arg)+Ba*S.X();
 By = Ba*S.Y();
 Bz = B0+Ba*S.Z();
 
 TVector3 B;
 B.SetXYZ(Bx,By,Bz);
 return B;
}


TVector3 rk4(TVector3 S0, Double_t R, Double_t Time, Double_t Freq, Double_t B0, Double_t B1, Double_t Bm, Double_t Arg_m, Double_t Ba, TVector3 S){

 TVector3 S1, S2, S3;
 TVector3 k1, k2, k3, k4;
 TVector3 B;
 B  = BField(Time, Freq, B0, B1, Bm, Arg_m, Ba, S);
 k1 = DT*TMath::TwoPi()*R*cross(S0,B);
 S1 = S0 + 0.5*k1;
 B  = BField(Time+0.5*DT, Freq, B0, B1, Bm, Arg_m, Ba, S);
 k2 = DT*TMath::TwoPi()*R*cross(S1,B);
 S2 = S0 + 0.5*k2;
 k3 = DT*TMath::TwoPi()*R*cross(S2,B);
 S3 = S0 + k3;
 B  = BField(Time+DT, Freq, B0, B1, Bm, Arg_m, Ba, S);
 k4 = DT*TMath::TwoPi()*R*cross(S3,B);
 return S0 + (1./6.)*(k1+2.*k2+2.*k3+k4);
 }

int main(int argc, char** argv)
{
  if(argc != 7) {
    cout << "Usage: " << argv[0] << " [Em] [Fm] [ra] [re] [Alpha] [Beta]" << endl;
    return 1;
  }
  Double_t Em = atoi(argv[1]);
  Double_t Fm = atoi(argv[2]);
  Double_t ra = atoi(argv[3]);
  Double_t re = atoi(argv[4]);
  Double_t Alpha = atoi(argv[5]);
  Double_t Beta = atoi(argv[6]);

 gSystem->Load("libPhysics.so");

 Int_t bins  = 1000000;
 Int_t step  = 10;

 TObjArray Hlist(0);


 TH1D *CosThetan13 = new TH1D("CosThetan13","cos #theta_{n_{1},3}",bins,0.,TMax);
 TH1D *CosThetan23 = new TH1D("CosThetan23","cos #theta_{n_{2},3}",bins,0.,TMax);
 Hlist.Add(CosThetan13);
 Hlist.Add(CosThetan23);

//B field parameters-------------------------------------------------------
 Double_t Bd = X*B0/Y;  			//dressing field
 Double_t Bm = Bd*Em; 			//modulation amplitude
 Double_t Wa = 0.000942478*ra;
 Double_t We = 0.000001*re;
 Double_t Fd = Rn*B0/Y;				//dressing frequency
 Double_t Wm = TMath::TwoPi()*Fm;
 Double_t Ba = Wa/2./TMath::Pi()/Rn*ra;
 Double_t Be = We/2./TMath::Pi()/Rn*re;
 Double_t Bnull = 0.;				//zero field
 Double_t Bz = B0;

 Double_t BCA = 0.;
 Double_t BCB = 0.;
 Double_t BC0 = Bd;

 TVector3 B;
 Double_t T=0;
 Double_t Arg;
 Double_t Argm = 0.;
 Double_t Tm;

 if(Fm == 0){
   Tm = TMax;
 }else{
   Tm = (1./Fm);
 }


 Int_t BinMod;
 BinMod = Tm/DT; 
 Int_t CycleMod;
 CycleMod = TMax/BinMod/DT;

 cout << "How many time step in one period? " << BinMod<<endl;
 cout << "How many cycles during the measurement period? " << CycleMod<<endl;


 TH1F *s1;
 TH1F *s2;
 TH1F *p1;
 TH1F *p2;
 TH1F *p1mc;

 TH1F *sigmax1 = new TH1F("sigmax1","#Sigma 1",CycleMod,0,TMax);
 TH1F *sigmax2 = new TH1F("sigmax2","#Sigma 2",CycleMod,0,TMax);
 TH1F *sigma = new TH1F("sigma","#Delta #Sigma",CycleMod,0,TMax);

 TH1F *bd= new TH1F("bd","B_{d}",CycleMod,0,TMax);

 TH1F *phi1  = new TH1F("phi1","#Phi",CycleMod,0,TMax);
 TH1F *phi1mc = new TH1F("phi1mc","#Phi_{mc}",CycleMod,0,TMax);
 TH1F *phi2  = new TH1F("phi2","#Phi_{2}",CycleMod,0,TMax);
 TH1F *phi2mc = new TH1F("phi2mc","#Phi_{2,mc}",CycleMod,0,TMax);

 TH1F *sig1  = new TH1F("sig1","cos#theta_{n3}",CycleMod,0,TMax);
 TH1F *sig2  = new TH1F("sig2","cos#theta_{n3}",CycleMod,0,TMax);

 Hlist.Add(sigmax1);
 Hlist.Add(sigmax2);
 Hlist.Add(sigma);

 Hlist.Add(bd);

 Hlist.Add(phi1);
 Hlist.Add(phi1mc);
 Hlist.Add(phi2);
 Hlist.Add(phi2mc);

 Hlist.Add(sig1);
 Hlist.Add(sig2);

 TVector3 Sn1;
 TVector3 SN1;
 TVector3 S0;
 TVector3 Sn2;
 TVector3 SN2;
 S0.SetXYZ(0.,0.,0.);
 Sn1.SetXYZ(1.,0.,0.);
 SN1.SetXYZ(1.,0.,0.);
 Sn2.SetXYZ(1.,0.,0.);
 SN2.SetXYZ(1.,0.,0.);
 TVector3 Bn1,Bn2,BN1, BN2;
 Double_t costheta1 = 1.; // Angle between UCN and 3He
 Double_t costheta2 = 1.; // Angle between UCN and 3He
 
 Double_t y1 = 0.;
 Double_t y2 = 0.;
 Double_t y3 = 0.;
 Double_t y4 = 0.;

 Double_t f1 = 0.;
 Double_t f2 = 0.;
 Double_t f3 = 0.;
 Double_t f4 = 0.;
  
 Double_t norm0 = 0.;
 Double_t N0 = neutronnumber;
 Double_t gamma = 1./taubeta+1./taucell+1./tau3;

 Char_t name1[1000];
 Char_t name2[1000];
 Char_t name3[1000];
 Char_t name4[1000];
 Char_t name5[1000];
 Int_t j = 0;
 Int_t k = 0;
 for(Int_t l=0; l < CycleMod;l++){

//   bd->SetBinContent(l+1,Bd);

//   sprintf(name1,"s1%d",l);
//   sprintf(name2,"s2%d",l);
//   sprintf(name3,"p1%d",l);
//   sprintf(name4,"p2%d",l);
//   sprintf(name5,"p1mc%d",l);


   s1 = new TH1F(name1,"Cos_{1}",BinMod,0,Tm);
   s2 = new TH1F(name2,"Cos_{2}",BinMod,0,Tm);
   p1 = new TH1F(name3,"#Phi_{1}",BinMod,0,Tm);
   p2 = new TH1F(name4,"#Phi_{2}",BinMod,0,Tm);

   p1mc = new TH1F(name5,"#Phi_{1,mc}",BinMod,0,Tm);


//   Hlist.Add(s1);
//   Hlist.Add(s2);
//   Hlist.Add(p1);
//   Hlist.Add(p1mc);
//   Hlist.Add(p2);

   f3 = 0.;
   f4 = 0.;
   Double_t sum1 = 0.;
   Double_t sum2 = 0.;
   Double_t sum3 = 0.;
   Double_t sum4 = 0.;

   for(Int_t i=0; i < BinMod; i++){

     k = i + l*CycleMod;
     T = i*DT + l*Tm;

     Argm = Wm*T;
     y1 = costheta1;
     y2 = costheta2;
     sum3 = sum3 + y1;
     sum4 = sum4 + y2;

     s1->SetBinContent(i+1,y1); 
     s2->SetBinContent(i+1,y2); 
 
     y3 = y3 + y1*DT;
     y4 = y4 + y2*DT;
    
     f1 = exp(-gamma*T + 1./tau3*y3)*(1./taubeta+1./tau3*(1.-y1));
     f2 = exp(-gamma*T + 1./tau3*y4)*(1./taubeta+1./tau3*(1.-y2));
     sum1 = sum1 + f1;
     sum2 = sum2 + f1;

     p1->SetBinContent(i+1,f1); 
     p2->SetBinContent(i+1,f2); 

     f3 = f3 + f1*DT;
     f4 = f4 + f2*DT;

     if (TMath::Cos(Argm)>=0.){
       Arg = 1.;
       Bz = B0;
       BN1 = BField(T,Fd, Bz, Bd, Bm, Arg, Bnull, S0);
       SN1 = rk4(SN1,R3,T,Fd,Bz,Bd, Bm, Arg, Bnull, S0);
       SN1 = SN1.Unit();
       Bz = B0 + Be;
       Bn1 = BField(T,Fd,Bz, Bd, Bm, Arg, Ba, SN1);
       Sn1 = rk4(Sn1,Rn,T,Fd,Bz,Bd, Bm, Arg, Ba, SN1);
       Sn1 = Sn1.Unit();
       costheta1 = Sn1.Dot(SN1);

       Bz = B0; 
       BN2 = BField(T,Fd, Bz, Bd, Bm, Arg, Bnull, S0);
       SN2 = rk4(SN2,R3,T,Fd,B0,Bd, Bm, Arg, Bnull, S0);
       SN2 = SN2.Unit();
       Bz = B0 - Be;
       Bn1 = BField(T,Fd,Bz, Bd, Bm, Arg, Ba, SN1);
       Sn2 = rk4(Sn2,Rn,T,Fd,Bz,Bd, Bm, Arg, Ba, SN2);
       Sn2 = Sn2.Unit();
       costheta2 = Sn2.Dot(SN2);
     } else
     {
       Arg = -1.; 
       Bz = B0;
       BN1 = BField(T,Fd, Bz, Bd, Bm, Arg, Bnull, S0);
       SN1 = rk4(SN1,R3,T,Fd,Bz,Bd, Bm, Arg, Bnull, S0);
       SN1 = SN1.Unit();
       Bz = B0 + Be;
       Bn1 = BField(T,Fd,Bz, Bd, Bm, Arg, Ba, SN1);
       Sn1 = rk4(Sn1,Rn,T,Fd,Bz,Bd, Bm, Arg, Ba, SN1);
       Sn1 = Sn1.Unit();
       costheta1 = Sn1.Dot(SN1);

       Bz = B0; 
       BN2 = BField(T,Fd, Bz, Bd, Bm, Argm, Bnull, S0);
       SN2 = rk4(SN2,R3,T,Fd,B0,Bd, Bm, Arg, Bnull, S0);
       SN2 = SN2.Unit();
       Bz = B0 - Be;
       Bn1 = BField(T,Fd,Bz, Bd, Bm, Arg, Ba, SN1);
       Sn2 = rk4(Sn2,Rn,T,Fd,Bz,Bd, Bm, Arg, Ba, SN2);
       Sn2 = Sn2.Unit();
       costheta2 = Sn2.Dot(SN2);
     }

    if(k%step == 1){
      j = j + 1;
      CosThetan13->SetBinContent(j,costheta1);
      CosThetan23->SetBinContent(j,costheta2);
    }
   }


   sum1 = sum1/BinMod;
   sum2 = sum2/BinMod;
   sum3 = sum3/BinMod;
   sum4 = sum4/BinMod;


   phi1->SetBinContent(l+1,sum1);
   phi2->SetBinContent(l+1,sum2);
   sig1->SetBinContent(l+1,sum3);
   sig2->SetBinContent(l+1,sum4);

   norm0 = N0*f3;

   for(Int_t j=0;j< norm0 ;j++){
     Double_t r1 = p1->GetRandom();
     p1mc->Fill(r1);
   }

   Double_t x5 = 0.;
   Double_t x6 = 0.;
   Double_t SigmaX1 = 0.;
   Double_t SigmaX2 = 0.;
   Double_t Sigma=0.;

   for(Int_t i=0; i < BinMod;i++){

     if ( i < BinMod/2){
       x5 = p1mc->GetBinContent(i+1);
       SigmaX1 = SigmaX1 + x5;
     } else
     {
       x6 = p1mc->GetBinContent(i+1);
       SigmaX2 = SigmaX2 + x6;
     }
   }


   Sigma = SigmaX1-SigmaX2;
   sigmax1->SetBinContent(l+1,SigmaX1);
   sigmax2->SetBinContent(l+1,SigmaX2);
   sigma->SetBinContent(l+1,Sigma);

   BCA = BC0 - Alpha*Sigma;
   BCB = -Beta*Sigma;
   Bd  = BCA+BCB;
   BC0 = BCA;

   SigmaX1 = 0.;
   SigmaX2 = 0.;
   Sigma = 0.;

    delete s1;
    delete s2;
    delete p1;
    delete p2;
    delete p1mc;

 }


 for(Int_t k=0; k< neutronnumber; k++){
   Double_t r1 = phi1->GetRandom();
   phi1mc->Fill(r1);
   Double_t r2 = phi2->GetRandom();
   phi2mc->Fill(r2);
 }


 TFile *f = new TFile("bloch.root","RECREATE");
 Hlist.Write();
 f->Close();
 return 0;
}
