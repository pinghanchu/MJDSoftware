// 2016.11.21 Pinghan Chu
// calibrate sum spectrum of all high gain, all natural detectors etc.
#include "GATAutoCal.hh"
#include "TMatrixD.h"
#include "TFile.h"

TGraphErrors* GetGraph(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr){
  const Int_t nPy = Py.size();

  Double_t A[nPy];
  Double_t AErr[nPy];
  Double_t B[nPy];
  Double_t BErr[nPy];

  for(Int_t j=0;j<nPy;j++){
    A[j]=Px.at(j);
    AErr[j] = PxErr.at(j);
    B[j]=Py.at(j);
    BErr[j] =PyErr.at(j);
  }
  TGraphErrors* fGraph = new TGraphErrors(nPy,A,B,AErr,BErr);
  fGraph->SetFillStyle(0);
  fGraph->SetFillColor(0);

  return fGraph;
}


int main(int argc, char** argv)
{
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [start run] [endrun] [energy name]" << endl;
    return 1;
  }
  
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  const char* EnergyName = Form("%s",fEName.c_str());


  ofstream fout(Form("./List/cal/calibrationlog_%d_%d.txt",fStartRun,fEndRun));
  ofstream fcal(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream fcov(Form("./List/cal/cov_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream freso(Form("./List/cal/reso_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream froi(Form("./List/cal/roi_%d_%d.txt",fStartRun,fEndRun),ios::app);
 
  fcal.precision(15);
  fcov.precision(15);
  freso.precision(15);
  froi.precision(15);
  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
 
  vector<string> DataName;
  DataName.push_back("00"); // Natural HG
  // DataName.push_back("01"); // Natural LG
  //DataName.push_back("10"); // Enriched HG
  //DataName.push_back("11"); // Enriched LG
  //DataName.push_back("HG");
  //DataName.push_back("LG");
  vector<string> TitleName;
  TitleName.push_back("Natural HG");
  //TitleName.push_back("Natural LG");
  //TitleName.push_back("Enriched HG");
  //TitleName.push_back("Enrichd LG");
  //TitleName.push_back("High Gain");
  //TitleName.push_back("Low Gain");

  const Int_t nDataName = DataName.size();

  TH1D *Energy[nDataName];
  Double_t Low, Up;
  Int_t rebin = 1;
  string FileName;
  vector<Double_t> EnergyROI;
  EnergyROI.push_back(60);
  EnergyROI.push_back(2039);
  EnergyROI.push_back(2614.5);
  Double_t LinearROI = 0;
  Double_t LinearROIErr = 0;
  Double_t ResoROI = 0;
  Double_t ResoROIErr = 0;
  Double_t fScale = 1;
  Double_t fOffset = 0;
  Double_t fScaleErr = 0;
  Double_t fOffsetErr = 0;
  Double_t P0=0;
  Double_t P1=0;
  Double_t P2=0;
  Double_t P0Err=0;
  Double_t P1Err=0;
  Double_t P2Err=0;
  Double_t PLinear[2];
  Double_t PReso[3];
  Double_t errsum = 0;
  Int_t Index = 0;
  vector<string> ParName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> Cov;
  vector<Double_t> cov;

  TFile fhist(Form("./Hist/hist_%d_%d.root",fStartRun,fEndRun),"read");
  ofstream fpar(Form("./Hist/par_%d_%d.txt",fStartRun,fEndRun));

  for(Int_t i = 0;i<nDataName;i++){
    fout << "=================="<< endl;
    fout << Form("%s",TitleName.at(i).c_str()) << endl;
    fout << "=================="<< endl;
 
    Energy[i] = (TH1D*)fhist.Get(Form("%s%s",fEName.c_str(),DataName.at(i).c_str()));
    FileName = Form("%s_%d_%d_%s",fEName.c_str(),fStartRun,fEndRun,DataName.at(i).c_str());
    ParName.clear();
    Par.clear();
    ParErr.clear();
    Cov.clear();
    Low = 2500;
    Up = 2700;
    rebin = 20;
    TH1F *h1 = (TH1F*)Energy[i]->Clone();
    TH1F *h2 = (TH1F*)Energy[i]->Clone();
    h1->Rebin(rebin);
    h2->Rebin(rebin);
    h1->GetXaxis()->SetRangeUser(Low,Up);
    Int_t maxbin = h1->GetMaximumBin();
    Double_t scaleenergy = h1->GetBinCenter(maxbin)/2614.5333;
    h1->GetXaxis()->SetRangeUser(2500*scaleenergy,2700*scaleenergy);
    Double_t scaleamps = h1->Integral()*h1->GetBinWidth(maxbin)/35000;
    fout << Low << " " << Up << " " << maxbin << " " <<  scaleenergy << " "<< scaleamps << endl;
    h1->GetXaxis()->SetRange();
    
    ds.MultiPeakFitterFull(h2,scaleenergy, scaleamps, FileName,&ParName,&Par,&ParErr, &Cov);

    for(size_t ip = 0;ip<ParName.size();ip++){
      fout << ip << " " << ParName.at(ip).c_str() << " " << Par.at(ip) << " " << ParErr.at(ip) << endl;
    }

    
    fOffset = -Par.at(25)/Par.at(26);
    fScale = 1./Par.at(26);
    fOffsetErr = ParErr.at(25)/Par.at(26);
    fScaleErr = ParErr.at(26)/pow(Par.at(26),2);
    P0 = Par.at(27);
    P1 = Par.at(28);
    P2 = Par.at(29);
    P0Err = ParErr.at(27);
    P1Err = ParErr.at(28);
    P2Err = ParErr.at(29);
   
    fcal << "000" << " " <<  DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	 << fOffset << " " << fOffsetErr << " " << fScale << " " << fScaleErr << endl;
    freso << "000" << " " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	  << P0 << " " << P0Err << " " 
	  << P1 << " " << P1Err << " "		\
	  << P2 << " " << P2Err <<endl;
    TMatrixD mat(2,2);
    for(Int_t i1 = 0;i1<2;i1++){
      for(Int_t i2 = 0;i2<2;i2++){
	Index = (i1+25)*65+(i2+25);
	mat(i1,i2) = Cov.at(Index);
	cov.push_back(Cov.at(Index));
      }
    }
    fcov << "000" << " " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	 << cov.at(0) << " " << cov.at(1) << " " << cov.at(2) << " " << cov.at(3) << endl;

    for(size_t ir = 0;ir<EnergyROI.size();ir++){
      Double_t energyROI = EnergyROI.at(ir);
      LinearROI = (energyROI-fOffset)/fScale;
      PLinear[0] = -1/fScale;
      PLinear[1] = -(energyROI-fOffset)/(fScale*fScale);
      Index = 0;
      errsum = 0;
      for(Int_t i1 = 0;i1<2;i1++){
	for(Int_t i2 = 0;i2<2;i2++){
	  Index = (i1+25)*65+(i2+25);
	  fout << i1 << " " << i2 << " " << Index << " " << Cov.at(Index) << " " << PLinear[i1] << " " << PLinear[i2] << endl;
	  if(i1==i2){
	    errsum = errsum + PLinear[i1]*TMath::Abs(Cov.at(Index))*PLinear[i2];
	  }else{
	    errsum = errsum + PLinear[i1]*Cov.at(Index)*PLinear[i2];
	  }
	}
      }

      LinearROIErr = TMath::Sqrt(errsum);
      
      ResoROI = TMath::Sqrt(P0*P0+P1*P1*energyROI+P2*P2*energyROI*energyROI);
      PReso[0] = 0.5/ResoROI*2*P0;
      PReso[1] = 0.5/ResoROI*2*P1*energyROI;
      PReso[2] = 0.5/ResoROI*2*P2*energyROI*energyROI;
      errsum = 0;
      for(Int_t i1 = 0;i1<3;i1++){
	for(Int_t i2 = 0;i2<3;i2++){
	  Index = (i1+27)*65+(i2+27);
	  fout << i1 << " " << i2 << " " << Index << " " << Cov.at(Index) << " " << PReso[i1] << " " << PReso[i2] << endl;
	  if(i1==i2){
	    errsum = errsum + PReso[i1]*TMath::Abs(Cov.at(Index))*PReso[i2];
	  }else{
	    errsum = errsum + PReso[i1]*Cov.at(Index)*PReso[i2];
	  }
	}
      }
      ResoROIErr = TMath::Sqrt(errsum);
      
      froi << "000 " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " 
	   << fEName.c_str() << " " << energyROI << " " 
	   << LinearROI << " " << LinearROIErr << " " << ResoROI << " " << ResoROIErr << endl;
    }
      
    vector<Double_t> Px;
    vector<Double_t> PxErr;
    vector<Double_t> Py;
    vector<Double_t> PyErr;
    
    for(size_t ii = 65;ii< 90;ii++){
      Int_t ij = ii+25;
      Px.push_back(Par.at(ii));
      PxErr.push_back(0);
      Py.push_back(Par.at(ij)-Par.at(ii));
      PyErr.push_back(ParErr.at(ij));
      fout <<  "000 " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " "
	   << fEName.c_str() << " " 
	   << Par.at(ii) << " 0 " << Par.at(ij) << " " << ParErr.at(ij) << endl;
    }

    TCanvas *c1 = new TCanvas("c1");
    TGraphErrors *gr = GetGraph(Px,PxErr,Py,PyErr);
    gr->SetTitle(";Energy (keV); #Delta E (keV)");
    gr->GetYaxis()->SetTitleOffset(1.3);
    gr->SetMarkerStyle(2);
    gr->SetMarkerSize(3);
    gr->SetLineColor(1);
    gr->SetMarkerColor(1);
    gr->Draw("AP");
    c1->Print(Form("./Plot/energydelta_%s_%d_%d_%s.pdf",fEName.c_str(),fStartRun,fEndRun,DataName.at(i).c_str()));

    delete h1;
    delete h2;
    delete gr;
    delete c1;    
  }
}


