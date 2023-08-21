#include "TMath.h"
#include "TF1.h"

#include "Math/IFunction.h"
#include "TSystem.h"
#include <cmath>
#include "TAxis.h"
#include "TPaveLabel.h"
#include <Riostream.h>


double P_NB(double N,double Nbar,double k)
{
 double z1 = (tgamma(N+k));
 double z2 = (tgamma(N+1)*tgamma(k));
 double z3 = pow((Nbar/k),N);
 double z4 = pow((1+(Nbar/k)),(N+k));
 
 return (z1*z3)/(z2*z4);
}

const int np=100;

void ep_t9_NBfit()
{

  TCanvas *c1 = new TCanvas();
  double n;
  double Pn,err_Pn,err_n;
  double PN[100],N[100],err_PN[100],err_N[100];
  
  fstream file("/home/soumya/Tsallis_Data/Pyhton_Code/ep_data/ep-H1-data/ep-H1-data/ep_py_data/t9-epH1_py.txt",ios::in);
  
  //Initialise all the elements of array to 0

 for(int i=0;i<100;i++)
  {
   PN[i] = 0;
   N[i] = 0;
   err_N[i]=0;
   err_PN[i] = 0;
  }
  
  int c=0;
  double sum =0;
  double avg = 0;
  
  while(!file.eof())
  {
   file>>n>>Pn>>err_n>>err_Pn;
   N[c] = n;
   PN[c] = Pn;
   err_N[c] = err_n;
   err_PN[c] = err_Pn;
   sum = sum+Pn;
   avg = avg+n*Pn;
   c++;
    if(n==0 && Pn==0) break;
  }
     
   TGraphErrors *gr = new TGraphErrors(c,N,PN,err_N,err_PN);
   gr->SetMarkerColor(4);
   gr->SetMarkerSize(1.5);
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlack);
   
   //For error bars
   gr->SetLineWidth(2);
   gStyle->SetEndErrorSize(7);
   
   gr->SetTitle(0);
   gr->GetXaxis()->SetTitle("n");
   gr->GetYaxis()->SetTitle("P(n)");
   gr->GetXaxis()->CenterTitle(true);
   gr->GetYaxis()->CenterTitle(true);
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->Draw("AP");
   gr->GetYaxis()->SetRangeUser(-0.02,0.3);

   //Drawing the function defined above with scaling to normalise it

   TF1 *func = new TF1("func","[2]*P_NB(x,[0],[1])",0,40);
   func->SetParameter(0,5);
   func->SetParameter(1,5);
   func->SetParameter(2,20);
   func->SetLineWidth(2);
   func->SetLineColor(kRed);
   
   //Fitting the Graph
   gr->Fit("func");
   func->Draw("SAME");
   cout<<"Chi2: "<<func->GetChisquare()<<endl;
   cout<<"N parameter: "<<func->GetParameter(0)<<" +- "<<func->GetParError(0)<<endl;
   cout<<"k parameter: "<<func->GetParameter(1)<<" +- "<<func->GetParError(1)<<endl;
   cout<<"Normalisation : "<<func->GetParameter(2)<<" +- "<<func->GetParError(2)<<endl;
   //cout<<"Degrees of freedom: "<<func->GetNdf()<<endl;
   
   TLegend *legend = new TLegend(0.6,0.4,0.8,0.7,NULL,"NDC");
   legend->SetTextSize(0.05);
   legend->SetFillStyle(10);
   legend->SetFillColor(10);
   legend->SetBorderSize(1);
   legend->SetLineColor(kWhite);	
   legend->AddEntry(func,"NB fitting function","l");
   legend->AddEntry(gr,"ep Data","ep");
   legend->Draw();
   
   gStyle->SetOptFit(1111);
   
   c1->SaveAs("/home/soumya/Tsallis_Data/Plots/NB_fit_Plots/ep_t9_NBfit.pdf");
   
}


