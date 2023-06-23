#include <iostream>
#include <vector>
#include "tools.h"
#include "Generator.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"


using namespace std;

int main(int argc, char *argv[])
{
  double X0 = 0, lambda = 6;
  double D0 = 1030, H = 6500;
  double Xmax = 250, Rmax=300;
  double c = 2.998e8;

  auto DepthToHeight = [H,D0](double X)
  {
    return H*log(D0/X);
  };

  auto HeightToDepth = [H,D0](double h)
  {
    return D0*exp(-h/H);
  };

  int N2 = 250;
  vector<double> R,T;
  //vector<double> xx,yy;
  std::vector<double> XX(N2,0);
  std::vector<double> curvature(N2,0);

  for(int j=0; j<N2; j++)
  {
    //h = (-1)*H*log(X_max[j-1]/X0);
    Generator *Gen = new Generator(X0,lambda,Xmax,Rmax);

  for(int i=0;i<100000;i++)
  {
    double h = DepthToHeight(Gen->GenerateDepth());

    //xx.push_back((double)0.1*i );
    //yy.push_back(Gen->GaisserHillas(HeightToDepth(0.1*i)));

    vector<double> d = Gen->GenerateDirection();
    double t = h/(abs(d[1])*c);
    double r = c*t*d[0];

    if(r<Rmax)
    {
      R.push_back(r);
      T.push_back(1e9*(t-h/c));
    }

    
  } 

    auto g = new TGraph(R.size(), R.data(), T.data());
    g->Fit("pol2");

    TF1* fit = g->GetFunction("pol2");

    curvature[j] = fit->GetParameter(2) * 1e-9;

    XX[j] = Xmax;
    Xmax++;

  }

  auto G = new TGraph(R.size(),R.data(),T.data());
  //auto G1 = new TGraph(xx.size(),xx.data(),yy.data());
  auto C = new TCanvas("c","c",1600,900);
  G->GetXaxis()->SetRangeUser(-5,Rmax+5);
  G->Draw("AP");

  C->SaveAs("Hmm.pdf");

  //Draw the shower plane curvature c as a function of X_max
  auto c2 = new TCanvas("c2"," c(X_max)",200,10,1500,1500);
  c2->cd();

  auto g2 = new TGraph(N2, XX.data(), curvature.data());
  g2->GetXaxis()->SetTitle("X_max[g/cm^2]");
  g2->GetYaxis()->SetTitle("c[s/m^2]");
  
  g2->Fit("pol1");
  g2->Draw("AP");

  c2->SaveAs("cXmax.pdf");

  return 0;
}
