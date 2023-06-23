#include <iostream>
#include <vector>
#include "tools.h"
#include "Generator.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLatex.h"


using namespace std;
int main(int argc, char *argv[])
{
  double X0 = 0, lambda = 1;
  double D0 = 1030, H = 6500;
  double Xmax = 150, Rmax=300;
  double c = 2.998e8;
  auto DepthToHeight = [H,D0](double X)
  {
    return H*log(D0/X);
  };

  auto HeightToDepth = [H,D0](double h)
  {
    return D0*exp(-h/H);
  };


  remove("Lambda.gif");
  for(lambda;lambda<=150;lambda++)
  {
    Generator *Gen = new Generator(X0,lambda,Xmax,Rmax);
    vector<double> R,T;

    for(int i=0;i<100000;i++)
    {
      double h = DepthToHeight(Gen->GenerateDepth());
      while(h<0){h = DepthToHeight(Gen->GenerateDepth());};

      vector<double> d = Gen->GenerateDirection();
      double t = h/(abs(d[1])*c);
      double r = c*t*d[0];

      if(r<Rmax && t>0 )
      {
        R.push_back(r);
        T.push_back(1e9*(t-h/c));
      }


    }
    auto G = new TGraph(R.size(),R.data(),T.data());
    G->Fit("pol2");
    G->SetTitle("#Delta t as a function of depth for #lambda = [1,150]");
    auto C = new TCanvas("c","c",1600,900);
    G->GetXaxis()->SetRangeUser(-5,Rmax+5);
    G->GetYaxis()->SetRangeUser(-1,25);
    G->GetXaxis()->SetTitle("Depth [gcm^{-2}]");
    G->GetYaxis()->SetTitle("#Delta t [ns]");
    G->Draw("AP");

    C->SaveAs("Lambda.gif+10");
  }



  return 0;
}
