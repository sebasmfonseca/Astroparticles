#include <iostream>
#include <vector>
#include "tools.h"
#include "Generator.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include "TLatex.h"


using namespace std;
int main(int argc, char *argv[])
{
  double X0 = 0, lambda = 30;
  double D0 = 1030, H = 6500, Rmax=300;
  double c = 2.998e8;

  auto DepthToHeight = [H,D0](double X)
  {
    return H*log(D0/X);
  };

  auto HeightToDepth = [H,D0](double h)
  {
    return D0*exp(-h/H);
  };


  remove("Boom.gif");

  vector<double> VV,L;

  for(lambda;lambda<=100;lambda+=5)
  {
    vector<double> C,XX;
    cout<<lambda<<endl;
    for(double Xmax=250;Xmax<500;Xmax+=5)
    {
      cout<<Xmax<<endl;
      Generator *Gen = new Generator(X0,lambda,Xmax,Rmax);
      vector<double> R,T;
      for(int i=0;i<100000;i++)
      {
        double h = DepthToHeight(Gen->GenerateDepth());
        while(h<0){h = DepthToHeight(Gen->GenerateDepth());};
        vector<double> d = Gen->GenerateDirection();
        double t = h/(abs(d[1])*c);
        double r = c*t*d[0];

        if(r<Rmax)
        {
          R.push_back(r);
          T.push_back(1e9*(t-h/c));
        }
      }
      auto G = new TGraph(R.size(),R.data(),T.data());
      G->Fit("pol2");
      TF1* fit = G->GetFunction("pol2");
      XX.push_back(Xmax);

      C.push_back(fit->GetParameter(2));
    }


    auto G = new TGraph(XX.size(),XX.data(),C.data());
    G->Fit("pol1");
    cout<<"FODASSE"<<endl;
    TF1* fit = G->GetFunction("pol1");
    VV.push_back(fit->GetParameter(1));
    L.push_back(lambda);
    auto cv = new TCanvas("c","c",1600,900);
    G->Draw("AP");

  }

  auto G = new TGraph(L.size(),L.data(),VV.data());
  auto cv = new TCanvas("c","c",1600,900);
  G->SetMarkerStyle(8);
  G->Draw("AP");
  cv->SaveAs("Beta.pdf");





  return 0;
}
