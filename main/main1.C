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

  Generator *Gen = new Generator(X0,lambda,Xmax,Rmax);
  double t0 = DepthToHeight(Xmax)/c;

  vector<double> R,T;
  vector<double> xx,yy;
  for(int i=0;i<100000;i++)
  {
    double h = DepthToHeight(Gen->GenerateDepth());

    // xx.push_back((double)0.1*i );
    // yy.push_back(Gen->GaisserHillas(HeightToDepth(0.1*i)));

    vector<double> d = Gen->GenerateDirection();
    double t = h/(abs(d[1])*c);
    double r = c*t*d[0];
    t-=h/c;
    if(r<Rmax && t>0 )
    {
      R.push_back(r);
      T.push_back(1e9*(t));
    }


  }
  auto G = new TGraph(R.size(),R.data(),T.data());
  G->Fit("pol2");
  TF1* fit = G->GetFunction("pol2");
  cout<<fit->GetParameter(2)<<endl;
  auto C = new TCanvas("c","c",1600,900);
  G->GetXaxis()->SetRangeUser(-5,Rmax+5);
  G->GetYaxis()->SetRangeUser(-1,25);
  G->Draw("AP");

  C->SaveAs("Hmm.pdf");

  return 0;
}
