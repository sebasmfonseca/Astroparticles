#include <iostream>
#include <vector>
#include "tools.h"
#include "Particle.h"
#include "Generator.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;
int main(int argc, char *argv[])
{

 Generator* Gen = new Generator();
 double dt = 1e-9, XMax = 5000;
 Particle Primary({0,0,XMax},{0,0,-1});
 Primary.PropagateSimple(dt);


 double t0 = Primary.GetTime();

/////////////////////////
vector<double> t1;
vector<double> R1;
 for(int i=1;i<1000;i++)
 {
   //cout<<i<<endl;
   Particle aux({0,0,XMax},Gen->UniformDirection());
   if(aux.GetDirection()[2]<-.17)
   {
     aux.PropagateSimple(dt);
     t1.push_back(aux.GetTime()-t0);
     R1.push_back(aux.GetRadius());
     //cout<<aux.GetTime()-t[0]<<endl;
   }


 }
 auto c = new TCanvas("c","c",1200,900);
 auto G1 = new TGraph(t1.size(),R1.data(),t1.data());
 G1->GetXaxis()->SetRangeUser(-10,1.1*XMax/.17);
 G1->GetYaxis()->SetRangeUser(-1e-5,1e-4);
 G1->Draw("AP");
 G1->SetMarkerStyle(20);
 G1->SetMarkerColor(kRed);
 c->SaveAs("Simple.pdf");
/////////////////////////

 double stddev = 500;
 vector<double> t2;
 vector<double> R2;
  for(int i=1;i<1000;i++)
  {
    //cout<<i<<endl;
    Particle aux({0,0,Gen->NormalHeight(XMax,stddev)},Gen->UniformDirection());
    if(aux.GetDirection()[2]<-.17)
    {
      aux.PropagateSimple(dt);
      t2.push_back(aux.GetTime()-t0);
      R2.push_back(aux.GetRadius());
      //cout<<aux.GetTime()-t[0]<<endl;
    }


  }
  auto G2 = new TGraph(t2.size(),R2.data(),t2.data());
  G2->GetXaxis()->SetRangeUser(-10,1.1*XMax/.17);
  G2->GetYaxis()->SetRangeUser(-1e-5,1e-4);
  G2->Draw("AP");
  G2->SetMarkerStyle(20);
  G2->SetMarkerColor(kRed);
  c->SaveAs("NormalDist.pdf");



 return 0;
}
