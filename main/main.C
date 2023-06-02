#include <iostream>
#include <vector>
#include "tools.h"
#include "Particle.h"
#include "Generator.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"

int main() {

  double photon_p = 1.; //photon momentum
  double theta = 0.; //photon angle
  const double c = 3*pow(10,8);
  double d = 0.; //distance travelled by the photon
  double h = 0.;

  int N1 = 3000;
  int N2 = 500;
  double X0 = 1030.; //(g cm^-2)
  double H = 6500; //(m)

  //std::vector<double> photon_pos {0,0,X_max};
  //Particle photon(22, photon_p, tools::Get_Dir_Spherical(1, 0, theta), photon_pos);
  //double v = photon.GetVelocity()*3*pow(10,8);

  double v = c;

  std::vector<double> r(N1,0);
  std::vector<double> Delta_t(N1,0);
  std::vector<double> X_max(N2,0);
  std::vector<double> curvature(N2,0);

  for(int j = 1; j < N2+1; j++){

    X_max[j-1] = j;
    h = (-1)*H*log(X_max[j-1]/X0);
    //std::cout << X_max[j] << '\n';
    
    double t0 = h/v;
    
    for(int i=0; i < N1; i++){
      r[i] = i*0.1;
      theta = atan(r[i]/h);
      d = h/cos(theta);
      Delta_t[i] = (d/v - t0)*pow(10,9);
    }

    auto g = new TGraph(N1, r.data(), Delta_t.data());
    g->Fit("pol2");

    TF1* fit = g->GetFunction("pol2");

    curvature[j] = fit->GetParameter(2);

    //std::cout << curv[j] << '\n';
    
    //delete fit;
    //delete g;
  }

  //std::cout << "t0: " << t << '\n';

  //TApplication app("app", nullptr, nullptr);
  

  //Draw Delta_t as a function of r
  auto c1 = new TCanvas("c1","Delta_t (r)",200,10,1500,1500);

  auto g1 = new TGraph(N1, r.data(), Delta_t.data());
  g1->GetXaxis()->SetTitle("r[m]");
  g1->GetYaxis()->SetTitle("t[ns]");
  g1->Draw("AC");

  c1->SaveAs("Deltatr.pdf");

  //Draw the shower plane curvature c as a function of X_max
  auto c2 = new TCanvas("c2"," c(X_max)",200,10,1500,1500);
  c2->cd();

  auto g2 = new TGraph(N2, X_max.data(), curvature.data());
  g2->GetXaxis()->SetTitle("X_max[g/cm^2]");
  g2->GetYaxis()->SetTitle("c[s/m^2]");

  g2->Fit("pol1");
  g2->Draw("AP");

  c2->SaveAs("cXmax.pdf");

  //app.Run();

  //app.Terminate();

  return 0;
}
