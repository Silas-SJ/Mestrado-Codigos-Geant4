#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
int readspecPion(){
string text;
double energia[200];
double fluxo[200];
int i=0;
std::ifstream datafile("spec-pion.dat",std::ifstream::in);
while(datafile >> energia[i] >> fluxo[i]){
 std::cout  << energia[i] << " " << fluxo[i] << endl;
 i++;
}
TCanvas * c1 = new TCanvas("Spectrum");
    c1->SetLogx();
    c1->SetLogy();
//TF1 * f1 = new TF1("f1","[0]*exp([1]*log([2]*x))",0,4000);
//TF1 * f1 = new TF1("f1","[0]*pow(10,[1]*log([2]*x))",0,4000);
TF1 * f1 = new TF1("f1","[0]*pow(x,[1])",100,2000);
f1->SetParameter(0,3e-5);
f1->SetParameter(1,-1.34);
//f1->SetParameter(2,1);
TGraph * gr = new TGraph(i,energia,fluxo);
gr->SetMarkerStyle(22);
gr->Draw("AP");
gr->Fit("f1");
    
double energy=0.0;
for (int i=0;i<50;i++){
  energy=energy+0.02;
  cout << "energy: " << energy << " flux [/cm2 sr s]: " << f1->Eval(energy) << endl;
}
cout << "Integral 0.01 -> 2 GeV:  " << f1->Integral(100,2000) << endl;
cout << "Integral x ano:  " << f1->Integral(100,2000)*3600*24*365 << endl;

return 0;
}
