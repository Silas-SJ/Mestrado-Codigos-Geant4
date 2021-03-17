#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
int readspec(){
string text;
double energia[200];
double fluxo[200];
int i=0;
std::ifstream datafile("spec-neutron.dat",std::ifstream::in);
while(datafile >> energia[i] >> fluxo[i]){
 std::cout  << energia[i] << " " << fluxo[i] << endl;
 i++;
}
TF1 * f1 = new TF1("f1","[0]*exp([1]*log([2]*x))",0,4000);
//TF1 * f1 = new TF1("f1","[0]*pow(10,[1]*log([2]*x))",0,4000);
f1->SetParameter(0,3000);
f1->SetParameter(1,-1.0);
f1->SetParameter(2,1);
TGraph * gr = new TGraph(i,energia,fluxo);
gr->SetMarkerStyle(22);
gr->Draw("AP");
gr->Fit("f1");
return 0;
}
