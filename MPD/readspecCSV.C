#include <iostream>
#include <string>
#include <stdio.h>
int readspecCSV(){
string text;
double energia[200];
double fluxo[200];
int i=0;
FILE * datafile;
datafile=fopen("spec-pion.csv","r");
while(fscanf(datafile,"%lf, %lf",&energia[i],&fluxo[i])>0){
 std::cout  << energia[i] << " " << fluxo[i] << endl;
 i++;
}
//TF1 * f1 = new TF1("f1","[0]*exp([1]*log10([2]*x))",0,4000);
TF1 * f1 = new TF1("f1","[0]*pow(10,[1]*log10([2]*x))",0,4000);
f1->SetParameter(0,3000);
f1->SetParameter(1,-1.0);
f1->SetParameter(2,1);
TGraph * gr = new TGraph(i,energia,fluxo);
gr->SetMarkerStyle(22);
gr->Draw("AP");
gr->Fit("f1");
return 0;
}
