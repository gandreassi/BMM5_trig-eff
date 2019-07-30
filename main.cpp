#include "fitter.cpp"

int main (int argc, char *argv[]) {
  
  if (argc > 1){
    return 0;
  }


	TChain* chain = new TChain("Events");
	chain->Add("/mnt/hadoop/scratch/gandreas/Charmonium_Dimuon0_Jpsi_NoVertexing/*.root");

	fitter f;
	f.makeDataSet(chain);
  	f.fit();
  	f.saveFitPdf("plots/fit.pdf");
  	cout<<f.getSignalYield()<<" +/- "<<f.getSignalYieldError()<<endl;

  return 0;
}
