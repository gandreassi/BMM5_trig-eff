#ifndef fitter_H
#define fitter_H

#include <iostream>
#include <cmath>
#include <string>
#include "RooAbsData.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooStats/SPlot.h"
#include "RooTreeDataStore.h"
#include "TFile.h"
#include "TAxis.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"
#include "TTreeFormula.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TText.h"
#include "TLine.h"
#include "TString.h"
#include "TApplication.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TChain.h"
#include "TObject.h"
#include "TParticle.h"
#include "TClass.h"
#include "TSystem.h"
#include "TEfficiency.h"

//  RooFit
#include "RooFit.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooBinning.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "TStyle.h"
#include "RooMinimizer.h"
#include "RooCustomizer.h"

class fitter {
	public:
		
		fitter();//constructor

		void makeDataSet(TChain*, float, float, std::string, std::string);
		void reduceDataSet(std::string, float, float);
		void resetDataSet();
		void fit();
		void preparePDF(bool);
		void saveFitPdf(std::string);
		float getSignalYield();
		float getSignalYieldError();

	private:
		RooWorkspace w;
		RooDataSet *data;
		RooDataSet *data_full;
		RooDataHist *binned_data;
		float M_min_def;
		float M_max_def;
		RooFitResult * result;
};

#endif