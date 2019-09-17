import ROOT as r
from uncertainties import ufloat
from sys import argv
from math import sqrt

year = argv[1]

lumi={	"2017" : {	"HLT_Dimuon0_Jpsi" : 75350805.341,
					"HLT_Dimuon0_Jpsi_NoVertexing" : 55714920.339,
					#"HLT_DoubleMu4_3_Jpsi_Displaced" : 14479239278.717,
					"HLT_DoubleMu4_Jpsi_Displaced" : 876993819.962,
					"HLT_DoubleMu4_Jpsi_NoVertexing" : 234517906.885},
		"2018" : {	"HLT_Dimuon0_Jpsi" : 2692017.753,
					"HLT_Dimuon0_Jpsi_NoVertexing" : 4036006.636,
					"HLT_DoubleMu4_3_Jpsi" : 1173487158.329,
					#"HLT_DoubleMu4_Jpsi_Displaced" : 212306014.543,
					"HLT_DoubleMu4_Jpsi_NoVertexing" : 52366197.343}}

chain = r.TChain("Events")
chain.Add("/eos/user/g/gandreas/jpsikmc/"+year+"/*.root")

f_map = r.TFile.Open("hists"+year+".root")
h_pass = f_map.Get("hpass")
h_pass.Sumw2()
h_pass.Scale(1./lumi[year]["HLT_DoubleMu4_Jpsi_Displaced"])
h_tot = f_map.Get("htot")
h_tot.Sumw2()
h_tot.Scale(1./lumi[year]["HLT_DoubleMu4_Jpsi_NoVertexing"])

eff_data = h_pass.Clone()
eff_data.Divide(h_tot)

binningX = eff_data.GetXaxis().GetXbins()
binningY = eff_data.GetYaxis().GetXbins()
eff_MC = eff_data.Clone("simulation")
eff_MC = r.TEfficiency("eff_MC","simulation;vertex probability;#epsilon", eff_data.GetNbinsX(), binningX.GetArray(), eff_data.GetNbinsY(), binningY.GetArray())
eff_MC.SetStatisticOption(r.TEfficiency.kBBayesian)
eff_MC.SetUseWeightedEvents()

#ge tpileup histogram for data and MC
data_PU_file = r.TFile.Open("DataPileup_"+year+".root")
data_PU_h = data_PU_file.Get("pileup")
data_PU_h.Sumw2()
data_PU_h.Scale(1./data_PU_h.Integral())
MC_PU_h = data_PU_h.Clone("MC_PU")
MC_PU_h.Reset()
chain.Draw("Pileup_nTrueInt>>MC_PU", "", "goff")
MC_PU_h.Sumw2()
MC_PU_h.Scale(1./MC_PU_h.Integral())
PU_weights = data_PU_h.Clone("PU_w")
PU_weights.Divide(MC_PU_h)

n_pass = ufloat(0.,0.)
n_tot = ufloat(0.,0.)

n_pass_MC = 0.

for event in chain:
	if len(event.mm_kin_vtx_prob)>0:
		if  abs(event.mm_gen_pdgId[0])==443 \
			and abs(event.mm_gen_mu1_pdgId[0])==13 and abs(event.mm_gen_mu2_pdgId[0])==13\
			and event.mm_gen_mu1_pdgId[0]*event.mm_gen_mu2_pdgId[0]<0:
			#and (event.L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 or event.L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4)\

			PU_weight = ufloat(PU_weights.GetBinContent(PU_weights.FindBin(chain.Pileup_nTrueInt)), PU_weights.GetBinError(PU_weights.FindBin(chain.Pileup_nTrueInt)))
			this_vertex_prob = event.mm_kin_vtx_prob[0]
			this_vertex_det = event.mm_kin_slxy[0]
			passed = False
			if event.HLT_DoubleMu4_Jpsi_Displaced:
				passed = True
				n_pass_MC += PU_weight
			if event.HLT_DoubleMu4_Jpsi_NoVertexing:
				this_bin = eff_data.FindFixBin(this_vertex_prob, this_vertex_det)
				n_pass += ufloat(eff_data.GetBinContent(this_bin), eff_data.GetBinError(this_bin))*PU_weight
				n_tot += PU_weight
				eff_MC.FillWeighted(passed, PU_weight.n, this_vertex_prob, this_vertex_det)


print "++++++++++++++++++++++++++++++++++++++++++++"
print "Global numbers:"
eps_data = n_pass/n_tot
print "Data-driven efficiency =", eps_data
eps_MC = n_pass_MC/n_tot
print "Simulated efficiency =", eps_MC
diff = eps_data-eps_MC
print "Relative difference (data-MC) =", diff.n/eps_data.n
print "Difference significance =", abs(diff.n)/diff.s
print "Ratio (data/MC): ", eps_MC/eps_data

eff_data.SetTitle("data; vertex probability ; vertex detachment significance;"); 

c = r.TCanvas("c")
eff_data.Draw("colz")
r.gStyle.SetOptStat(0)
c.SaveAs("plots/data_eff"+year+".pdf")
c.SaveAs("plots/data_eff"+year+".root")

c2 = r.TCanvas("c2")
eff_MC.Draw("colz")
r.gPad.Update()
eff_MC.SetLineColor(r.kRed)
r.gPad.Update()
c2.SaveAs("plots/MC_eff"+year+".pdf")
c2.SaveAs("plots/MC_eff"+year+".root")

#make ratio between data and MC efficiencies
eff_ratio = eff_data.Clone()
eff_ratio.Divide(eff_MC.GetPaintedHistogram())
c3 = r.TCanvas("c3")
eff_ratio.Draw("colz")
c3.SaveAs("plots/eff_ratio"+year+".pdf")
c3.SaveAs("plots/eff_ratio"+year+".root")
#make difference significance histogram
eff_diffs = eff_data.Clone()
for b in range(eff_diffs.GetNbinsX()*eff_diffs.GetNbinsY()):
	eff_MC_bin_avg_error = 0.5*(eff_MC.GetEfficiencyErrorUp(b)+eff_MC.GetEfficiencyErrorLow(b))
	if not eff_MC_bin_avg_error==0:	
		diffs = (eff_data.GetBinContent(b)-eff_MC.GetEfficiency(b))/sqrt((eff_data.GetBinError(b))**2+(eff_MC_bin_avg_error)**2)
	else:
		diffs = 0
	eff_diffs.SetBinContent(b, diffs)
c4 = r.TCanvas("c4")
#eff_diffs.GetZaxis().SetRangeUser(0,5)
eff_diffs.Draw("colz")
c4.SaveAs("plots/eff_diffs"+year+".pdf")
c4.SaveAs("plots/eff_diffs"+year+".root")

#raw_input()
