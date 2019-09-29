import ROOT as r
from uncertainties import ufloat
from sys import argv
from math import sqrt

year = argv[1]
SN = argv[2]

lumi={	"2017" : {	"HLT_Dimuon0_Jpsi" : 42974747.938,
					"HLT_Dimuon0_Jpsi_NoVertexing" : 37468393.496,
					"HLT_DoubleMu4_3_Jpsi_Displaced" : 7376118661.867,
					"HLT_DoubleMu4_Jpsi_Displaced" :   876993819.962,
					"HLT_DoubleMu4_Jpsi_NoVertexing" : 234517906.885},
		"2018" : {	"HLT_Dimuon0_Jpsi" : 2692017.753,
					"HLT_Dimuon0_Jpsi_NoVertexing" : 4036006.636,
					"HLT_DoubleMu4_3_Jpsi" : 1173487158.329,
					"HLT_DoubleMu4_Jpsi_Displaced" : 212306014.543,
					"HLT_DoubleMu4_Jpsi_NoVertexing" : 52366197.343}}
# lumi={	"2017" : {	"HLT_Dimuon0_Jpsi" : 75350805.341,
# 					"HLT_Dimuon0_Jpsi_NoVertexing" : 55714920.339,
# 					"HLT_DoubleMu4_3_Jpsi_Displaced" : 14479239278.717,
# 					"HLT_DoubleMu4_Jpsi_Displaced" :   876993819.962,
# 					"HLT_DoubleMu4_Jpsi_NoVertexing" : 234517906.885},
# 		"2018" : {	"HLT_Dimuon0_Jpsi" : 2692017.753,
# 					"HLT_Dimuon0_Jpsi_NoVertexing" : 4036006.636,
# 					"HLT_DoubleMu4_3_Jpsi" : 5889941595.296, 908127530.647321224
# 					"HLT_DoubleMu4_Jpsi_Displaced" : 1038251431.797, 168079409.930
# 					"HLT_DoubleMu4_Jpsi_NoVertexing" : 258113597.156}} 41309546.189

if (year=="2017" and SN=="S"):
	trig_A = "HLT_Dimuon0_Jpsi"
	trig_B = "HLT_Dimuon0_Jpsi_NoVertexing"
	L1_seed = "1" #take everything
elif (year=="2017" and SN=="N"):
	trig_A = "HLT_DoubleMu4_Jpsi_Displaced"
	trig_B = "HLT_DoubleMu4_Jpsi_NoVertexing"
	L1_seed = "1"#"(L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4)"
elif (year=="2018" and SN=="S"):
	trig_A = "HLT_Dimuon0_Jpsi"
	trig_B = "HLT_Dimuon0_Jpsi_NoVertexing"
	L1_seed = "1" #take everything
elif (year=="2018" and SN=="N"):
	trig_A = "HLT_DoubleMu4_Jpsi_Displaced"
	trig_B = "HLT_DoubleMu4_Jpsi_NoVertexing"
	L1_seed = "1"#"(L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4)"
else: raise Exception("ERROR: unexpected combination of year and S/N flag (signal/normalisation).")

chain = r.TChain("Events")
chain.Add("/eos/user/g/gandreas/jpsikmc/"+year+"/*.root")
#deactivate useless branches
chain.SetBranchStatus("*", 0)
chain.SetBranchStatus("HLT*", 1)
chain.SetBranchStatus("L1*", 1)
chain.SetBranchStatus("mm_*", 1)
chain.SetBranchStatus("bkmm_*", 1)
chain.SetBranchStatus("Pileup_nTrueInt", 1)

f_map = r.TFile.Open("hists"+year+SN+".root")
h_pass = f_map.Get("hpass")
h_pass.Sumw2()
h_pass.Scale(1./lumi[year][trig_A])
h_tot = f_map.Get("htot")
h_tot.Sumw2()
h_tot.Scale(1./lumi[year][trig_B])

eff_data = h_pass.Clone()
eff_data.Divide(h_tot)

binningX = eff_data.GetXaxis().GetXbins()
binningY = eff_data.GetYaxis().GetXbins()
eff_MC = eff_data.Clone("simulation")
eff_MC = r.TEfficiency("eff_MC","simulation ; vertex probability ; vertex detachment significance;", eff_data.GetNbinsX(), binningX.GetArray(), eff_data.GetNbinsY(), binningY.GetArray())
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

f_buffer = r.TFile.Open("f_buffer.root", "recreate")
f_buffer.cd()
matching_cut_jpsi = "abs(mm_gen_pdgId[bkmm_mm_index[0]])==443 && abs(mm_gen_mu1_pdgId[bkmm_mm_index[0]])==13 && "\
					 "abs(mm_gen_mu2_pdgId[bkmm_mm_index[0]])==13 &&" \
					 "(mm_gen_mu1_pdgId[bkmm_mm_index[0]]*mm_gen_mu2_pdgId[bkmm_mm_index[0]])<0 &&" \
					 "abs(bkmm_gen_pdgId[0])==521"
additional_cut = matching_cut_jpsi+" && mm_kin_cosAlpha[bkmm_mm_index[0]]>0.9"
chain_matched = chain.CopyTree(matching_cut_jpsi)
del chain

### L1 seed selection
#L1_eff = float(chain_matched.GetEntries(matching_cut_jpsi+" && {0} && {1}".format(trig_B, L1_seed)))/chain_matched.GetEntries(matching_cut_jpsi+" && {0}".format(trig_B))
#eff_data.Scale(1./L1_eff)
#print L1_eff, "<<< this number should be 1 if you are running on S (signal triggers)."
#if SN=="S" and L1_eff!=1:
#	raise Exception("45.961216, 12.976180")

if SN=="N": chain_matched_L1 = chain_matched.CopyTree(additional_cut)
else: chain_matched_L1 = chain_matched.Clone()
del chain_matched

n_pass = ufloat(0.,0.)
n_tot = ufloat(0.,0.)

n_pass_MC = 0.
L1_tot=0
L1_pass = 0.

h_PU_weights = r.TH1F("PU_weights", "PU weights", 100, 0, 4)
#make trigger map on MC
for event in chain_matched_L1:
		PU_weight = ufloat(PU_weights.GetBinContent(PU_weights.FindBin(chain_matched_L1.Pileup_nTrueInt)), PU_weights.GetBinError(PU_weights.FindBin(chain_matched_L1.Pileup_nTrueInt)))
		h_PU_weights.Fill(PU_weight.n)
		this_vertex_prob = event.mm_kin_vtx_prob[0]
		this_vertex_det = event.mm_kin_slxy[0]
		if getattr(event, trig_B):
			eff_MC.FillWeighted(getattr(event, trig_A), PU_weight.n, this_vertex_prob, this_vertex_det)
		mm_index = event.bkmm_mm_index[0]
		this_vertex_prob = event.mm_kin_vtx_prob[mm_index]
		this_vertex_det = event.mm_kin_slxy[mm_index]
		passed = False
		if getattr(event, trig_A):
			passed = True
			n_pass_MC += PU_weight
		if getattr(event, trig_B):
			this_bin = eff_data.FindBin(this_vertex_prob, this_vertex_det)
			n_pass += ufloat(eff_data.GetBinContent(this_bin), eff_data.GetBinError(this_bin))*PU_weight
			n_tot += PU_weight

c_PU = r.TCanvas("c_PU")
h_PU_weights.Draw()
c_PU.SaveAs("plots/PU_weights.pdf")
				
f_buffer.Close()

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
eff_data.GetZaxis().SetRangeUser(0,1)
r.gStyle.SetOptStat(0)
c.SaveAs("plots/data_eff"+year+SN+".pdf")
c.SaveAs("plots/data_eff"+year+SN+".root")

c2 = r.TCanvas("c2")
eff_MC.Draw("colz")
r.gPad.Update()
eff_MC.GetPaintedHistogram().GetZaxis().SetRangeUser(0,1)
eff_MC.SetLineColor(r.kRed)
r.gPad.Update()
c2.SaveAs("plots/MC_eff"+year+SN+".pdf")
c2.SaveAs("plots/MC_eff"+year+SN+".root")

#make ratio between data and MC efficiencies
eff_ratio = eff_data.Clone()
eff_ratio.GetZaxis().UnZoom()
eff_ratio.Divide(eff_MC.GetPaintedHistogram())
c3 = r.TCanvas("c3")
eff_ratio.Draw("colz")
c3.SaveAs("plots/eff_ratio"+year+SN+".pdf")
c3.SaveAs("plots/eff_ratio"+year+SN+".root")

#make difference significance histogram
eff_diffs = eff_data.Clone()
eff_diffs.Reset()
eff_diffs.GetZaxis().UnZoom()
for b in range((eff_diffs.GetNbinsX()+2)*(eff_diffs.GetNbinsY()+2)):
	eff_MC_bin_avg_error = 0.5*(eff_MC.GetEfficiencyErrorUp(b)+eff_MC.GetEfficiencyErrorLow(b))
	if not eff_MC_bin_avg_error==0:	
		diffs = (eff_data.GetBinContent(b)-eff_MC.GetEfficiency(b))/sqrt((eff_data.GetBinError(b))**2+(eff_MC_bin_avg_error)**2)
	else:
		diffs = 0
	eff_diffs.SetBinContent(b, diffs)
c4 = r.TCanvas("c4")
#eff_diffs.GetZaxis().SetRangeUser(0,5)
eff_diffs.Draw("colz")
c4.SaveAs("plots/eff_diffs"+year+SN+".pdf")
c4.SaveAs("plots/eff_diffs"+year+SN+".root")


#raw_input()
