import ROOT as r
from uncertainties import ufloat

chain = r.TChain("Events")
chain.Add("/eos/user/g/gandreas/jpsikmc/2017/*.root")

f_map = r.TFile.Open("hists2017.root")
h_pass = f_map.Get("hpass")
h_pass.Sumw2()
h_pass.Scale(1./0.001694)
h_tot = f_map.Get("htot")
h_tot.Sumw2()
h_tot.Scale(1./0.001254)

eff_data = r.TEfficiency(h_pass, h_tot)
eff_MC = eff_data.Clone("simulation")
eff_MC.GetPassedHistogram().Reset()
eff_MC.GetTotalHistogram().Reset()

#ge tpileup histogram for data and MC
data_PU_file = r.TFile.Open("DataPileup_2017.root")
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
		if event.HLT_Dimuon0_Jpsi_NoVertexing\
			and abs(event.mm_gen_pdgId[0])==443 \
			and abs(event.mm_gen_mu1_pdgId[0])==13 and abs(event.mm_gen_mu2_pdgId[0])==13\
			and event.mm_gen_mu1_pdgId[0]*event.mm_gen_mu2_pdgId[0]<0:

			PU_weight = ufloat(PU_weights.GetBinContent(PU_weights.FindBin(chain.Pileup_nTrueInt)), PU_weights.GetBinError(PU_weights.FindBin(chain.Pileup_nTrueInt)))

			this_vertex_prob = event.mm_kin_vtx_prob[0]
			this_bin = eff_data.FindFixBin(this_vertex_prob)
			n_pass += ufloat(eff_data.GetEfficiency(this_bin), 0.5*(eff_data.GetEfficiencyErrorLow(this_bin)+eff_data.GetEfficiencyErrorUp(this_bin)))*PU_weight
			n_tot += PU_weight
			passed = False
			if event.HLT_Dimuon0_Jpsi:
				passed = True
				n_pass_MC += PU_weight
			eff_MC.FillWeighted(passed, PU_weight.n, event.mm_kin_vtx_prob[0])


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

eff_data.SetTitle("data; vertex probability ; #epsilon"); 

c = r.TCanvas("c")
eff_data.Draw("AP")
r.gPad.Update()
eff_data.GetPaintedGraph().SetMinimum(0)
eff_data.GetPaintedGraph().SetMaximum(1.05)
eff_MC.Draw("Psame")
r.gPad.Update()
eff_MC.GetPaintedGraph().SetMinimum(eff_data.GetPaintedGraph().GetMinimum())
eff_MC.GetPaintedGraph().SetMaximum(eff_data.GetPaintedGraph().GetMaximum())
eff_MC.SetLineColor(r.kRed)
r.gPad.Update()
c.BuildLegend()
# eff_MC.GetPaintedGraph().SetMarkerStyle(4)
# eff_MC.GetPaintedGraph().SetMarkerColor(r.kRed)
# eff_MC.GetPaintedGraph().SetLineColor(r.kRed)
c.SaveAs("plots/compare_effs.pdf")
c.SaveAs("plots/compare_effs.root")

#raw_input()
