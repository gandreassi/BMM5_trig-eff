import ROOT as r
from uncertainties import ufloat

chain = r.TChain("Events")
#chain.Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3+MINIAODSIM/70F7ABE3-69AA-E811-882F-B499BAAC0694.root")
chain.Add("/eos/user/g/gandreas/jpsikmc/2017/12F836CC-DA41-E811-AB46-0CC47A57D036.root")

f_map = r.TFile.Open("hists2017.root")
h_pass = f_map.Get("hpass")
h_pass.Sumw2()
h_pass.Scale(1./0.001694)
h_tot = f_map.Get("htot")
h_tot.Sumw2()
h_tot.Scale(1./0.001254)

eff_data = h_pass.Clone("h_eff_data")
eff_data.Divide(h_tot)

#eff_data = r.TEfficiency(h_pass, h_tot)
#eff_MC = r.TEfficiency("eff_MC", "eff_MC", eff_data.GetNbinsX(), eff_data.GetXaxis().GetXmin(), eff_data.GetXaxis().GetXmax())
h_tot_MC  = eff_data.Clone("h_tot_MC")
h_pass_MC = eff_data.Clone("h_pass_MC")
h_tot_MC.Reset()
h_pass_MC.Reset()

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
			n_pass += ufloat(eff_data.GetBinContent(this_bin), eff_data.GetBinError(this_bin))*PU_weight
			n_tot += PU_weight
			passed = False
			if event.HLT_Dimuon0_Jpsi:
				passed = True
				n_pass_MC += PU_weight
				h_pass_MC.Fill(event.mm_kin_vtx_prob[0], PU_weight.n)
			h_tot_MC.Fill(event.mm_kin_vtx_prob[0], PU_weight.n)

print n_pass_MC
print n
h_pass_MC.Draw()
raw_input()
h_pass_MC.Sumw2()
h_tot_MC.Sumw2()
eff_MC = h_pass_MC.Clone("h_eff_MC")
eff_MC.Divide(h_tot_MC)

print "++++++++++++++++++++++++++++++++++++++++++++"
print "Global numbers:"
eps_data = n_pass/n_tot
print "Data-driven efficiency =", eps_data
eps_MC = n_pass_MC/n_tot
print "Simulated efficiency =", eps_MC
print "Ratio (MC/data): ", eps_MC/eps_data

eff_data.SetTitle("Vertex trigger efficiency; vertex probability ; #epsilon"); 

c = r.TCanvas("c")
eff_data.Draw("E")
r.gPad.Update()
eff_data.GetYaxis().SetRangeUser(0,1.05)
eff_MC.Draw("Esame")
eff_MC.SetLineColor(r.kRed)
r.gPad.Update()
# eff_MC.GetPaintedGraph().SetMarkerStyle(4)
# eff_MC.GetPaintedGraph().SetMarkerColor(r.kRed)
# eff_MC.GetPaintedGraph().SetLineColor(r.kRed)
c.SaveAs("plots/compare_effs.pdf")
c.SaveAs("plots/compare_effs.root")

#raw_input()
