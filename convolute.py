import ROOT as r 

chain = r.TChain("Events")
#chain.Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3+MINIAODSIM/70F7ABE3-69AA-E811-882F-B499BAAC0694.root")
chain.Add("/eos/user/g/gandreas/jpsikmc/*.root")

f_map = r.TFile.Open("hists.root")
h_pass = f_map.Get("hpass")
h_pass.Scale(1./0.001694)
h_tot = f_map.Get("htot")
h_tot.Scale(1./0.001254)

eff_data = r.TEfficiency(h_pass, h_tot)
eff_MC = eff_data.Clone("MC")
eff_MC.GetPassedHistogram().Reset()
eff_MC.GetTotalHistogram().Reset()

n_pass = 0.
n_tot = 0

n_pass_MC = 0.
n_tot_MC = 0

for event in chain:
	if len(event.mm_kin_vtx_prob)>0:
		if event.HLT_Dimuon0_Jpsi_NoVertexing\
			and abs(event.mm_gen_pdgId[0])==443 \
			and abs(event.mm_gen_mu1_pdgId[0])==13 and abs(event.mm_gen_mu2_pdgId[0])==13\
			and event.mm_gen_mu1_pdgId[0]*event.mm_gen_mu2_pdgId[0]<0:

			this_vertex_prob = event.mm_kin_vtx_prob[0]
			this_bin = eff_data.FindFixBin(this_vertex_prob)
			n_pass += eff_data.GetEfficiency(this_bin)
			n_tot += 1
			passed = False
			if event.HLT_Dimuon0_Jpsi:
				passed = True
				n_pass_MC += 1
			eff_MC.Fill(passed, event.mm_kin_vtx_prob[0])



print "++++++++++++++++++++++++++++++++++++++++++++"
print "Global numbers:"
print "Data-driven efficiency =", n_pass/n_tot
print "Simulated efficiency =", n_pass_MC/n_tot

print n_pass_MC

eff_data.SetTitle("Vertex trigger efficiency; vertex probability ; #epsilon"); 

c = r.TCanvas("c")
eff_data.Draw("AP")
r.gPad.Update()
eff_data.GetPaintedGraph().SetMinimum(0)
eff_data.GetPaintedGraph().SetMaximum(1.05)
eff_MC.Draw("Psame")
r.gPad.Update()
eff_MC.GetPaintedGraph().SetMinimum(eff_data.GetPaintedGraph().GetMinimum())
eff_MC.GetPaintedGraph().SetMaximum(eff_data.GetPaintedGraph().GetMaximum())
# eff_MC.GetPaintedGraph().SetMarkerStyle(4)
# eff_MC.GetPaintedGraph().SetMarkerColor(r.kRed)
# eff_MC.GetPaintedGraph().SetLineColor(r.kRed)
c.SaveAs("plots/compare_effs.pdf")
c.SaveAs("plots/compare_effs.root")

#raw_input()
