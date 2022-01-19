TH2F* Divide2DHisto(const TH2F* h1, const TH2F* h2, const TString& name) {

  if(h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()) {
    throw std::runtime_error("Histograms should have the same number of bins!");
  }

  TH2F* res = dynamic_cast<TH2F*>(h1->Clone(name));
  res->Divide(h2);

  for(int i = 1; i <= h1->GetNbinsX(); ++i){
    for(int j = 1; j <= h1->GetNbinsY(); ++j) {
      const auto a = h1->GetBinContent(i, j);
      const auto b = h2->GetBinContent(i, j);
      const auto da = h1->GetBinError(i, j);
      const auto db = h2->GetBinError(i, j);
      if(a == 0 || b == 0) {
        continue;
      }
      double error = res->GetBinContent(i, j) * sqrt(da * da / a / a + db * db / b / b);
      res->SetBinError(i, j, error);
    }
  }
  return res;
}


void PlotHistos(TH1* mc, TH1* reco){

  reco->SetMarkerStyle(21);
  reco->SetMarkerColor(kRed);
  reco->SetLineColor(kRed);
  reco->SetMarkerSize(2);
  reco->SetStats(0);
  reco->SetTitleSize(0.06);
  reco->GetXaxis()->SetRangeUser(0.2, 3);
  reco->GetXaxis()->SetTitleOffset(0.9);
  reco->GetXaxis()->SetLabelSize(0.05);
  reco->GetYaxis()->SetTitleOffset(0.7);
  reco->GetYaxis()->SetLabelSize(0.05);
  reco->GetYaxis()->SetTitleSize(0.06);
  reco->GetXaxis()->SetTitleSize(0.06);
//reco->GetYaxis()->SetRangeUser(0, 3000);
  reco->Draw("pe1");
  
  mc->SetMarkerStyle(22);
  mc->SetMarkerColor(kBlue);
  mc->SetLineColor(kBlue);
  mc->SetMarkerSize(2);
  mc->SetStats(0);
  mc->SetTitleSize(0.06);
  mc->GetXaxis()->SetRangeUser(0.2, 3);
  mc->GetXaxis()->SetTitleOffset(0.9);
  mc->GetXaxis()->SetLabelSize(0.05);
  mc->GetYaxis()->SetTitleOffset(0.7);
  mc->GetYaxis()->SetLabelSize(0.05);
  mc->GetYaxis()->SetTitleSize(0.06);
  mc->GetXaxis()->SetTitleSize(0.06);
//  mc->GetYaxis()->SetRangeUser(0, 3000);
  mc->Draw("pe1,same");
  
  auto legend = new TLegend(0.13,0.3,0.65,0.2);
  legend->AddEntry(mc,"Simulated","ep");
  legend->AddEntry(reco,"Efficiency corrected","ep");
  legend -> SetLineWidth (0);
  legend->SetTextSize(0.039);
  legend->Draw();
  

}


void plot_spectra1() {
    gStyle->SetFillColor(-1);
    //gStyle->SetFillStyle(4000);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetFrameFillColor(-1);

  auto file_dcm = TFile::Open("new_dcm_100_efficiency_pt_y_yield_bdt_cut_0.9.root");
  auto file_urqmd = TFile::Open("new_urqmd_efficiency_pt_y_yield_bdt_cut_0.9.root");

  auto rec_urqmd = file_urqmd->Get<TH2F>("recons_urqmd");
  auto sim_urqmd = file_urqmd->Get<TH2F>("Mc_urqmd");
  auto rec_dcm = file_dcm->Get<TH2F>("recons");
  auto sim_dcm = file_dcm->Get<TH2F>("Mc");

  auto eff_urqm = Divide2DHisto(rec_urqmd, sim_urqmd, "urqmd_eff");
  auto eff_dcm = Divide2DHisto(rec_dcm, sim_dcm, "dcm_eff");

  auto corr_urqmd = Divide2DHisto(rec_urqmd, eff_dcm, "corr_urqmd");
  
  int bin1 = 5;
  int bin2 = 5;
  auto corr_urqmd_y = corr_urqmd->ProjectionX("corr_urqmd_y", bin1,bin2);
  auto mc_urqmd_y = sim_urqmd->ProjectionX("mc_urqmd_y", bin1,bin2);
  
  int bin3 = 7;
  int bin4 = 7;
  auto corr_urqmd_pT = corr_urqmd->ProjectionY("corr_urqmd_pT", bin3,bin4);
  auto mc_urqmd_pT = sim_urqmd->ProjectionY("mc_urqmd_pT", bin3,bin4);
  
  
  corr_urqmd_y->SetTitle(Form("p_{T} = [%0.1f,%0.1f]",corr_urqmd->GetYaxis()->GetBinLowEdge(bin1), corr_urqmd->GetYaxis()->GetBinUpEdge(bin2)));
  corr_urqmd_y->GetXaxis()->SetTitle("#it{y}_{Lab}");
  corr_urqmd_y->GetYaxis()->SetTitle("Yield");
  mc_urqmd_y->SetTitle("X Projection of Corrected URQMD and MC URQMD");
  mc_urqmd_y->GetXaxis()->SetTitle("#it{y}_{Lab}");
  mc_urqmd_y->GetYaxis()->SetTitle("Yield");
  
  corr_urqmd_pT->SetTitle(Form("#it{y}_{LAB} =[%0.1f,%0.1f]",corr_urqmd->GetYaxis()->GetBinLowEdge(bin3), corr_urqmd->GetYaxis()->GetBinUpEdge(bin4)));
  corr_urqmd_pT->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  corr_urqmd_pT->GetYaxis()->SetTitle("Yield");
  mc_urqmd_pT->SetTitle("Y Projection of Corrected URQMD and MC URQMD");
  mc_urqmd_pT->GetXaxis()->SetTitle("(GeV/#it{c})");
  mc_urqmd_pT->GetYaxis()->SetTitle("Yield");

  auto c = new TCanvas (" canvas ","", 800,400);
  c->Draw();
  gStyle->SetErrorX(0);
  c->Divide(2,1);
  TPad *pad1 = (TPad*)c->cd(1);
  pad1->SetTopMargin (0.01);
  pad1->SetBottomMargin (0.15);
  pad1->SetLogy();
  pad1->Draw ();
  pad1->cd();
  PlotHistos(mc_urqmd_pT, corr_urqmd_pT);
  auto latex = new TLatex ();
  latex->SetNDC ();
  latex->SetTextSize (0.0449);
  latex->DrawLatex (0.5 ,0.83, "CBM performance");
  latex->DrawLatex (0.5 ,0.73, "URQMD, Au+Au" );
  latex->DrawLatex (0.5 ,0.66, "@ 12#it{A} GeV/#it{c}" );
  latex->DrawLatex (0.5 ,0.94, "#Lambda hyperons");

  latex->Draw();
  
  
  TPad *pad2 = (TPad*)c->cd(2);
  pad2->SetTopMargin (0.01);
  pad2->SetBottomMargin (0.15);
  pad2->SetLeftMargin (0.15);
  pad2->Draw ();
  pad2->cd ();
  PlotHistos( mc_urqmd_y, corr_urqmd_y);

  
}
