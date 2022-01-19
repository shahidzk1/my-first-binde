#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

void plot_histo(){

  TFile* file = TFile::Open("data/pt_y_yield_bdt_cut_0.9.root", "r");
  auto tree = file->Get<TTree>("t1");

  float y, pT, yield;
  tree->SetBranchAddress("y_min", &y);
  tree->SetBranchAddress("pT_min", &pT);
  tree->SetBranchAddress("yield", &yield);

  TH2F* h2 = new TH2F("h2", "Lambda yield", 15, 0., 3., 15, 0., 3.);

  const auto n_enties = tree->GetEntries();
  for(int i=0; i<n_enties; ++i){
    tree->GetEntry(i);
    const int y_bin = (y+0.1)/0.2 + 1;
    const int pT_bin = (pT+0.1)/0.2 + 1;
    h2->SetBinContent(y_bin, pT_bin, yield);
  }

  h2->Draw("colz");
}

int main(){
  plot_histo();
  return 0;
}
