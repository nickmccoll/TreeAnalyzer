//MAKE TOTAL BACKGROUND
{
std::vector<TString> bkgs = {"qg","losttw","mt","mw"};
TString reg = "mu_L_LP_full";

TH2 * outH=0;

for(unsigned int iB = 0; iB < bkgs.size();++iB){
    TFile * fmc = new TFile("HHlnujj_"+bkgs[iB]+"_distributions.root");
    TH2 * hmc = 0;
    fmc->GetObject(bkgs[iB]+"_"+reg+"_hbbMass_hhMass",hmc);
    if(hmc==0){
        std::cout << bkgs[iB]+"_"+reg+"_hbbMass_hhMass"<<std::endl;
        continue;
    }
    TFile * ftemp = new TFile("HHlnujj_"+bkgs[iB]+"_2D_template_debug.root");
    TH2 * hmtemp = 0;
    ftemp->GetObject(bkgs[iB]+"_"+reg,hmtemp);
    if(hmtemp==0)continue;
    hmtemp->Scale(hmc->Integral()/hmtemp->Integral());
    if(outH == 0) outH=(TH2*)hmtemp->Clone();
    else outH->Add(hmtemp);
}

outH->Draw("COLZ");



}


{
  TFile * f = new TFile("bkgInputs/HHlnujj_mt_emu_M_I_none_MJJ_fit.root");
  TCanvas* c=0;
  f->GetObject("can_1000to1100",c);
auto l  = c->GetListOfPrimitives();
for(unsigned int iP = 0; iP < l->GetSize(); ++iP){
    auto * p  =   l->At(iP);
    std::string nm = p->ClassName();
    if(nm.find("TH1D") == std::string::npos) continue;
    ((TH1*)(p))->GetXaxis()->SetTitle("#it{m}_{H#rightarrowbb} [GeV]");
    ((TH1*)(p))->GetYaxis()->SetTitle("N. of events");    
}
c->cd();
TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.03);
latex.DrawLatex(0.15,0.88,"e#mu, bM, WI, -ExB -#it{m}_{D} -#it{p/m} -#tau_{0.75} 1<#it{m}_{HH}<1.1 TeV]");
c->Draw();
c->Print("baseline/mt_examplefit.pdf");
}  


{
  TFile * f = new TFile("bkgInputs/HHlnujj_mw_emu_LMT_I_none_MJJ_fit.root");
  TCanvas* c=0;
  f->GetObject("can_900to1000",c);
auto l  = c->GetListOfPrimitives();
for(unsigned int iP = 0; iP < l->GetSize(); ++iP){
    auto * p  =   l->At(iP);
    std::string nm = p->ClassName();
    if(nm.find("TH1D") == std::string::npos) continue;
    ((TH1*)(p))->GetXaxis()->SetTitle("#it{m}_{H#rightarrowbb} [GeV]");
    ((TH1*)(p))->GetYaxis()->SetTitle("N. of events");    
}
c->cd();
TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.03);
latex.DrawLatex(0.15,0.88,"All categories, -ExB -#it{m}_{D} -#it{p/m} -#tau_{0.75} [1<#it{m}_{HH}<1.1 TeV]");
c->Draw();
c->Print("baseline/mw_examplefit.pdf");
}  

  
{
  TFile * f = new TFile("signalInputs/HHlnujj_radHH_emu_M_I_ltmb_MJJ_fit.root");
  TCanvas* c=0;
  f->GetObject("can_m1600",c);
auto l  = c->GetListOfPrimitives();
for(unsigned int iP = 0; iP < l->GetSize(); ++iP){
    auto * p  =   l->At(iP);
    std::string nm = p->ClassName();
    if(nm.find("TH1D") == std::string::npos) continue;
    ((TH1*)(p))->GetXaxis()->SetTitle("#it{m}_{H#rightarrowbb} [GeV]");
    ((TH1*)(p))->GetYaxis()->SetTitle("arbitrary units");    
}
c->cd();
TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.03);
latex.DrawLatex(0.15,0.88,"e#mu, bM, WI, -ExB -#tau_{0.75} [1.6 TeV X]");
c->Draw();
c->Print("baseline/mhbbFitExample.pdf");
}  

///must uncomment line in makeSignalInputs to get "forAN" vers
{
  TFile * f = new TFile("baseline/radHH_2D_fit_forAN.root");
  TCanvas* c=0;
  f->GetObject("mu_M_LP_full_COMP",c);
auto l  = c->GetListOfPrimitives();
for(unsigned int iP = 0; iP < l->GetSize(); ++iP){
    auto * p  =   l->At(iP);
    std::string nm = p->ClassName();
    if(nm.find("TH1D") == std::string::npos) continue;
    ((TH1*)(p))->GetXaxis()->SetTitle("#it{m}_{H#rightarrowbb} [GeV]");
    ((TH1*)(p))->GetYaxis()->SetTitle("arbitrary units");    
}
c->cd();
TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.03);
latex.DrawLatex(0.15,0.88,"#mu, bM, WL [1.6 TeV X]");
c->Draw();
c->Print("baseline/mhhFitExample.pdf");
}  

{
  TFile * f = new TFile("signalInputs/HHlnujj_radHH_mu_M_LP_full_2D_fit.root");
  TH2* c=0;
  f->GetObject("pdf_m1600__MJ_MR",c);
  c->GetXaxis()->SetTitle("#it{m}_{H#rightarrowbb} [GeV]");
  c->GetYaxis()->SetTitle("#it{m}_{HH} [GeV]");    
  TCanvas * can = new TCanvas();
  c->Draw("COLZ");
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  latex.SetTextColor(kWhite);
  latex.DrawLatex(0.15,0.88,"#mu, bM, WL [1.6 TeV X]");
  can->Draw();
  can->Print("baseline/sig_corr2DEx.pdf");
}
  