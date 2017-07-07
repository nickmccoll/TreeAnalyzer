{
  TFile * f = new TFile("zMuMu_plots.root");  
auto plot= [&] (TString* num, TString den, TString name ){
  Plotter * p = new Plotter();
  TH1 * hd = 0;
  f->GetObject(den,hd);
  if(!hd) return;
  p->addHistLine(hd,den);
  for(unsigned int iN = 0; num[iN][0]; ++iN){
    TH1 * hn = 0;
    f->GetObject(num[iN],hn);
    if(!hn) continue;
    p->addHistLine(hn,num[iN]);
  }
  p->drawSplitRatio(0,"stack",true,false,name);
};

TString enum1[] = {
"electron_OS_passVetoID_pt"       ,
"electron_OS_passLooseID_pt"      ,
"electron_OS_passMedID_pt"        ,
"electron_OS_passTightID_pt"      ,
"electron_OS_passHEEPID_pt"       ,""};
TString eden = "electron_OS_pt";
plot(enum1,eden,"ele");


TString enum2[] = {
"electron_OS_passVetoID_noISO_pt" ,
"electron_OS_passLooseID_noISO_pt",
"electron_OS_passMedID_noISO_pt"  ,
"electron_OS_passTightID_noISO_pt",
"electron_OS_passHEEPID_noISO_pt" , ""};
plot(enum2,eden,"ele_noIso");


TString mnum[] = {
  "muon_OS_passSoftID_pt" ,
  "muon_OS_passLooseID_pt",
  "muon_OS_passMedID_pt"  ,
  "muon_OS_passTightID_pt",
  "muon_OS_passMed16ID_pt",
  "muon_OS_passHighPT_pt" ,""};
  TString mden = "muon_OS_pt";
  plot(mnum,mden,"muon");


}