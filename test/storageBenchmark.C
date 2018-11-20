
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TStopwatch.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

TStopwatch timer;
UInt_t nEvents = 100000;
int avg = 12;
int stdDev = 4;
UInt_t minJ = 3;
UInt_t maxJ = 100;
TRandom3 * rnd = new TRandom3(1234);


class ReduceMantissaToNbitsRounding {
public:
    ReduceMantissaToNbitsRounding(int bits) :
        shift(23-bits), mask((0xFFFFFFFF >> (shift)) << (shift)),
        test(1 << (shift-1)), maxn((1<<bits)-2) {
        assert(bits <= 23); // "max mantissa size is 23 bits"
    }
    float operator()(float f) const {
        constexpr uint32_t low23 = (0x007FFFFF); // mask to keep lowest 23 bits = mantissa
        constexpr uint32_t  hi9  = (0xFF800000); // mask to keep highest 9 bits = the rest
        union { float flt; uint32_t i32; } conv;
        conv.flt=f;
        if (conv.i32 & test) { // need to round
            uint32_t mantissa = (conv.i32 & low23) >> shift;
            if (mantissa < maxn) mantissa++;
            conv.i32 = (conv.i32 & hi9) | (mantissa << shift);
        } else {
            conv.i32 &= mask;
        }
        return conv.flt;
    }
private:
    const int shift;
    const uint32_t mask, test, maxn;
};

void storeArray(int bits_=-1){
    float col1 [maxJ];
    float col2 [maxJ];
    float col3 [maxJ];
    float col4 [maxJ];
    float col6 [maxJ];
    float col7 [maxJ];
    float col8 [maxJ];
    float col9 [maxJ];
    float col10[maxJ];
    float col11[maxJ];
    float col12[maxJ];
    float col13[maxJ];
    float col14[maxJ];
    float col15[maxJ];
    UInt_t   nJets;

    TFile *f = new TFile("testArray.root","recreate","");
    TTree *t = new TTree("t","Reconst ntuple");
    t->Branch("nJets",&nJets,"nJets/i");
    t->Branch("col1", col1 ,"col1[nJets]/F");
    t->Branch("col2", col2 ,"col2[nJets]/F");
    t->Branch("col3", col3 ,"col3[nJets]/F");
    t->Branch("col4", col4 ,"col4[nJets]/F");
    t->Branch("col6", col6 ,"col6[nJets]/F");
    t->Branch("col7", col7 ,"col7[nJets]/F");
    t->Branch("col8", col8 ,"col8[nJets]/F");
    t->Branch("col9", col9 ,"col9[nJets]/F");
    t->Branch("col10",col10,"col10[nJets]/F");
    t->Branch("col11",col11,"col11[nJets]/F");
    t->Branch("col12",col12,"col12[nJets]/F");
    t->Branch("col13",col13,"col13[nJets]/F");
    t->Branch("col14",col14,"col14[nJets]/F");
    t->Branch("col15",col15,"col15[nJets]/F");
    for (UInt_t i=0;i<nEvents;i++) {
        UInt_t nt = rnd->Gaus(avg,stdDev);
        nt = std::min(std::max(nt,minJ),maxJ);
        nJets=nt;
        for (UInt_t n=0;n<nt;n++) {
            col1 [n] = rnd->Gaus(100,30)*rnd->Gaus(1 ,1./ 5.);
            col2 [n] = rnd->Gaus(100,30)*rnd->Gaus(2 ,2./ 5.);
            col3 [n] = rnd->Gaus(100,30)*rnd->Gaus(3 ,3./ 5.);
            col4 [n] = rnd->Gaus(100,30)*rnd->Gaus(4 ,4./ 5.);
            col6 [n] = rnd->Gaus(100,30)*rnd->Gaus(6 ,6./ 5.);
            col7 [n] = rnd->Gaus(100,30)*rnd->Gaus(7 ,7./ 5.);
            col8 [n] = rnd->Gaus(100,30)*rnd->Gaus(8 ,8./ 5.);
            col9 [n] = rnd->Gaus(100,30)*rnd->Gaus(9 ,9./ 5.);
            col10[n] = rnd->Gaus(100,30)*rnd->Gaus(10,10./ 5.);
            col11[n] = rnd->Gaus(100,30)*rnd->Gaus(11,11./ 5.);
            col12[n] = rnd->Gaus(100,30)*rnd->Gaus(12,12./ 5.);
            col13[n] = rnd->Gaus(100,30)*rnd->Gaus(13,13./ 5.);
            col14[n] = rnd->Gaus(100,30)*rnd->Gaus(14,14./ 5.);
            col15[n] = rnd->Gaus(100,30)*rnd->Gaus(15,15./ 5.);
        }
        t->Fill();


    }
    t->Write();
    f->Close();

}

void storeVector(int bits_=-1){
    std::vector<float>* col1 =new std::vector<float>;
    std::vector<float>* col2 =new std::vector<float>;
    std::vector<float>* col3 =new std::vector<float>;
    std::vector<float>* col4 =new std::vector<float>;
    std::vector<float>* col6 =new std::vector<float>;
    std::vector<float>* col7 =new std::vector<float>;
    std::vector<float>* col8 =new std::vector<float>;
    std::vector<float>* col9 =new std::vector<float>;
    std::vector<float>* col10=new std::vector<float>;
    std::vector<float>* col11=new std::vector<float>;
    std::vector<float>* col12=new std::vector<float>;
    std::vector<float>* col13=new std::vector<float>;
    std::vector<float>* col14=new std::vector<float>;
    std::vector<float>* col15=new std::vector<float>;

    TFile *f = new TFile("storeVector.root","recreate","");
    TTree *t = new TTree("t","Reconst ntuple");
    t->Branch("col1", &col1 );
    t->Branch("col2", &col2 );
    t->Branch("col3", &col3 );
    t->Branch("col4", &col4 );
    t->Branch("col6", &col6 );
    t->Branch("col7", &col7 );
    t->Branch("col8", &col8 );
    t->Branch("col9", &col9 );
    t->Branch("col10",&col10);
    t->Branch("col11",&col11);
    t->Branch("col12",&col12);
    t->Branch("col13",&col13);
    t->Branch("col14",&col14);
    t->Branch("col15",&col15);
    for (UInt_t i=0;i<nEvents;i++) {

        col1   ->clear();
        col2   ->clear();
        col3   ->clear();
        col4   ->clear();
        col6   ->clear();
        col7   ->clear();
        col8   ->clear();
        col9   ->clear();
        col10  ->clear();
        col11  ->clear();
        col12  ->clear();
        col13  ->clear();
        col14  ->clear();
        col15  ->clear();


        UInt_t nt = rnd->Gaus(avg,stdDev);
        nt = std::min(std::max(nt,minJ),maxJ);
        for (UInt_t n=0;n<nt;n++) {
         col1 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(1 ,1./ 5.) );
         col2 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(2 ,2./ 5.) );
         col3 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(3 ,3./ 5.) );
         col4 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(4 ,4./ 5.) );
         col6 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(6 ,6./ 5.) );
         col7 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(7 ,7./ 5.) );
         col8 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(8 ,8./ 5.) );
         col9 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(9 ,9./ 5.) );
         col10->push_back(rnd->Gaus(100,30)*rnd->Gaus(10,10./ 5.));
         col11->push_back(rnd->Gaus(100,30)*rnd->Gaus(11,11./ 5.));
         col12->push_back(rnd->Gaus(100,30)*rnd->Gaus(12,12./ 5.));
         col13->push_back(rnd->Gaus(100,30)*rnd->Gaus(13,13./ 5.));
         col14->push_back(rnd->Gaus(100,30)*rnd->Gaus(14,14./ 5.));
         col15->push_back(rnd->Gaus(100,30)*rnd->Gaus(15,15./ 5.));
        }
        t->Fill();


    }
    t->Write();
    f->Close();

}

void storeArrayWVector(){
    std::vector<float>* col1 =new std::vector<float>;
    std::vector<float>* col2 =new std::vector<float>;
    std::vector<float>* col3 =new std::vector<float>;
    std::vector<float>* col4 =new std::vector<float>;
    std::vector<float>* col6 =new std::vector<float>;
    std::vector<float>* col7 =new std::vector<float>;
    std::vector<float>* col8 =new std::vector<float>;
    std::vector<float>* col9 =new std::vector<float>;
    std::vector<float>* col10=new std::vector<float>;
    std::vector<float>* col11=new std::vector<float>;
    std::vector<float>* col12=new std::vector<float>;
    std::vector<float>* col13=new std::vector<float>;
    std::vector<float>* col14=new std::vector<float>;
    std::vector<float>* col15=new std::vector<float>;
    UInt_t   nJets;


    TFile *f = new TFile("storeArrayWVector.root","recreate","");
    TTree *t = new TTree("t","Reconst ntuple");
    t->Branch("nJets",&nJets,"nJets/i");
      t->Branch("col1", &((*col1 )[0]),"col1[nJets]/F");
      t->Branch("col2", &((*col2 )[0]),"col2[nJets]/F");
      t->Branch("col3", &((*col3 )[0]),"col3[nJets]/F");
      t->Branch("col4", &((*col4 )[0]),"col4[nJets]/F");
      t->Branch("col6", &((*col6 )[0]),"col6[nJets]/F");
      t->Branch("col7", &((*col7 )[0]),"col7[nJets]/F");
      t->Branch("col8", &((*col8 )[0]),"col8[nJets]/F");
      t->Branch("col9", &((*col9 )[0]),"col9[nJets]/F");
      t->Branch("col10",&((*col10)[0]),"col10[nJets]/F");
      t->Branch("col11",&((*col11)[0]),"col11[nJets]/F");
      t->Branch("col12",&((*col12)[0]),"col12[nJets]/F");
      t->Branch("col13",&((*col13)[0]),"col13[nJets]/F");
      t->Branch("col14",&((*col14)[0]),"col14[nJets]/F");
      t->Branch("col15",&((*col15)[0]),"col15[nJets]/F");
    for (UInt_t i=0;i<nEvents;i++) {

        col1   ->clear();
        col2   ->clear();
        col3   ->clear();
        col4   ->clear();
        col6   ->clear();
        col7   ->clear();
        col8   ->clear();
        col9   ->clear();
        col10  ->clear();
        col11  ->clear();
        col12  ->clear();
        col13  ->clear();
        col14  ->clear();
        col15  ->clear();


        UInt_t nt = rnd->Gaus(avg,stdDev);
        nt = std::min(std::max(nt,minJ),maxJ);
        nJets = nt;
        for (UInt_t n=0;n<nt;n++) {
         col1 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(1 ,1./ 5.) );
         col2 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(2 ,2./ 5.) );
         col3 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(3 ,3./ 5.) );
         col4 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(4 ,4./ 5.) );
         col6 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(6 ,6./ 5.) );
         col7 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(7 ,7./ 5.) );
         col8 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(8 ,8./ 5.) );
         col9 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(9 ,9./ 5.) );
         col10->push_back(rnd->Gaus(100,30)*rnd->Gaus(10,10./ 5.));
         col11->push_back(rnd->Gaus(100,30)*rnd->Gaus(11,11./ 5.));
         col12->push_back(rnd->Gaus(100,30)*rnd->Gaus(12,12./ 5.));
         col13->push_back(rnd->Gaus(100,30)*rnd->Gaus(13,13./ 5.));
         col14->push_back(rnd->Gaus(100,30)*rnd->Gaus(14,14./ 5.));
         col15->push_back(rnd->Gaus(100,30)*rnd->Gaus(15,15./ 5.));
        }

        t->SetBranchAddress("col1" ,&((*col1 )[0]));
        t->SetBranchAddress("col2" ,&((*col2 )[0]));
        t->SetBranchAddress("col3" ,&((*col3 )[0]));
        t->SetBranchAddress("col4" ,&((*col4 )[0]));
        t->SetBranchAddress("col6" ,&((*col6 )[0]));
        t->SetBranchAddress("col7" ,&((*col7 )[0]));
        t->SetBranchAddress("col8" ,&((*col8 )[0]));
        t->SetBranchAddress("col9" ,&((*col9 )[0]));
        t->SetBranchAddress("col10",&((*col10)[0]));
        t->SetBranchAddress("col11",&((*col11)[0]));
        t->SetBranchAddress("col12",&((*col12)[0]));
        t->SetBranchAddress("col13",&((*col13)[0]));
        t->SetBranchAddress("col14",&((*col14)[0]));
        t->SetBranchAddress("col15",&((*col15)[0]));

        t->Fill();


    }
    t->Write();
    f->Close();

}

void storeArrayWVectorBr(){
    std::vector<float>* col1 =new std::vector<float>;
    std::vector<float>* col2 =new std::vector<float>;
    std::vector<float>* col3 =new std::vector<float>;
    std::vector<float>* col4 =new std::vector<float>;
    std::vector<float>* col6 =new std::vector<float>;
    std::vector<float>* col7 =new std::vector<float>;
    std::vector<float>* col8 =new std::vector<float>;
    std::vector<float>* col9 =new std::vector<float>;
    std::vector<float>* col10=new std::vector<float>;
    std::vector<float>* col11=new std::vector<float>;
    std::vector<float>* col12=new std::vector<float>;
    std::vector<float>* col13=new std::vector<float>;
    std::vector<float>* col14=new std::vector<float>;
    std::vector<float>* col15=new std::vector<float>;
    UInt_t   nJets;


    TFile *f = new TFile("storeArrayWVectorBr.root","recreate","");
    TTree *t = new TTree("t","Reconst ntuple");
      t->Branch("nJets",&nJets,"nJets/i");
      TBranch * brcol1  = t->Branch("col1", 0,"col1[nJets]/F");
      TBranch * brcol2  = t->Branch("col2", 0,"col2[nJets]/F");
      TBranch * brcol3  = t->Branch("col3", 0,"col3[nJets]/F");
      TBranch * brcol4  = t->Branch("col4", 0,"col4[nJets]/F");
      TBranch * brcol6  = t->Branch("col6", 0,"col6[nJets]/F");
      TBranch * brcol7  = t->Branch("col7", 0,"col7[nJets]/F");
      TBranch * brcol8  = t->Branch("col8", 0,"col8[nJets]/F");
      TBranch * brcol9  = t->Branch("col9", 0,"col9[nJets]/F");
      TBranch * brcol10 = t->Branch("col10",0,"col10[nJets]/F");
      TBranch * brcol11 = t->Branch("col11",0,"col11[nJets]/F");
      TBranch * brcol12 = t->Branch("col12",0,"col12[nJets]/F");
      TBranch * brcol13 = t->Branch("col13",0,"col13[nJets]/F");
      TBranch * brcol14 = t->Branch("col14",0,"col14[nJets]/F");
      TBranch * brcol15 = t->Branch("col15",0,"col15[nJets]/F");
    for (UInt_t i=0;i<nEvents;i++) {

        col1   ->clear();
        col2   ->clear();
        col3   ->clear();
        col4   ->clear();
        col6   ->clear();
        col7   ->clear();
        col8   ->clear();
        col9   ->clear();
        col10  ->clear();
        col11  ->clear();
        col12  ->clear();
        col13  ->clear();
        col14  ->clear();
        col15  ->clear();


        UInt_t nt = rnd->Gaus(avg,stdDev);
        nt = std::min(std::max(nt,minJ),maxJ);
        nJets = nt;
        for (UInt_t n=0;n<nt;n++) {
         col1 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(1 ,1./ 5.) );
         col2 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(2 ,2./ 5.) );
         col3 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(3 ,3./ 5.) );
         col4 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(4 ,4./ 5.) );
         col6 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(6 ,6./ 5.) );
         col7 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(7 ,7./ 5.) );
         col8 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(8 ,8./ 5.) );
         col9 ->push_back(rnd->Gaus(100,30)*rnd->Gaus(9 ,9./ 5.) );
         col10->push_back(rnd->Gaus(100,30)*rnd->Gaus(10,10./ 5.));
         col11->push_back(rnd->Gaus(100,30)*rnd->Gaus(11,11./ 5.));
         col12->push_back(rnd->Gaus(100,30)*rnd->Gaus(12,12./ 5.));
         col13->push_back(rnd->Gaus(100,30)*rnd->Gaus(13,13./ 5.));
         col14->push_back(rnd->Gaus(100,30)*rnd->Gaus(14,14./ 5.));
         col15->push_back(rnd->Gaus(100,30)*rnd->Gaus(15,15./ 5.));
        }
        if(nt==0){
            brcol1 ->SetAddress(0);
            brcol2 ->SetAddress(0);
            brcol3 ->SetAddress(0);
            brcol4 ->SetAddress(0);
            brcol6 ->SetAddress(0);
            brcol7 ->SetAddress(0);
            brcol8 ->SetAddress(0);
            brcol9 ->SetAddress(0);
            brcol10->SetAddress(0);
            brcol11->SetAddress(0);
            brcol12->SetAddress(0);
            brcol13->SetAddress(0);
            brcol14->SetAddress(0);
            brcol15->SetAddress(0);
        } else {

        brcol1 ->SetAddress(&((*col1 )[0]));
        brcol2 ->SetAddress(&((*col2 )[0]));
        brcol3 ->SetAddress(&((*col3 )[0]));
        brcol4 ->SetAddress(&((*col4 )[0]));
        brcol6 ->SetAddress(&((*col6 )[0]));
        brcol7 ->SetAddress(&((*col7 )[0]));
        brcol8 ->SetAddress(&((*col8 )[0]));
        brcol9 ->SetAddress(&((*col9 )[0]));
        brcol10->SetAddress(&((*col10)[0]));
        brcol11->SetAddress(&((*col11)[0]));
        brcol12->SetAddress(&((*col12)[0]));
        brcol13->SetAddress(&((*col13)[0]));
        brcol14->SetAddress(&((*col14)[0]));
        brcol15->SetAddress(&((*col15)[0]));
        }
        t->Fill();


    }
    t->Write();
    f->Close();

}

void readVectorAsArray(){
    UInt_t   nJets;
    float col1  [500];
    float col2  [500];
    float col3  [500];
    float col4  [500];
    float col6  [500];
    float col7  [500];
    float col8  [500];
    float col9  [500];
    float col10 [500];
    float col11 [500];
    float col12 [500];
    float col13 [500];
    float col14 [500];
    float col15 [500];



    TFile *f = new TFile("storeArrayWVectorBr.root","read","");
    TTree *t = 0;
    f->GetObject("t",t);

    t->SetBranchAddress("nJets",&nJets);
    t->SetBranchAddress("col1" ,col1  );
    t->SetBranchAddress("col2" ,col2  );
    t->SetBranchAddress("col3" ,col3  );
    t->SetBranchAddress("col4" ,col4  );
    t->SetBranchAddress("col6" ,col6  );
    t->SetBranchAddress("col7" ,col7  );
    t->SetBranchAddress("col8" ,col8  );
    t->SetBranchAddress("col9" ,col9  );
    t->SetBranchAddress("col10",col10 );
    t->SetBranchAddress("col11",col11 );
    t->SetBranchAddress("col12",col12 );
    t->SetBranchAddress("col13",col13 );
    t->SetBranchAddress("col14",col14 );
    t->SetBranchAddress("col15",col15 );

    Long64_t ievent = 0;


    while(t->GetEvent(ievent)){

        std::cout<< ievent <<" : ";
        for(unsigned int iJ = 0; iJ < nJets;++iJ)
            std::cout <<col1[iJ]<<",";
        std::cout<<std::endl;

        ievent++;

    }
    f->Close();

}


void readVectorAsTTreeReader(){


    TFile *f = new TFile("storeArrayWVectorBr.root","read","");
    TTreeReader t("t", f);

    TTreeReaderValue<UInt_t> nJets(t,"nJets");
    TTreeReaderArray<float> col1 (t,"col1" );
    TTreeReaderArray<float> col2 (t,"col2" );
    TTreeReaderArray<float> col3 (t,"col3" );
    TTreeReaderArray<float> col4 (t,"col4" );
    TTreeReaderArray<float> col6 (t,"col6" );
    TTreeReaderArray<float> col7 (t,"col7" );
    TTreeReaderArray<float> col8 (t,"col8" );
    TTreeReaderArray<float> col9 (t,"col9" );
    TTreeReaderArray<float> col10(t,"col10");
    TTreeReaderArray<float> col11(t,"col11");
    TTreeReaderArray<float> col12(t,"col12");
    TTreeReaderArray<float> col13(t,"col13");
    TTreeReaderArray<float> col14(t,"col14");
    TTreeReaderArray<float> col15(t,"col15");

    Long64_t ievent = 0;


    while (t.Next()) {
        std::cout<< ievent <<" : ";
        for(unsigned int iJ = 0; iJ < *nJets;++iJ)
            std::cout <<col1[iJ]<<",";
        std::cout<<std::endl;

        ievent++;
    }


    f->Close();

}

template<typename T>
class ReaderData{
public:
    void set(TTreeReader& reader,const char* branchname){
        if(data){
            delete data;
            data=0;
        }
        data = new TTreeReaderValue<T> (reader,branchname);
    }
    T* operator->() { return data->operator(); }
    T& operator*() {  return data->operator*(); }
    TTreeReaderValue<T> * data =0;
};

template<typename T>
class ReaderDataArray{
public:
    void set(TTreeReader& reader,const char* branchname){
        if(data){
            delete data;
            data=0;
        }
        data = new TTreeReaderArray<T> (reader,branchname);
    }

    T* operator->() {  return data->operator(); }
    T& operator*() {return data->operator*(); }


    T &At(std::size_t idx) { return data->At(idx); }
    const T &At(std::size_t idx) const {   return data->At(idx); }
    T &operator[](std::size_t idx) { return (*data)[idx]; }
    const T &operator[](std::size_t idx) const { return (*data)[idx];  }
    std::size_t size() const {
        return data->GetSize();
    }


    using iterator = typename TTreeReaderArray<T>::template Iterator_t<TTreeReaderArray<T>>;
    using const_iterator = typename TTreeReaderArray<T>::template Iterator_t<const TTreeReaderArray<T>>;

    iterator begin() { return data->begin(); }
    iterator end() {  return data->end(); }
    const_iterator begin() const {  return data->begin(); }
    const_iterator end() const { return data->end(); }
    const_iterator cbegin() const { return data->cbegin(); }
    const_iterator cend() const  { return data->cend(); }

    TTreeReaderArray<T> * data =0;
};

void readVectorAsTTreeReaderWithCopies(){


    TFile *f = new TFile("storeArrayWVectorBr.root","read","");
    TTreeReader t("t", f);

    ReaderData<UInt_t> nJets;
    ReaderDataArray<float> col1  ;
    ReaderDataArray<float> col2  ;
    ReaderDataArray<float> col3  ;
    ReaderDataArray<float> col4  ;
    ReaderDataArray<float> col6  ;
    ReaderDataArray<float> col7  ;
    ReaderDataArray<float> col8  ;
    ReaderDataArray<float> col9  ;
    ReaderDataArray<float> col10 ;
    ReaderDataArray<float> col11 ;
    ReaderDataArray<float> col12 ;
    ReaderDataArray<float> col13 ;
    ReaderDataArray<float> col14 ;
    ReaderDataArray<float> col15 ;
//    t.GetTree()->SetBranchStatus("*",0);
    if(true){
        nJets .set(t,"nJets");
        col1  .set(t,"col1" );
        col2  .set(t,"col2" );
        col3  .set(t,"col3" );
        col4  .set(t,"col4" );
        col6  .set(t,"col6" );
        col7  .set(t,"col7" );
        col8  .set(t,"col8" );
        col9  .set(t,"col9" );
        col10 .set(t,"col10");
        col11 .set(t,"col11");
        col12 .set(t,"col12");
        col13 .set(t,"col13");
        col14 .set(t,"col14");
        col15 .set(t,"col15");
    }


    Long64_t ievent = 0;
//    std::cout << t.GetTree()->GetBranchStatus("nJets")<<" "<< t.GetTree()->GetBranchStatus("col1")<<std::endl;
//    t.GetTree()->SetBranchStatus("nJets",1);
//    t.GetTree()->SetBranchStatus("col1",1);

    auto setEventRange=[&](unsigned int startEvent = 0, unsigned int numEvents = -1){
        if(numEvents >0 )
            t.SetEntriesRange(startEvent, startEvent  + numEvents );
        else
            t.SetEntriesRange(startEvent, startEvent-1 );
    };

    setEventRange(0,-1);
    while (t.Next()) {
//        std::cout << t.GetCurrentEntry()<<std::endl;
        std::cout<< ievent <<" : ";
        for(const auto& c : col1 ) std::cout <<c<<",";
        std::cout<<std::endl;


        ievent++;
    }


    f->Close();

}


void go(int test){
    switch (test){
    case 0:
        timer.Start();
        storeArray();
        timer.Stop();
        break;
    case 1:
        timer.Start();
        storeVector();
        timer.Stop();
        break;
    case 2:
        timer.Start();
        storeArrayWVector();
        timer.Stop();
        break;

    case 3:
        timer.Start();
        storeArrayWVectorBr();
        timer.Stop();
        break;
    case 4:
        timer.Start();
        readVectorAsArray();
        timer.Stop();
        break;
    case 5:
        timer.Start();
        readVectorAsTTreeReader();
        timer.Stop();
        break;
    case 6:
        timer.Start();
        readVectorAsTTreeReaderWithCopies();
        timer.Stop();
        break;

    }

    printf("RT=%6.2f s  Cpu=%6.2f s\n",timer.RealTime(),timer.CpuTime());



};


#endif

void storageBenchmark(int test){
    go(test);
}
