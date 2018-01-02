#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <typeinfo>
#include <stdexcept>

#include <TIterator.h>
#include <TString.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooStats/ModelConfig.h>
#include "../interface/utils.h"


using namespace std;



void utils::printRDH(RooAbsData *data) {
  std::vector<std::string> varnames, catnames;
  const RooArgSet *b0 = data->get();
  TIterator *iter = b0->createIterator();
  for (RooAbsArg *a = 0; (a = (RooAbsArg *)iter->Next()) != 0; ) {
    if (a->InheritsFrom("RooRealVar")) {
      varnames.push_back(a->GetName());
    } else if (a->InheritsFrom("RooCategory")) {
      catnames.push_back(a->GetName());
    }
  }
  delete iter;
  size_t nv = varnames.size(), nc = catnames.size();
  printf(" bin  ");
  for (size_t j = 0; j < nv; ++j) { printf("%16.16s  ", varnames[j].c_str()); }
  for (size_t j = 0; j < nc; ++j) { printf("%16.16s  ", catnames[j].c_str()); }
  printf("  weight\n");
  for (int i = 0, nb = data->numEntries(); i < nb; ++i) {
    const RooArgSet *bin = data->get(i);
    printf("%4d  ",i);
    for (size_t j = 0; j < nv; ++j) { printf("%16g  ",    bin->getRealValue(varnames[j].c_str())); }
    for (size_t j = 0; j < nc; ++j) { printf("%16.16s  ", bin->getCatLabel(catnames[j].c_str())); }
    printf("%12.7f\n", data->weight());
  }
}
