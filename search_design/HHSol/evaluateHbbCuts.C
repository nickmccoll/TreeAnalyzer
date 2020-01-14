#include <AnalysisSupport/Utilities/interface/Types.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdio.h>

#include "../background_estimation/predTools/DataCardMaker.h"

#if !defined(__CINT__) || defined(__MAKECINT__)

using namespace std;

    struct BTagCutInfo {
        static BTagCutInfo create(const string& variableLabel, const vector<double>& cutValues ) {
            BTagCutInfo out;
            out.variableLabel = variableLabel;
            out.cutValues = cutValues;
            return out;
        }

        vector<string> cutLabels = {"L","M","T","VT"};
        string variableLabel = "";
        vector<double> cutValues;
    };

    struct SignalInfo {
        static SignalInfo create(const string& modelName, const string& targetMassLabel,
                const double normalizationScaleFactor ) {
            SignalInfo out;
            out.modelName = modelName;
            out.targetMassLabel = targetMassLabel;
            out.normalizationScaleFactor = normalizationScaleFactor;
            return out;
        }
        string modelName = "";
        string targetMassLabel = "";
        double normalizationScaleFactor = 1.0;
    };

    struct FileData {
        vector<string> categoryLabels;
        vector<int> categoryCutIndicies;

        vector<string> categorySignalHistogramNames;
        vector<string> categoryBackgroundHistogramNames;

        vector<TH1*> categorySignalHistograms;
        vector<TH1*> categoryBackgroundHistograms;
    };

    struct CardCategoryData {
        string label = "";
        double signalCounts = 0;
        double backgroundCounts = 0;
    };

    struct CardData {
        string dataCardLabel = "";
        vector<CardCategoryData> categoryData;
    };

    string getCategoryCardFileName(const string& runLabel, const string& categoryLabel) {
        string cardFileName = "dataCard_";
        cardFileName += categoryLabel;
        cardFileName += "_" + runLabel;
        cardFileName += ".txt";
        return cardFileName;
    }

    string getCombinedCardFileName(const string& runLabel) {
        return getCategoryCardFileName("combined",runLabel);
    }

    class HistogramFileProcessor {
    public:

        void loadInputHistogram(const string& fileName) {
            clearInputFileIfNecessary();
            inputFile = new TFile(fileName.c_str(),"read");
        }

        void setRunData(const SignalInfo& signalInfo, const BTagCutInfo& btagCutInfo) {
            this->signalInfo = signalInfo;
            this->btagCutInfo = btagCutInfo;
        }

        FileData makeFileData() {
            FileData output = initializeDataWithCategoryLabels();

            loadBackgroundHistogramsFromFile(output);
            loadSignalHistogramsFromFile(output);

            return output;
        }

    private:
        FileData initializeDataWithCategoryLabels() {
            FileData output;
            for(unsigned int iP = 0; iP < purityCutLables.size(); iP++) {
                for(unsigned int iL = 0; iL < leptonCutLabels.size(); iL++) {
                    for(unsigned int iB = 0; iB < btagCutInfo.cutValues.size(); iB++) {
                        const string builtInCutLabel =
                                purityCutLables[iP] + "_" + leptonCutLabels[iL];
                        const string histogramNamePostFix =
                                builtInCutLabel + "_" + signalInfo.targetMassLabel
                                + "_" +  btagCutInfo.variableLabel;

                        const string categoryLabel =
                                builtInCutLabel + "_" + btagCutInfo.cutLabels[iB];
                        output.categoryLabels.push_back(categoryLabel);

                        output.categoryCutIndicies.push_back(iB);

                        const string signalHistogramName =
                                signalInfo.modelName + "_" + signalInfo.targetMassLabel + "_" + histogramNamePostFix;
                        output.categorySignalHistogramNames.push_back(signalHistogramName);

                        const string backgroundHistogramName = backgroundLabel
                                + "_" + histogramNamePostFix;
                        output.categoryBackgroundHistogramNames.push_back(backgroundHistogramName);
                    }
                }
            }
            return output;
        }

        void loadBackgroundHistogramsFromFile(FileData& data) {
            for(unsigned int iC = 0; iC < data.categoryBackgroundHistogramNames.size(); ++iC) {
                TH1* histogram = histogramFromFile(data.categoryBackgroundHistogramNames[iC]);
                data.categoryBackgroundHistograms.push_back(histogram);
            }
        }

        void loadSignalHistogramsFromFile(FileData& data) {
            for(unsigned int iC = 0; iC < data.categorySignalHistogramNames.size(); ++iC) {
                TH1* histogram = histogramFromFile(data.categorySignalHistogramNames[iC]);
                data.categorySignalHistograms.push_back(histogram);
            }
        }

        TH1 * histogramFromFile(const string& histogramName) {
            TH1 * histogram = 0;

            inputFile->GetObject(histogramName.c_str(),histogram);
            if(histogram == 0) {
                string errorString = "histogramFromFile ->  " + histogramName + " Not found.";
                throw std::invalid_argument(errorString);
            }
            histogram = (TH1*)histogram->Clone();
            histogram->SetDirectory(0);

            return histogram;
        }


    public:

        CardData makeCardData(const FileData& fileData) {
            CardData cardData;
            cardData.dataCardLabel = makeDataCardLabel();

            for(unsigned int iC = 0; iC < fileData.categoryLabels.size(); ++iC) {
                cardData.categoryData.push_back(getCategoryData(fileData,iC));
            }

            return cardData;
        }

    private:

        string makeDataCardLabel() {
            string out(signalInfo.modelName,0,1);

            out += "_" + signalInfo.targetMassLabel + "_" + btagCutInfo.variableLabel[0];

            for(const auto c : btagCutInfo.cutValues)
                out += "_" + makeBTagValueLabel(c);

            return out;
        }

        string makeBTagValueLabel(const float btagCutValue) {
            string cutValueAsString = ASTypes::flt2Str(btagCutValue);
            replace(cutValueAsString.begin(), cutValueAsString.end(), '.', 'p');
            return cutValueAsString;
        }

        CardCategoryData getCategoryData(const FileData& fileData, const int categoryIndex) {
            CardCategoryData categoryData;

            const int variableCutIndex = fileData.categoryCutIndicies[categoryIndex];

            const TH1* categorySignalHistogram = fileData.categorySignalHistograms[categoryIndex];
            const TH1* categoryBackgroundHistogram = fileData.categoryBackgroundHistograms[categoryIndex];

            const double categorySignalCount = getCategoryCount(categorySignalHistogram, variableCutIndex);
            double rawCategoryBackgroundCount = getCategoryCount(categoryBackgroundHistogram, variableCutIndex);
            const double categoryBackgroundCount = max(minimumBackgroundCount,rawCategoryBackgroundCount);

            categoryData.label = fileData.categoryLabels[categoryIndex];
            categoryData.backgroundCounts = categoryBackgroundCount;
            categoryData.signalCounts = categorySignalCount * signalInfo.normalizationScaleFactor;

            return categoryData;
        }

        double getCategoryCount(const TH1* histogram, const unsigned int variableCutIndex) {
            if(variableCutIndex >= btagCutInfo.cutValues.size()) {
                throw std::invalid_argument("getCategoryCount ->  Invalid cut index.");
            }

            double cutValue = btagCutInfo.cutValues[variableCutIndex];

            double categoryCount = getHistogramIntegralGreaterThan(histogram, cutValue);

            if(variableCutIndex + 1 != btagCutInfo.cutValues.size()){
                double nextCutValue = btagCutInfo.cutValues[variableCutIndex + 1];
                double nextCategoryCount = getHistogramIntegralGreaterThan(histogram, nextCutValue);
                categoryCount -= nextCategoryCount;
            }

            return categoryCount;
        }

        double getHistogramIntegralGreaterThan(const TH1* histogram, const double cutValue) {
            int cutBinIndex = histogram->FindFixBin(cutValue);
            return histogram->Integral(cutBinIndex,-1);
        }

    public:

        void clearFileData(FileData& input) {
            for(auto* histogram : input.categorySignalHistograms)
                delete histogram;
            for(auto* histogram : input.categoryBackgroundHistograms)
                delete histogram;
        }


        ~HistogramFileProcessor() {
            clearInputFileIfNecessary();
        }

    private:

        void clearInputFileIfNecessary(){
            if(inputFile){
                inputFile->Close();
                delete inputFile;
            }
        }

        TFile * inputFile = 0;

        BTagCutInfo btagCutInfo;
        SignalInfo signalInfo;

        const vector<string> leptonCutLabels = {"e","mu"};
        const vector<string> purityCutLables = {"LP","HP"};
        const string backgroundLabel = "bkg";
        const double minimumBackgroundCount = 0.01;
    };


    class DataCardCreator {
    public:
        void createCards(const CardData& cardData){
            for(const auto& category : cardData.categoryData) {
                makeDataCard(cardData.dataCardLabel,category);
            }
            combineDataCards(cardData);
        }

    private:

        void makeDataCard(const string& runLabel, const CardCategoryData& categoryData) {

            string cardFileName = getCategoryCardFileName(runLabel, categoryData.label);
            std::ofstream f (cardFileName,std::ios::out|std::ios::trunc);

            int observedDataValue = std::round(categoryData.backgroundCounts);

            f << "imax 1\n";
            f << "jmax  1\n";
            f << "kmax *\n";
            f << "bin " << categoryData.label << "\n";
            f << "observation  " << ASTypes::int2Str(observedDataValue) << "\n";
            f << "-------------------------\n";
            f << "bin\t";
            f << categoryData.label << "\t" << categoryData.label << "\n";
            f << "process" << "\t" << "signal" << "\t" << "bkg" << "\n";
            f << "process" << "\t" << "0" << "\t" << "1" << "\n";
            f << "rate" << "\t" << ASTypes::flt2Str(categoryData.signalCounts) << "\t"
                    << ASTypes::flt2Str(categoryData.backgroundCounts) << "\n";
            f << "signalN" << "\t" << "lnN " << "1.1" << "\t" << "-" << "\n";

            f.close();
        }

        void combineDataCards(const CardData& cardData) {
            string combinedFileNames =  getCombinedCardFileName(cardData.dataCardLabel);

            string shellExpression = "combineCards.py";

            for(auto& category : cardData.categoryData) {
                shellExpression += " " + getCategoryCardFileName(cardData.dataCardLabel, category.label);
            }

            shellExpression += " > " + getCombinedCardFileName(cardData.dataCardLabel);

            gSystem->Exec(shellExpression.c_str());
        }



    };


    class LimitEvaluator {
    public:

        double getLimit(const string& cardLabel) {
            runLimit(cardLabel);
            return extractExpectedLimitFromFile(cardLabel);
        }

    private:

        void runLimit(const string& cardLabel, float rMax = 2) {
            string limitString =  "combine -M AsymptoticLimits --run expected ";
            limitString +=  " --rMax " + ASTypes::flt2Str(rMax);
            limitString += " -n _" + cardLabel;
            limitString += " " + getCombinedCardFileName(cardLabel);

            string shellScript = "#/bin/bash \n";
            shellScript += ". ~/.zprofile \n";
            shellScript += limitString + " \n";

            gSystem->Exec(shellScript.c_str());
        }

        double extractExpectedLimitFromFile(const string& cardLabel){
            string fileName = "higgsCombine_" + cardLabel + ".AsymptoticLimits.mH120.root";

            TFile * file = new TFile(fileName.c_str(),"read");
            TTree * tree = 0;
            file->GetObject("limit",tree);

            const double limit = extractExpectedLimitFromTree(tree);

            file->Close();
            delete file;

            return limit;
        }

        double extractExpectedLimitFromTree(TTree * tree) {
            double output = -1;

            double limit;
            float quantileExpected;
            tree->SetBranchAddress("limit",&limit);
            tree->SetBranchAddress("quantileExpected",&quantileExpected);

            int eventNumber = 0;
            while(tree->GetEntry(eventNumber)){
                if(quantileExpected>0.49 && quantileExpected<0.51){
                    output = limit;
                    break;
                }
                ++eventNumber;
            }

            return output;
        }


    };

    class ResultsWriter {
    public:
        void addResult(const string& label, const double result) {
            labels.push_back(label);
            results.push_back(result);
        }

        void writeResultsToFile(const string& fileName) {

            std::ofstream f (fileName,std::ios::out|std::ios::trunc);

            for(unsigned int iR = 0; iR < results.size(); ++iR) {
                f << labels[iR] << "\t" << ASTypes::flt2Str(results[iR]) << "\n";
            }

            f.close();
        }

        void reset() {
            results.clear();
            labels.clear();
        }

    private:
        vector<double> results;
        vector<string> labels;
    };


    class HbbCutEvaluator {
    public:

        void setFileProcessor(HistogramFileProcessor* histogramProcessor) {
            this->histogramProcessor = histogramProcessor;
        }

        void evaluateHbbCut(const SignalInfo& signalInfo, const BTagCutInfo& btagCutInfo){
            histogramProcessor->setRunData(signalInfo,btagCutInfo);

            auto fileData = histogramProcessor->makeFileData();
            auto cardData = histogramProcessor->makeCardData(fileData);
            histogramProcessor->clearFileData(fileData);

            cardCreator.createCards(cardData);

            const double limit = limitEvaluator.getLimit(cardData.dataCardLabel);

            writer.addResult(cardData.dataCardLabel,limit);

        }

        void writeAllResultsToFile(const string fileLabel) {
            const string resultsFileName = "allresults_" + fileLabel + ".txt";
            writer.writeResultsToFile(resultsFileName);
            writer.reset();
        }

    private:
        HistogramFileProcessor * histogramProcessor =0;
        DataCardCreator cardCreator;
        LimitEvaluator limitEvaluator;
        ResultsWriter writer;

    };



#endif

void evaluateHbbCuts(string inputFileName = "../hSolTrees_hbbTaggingCutHistograms_2018.root",
        string label = "2018"){

//    BTagCutInfo deepCSVCutInfo;
//    deepCSVCutInfo.cutValues = {4,5,6};
//    deepCSVCutInfo.variableLabel = "dCSVCat";
//
//    BTagCutInfo deepAK8CutInfo;
//    deepAK8CutInfo.cutValues = {0.6,0.95};
//    deepAK8CutInfo.variableLabel = "ak8Tag";

    vector<SignalInfo> signalInfos;
    signalInfos.push_back(SignalInfo::create("graviton","m800",0.04));
    signalInfos.push_back(SignalInfo::create("graviton","m1000",0.02));
    signalInfos.push_back(SignalInfo::create("graviton","m1200",0.02));
    signalInfos.push_back(SignalInfo::create("graviton","m1400",0.015));
    signalInfos.push_back(SignalInfo::create("graviton","m2000",0.005));
    signalInfos.push_back(SignalInfo::create("graviton","m2500",0.006));
    signalInfos.push_back(SignalInfo::create("graviton","m3000",0.006));

    vector<BTagCutInfo> cutInfos;
//    cutInfos.push_back(BTagCutInfo::create("dCSVCat",{4,5,6}));
    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.9,0.98}));
    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.9,0.98}));
    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.98}));
//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.87,0.97}));

//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.6,0.9,0.98}));
//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.4,0.9,0.98}));
//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.92,0.98}));
//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.89,0.98}));
//
//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.9,0.96}));
//    cutInfos.push_back(BTagCutInfo::create("ak8Tag",{0.5,0.9,0.99}));

//    vector<double> highValues = {0.97,0.975,0.98};
//    vector<double> lowValues = {0.88,0.89,0.9,0.91,0.92};
//
//    for(const auto iL : lowValues) {
//        for(const auto iH : highValues) {
//            cutInfos.push_back(BTagCutInfo::create("ak8Tag",{iL,iH}));
//        }
//    }

//    vector<double> lowValues ={0.50,0.56,0.60,0.65,0.70,0.75,0.80,0.85};
//
//    for(const auto iL : lowValues) {
//            cutInfos.push_back(BTagCutInfo::create("ak8Tag",{iL,0.9,0.975}));
//    }



    HistogramFileProcessor processor;
    processor.loadInputHistogram(inputFileName);

    HbbCutEvaluator evaluator;
    evaluator.setFileProcessor(&processor);

    for(const auto& signal : signalInfos) {
        for(const auto& cut : cutInfos) {
            evaluator.evaluateHbbCut(signal,cut);
        }
    }

//    evaluator.evaluateHbbCut("m1000",deepCSVCutInfo);

    evaluator.writeAllResultsToFile(label);
}
