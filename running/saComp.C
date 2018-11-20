#include <sstream>
#include "TFormula.h"
int main(int argc, char* argv[])
{
    int treeInt = 0;
    int randSeed = 0;
    float xSec = 0;
    float numEvent = 0;


    if(argc > 3){
        std::stringstream convertTreeInt(argv[2]);
        convertTreeInt >> treeInt;
        std::stringstream convertRandSeed(argv[3]);
        convertRandSeed >> randSeed;
    }
    if(argc > 6){
        TFormula xsec("xsec",argv[5]);
        xSec = xsec.Eval(0);
//        std::stringstream convertxSec(argv[4]);
//        convertxSec >> xSec;
        std::stringstream convertnumEvent(argv[6]);
        convertnumEvent >> numEvent;
    }
    switch(argc) {
        case 7:
            __MACRO__NAME__(argv[1], treeInt,randSeed, argv[4],xSec,numEvent);
            return 0;
        default:
            std::cout << "Cannot take: "<< (argc -1)<<" arguments!"<<std::endl;
    }
    return 0;
}
