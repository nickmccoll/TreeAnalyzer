#include <sstream>
int main(int argc, char* argv[])
{
    int treeInt = 0;
    float xSec = 0;
    float numEvent = 0;
    if(argc > 2){
        std::stringstream convertTreeInt(argv[2]);
        convertTreeInt >> treeInt;
    }
    if(argc > 5){
        std::stringstream convertxSec(argv[2]);
        convertxSec >> xSec;
        std::stringstream convertnumEvent(argv[2]);
        convertnumEvent >> numEvent;
    }
    switch(argc) {
        case 4:
            __MACRO__NAME__(argv[1], treeInt, argv[3]);
            return 0;
        case 6:
            __MACRO__NAME__(argv[1], treeInt, argv[3],xSec,numEvent);
            return 0;
        default:
            std::cout << "Cannot take: "<< (argc -1)<<" arguments!"<<std::endl;
    }
    return 0;
}
