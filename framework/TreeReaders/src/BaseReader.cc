
#include "TreeReaders/interface/BaseReader.h"



namespace TAna{

void BaseReader::initialize(TreeReaderWrapper * wrapper) {
    std::cout << " ++  Attempting to initialize : " << readerName << " ("<<branchName<<")"                     << std::endl;
    std::cout << " ++  ";
    setup(wrapper);
    std::cout << std::endl << " ++  success"<< std::endl;
}
}
