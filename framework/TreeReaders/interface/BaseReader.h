#ifndef TREEANALYZER_TREEREADERS_BASEREADER_H
#define TREEANALYZER_TREEREADERS_BASEREADER_H
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReaderWrapper.h"
namespace TAna{

typedef   ReaderData<double>           rd_double  ;
typedef   ReaderData<float >           rd_float   ;
typedef   ReaderData<int   >           rd_int     ;
typedef   ReaderData<ASTypes::int8  >  rd_int8    ;
typedef   ReaderData<ASTypes::int16 >  rd_int16   ;
typedef   ReaderData<ASTypes::size8 >  rd_size8   ;
typedef   ReaderData<ASTypes::size16>  rd_size16  ;
typedef   ReaderData<ASTypes::size  >  rd_size    ;
typedef   ReaderData<ASTypes::size64>  rd_size64  ;
typedef   ReaderDataArray<double>           ra_double  ;
typedef   ReaderDataArray<float >           ra_float   ;
typedef   ReaderDataArray<int   >           ra_int     ;
typedef   ReaderDataArray<ASTypes::int8  >  ra_int8    ;
typedef   ReaderDataArray<ASTypes::int16 >  ra_int16   ;
typedef   ReaderDataArray<ASTypes::size8 >  ra_size8   ;
typedef   ReaderDataArray<ASTypes::size16>  ra_size16  ;
typedef   ReaderDataArray<ASTypes::size  >  ra_size    ;
typedef   ReaderDataArray<ASTypes::size64>  ra_size64  ;


class BaseReader {
public:
	BaseReader(std::string readerName,std::string branchName) : readerName(readerName), branchName(branchName) {};
	virtual ~BaseReader() {}
	//external functions
	virtual void initialize(TreeReaderWrapper * wrapper) final;

	virtual void setup(TreeReaderWrapper * wrapper) = 0; //set active branches
	virtual void processVars() = 0; //build objects

	const std::string readerName ="";
	const std::string branchName ="";

};
}

#endif
