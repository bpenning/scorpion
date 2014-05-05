#include "reader_base.hh"
#include <iostream>
cParser::cParser(){};
cParser::~cParser(){};
cParser::cParser(readerType type){
    rtype=type;
}
void cParser::Print()
{
    switch(rtype)
    {
	case d2type:
	    std::cout << "D2 Reader" << std::endl;
	    break;
    }
    return;
}
