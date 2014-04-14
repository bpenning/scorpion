#include "filepair_class.hh"

FilePair::FilePair() {}

FilePair::FilePair(const double & xsec, const std::vector<std::string> & filelist,const std::vector<int> & filereader) :
    mCrossSection(xsec),
    mFileList(filelist),
    mReaderList(filereader)
{
  /*  for(std::vector<std::int>::const_iterator read=reader.begin(); read!=reader.end(); read++) 
    {
	switch (*read)
	{
	    case 0:  
                cReader * read = &D2Parser();
		mReaderList.push_back(read);
	    default:
                cReader * read = &cParser();
		mReaderList.push_back(read);
	} 

    }
    //std::cout << "Initialised File Pair" << std::endl;
*/
}

FilePair::FilePair(const double & xsec, const std::vector<std::string> & filelist) :
    mCrossSection(xsec),
    mFileList(filelist),
    mReaderList()
{
    //std::cout << "Initialised File Pair" << std::endl;
}

FilePair::~FilePair() {}

const double FilePair::GetCrossSection() const {
    return mCrossSection;
}

const std::vector<std::string> FilePair::GetFileList() const {
    return mFileList;
}
const std::vector<int> FilePair::GetReaderList() const {
    return mReaderList;
}

void FilePair::Print() {
    std::cout << "  cross-section: " << mCrossSection << std::endl;
    std::cout << "  file list: " << std::endl;
    for(std::vector<std::string>::const_iterator it = mFileList.begin(); it != mFileList.end(); it++) {
	std::cout << "    " << (*it) << std::endl;
    }
    std::cout << "  file reader: " << std::endl;
    for(std::vector<int>::const_iterator it = mReaderList.begin(); it != mReaderList.end(); it++) {
	std::cout << "    " << (*it) << std::endl;
    }

}
