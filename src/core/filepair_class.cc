#include "filepair_class.hh"

FilePair::FilePair() {}

FilePair::FilePair(const double & xsec, const std::vector<std::string> & filelist) :
  mCrossSection(xsec),
  mFileList(filelist)
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

void FilePair::Print() {
  std::cout << "  cross-section: " << mCrossSection << std::endl;
  std::cout << "  file list: " << std::endl;
  for(std::vector<std::string>::const_iterator it = mFileList.begin(); it != mFileList.end(); it++) {
    std::cout << "    " << (*it) << std::endl;
  }
}
