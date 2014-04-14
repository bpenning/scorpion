#ifndef __FILEPAIRCLASS__
#define __FILEPAIRCLASS__

#include <iostream>
#include <string>
#include <vector>

class FilePair {

public:
  FilePair();
  FilePair(const double & xsec, const std::vector<std::string> & filelist,const std::vector<int> & filereader);
  FilePair(const double & xsec, const std::vector<std::string> & filelist);
  ~FilePair();
  
  const double GetCrossSection() const;
  const std::vector<std::string> GetFileList() const;
  const std::vector<int> GetReaderList() const;
  void Print();
  

private:
  double mCrossSection;
  std::vector<std::string> mFileList;
  std::vector<int> mReaderList;

};

#endif
