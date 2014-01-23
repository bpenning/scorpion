#ifndef __FILEOBJECTCLASS__
#define __FILEOBJECTCLASS__

#include <iostream>
#include <string>
#include <vector>
#include <map>

class FileObject {

public:
  //FileObject(const std::string & name, const double & crosssection, const std::vector<std::string> & filelist, const std::string & misc);
  FileObject(const std::string & name, const double & crosssection, const std::map<std::string, std::vector<std::string> > & filelistmap, const std::string & misc);
  ~FileObject();
  FileObject(const FileObject & fileobj);

  std::string GetName() const;
  std::string GetInfo() const;
  double GetCrossSection() const;
  const std::vector<std::string> & GetFileList(const std::string & experiment) const;
  std::map<std::string, std::vector<std::string> > GetFileListMap() const;

  void SetName(const std::string & newname);
  void SetCrossSection(const double & newxsec);
  void SetInfo(const std::string & newinfo);

  bool isExpInMap(const std::string & exp) const;
  void PrintMap();

  FileObject & operator=(const FileObject & fileobj);
  bool operator==(const FileObject & fileobj);
  bool operator!=(const FileObject & fileobj);

private:
  std::string mName;
  double mCrossSection;
  //std::vector<std::string> mFileList;
  std::map<std::string, std::vector<std::string> > mFileListMap;
  std::string mMisc;
  
};


#endif
