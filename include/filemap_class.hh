#ifndef __FILEMAPCLASS__
#define __FILEMAPCLASS__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "filepair_class.hh"

class FileMap {

public:
  FileMap(const std::string & name, const std::map<std::string, FilePair > & filelistmap);
  ~FileMap();
  FileMap(const FileMap & filemap);

  std::string GetName() const;
  const double GetCrossSection(const std::string & experiment) const;
  const std::vector<std::string> GetFileList(const std::string & experiment) const;
  std::map<std::string, FilePair> GetFileListMap() const;

  void SetName(const std::string & newname);

  bool isExpInMap(const std::string & exp) const;
  void Print();

  FileMap & operator=(const FileMap & filemap);
  bool operator==(const FileMap & filemap);
  bool operator!=(const FileMap & filemap);

private:
  std::string mName;
  std::map<std::string, FilePair > mFileListMap;
  
};


#endif
