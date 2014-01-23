#include "filemap_class.hh"

FileMap::FileMap(const std::string & name, const std::map<std::string, FilePair > & filelistmap) : 
  mName(name),
  mFileListMap(filelistmap)
{ 
  //std::cout << "Initialised File Map" << std::endl;
}

FileMap::~FileMap() {}

FileMap::FileMap(const FileMap & fileobj) {
  mName = fileobj.GetName();
  mFileListMap = fileobj.GetFileListMap();
}

std::string FileMap::GetName() const {
  return mName;
}

const double FileMap::GetCrossSection(const std::string & experiment) const {
  for( std::map<std::string, FilePair >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    if(ii->first == experiment) {
      return (ii->second).GetCrossSection();
    }
  }
  std::cout << "Experiment not found. No cross-section returned" << std::endl;

}

std::map<std::string, FilePair > FileMap::GetFileListMap() const {
  return mFileListMap;
}

const std::vector<std::string> FileMap::GetFileList(const std::string & experiment) const {

  for( std::map<std::string, FilePair >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    if(ii->first == experiment) {
      return (ii->second).GetFileList();
    }
  }
  std::cout << "Experiment not found. No files returned" << std::endl;
  
}

void FileMap::SetName(const std::string & newname) {
  mName = newname;
  return;
}

bool FileMap::isExpInMap(const std::string & exp) const {

  for( std::map<std::string, FilePair >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    if(ii->first == exp) {
      return true; 
    }
  }

  return false;

}

void FileMap::Print() {
  
  std::cout << "Process: " << mName << std::endl;
  for( std::map<std::string, FilePair >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    std::cout << ii->first << " : " << std::endl;
    FilePair temp(ii->second);
    temp.Print();
  }

  return;
}


bool FileMap::operator==(const FileMap & fileobj) {
  return(this->GetName() == fileobj.GetName());
}

bool FileMap::operator!=(const FileMap & fileobj) {
  return(this->GetName() != fileobj.GetName());
}

FileMap & FileMap::operator=(const FileMap & fileobj) {
  if(this != &fileobj) {
    this->SetName(fileobj.GetName());
  }

  return *this;
}
