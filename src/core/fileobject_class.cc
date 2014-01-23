#include "fileobject_class.hh"

//FileObject::FileObject(const std::string & name, const double & crosssection, const std::vector<std::string> & filelist, const std::string & misc) : 
//  mName(name), 
//  mCrossSection(crosssection), 
//  mFileList(filelist), 
//  mMisc(misc)
//{ 
//  std::cout << "Initialised Object" << std::endl;
//}

FileObject::FileObject(const std::string & name, const double & crosssection, const std::map<std::string, std::vector<std::string> > & filelistmap, const std::string & misc) : 
  mName(name), 
  mCrossSection(crosssection), 
  mFileListMap(filelistmap), 
  mMisc(misc)
{ 
  std::cout << "Initialised Object Map" << std::endl;
}

FileObject::~FileObject() {}

FileObject::FileObject(const FileObject & fileobj) {
  mName = fileobj.GetName();
  mCrossSection = fileobj.GetCrossSection();
  mMisc = fileobj.GetInfo();
  mFileListMap = fileobj.GetFileListMap();
}

std::string FileObject::GetName() const {
  return mName;
}

std::string FileObject::GetInfo() const {
  return mMisc;
}

double FileObject::GetCrossSection() const {
  return mCrossSection;
}

std::map<std::string, std::vector<std::string> > FileObject::GetFileListMap() const {
  return mFileListMap;
}

const std::vector<std::string> & FileObject::GetFileList(const std::string & experiment) const {

  for( std::map<std::string, std::vector<std::string> >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    if(ii->first == experiment) {
      return ii->second;
    }
  }
  
}

void FileObject::SetName(const std::string & newname) {
  mName = newname;
  return;
}

void FileObject::SetCrossSection(const double & newxsec) {
  mCrossSection = newxsec;
  return;
}

void FileObject::SetInfo(const std::string & newinfo) {
  mMisc = newinfo;
  return;
}

bool FileObject::isExpInMap(const std::string & exp) const {

  for( std::map<std::string, std::vector<std::string> >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    if(ii->first == exp) {
      return true; 
    }
  }

  return false;

}

void FileObject::PrintMap() {
  
  for( std::map<std::string, std::vector<std::string> >::const_iterator ii=mFileListMap.begin(); ii!=mFileListMap.end(); ii++) {
    std::cout << ii->first << " : " << std::endl;
    std::vector<std::string> tmp = ii->second;
    for(int i=0;i<tmp.size();i++) {
      std::cout << tmp.at(i) << std::endl;
    }
  }

  return;
}

bool FileObject::operator==(const FileObject & fileobj) {
  return(this->GetName() == fileobj.GetName());
}

bool FileObject::operator!=(const FileObject & fileobj) {
  return(this->GetName() != fileobj.GetName());
}

FileObject & FileObject::operator=(const FileObject & fileobj) {
  if(this != &fileobj) {
    this->SetName(fileobj.GetName());
  }

  return *this;
}
