#include <iostream>
#include <string>
#include <vector>
#include "analysis_alphatb.hh"
#include "analysis_atlas5.hh"
#include "fileobject_class.hh"
#include "filemap_class.hh"
#include "filepair_class.hh"
#include "analysis_manager.hh"
#include "analysis_base.hh"
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "boost/python/suite/indexing/map_indexing_suite.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(libjad_DelphesAnalysis)
{

    class_<AnalysisManager>("AnalysisManager")
	.def(init<const std::string &, const bool &>())
	//.def("run", &AnalysisManager::Run)
	.def("run", &AnalysisManager::RunMany)
	.def("add", &AnalysisManager::Add)
	.def("limit", &AnalysisManager::Limit)
	.def("write", &AnalysisManager::Write);

    //wrap the abstract/pure virtual base class
    class_<AnalysisBase, boost::noncopyable>("AnalysisBase", no_init)
	.def("getname", &AnalysisBase::GetName)
	.def("setname", &AnalysisBase::SetName)
	.def("setexp", &AnalysisBase::SetExperiment)
	.def("getexp", &AnalysisBase::GetExperiment);


    class_<AlphaTb, bases<AnalysisBase> >("AlphaTb", init<const std::string &, const std::string &, const unsigned int &>())
	.def(init<
		const std::string &, 
		const std::string &, 
		const unsigned int &, 
		const double &, 
		//const std::vector<int> &, 
		const std::vector<double> & 
		>())
	.def(init<
		const std::string &, 
		const std::string &, 
		const unsigned int &, 
		const double &, 
		const std::vector<double> &,
		const std::vector<double> &,
		const std::vector<int> &,
		const std::string &,
		const bool &
		>());
    ;
    class_<ATLAS5, bases<AnalysisBase> >("ATLAS5", init<const std::string &, const std::string &, const unsigned int &>())
	.def(init<
		const std::string &, 
		const std::string &, 
		const unsigned int &, 
		const double &, 
		//const std::vector<int> &, 
		const std::vector<double> & 
		>())
	.def(init<
		const std::string &, 
		const std::string &, 
		const unsigned int &, 
		const double &, 
		const std::vector<double> &,
		const std::vector<double> &,
		const std::vector<int> &,
		const std::string &,
		const bool &
		>());
    ;

    class_<FileObject>("FileObject", 
	    init<const std::string &, 
	    const double &, 
	    const std::map<std::string, std::vector<std::string> > &, 
	    const std::string &
	    > ())
	.def("printmap", &FileObject::PrintMap)
	.def("getname", &FileObject::GetName);

    class_<FilePair>("FilePair",
	    init<const double &,
	    const std::vector<std::string> &,
	    const int &
	    > ())
	.def("Print", &FilePair::Print);

    class_<FileMap>("FileMap",
	    init<const std::string &,
	    const std::map<std::string, FilePair> &
	    > ())
	.def("Print", &FileMap::Print);

    class_<std::map<std::string, FilePair> >("jad_FilePairMap")
	.def(map_indexing_suite<std::map<std::string, FilePair> >());

    class_<std::map<std::string, std::vector<std::string> > >("jad_FileObjectMap")
	.def(map_indexing_suite<std::map<std::string, std::vector<std::string> > >());

    class_<std::vector<std::string> >("jad_StringVector")
	.def(vector_indexing_suite<std::vector<std::string> >());

    class_<std::vector<double> >("jad_DoubleVector")
	.def(vector_indexing_suite<std::vector<double> >());

    class_<std::vector<int> >("jad_IntVector")
	.def(vector_indexing_suite<std::vector<int> >());

    class_<std::vector<FileObject> >("jad_FileObject")
	.def(vector_indexing_suite<std::vector<FileObject> >());

    class_<std::vector<FileMap> >("jad_FileMap")
	.def(vector_indexing_suite<std::vector<FileMap> >());

}
