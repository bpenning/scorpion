#! /bin/bash
/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/bin/g++ \
-Iinclude \
-I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/include/python2.6 \
-I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/boost/1.47.0/include \
-I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5//include \
-I/home/hep/jm1103/LandS/include \
-I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include \
-L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/boost/1.47.0/lib \
-L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/lib \
-L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5//lib \
-L/home/hep/jm1103/LandS -L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/lib \
-g -o sarah.x  sarah.cpp -lboost_python -lpython2.6 -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMathMore -lMinuit -lRooFit -lRooFitCore -pthread -lm -ldl -llimitcode -rdynamic -pthread -m64  obj/alphat_functions.o obj/analysis_alphatb8.o obj/analysis_alphatb.o obj/analysis_alphat.o obj/analysis_atlas5.o obj/analysis_base.o obj/analysis_cms_singel_lepton_20fb.o obj/analysis_lp5.o obj/analysis_manager.o obj/analysis_monojet5.o obj/analysis_os5.o obj/analysis_ss5.o obj/analysis_ss8high.o obj/analysis_ss8low.o obj/analysis_ssb8.o obj/analysis_zerolepdoj.o obj/delphesdictionary.o obj/doj_functions.o obj/filemap_class.o obj/fileobject_class.o obj/filepair_class.o obj/jad_jet_class.o obj/jad_lepton_class.o obj/jad_particle_class.o obj/jBlockClasses.o obj/pair_info.o obj/TreeReader.o
