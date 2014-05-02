INCDIR=include

g++ -I ${INCDIR} -c `root-config --cflags` main_old.cc  
g++ -I ${INCDIR} -c `root-config --cflags` src/core/jad_jet_class.cc
g++ -I ${INCDIR} -c `root-config --cflags` src/core/jad_lepton_class.cc
g++ -I ${INCDIR} -c `root-config --cflags` src/core/jad_particle_class.cc
g++ -I ${INCDIR} -c `root-config --cflags` src/core/jad_ob_class.cc
g++ -I ${INCDIR} -c `root-config --cflags` testbox/DelphesClasses.cc 
g++ -I ${INCDIR} -c `root-config --cflags` testbox/ExRootTreeBranch.cc 
g++ -I ${INCDIR} -c `root-config --cflags` src/core/D3Reader.cc 
g++ -I ${INCDIR} -c `root-config --cflags` testbox/DelphesFactory.cc 

g++ -I ${INCDIR} -o test_main `root-config --libs` main_old.o jad_jet_class.o jad_ob_class.o jad_particle_class.o jad_lepton_class.o D3Reader.o libDelphes.so
#g++ -I ${INCDIR} -o test_main `root-config --libs` main_old.o jad_jet_class.o jad_ob_class.o jad_particle_class.o jad_lepton_class.o DelphesClasses.o D3Reader.o DelphesFactory.o ExRootTreeBranch.o 
rm *.o # DelphesFactory.o ExRootTreeBranch.o 
