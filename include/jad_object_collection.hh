#ifndef __JOBJECTCOLL__
#define __JOBJECTCOLL__

#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jad_particle_class.hh"

class jobjectCollection{
    public:
	jobjectCollection();
	jobjectCollection(std::vector <jjet> jjetvec);
	jobjectCollection(std::vector <jlepton> jleptonvec);
	jobjectCollection(std::vector <jparticle> jpartvec);

	jobjectCollection(std::vector <jjet> jjetvec,double ptcut,double etacut,bool upeta);
	jobjectCollection(std::vector <jlepton> jleptonvec,double ptcut, double etacut, bool upeta);
	jlepton GetLepton(int num) const;
	jjet GetJet(int num) const;
	jparticle GetParticle(int num) const;
	int GetEntries() const;
	~jobjectCollection();
    private:
	std::vector <jjet> mjetvec;
	std::vector <jlepton> mleptonvec;
	std::vector <jparticle> mparticlevec;
};
#endif
