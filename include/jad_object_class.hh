#ifndef __JOBJECTCLASS__
#define __JOBJECTCLASS__

#include "TLorentzVector.h"

class jobject : public TLorentzVector{
    public:
	enum jobjectType {
	    mlepton,
	    mjet,
	    mparticle,
	    mphoton,
	    mtrack
	};

	jobject();
	jobject(double Px, double Py, double Pz, double E, jobjectType jtype);
	~jobject();
	jobjectType gettype();
	bool operator<(const jobject & rhs) const;
	bool operator>(const jobject & rhs) const;
    private:
	jobjectType jtype;
};

#endif
