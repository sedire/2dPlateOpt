#ifndef _PLATE_2D_VARVECT_
#define _PLATE_2D_VARVECT_ 1

#include <vector>
#include "plate_var_types.h"
using std::vector;

template<class PL_NUM>
class VarVect
{
public:
	VarVect();
	VarVect( int _varNum );
	~VarVect();
	void setup( int _varNum );

	vector<PL_NUM> Nk;
	vector<PL_NUM> Nk1;
	vector<PL_NUM> d1N;
	vector<PL_NUM> d2N;

	vector<PL_NUM> Nk0;			//don't really know why we need these. some computational tricks, probably
	vector<PL_NUM> d1N0;
	vector<PL_NUM> d2N0;
};

template<class PL_NUM>
VarVect<PL_NUM>::VarVect() {}

template<class PL_NUM>
VarVect<PL_NUM>::~VarVect() {}

template<class PL_NUM>
VarVect<PL_NUM>::VarVect( int _varNum )
{
	Nk.resize( _varNum, 0.0 );
	Nk1.resize( _varNum, 0.0 );
	d1N.resize( _varNum, 0.0 );
	d2N.resize( _varNum, 0.0 );

	Nk0.resize( _varNum, 0.0 );
	d1N0.resize( _varNum, 0.0 );
	d2N0.resize( _varNum, 0.0 );
}

template<class PL_NUM>
void VarVect<PL_NUM>::setup( int _varNum )
{
	Nk.resize( _varNum, 0.0 );
	Nk1.resize( _varNum, 0.0 );
	d1N.resize( _varNum, 0.0 );
	d2N.resize( _varNum, 0.0 );

	Nk0.resize( _varNum, 0.0 );
	d1N0.resize( _varNum, 0.0 );
	d2N0.resize( _varNum, 0.0 );
}

#endif