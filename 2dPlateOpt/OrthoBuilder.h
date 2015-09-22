#ifndef _PLATE_2D_ORTHOBUILDER_
#define _PLATE_2D_ORTHOBUILDER_ 1

#include "plate_var_types.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <fstream>
#include "omp.h"
#include <Eigen/Eigen>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::bad_alloc;
using std::complex;
using namespace Eigen;

//class SolInfo
//{
//public:
//	SolInfo();
//	~SolInfo();
//
//	vector<PL_NUM> o;		//omega matrix to restore the solution
//
//	void setup( int _varNum );
//	void flushO();
//};

class OrthoBuilder
{
public:
	//vector<SolInfo> solInfoMap;
	OrthoBuilder( int _varNum, int Km );
	virtual ~OrthoBuilder();
	virtual void setParams();
	virtual void orthonorm( int y, PL_NUM NtoOrt[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] ) {};	//this version should orthonorm all 
	virtual void orthonorm( int y, Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>* NtoOrt ) {};	//this version should orthonorm all 
	virtual void buildSolution( vector<VarVect>* _mesh ) {};
	//virtual void flushO( int x );
	virtual void setOmegasZero();
	virtual void setInitVects( const vector<PL_NUM>& N1, const vector<PL_NUM>& N2, const vector<PL_NUM>& N3, const vector<PL_NUM>& N4, const vector<PL_NUM>& N5 );
	virtual void setOrthoDoneInfo( int y );
	virtual void resetOrthoDoneInfo();
	inline virtual void setNextSolVects( int n, const PL_NUM decompVect[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] ) {};
	inline virtual int checkOrtho( int n, 
									PL_NUM vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES], 
									//const Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>& vectSetOrtho,
									const Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>& vectSetOrig ) { return 1; };

	inline PL_NUM getInfNorm( PL_NUM* vect, int vectSize );
	Matrix<PL_NUM, Dynamic, Dynamic>* zi;
protected:
	const int varNum;
	const int Km;

	vector<bool> orthoDone;
	vector<vector<PL_NUM> > omega;
	PL_NUM omega2[EQ_NUM * NUMBER_OF_LINES / 2];
	PL_NUM Cx[EQ_NUM * NUMBER_OF_LINES / 2];
	PL_NUM Cx1[EQ_NUM * NUMBER_OF_LINES / 2];
private:
	OrthoBuilder();
};

class OrthoBuilderGSh : public OrthoBuilder			//use this
{
public:
	OrthoBuilderGSh( int _varNum, int Km );
	~OrthoBuilderGSh() {};

	void orthonorm( int y, PL_NUM NtoOrt[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] );	//this version should orthonorm all 
	void orthonorm( int y, Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>* NtoOrt );	//this version should orthonorm all 
	void buildSolution( vector<VarVect>* _mesh );

	inline void setNextSolVects( int n, const PL_NUM decompVect[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] );
	inline int checkOrtho( int n, 
							PL_NUM vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES], 
							//const Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>& vectSetOrtho,
							const Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>& vectSetOrig  );
private:
	OrthoBuilderGSh();

	HouseholderQR<Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> > qr;
	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES> qrQ;
	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> qrR;
};

#endif