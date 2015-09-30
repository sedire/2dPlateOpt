#ifndef _PLATE_2D_SOLVER_
#define _PLATE_2D_SOLVER_ 1

#include <time.h>
#include "plate_var_types.h"
#include "OrthoBuilder.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include "hyperDual.h"

using std::cout;
using std::vector;
using std::ofstream;
using std::string;
using std::stringstream;
using std::setprecision;

template<class PL_NUM>
class Solver
{
public:
	Solver( N_PRES _E1, N_PRES _E2, PL_NUM _By0, N_PRES _ap, N_PRES _bp );
	~Solver();

	void setTask();
	void setCurrentParms( int _currentType, N_PRES _J0scale, const vector<PL_NUM>& currentParms );
	void setStressParms( int _stressType, const vector<N_PRES>& stressParms );
	void reset();

	void calc_nonlin_system_run_test( long  _x, long _t );

	void pre_step();
	PL_NUM do_step();
	void dump_sol();
	void dump_whole_sol( int var );
	void dump_check_sol();
	void dump_check_sol2D();

	N_PRES increaseTime();
	N_PRES getCurTime();

private:
	N_PRES E1;				//Young's modulus
	N_PRES E2;				//Young's modulus
	const N_PRES nu21;			//Poisson's ratio	
	const N_PRES nu23;			//Poisson's ratio	
	const N_PRES rho;				//composite's density
	const N_PRES G23;				//shear modulus

	const N_PRES B11;
	const N_PRES B22;
	const N_PRES B12;
	const N_PRES B66;
	PL_NUM By0;
	const PL_NUM By1;				// in considered boundary-value problem
	const PL_NUM By2;   
	const N_PRES beta;		//parameter at Newmark's time integration scheme

	const N_PRES mu;
	const N_PRES sigma_x;			//electric conductivity
	const N_PRES sigma_x_mu;
	const N_PRES sigma_y;
	const N_PRES sigma_y_mu;
	const N_PRES sigma_z;

	int currentType;
	PL_NUM J0;
	PL_NUM tauS;
	PL_NUM tauE;
	N_PRES J0scale;
	int stressType;
	N_PRES p0;		//constant mechanical load
	N_PRES tauP;
	N_PRES impRadFactor;
	N_PRES impRadSq;

	const N_PRES eps_0;
	const N_PRES eps_x;
	const N_PRES eps_x_0;

	const N_PRES hp;				//thickness of the plate
	N_PRES ap;				//width of the plate
	N_PRES bp;				//length of the plate

	const int Km;				//number of steps by space
	const int nx;				//number of lines in Method of Lines
	const int eq_num;				//number of equations in main system
	const int varNum;		//number of variables on the whole line. equals eq_num * number_of_lines

	const int maxNewtonIt;
	int newtonIt;

	const N_PRES dx;
	const N_PRES dy;

	const N_PRES dt;			//time step
	N_PRES cur_t;
	int curTimeStep;

	const N_PRES al;			//some weird var for normalization. It is said that this will improve the numerical scheme. must be equal to density

	vector<VarVect<PL_NUM> > mesh;		//2d mesh for solution.

	SparseMatrix<PL_NUM> matrA;
	SparseMatrix<PL_NUM> matrAprev;
	Matrix<PL_NUM, Dynamic, 1> vectF;
	Matrix<PL_NUM, Dynamic, 1> vectFprev;

	PL_NUM newmark_A[EQ_NUM * NUMBER_OF_LINES];
	PL_NUM newmark_B[EQ_NUM * NUMBER_OF_LINES];
	//PL_NUM decompVect[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES];
	PL_NUM decompVectOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES];

	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> decompVect;
	//Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> decompVectOrtho;

	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> tempPhi;
	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> tempPhi2;
	Matrix<PL_NUM, Dynamic, Dynamic>* Phi;

	//Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES, RowMajor> Ma;
	//bounds on the iterators of Matrix A
	//int lb[EQ_NUM * NUMBER_OF_LINES];
	//int rb[EQ_NUM * NUMBER_OF_LINES];

	const N_PRES rgkU;
	const N_PRES rgkV;
	N_PRES rgkC1, rgkC2, rgkC3, rgkC4;
	N_PRES rgkD21, rgkD32, rgkD31, rgkD43, rgkD42, rgkD41;

	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> rgkF1;
	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> rgkF2;
	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> rgkF3;
	Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1> rgkF4;

	OrthoBuilder<PL_NUM>* orthoBuilder;

	void rgkCalc( int thrNum, int hom, const vector<PL_NUM>& x, PL_NUM* x1 );
	void rgkCalc( const Matrix<PL_NUM, Dynamic, Dynamic>& x, Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>* x1 );	

	void calcNewmarkAB( int _x, int mode );		//don't really know why we need this stuff with mode need to FIX
	void calc_system( int _x );
	void walkthrough( int mode );
	void updateDerivs();

	int checkConv();
};

template<class PL_NUM>
Solver<PL_NUM>::Solver( N_PRES _E1, N_PRES _E2, PL_NUM _By0, N_PRES _ap, N_PRES _bp ):
	E1( _E1 ),
	E2( _E2 ),
	//E1( 102970000000 ),
	//E2( 7550000000 ),
	//E2( 102970000000 ),
	nu21( 0.3 ),
	nu23( 0.3 ),
	rho( 1594 ),
	G23( E2 / 2.0 / ( 1 + nu23 ) ),

	B11( E1 * E1 / ( E1 - nu21 * nu21 * E2 ) ),
	B22( E2 / ( 1 - nu21 * nu21 * E2 / E1 ) ),
	B12( nu21 * E2 * E1 / ( E1 - nu21 * nu21 * E2 ) ),
	B66( G23 ),

	By0( _By0 ),
	//By0( 1.0 ),
	By1( 2.0 * By0 ),
	By2( 0.0 ),

	beta( 0.25 ),

	mu( 0.00000125664 ),
	sigma_x( 39000 ),
	sigma_x_mu( sigma_x * mu ),
	sigma_y( sigma_x * 0.0001 ),
	//sigma_y( sigma_x ),
	sigma_y_mu( sigma_y * mu ),
	sigma_z( sigma_y ),

	currentType( currentExpSin ),
	J0( 0.0 ),
	tauS( 1.0 ),
	tauE( 1.0 ),
	stressType( stressCentered ),
	p0( 0.0 ),
	tauP( 1.0 ),
	impRadFactor( 1.0 ),
	impRadSq( 0.0 ),

	eps_0( 0.000000000008854 ),
	eps_x( 0.0000000002501502912 ),
	eps_x_0( eps_x - eps_0 ),

	hp( 0.0021 ),
	//ap( 0.1524 * 2.0 ),		//len in x-dir
	//ap( 0.1524 ),		//len in x-dir
	//bp( 0.1524 ),		//width in y-dir
	ap( _ap ),
	bp( _bp ),

	Km( NODES_ON_Y ),
	nx( NUMBER_OF_LINES ),
	eq_num( EQ_NUM ),
	varNum( NUMBER_OF_LINES * EQ_NUM ),

	maxNewtonIt( MAX_NEWTON_IT ),
	newtonIt( 0 ),

	dx( ap / ( nx + 1 ) ),
	dy( bp / ( Km - 1 ) ),

	dt( 0.0001 ),
	cur_t( 0.0001 ),
	curTimeStep( 1 ),

	al( 1.0 ),

	rgkU( 0.3 ),
	rgkV( 0.6 ),

	Phi( 0 ),

	orthoBuilder( 0 )
{
	setTask();
}

template<class PL_NUM>
void Solver<PL_NUM>::setCurrentParms( int _currentType, N_PRES _J0scale, const vector<PL_NUM>& currentParms )
{
	J0scale = _J0scale;
	currentType = _currentType;

	if( _currentType == currentConst && currentParms.size() == 1 )
	{
		J0 = currentParms[0] * J0scale;
	}
	else if( _currentType == currentSin && currentParms.size() == 2 )
	{
		J0 = currentParms[0] * J0scale;
		tauS = currentParms[1];
	}
	else if( _currentType == currentExpSin && currentParms.size() == 3 )
	{
		J0 = currentParms[0] * J0scale;
		tauS = currentParms[1];
		tauE = currentParms[2];
	}
	else
	{
		cout << "ERROR: setCurrentParms: current type is unsupported or wrong parms vector size\n";
		return;
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::setStressParms( int _stressType, const vector<N_PRES>& stressParms )
{
	stressType = _stressType;

	if( _stressType == stressWhole && stressParms.size() == 1 )
	{
		p0 = stressParms[0];
	}
	else if( _stressType == stressSin && stressParms.size() == 2 )
	{
		p0 = stressParms[0];
		tauP = stressParms[1];

		//cout << " --- " << p0 << " " << tauP << endl;
	}
	else if( _stressType == stressCentered && stressParms.size() == 3 )
	{
		p0 = stressParms[0];
		tauP = stressParms[1];
		impRadFactor = stressParms[2];
		impRadSq = ap * ap / impRadFactor / impRadFactor;
	}
	else
	{
		cout << "ERROR: setStressParms: stress type is unsupported or wrong parms vector size\n";
		return;
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::reset()		//TODO: UNTESTED!
{
	ofstream of1( "test_sol.txt" );
	of1.close();

	cur_t = dt;
	curTimeStep = 1;

	if( orthoBuilder != 0 )
	{
		delete orthoBuilder;
	}
	orthoBuilder = new OrthoBuilderGSh<PL_NUM>( varNum, Km );

	mesh.resize( 0 );
	mesh.resize( Km );
	for( int i = 0; i < mesh.size(); ++i ){
		mesh[i].setup( varNum );
	}

	for( int i = 0; i < ABM_STAGE_NUM; ++i )
	{
		Phi[i] = Matrix<PL_NUM, Dynamic, Dynamic>::Zero( EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1 );
	}
}

template<class PL_NUM>
Solver<PL_NUM>::~Solver()
{
	if( orthoBuilder != 0 )
	{
		delete( orthoBuilder );
	}
	if( Phi != 0 )
	{
		delete[] Phi;
	}
}

template<class PL_NUM>
N_PRES Solver<PL_NUM>::increaseTime()
{
	cur_t += dt;
	++curTimeStep;
	return cur_t;
}

template<class PL_NUM>
N_PRES Solver<PL_NUM>::getCurTime()
{
	return cur_t;
}

template<class PL_NUM>
void Solver<PL_NUM>::setTask()
{
	ofstream of1( "test_sol.txt" );
	of1.close();

	rgkC2 = ( 2.0 * rgkV - 1.0 ) / 12.0 / rgkU / ( rgkV - rgkU ) / ( 1.0 - rgkU );
	rgkC3 = ( 1.0 - 2.0 * rgkU ) / 12.0 / rgkV / ( rgkV - rgkU ) / ( 1.0 - rgkV );
	rgkC4 = ( 6.0 * rgkU * rgkV - 4.0 * rgkU - 4.0 * rgkV + 3.0 ) / 12.0 / ( 1.0 - rgkU ) / ( 1.0 - rgkV );
	rgkC1 = 1.0 - rgkC2 - rgkC3 - rgkC4;
	rgkD21 = rgkU;
	rgkD32 = 1.0 / 24.0 / rgkC3 / rgkU / ( 1.0 - rgkV );
	rgkD31 = rgkV - rgkD32;
	rgkD43 = ( 1.0 - 2.0 * rgkU ) / 12.0 / rgkC4 / rgkV / ( rgkV - rgkU );
	rgkD42 = -( rgkV * ( 4.0 * rgkV - 5.0 ) - rgkU + 2.0 ) / 24.0 / rgkC4 / rgkU / ( rgkV - rgkU ) / ( 1.0 - rgkV );
	rgkD41 = 1.0 - rgkD42 - rgkD43;

	orthoBuilder = new OrthoBuilderGSh<PL_NUM>( varNum, Km );

	mesh.resize( Km );
	for( int i = 0; i < mesh.size(); ++i ){
		mesh[i].setup( varNum );
	}

	matrA = SparseMatrix<PL_NUM>( EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES );
	matrAprev = SparseMatrix<PL_NUM>( EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES );
	vectF = Matrix<PL_NUM, Dynamic, 1>::Zero( EQ_NUM * NUMBER_OF_LINES, 1 );
	vectFprev = Matrix<PL_NUM, Dynamic, 1>::Zero( EQ_NUM * NUMBER_OF_LINES, 1 );

	//filling A
	int i = 0;
	int r = i + 1;
	int rr = i + 2;
	int j = i - 1;
	int jj = i - 2;
	
	//for the left line:
	matrA.insert( 0 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrA.insert( 0 + i * eq_num, 2 + i * eq_num ) = 0.0;

	matrA.insert( 1 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrA.insert( 1 + i * eq_num, 3 + i * eq_num ) = 0.0;

	matrA.insert( 2 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 3 + r * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 9 + r * eq_num ) = 0.0;

	matrA.insert( 3 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 2 + r * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 3 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrA.insert( 4 + i * eq_num, 5 + i * eq_num ) = 0.0;

	matrA.insert( 5 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrA.insert( 5 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 5 + i * eq_num, 6 + i * eq_num ) = 0.0;

	matrA.insert( 6 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 5 + r * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 7 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 9 + i * eq_num ) = 0.0;
	
	matrA.insert( 7 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + rr * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + ( rr + 1 ) * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 5 + r * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 6 + r * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 8 + r * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 9 + r * eq_num ) = 0.0;

	matrA.insert( 8 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrA.insert( 8 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrA.insert( 8 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrA.insert( 8 + i * eq_num, 9 + r * eq_num ) = 0.0;

	matrA.insert( 9 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrA.insert( 9 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 9 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 9 + i * eq_num, 9 + i * eq_num ) = 0.0;

	i = nx - 1;
	r = i + 1;
	rr = i + 2;
	j = i - 1;
	jj = i - 2;
	
	//for the right line:
	matrA.insert( 0 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrA.insert( 0 + i * eq_num, 2 + i * eq_num ) = 0.0;

	matrA.insert( 1 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrA.insert( 1 + i * eq_num, 3 + i * eq_num ) = 0.0;

	matrA.insert( 2 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 3 + j * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrA.insert( 2 + i * eq_num, 9 + j * eq_num ) = 0.0;

	matrA.insert( 3 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 2 + j * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 3 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 3 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrA.insert( 4 + i * eq_num, 5 + i * eq_num ) = 0.0;

	matrA.insert( 5 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrA.insert( 5 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 5 + i * eq_num, 6 + i * eq_num ) = 0.0;

	matrA.insert( 6 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 5 + j * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 7 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 6 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrA.insert( 7 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + ( jj - 1 ) * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + jj * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 5 + j * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 6 + j * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 8 + j * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 9 + j * eq_num ) = 0.0;
	matrA.insert( 7 + i * eq_num, 9 + i * eq_num ) = 0.0;
	
	matrA.insert( 8 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrA.insert( 8 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrA.insert( 8 + i * eq_num, 9 + j * eq_num ) = 0.0;
	matrA.insert( 8 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrA.insert( 9 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrA.insert( 9 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrA.insert( 9 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrA.insert( 9 + i * eq_num, 9 + i * eq_num ) = 0.0;

	for( int i = 1; i < nx - 1; ++i )
	{
		r = i + 1;
		rr = i + 2;
		j = i - 1;
		jj = i - 2;

		matrA.insert( 0 + i * eq_num, 1 + r * eq_num ) = 0.0;
		matrA.insert( 0 + i * eq_num, 1 + j * eq_num ) = 0.0;
		matrA.insert( 0 + i * eq_num, 2 + i * eq_num ) = 0.0;

		matrA.insert( 1 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrA.insert( 1 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrA.insert( 1 + i * eq_num, 3 + i * eq_num ) = 0.0;

		matrA.insert( 2 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 0 + i * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 1 + r * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 1 + j * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 3 + r * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 3 + j * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 4 + j * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 9 + i * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 9 + r * eq_num ) = 0.0;
		matrA.insert( 2 + i * eq_num, 9 + j * eq_num ) = 0.0;

		matrA.insert( 3 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 1 + i * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 2 + r * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 2 + j * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 3 + i * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 5 + i * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrA.insert( 3 + i * eq_num, 9 + i * eq_num ) = 0.0;

		matrA.insert( 4 + i * eq_num, 5 + i * eq_num ) = 0.0;

		matrA.insert( 5 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrA.insert( 5 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrA.insert( 5 + i * eq_num, 4 + j * eq_num ) = 0.0;
		matrA.insert( 5 + i * eq_num, 6 + i * eq_num ) = 0.0;

		matrA.insert( 6 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 4 + j * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 5 + i * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 5 + r * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 5 + j * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 6 + i * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 7 + i * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrA.insert( 6 + i * eq_num, 9 + i * eq_num ) = 0.0;

		matrA.insert( 7 + i * eq_num, 1 + i * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 4 + j * eq_num ) = 0.0;
		if( i != nx - 2 )
		{
			matrA.insert( 7 + i * eq_num, 4 + rr * eq_num ) = 0.0;
		}
		if( i != 1 )
		{
			matrA.insert( 7 + i * eq_num, 4 + jj * eq_num ) = 0.0;
		}
		matrA.insert( 7 + i * eq_num, 5 + i * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 5 + r * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 5 + j * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 6 + i * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 6 + r * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 6 + j * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 8 + r * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 8 + j * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 9 + i * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 9 + r * eq_num ) = 0.0;
		matrA.insert( 7 + i * eq_num, 9 + j * eq_num ) = 0.0;

		matrA.insert( 8 + i * eq_num, 0 + i * eq_num ) = 0.0;
		matrA.insert( 8 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrA.insert( 8 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrA.insert( 8 + i * eq_num, 9 + i * eq_num ) = 0.0;
		matrA.insert( 8 + i * eq_num, 9 + r * eq_num ) = 0.0;
		matrA.insert( 8 + i * eq_num, 9 + j * eq_num ) = 0.0;

		matrA.insert( 9 + i * eq_num, 1 + i * eq_num ) = 0.0;
		matrA.insert( 9 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrA.insert( 9 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrA.insert( 9 + i * eq_num, 9 + i * eq_num ) = 0.0;
	}

	//filling A
	i = 0;
	r = i + 1;
	rr = i + 2;
	j = i - 1;
	jj = i - 2;
	
	//for the left line:
	matrAprev.insert( 0 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrAprev.insert( 0 + i * eq_num, 2 + i * eq_num ) = 0.0;

	matrAprev.insert( 1 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrAprev.insert( 1 + i * eq_num, 3 + i * eq_num ) = 0.0;

	matrAprev.insert( 2 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 3 + r * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 9 + r * eq_num ) = 0.0;

	matrAprev.insert( 3 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 2 + r * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 3 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrAprev.insert( 4 + i * eq_num, 5 + i * eq_num ) = 0.0;

	matrAprev.insert( 5 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrAprev.insert( 5 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 5 + i * eq_num, 6 + i * eq_num ) = 0.0;

	matrAprev.insert( 6 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 5 + r * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 7 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 9 + i * eq_num ) = 0.0;
	
	matrAprev.insert( 7 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 1 + r * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + r * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + rr * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + ( rr + 1 ) * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 5 + r * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 6 + r * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 8 + r * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 9 + r * eq_num ) = 0.0;

	matrAprev.insert( 8 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrAprev.insert( 8 + i * eq_num, 0 + r * eq_num ) = 0.0;
	matrAprev.insert( 8 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrAprev.insert( 8 + i * eq_num, 9 + r * eq_num ) = 0.0;

	matrAprev.insert( 9 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrAprev.insert( 9 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 9 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 9 + i * eq_num, 9 + i * eq_num ) = 0.0;

	i = nx - 1;
	r = i + 1;
	rr = i + 2;
	j = i - 1;
	jj = i - 2;
	
	//for the right line:
	matrAprev.insert( 0 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrAprev.insert( 0 + i * eq_num, 2 + i * eq_num ) = 0.0;

	matrAprev.insert( 1 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrAprev.insert( 1 + i * eq_num, 3 + i * eq_num ) = 0.0;

	matrAprev.insert( 2 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 3 + j * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 9 + i * eq_num ) = 0.0;
	matrAprev.insert( 2 + i * eq_num, 9 + j * eq_num ) = 0.0;

	matrAprev.insert( 3 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 2 + j * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 3 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 3 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrAprev.insert( 4 + i * eq_num, 5 + i * eq_num ) = 0.0;

	matrAprev.insert( 5 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrAprev.insert( 5 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 5 + i * eq_num, 6 + i * eq_num ) = 0.0;

	matrAprev.insert( 6 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 5 + j * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 7 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 6 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrAprev.insert( 7 + i * eq_num, 1 + j * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + ( jj - 1 ) * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + jj * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + j * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 5 + j * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 5 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 6 + j * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 6 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 8 + j * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 9 + j * eq_num ) = 0.0;
	matrAprev.insert( 7 + i * eq_num, 9 + i * eq_num ) = 0.0;
	
	matrAprev.insert( 8 + i * eq_num, 0 + j * eq_num ) = 0.0;
	matrAprev.insert( 8 + i * eq_num, 0 + i * eq_num ) = 0.0;
	matrAprev.insert( 8 + i * eq_num, 9 + j * eq_num ) = 0.0;
	matrAprev.insert( 8 + i * eq_num, 9 + i * eq_num ) = 0.0;

	matrAprev.insert( 9 + i * eq_num, 1 + i * eq_num ) = 0.0;
	matrAprev.insert( 9 + i * eq_num, 4 + i * eq_num ) = 0.0;
	matrAprev.insert( 9 + i * eq_num, 8 + i * eq_num ) = 0.0;
	matrAprev.insert( 9 + i * eq_num, 9 + i * eq_num ) = 0.0;

	for( int i = 1; i < nx - 1; ++i )
	{
		r = i + 1;
		rr = i + 2;
		j = i - 1;
		jj = i - 2;

		matrAprev.insert( 0 + i * eq_num, 1 + r * eq_num ) = 0.0;
		matrAprev.insert( 0 + i * eq_num, 1 + j * eq_num ) = 0.0;
		matrAprev.insert( 0 + i * eq_num, 2 + i * eq_num ) = 0.0;

		matrAprev.insert( 1 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrAprev.insert( 1 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrAprev.insert( 1 + i * eq_num, 3 + i * eq_num ) = 0.0;

		matrAprev.insert( 2 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 0 + i * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 1 + r * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 1 + j * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 3 + r * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 3 + j * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 4 + j * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 9 + i * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 9 + r * eq_num ) = 0.0;
		matrAprev.insert( 2 + i * eq_num, 9 + j * eq_num ) = 0.0;

		matrAprev.insert( 3 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 1 + i * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 2 + r * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 2 + j * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 3 + i * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 5 + i * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrAprev.insert( 3 + i * eq_num, 9 + i * eq_num ) = 0.0;

		matrAprev.insert( 4 + i * eq_num, 5 + i * eq_num ) = 0.0;

		matrAprev.insert( 5 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrAprev.insert( 5 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrAprev.insert( 5 + i * eq_num, 4 + j * eq_num ) = 0.0;
		matrAprev.insert( 5 + i * eq_num, 6 + i * eq_num ) = 0.0;

		matrAprev.insert( 6 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 4 + j * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 5 + i * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 5 + r * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 5 + j * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 6 + i * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 7 + i * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrAprev.insert( 6 + i * eq_num, 9 + i * eq_num ) = 0.0;

		matrAprev.insert( 7 + i * eq_num, 1 + i * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 4 + r * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 4 + j * eq_num ) = 0.0;
		if( i != nx - 2 )
		{
			matrAprev.insert( 7 + i * eq_num, 4 + rr * eq_num ) = 0.0;
		}
		if( i != 1 )
		{
			matrAprev.insert( 7 + i * eq_num, 4 + jj * eq_num ) = 0.0;
		}
		matrAprev.insert( 7 + i * eq_num, 5 + i * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 5 + r * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 5 + j * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 6 + i * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 6 + r * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 6 + j * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 8 + r * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 8 + j * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 9 + i * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 9 + r * eq_num ) = 0.0;
		matrAprev.insert( 7 + i * eq_num, 9 + j * eq_num ) = 0.0;

		matrAprev.insert( 8 + i * eq_num, 0 + i * eq_num ) = 0.0;
		matrAprev.insert( 8 + i * eq_num, 0 + r * eq_num ) = 0.0;
		matrAprev.insert( 8 + i * eq_num, 0 + j * eq_num ) = 0.0;
		matrAprev.insert( 8 + i * eq_num, 9 + i * eq_num ) = 0.0;
		matrAprev.insert( 8 + i * eq_num, 9 + r * eq_num ) = 0.0;
		matrAprev.insert( 8 + i * eq_num, 9 + j * eq_num ) = 0.0;

		matrAprev.insert( 9 + i * eq_num, 1 + i * eq_num ) = 0.0;
		matrAprev.insert( 9 + i * eq_num, 4 + i * eq_num ) = 0.0;
		matrAprev.insert( 9 + i * eq_num, 8 + i * eq_num ) = 0.0;
		matrAprev.insert( 9 + i * eq_num, 9 + i * eq_num ) = 0.0;
	}

	Phi = new Matrix<PL_NUM, Dynamic, Dynamic>[ABM_STAGE_NUM];
	for( int i = 0; i < ABM_STAGE_NUM; ++i )
	{
		Phi[i] = Matrix<PL_NUM, Dynamic, Dynamic>::Zero( EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1 );
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::calcNewmarkAB( int _x, int mode )
{
	if( mode == 0 )
	{
		for( int i = 0; i < varNum; ++i)
		{
			newmark_A[i] = -mesh[_x].Nk0[i] / beta / dt / dt - mesh[_x].d1N0[i] / beta / dt
							- ( 0.5 - beta ) / beta * mesh[_x].d2N0[i];
			newmark_B[i] = -0.5 * mesh[_x].Nk0[i] / beta / dt + ( 1.0 - 0.5 / beta ) * mesh[_x].d1N0[i]
							- 0.5 * dt * ( 1.0 - 1.0 / beta * ( 0.5 - beta ) ) * mesh[_x].d2N0[i];
		}
	}
	else
	{
		for( int i = 0; i < varNum; ++i)
		{
			newmark_A[i] = -mesh[_x].Nk1[i] / beta / dt / dt - mesh[_x].d1N[i] / beta / dt
							- ( 0.5 - beta ) / beta * mesh[_x].d2N[i];
			newmark_B[i] = -0.5 * mesh[_x].Nk1[i] / beta / dt + ( 1.0 - 0.5 / beta ) * mesh[_x].d1N[i]
							- 0.5 * dt * ( 1.0 - 1.0 / beta * ( 0.5 - beta ) ) * mesh[_x].d2N[i];
		}
	}
}

template<>
void Solver<HPD<N_PRES, 2> >::calc_system( int _x );

template<class PL_NUM>
void Solver<PL_NUM>::calc_system( int _x )
{
	N_PRES h = hp;
	N_PRES Btdt = 2 * dt * beta;

	PL_NUM Jx = 0.0;
	if( currentType == currentConst )
	{
		Jx = J0;
	}
	else if( currentType == currentSin )
	{
		Jx = J0 * sin( M_PI / tauS * cur_t );
	}
	else if( currentType == currentExpSin )
	{
		Jx = J0 * exp( -cur_t / tauE ) * sin( M_PI / tauS * cur_t );
	}

	N_PRES Pimp = 0.0;
	if( stressType == stressWhole )
	{
		Pimp = p0;
	}
	else if( stressType == stressSin )
	{
		Pimp = p0 * sin( M_PI / tauP * cur_t );
	}

	//strip load
	//PL_NUM cur_X = _x * dy - bp / 2.0;
	//if( fabs( cur_X ) <= h / 10.0 && cur_t + dt <= tauP )
	//{
	//	Pimp = p0 * sqrt( 1.0 - cur_X / h * 10.0 * cur_X / h * 10.0 ) * sin( _MMM_PI * ( cur_t + dt ) / tauP );
	//}

	int i = 0;
	int r = i + 1;
	int rr = i + 2;
	int j = i - 1;
	int jj = i - 2;
	
	//for the left line:
	matrA.coeffRef( 0 + i * eq_num, 1 + r * eq_num ) = -1.0 / ( 2.0 * dx );
	matrA.coeffRef( 0 + i * eq_num, 2 + i * eq_num ) = 1.0 / ( h  *B66 );

	matrA.coeffRef( 1 + i * eq_num, 0 + r * eq_num ) = -B12 / ( B22 * 2.0 * dx );
	matrA.coeffRef( 1 + i * eq_num, 3 + i * eq_num ) = 1.0 / ( h * B22 );

	matrA.coeffRef( 2 + i * eq_num, 0 + r * eq_num ) = ( ( 3.0 * B12 * B12 - 4.0 * B11 * B22 ) * h ) / ( 4.0 * B22 * dx * dx );
	matrA.coeffRef( 2 + i * eq_num, 0 + i * eq_num ) = ( 2 * B11 * h ) / dx / dx - ( B12 * B12 * h ) / ( B22 * dx * dx ) + ( h * ( 8.0 * rho + By1 * By1 * dt * sigma_z ) ) 
				/ ( 8.0 * beta * dt * dt );
	matrA.coeffRef( 2 + i * eq_num, 1 + r * eq_num ) = ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt * dx );
	matrA.coeffRef( 2 + i * eq_num, 3 + r * eq_num ) = -B12 / ( B22 * 2.0 * dx );
	matrA.coeffRef( 2 + i * eq_num, 4 + r * eq_num ) = -( ( By1 * eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] ) / ( 8.0 * beta * dt * dx ) );
	matrA.coeffRef( 2 + i * eq_num, 8 + i * eq_num ) = -( ( eps_x_0 * h * ( 2.0 * beta * dt * ( newmark_B[4 + r * eq_num] * By1 - 2.0 * newmark_B[1 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] )
				- 2.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + r * eq_num] + By1 * mesh[_x].Nk[4 + r * eq_num] ) ) / ( 8.0 * beta * dt * dx ) );
	matrA.coeffRef( 2 + i * eq_num, 9 + i * eq_num ) = ( 8.0 * beta * dt * h * ( -2.0 * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] ) + 3.0 * eps_x_0 * h * mu 
				* ( 2.0 * beta * newmark_B[1 + r * eq_num] * dt + mesh[_x].Nk[1 + r * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 12.0 * beta * dt * dx * mu );
	matrA.coeffRef( 2 + i * eq_num, 9 + r * eq_num ) = ( 2.0 * h * mesh[_x].Nk[9 + i * eq_num] ) / ( 3.0 * dx * mu );

	matrA.coeffRef( 3 + i * eq_num, 0 + r * eq_num ) = -( ( B12 * eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * B22 * dt * dx ) );
	matrA.coeffRef( 3 + i * eq_num, 1 + r * eq_num ) = -( ( B66 * h ) / ( 4 * dx * dx ) );
	matrA.coeffRef( 3 + i * eq_num, 1 + i * eq_num ) = ( B66 * h ) / dx / dx + ( h * ( 2.0 * rho + dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) )
				/ ( 2.0 * beta * dt * dt);
	matrA.coeffRef( 3 + i * eq_num, 2 + r * eq_num ) = -1.0 / ( 2.0 * dx );
	matrA.coeffRef( 3 + i * eq_num, 3 + i * eq_num ) = ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * beta * B22 * dt );
	matrA.coeffRef( 3 + i * eq_num, 4 + i * eq_num ) = -( ( By1 * h * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 4.0 * beta * dt ) );
	matrA.coeffRef( 3 + i * eq_num, 5 + i * eq_num ) = -( ( By1 * eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt ) );
	matrA.coeffRef( 3 + i * eq_num, 8 + i * eq_num ) = -( ( eps_x_0 * ( 2.0 * beta * dt * ( B22 * newmark_B[5 + i * eq_num] * By1 * dx * h - 2.0 * newmark_B[3 + i * eq_num] * dx
				* mesh[_x].Nk[9 + i * eq_num]
				+ B12 * newmark_B[0 + r * eq_num] * h * mesh[_x].Nk[9 + i * eq_num] ) + B12 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[0 + r * eq_num] - 2.0 * dx
				* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[3 + i * eq_num] 
				+ B22 * By1 * dx * h * mesh[_x].Nk[5 + i * eq_num] ) ) / ( 4.0 * beta * B22 * dt * dx ) ) + h * mesh[_x].Nk[9 + i * eq_num] * sigma_x;
	matrA.coeffRef( 3 + i * eq_num, 9 + i * eq_num ) = ( 1.0 / ( 4.0 * beta * B22 * dt * dx ) ) * ( ( -B12 ) * eps_x_0 * h * mesh[_x].Nk[0 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num] 
				+ 2.0 * dx * eps_x_0 * mesh[_x].Nk[3 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + B22 * dx * h * ( 4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num]
				- By1 * mesh[_x].Nk[4 + i * eq_num] ) * sigma_x 
				+ 2.0 * beta * dt * ( eps_x_0 * ( 2.0 * newmark_B[3 + i * eq_num] * dx - B12 * newmark_B[0 + r * eq_num] * h ) * mesh[_x].Nk[8 + i * eq_num]
				+ B22 * dx * h * ( 2.0 * Jx - newmark_B[4 + i * eq_num] * By1 * sigma_x + 4.0 * newmark_B[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x
				+ 2.0 * mesh[_x].Nk[8 + i * eq_num] * sigma_x ) ) );

	matrA.coeffRef( 4 + i * eq_num, 5 + i * eq_num ) = 1.0;

	matrA.coeffRef( 5 + i * eq_num, 4 + r * eq_num ) = -B12 / ( B22 * dx * dx );
	matrA.coeffRef( 5 + i * eq_num, 4 + i * eq_num ) = 2.0 * B12 / ( B22 * dx * dx );
	matrA.coeffRef( 5 + i * eq_num, 6 + i * eq_num ) = -12.0 / ( B22 * h * h * h );

	matrA.coeffRef( 6 + i * eq_num, 4 + i * eq_num ) = -( ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 12.0 * beta * B22 * dt * dx * dx ) );
	matrA.coeffRef( 6 + i * eq_num, 4 + r * eq_num ) = ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 24.0 * beta * B22 * dt * dx * dx );
	matrA.coeffRef( 6 + i * eq_num, 5 + i * eq_num ) = -( ( B66 * h * h * h ) / ( 3.0 * dx * dx ) ) - ( h * h * h * rho ) / ( 12.0 * beta * dt * dt ) 
				 - ( h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 24.0 * beta * dt );
	matrA.coeffRef( 6 + i * eq_num, 5 + r * eq_num ) = ( B66 * h * h * h ) / ( 6.0 * dx * dx );
	matrA.coeffRef( 6 + i * eq_num, 6 + i * eq_num ) = ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * beta * B22 * dt );
	matrA.coeffRef( 6 + i * eq_num, 7 + i * eq_num ) = 1.0;
	matrA.coeffRef( 6 + i * eq_num, 8 + i * eq_num ) = ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num]
				+ ( -2.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] ) / ( 2.0 * beta * dt ) ) ) / ( 12.0 * B22 * dx * dx )
				+ ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num] / ( 2.0 * beta * dt ) ) ) / B22;
	matrA.coeffRef( 6 + i * eq_num, 9 + i * eq_num ) = ( B12 * eps_x_0 * h * h * h * ( -2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num] + ( -2.0 * mesh[_x].Nk[4 + i * eq_num]
				+ mesh[_x].Nk[4 + r * eq_num] ) / ( 2.0 * beta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] )
				/ ( 12.0 * B22 * dx * dx ) + ( eps_x_0 * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num] / ( 2.0 * beta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] )
				/ B22 - ( 1.0 / 6.0 ) * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2.0 * beta * dt ) ) * sigma_x;
	
	matrA.coeffRef( 7 + i * eq_num, 1 + i * eq_num ) = ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 4.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + r * eq_num] )
				 * ( 2.0 * beta * newmark_B[5 + r * eq_num] * dt + mesh[_x].Nk[5 + r * eq_num] ) 
				 - 108.0 * beta * By1 * dt * dx * dx * h * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 432.0 * beta * beta * dt * dt * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 1 + r * eq_num ) = -( ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 4.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + r * eq_num] )
				 * ( 2.0 * beta * newmark_B[5 + r * eq_num] * dt + mesh[_x].Nk[5 + r * eq_num] ) ) 
				 / ( 1728.0 * beta * beta * dt * dt * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 4 + i * eq_num ) = -( ( 1.0 / ( 288.0 * beta * beta * B12 * B22 * dt * dt * dx * dx * dx * dx ) ) * ( 24.0
				 * beta * beta * ( B12 * B12 * B12 - B11 * B12 * B22 + 2.0 * B12 * B12 * B66 - 10.0 * B11 * B22 * B66 ) * dt * dt
				 * h * h * h + B12 * B22 * By1 * dx * dx * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] - 2.0 * beta * B12 * B22 * dx * dx * h
				 * ( ( -newmark_B[5 + r * eq_num] ) * By1 * dt * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] + 3.0 * ( 6.0 * dx * dx * ( 8.0 * rho + By1 * By1 * dt * sigma_x )
				 + h * h * ( 8.0 * rho + 4.0 * dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_y + By1 * By1 * dt * sigma_z ) ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 4 + r * eq_num ) = -( ( 1.0 / ( 1152.0 * beta * beta * B12 * B22 * dt * dt * dx * dx * dx * dx ) ) * ( h * h * h * ( 192.0
				 * beta * beta * ( B12 * B12 * B12 - B11 * B12 * B22 + 2.0 * B12 * B12 * B66 + 4.0 * B11 * B22 * B66 ) * dt * dt 
				 - B12 * B22 * By1 * dx * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] + 2.0 * beta * B12 * B22 * dx * dx
				 * ( ( -newmark_B[5 + r * eq_num] ) * By1 * dt * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] + 48.0 * rho
				 + 8.0 * dt * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + i * eq_num] + 2.0 * mesh[_x].Nk[9 + r * eq_num] ) * sigma_y + 6.0 * By1 * By1 * dt * sigma_z ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 4 + rr * eq_num ) = ( ( 3.0 * B12 * B12 * B12 - 3.0 * B11 * B12 * B22 + 6.0 * B12 * B12 * B66 + 2.0 * B11 * B22 * B66 ) * h * h * h)
				 / ( 12.0 * B12 * B22 * dx * dx * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 4 + ( rr + 1 ) * eq_num ) = -( ( ( B12 * B12 - B11 * B22 + 2.0 * B12 * B66 ) * h * h * h ) / ( 12.0 * B22 * dx * dx * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 5 + i * eq_num ) = ( eps_x_0 * ( -6.0 * dx * dx * h + h * h * h ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 12.0 * beta * dt
				 * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 5 + r * eq_num ) = -( ( 1.0 / ( 3456.0 * beta * beta * dt * dt * dx * dx ) ) * ( eps_x_0 * h * h * h * ( mesh[_x].Nk[9 + i * eq_num]
				 * ( 8.0 * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[1 + i * eq_num] - 2.0 * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[1 + r * eq_num]
				 + 8.0 * mesh[_x].Nk[9 + i * eq_num] * ( -4.0 * mesh[_x].Nk[1 + i * eq_num] + mesh[_x].Nk[1 + r * eq_num] ) + 12.0 * By1 * mesh[_x].Nk[4 + i * eq_num]
				 - 3.0 * By1 * mesh[_x].Nk[4 + r * eq_num] )
				 + 2.0 * beta * dt * ( 12.0 * newmark_B[4 + i * eq_num] * By1 * mesh[_x].Nk[9 + i * eq_num] - 3.0 * newmark_B[4 + r * eq_num] * By1
				 * mesh[_x].Nk[9 + i * eq_num] - 32.0 * newmark_B[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] 
				 + 8.0 * newmark_B[1 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + 8.0 * newmark_B[1 + i * eq_num]
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num] - 2.0 * newmark_B[1 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num]
				 + 24.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 24.0 * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num]
				 + 24.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + r * eq_num] ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 6 + i * eq_num ) = ( 2.0 * ( B12 + 2.0 * B66 ) ) / ( B22 * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 6 + r * eq_num ) = -( ( B12 + 2.0 * B66 ) / ( B22 * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 8 + i * eq_num ) = -( ( 1.0 / ( 72.0 * beta * dt * dx * dx ) ) * ( h * ( 36.0 * dx * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[5 + i * eq_num]
				 + eps_x_0 * h * h * ( -6.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + i * eq_num] + ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] )
				 * mesh[_x].Nk[5 + r * eq_num] ) + 2.0 * beta * dt * ( 6.0 * newmark_B[5 + i * eq_num] * eps_x_0 * ( 6.0 * dx * dx - h * h ) * mesh[_x].Nk[9 + i * eq_num]
				 + newmark_B[5 + r * eq_num] * eps_x_0 * h * h * ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] ) + 18.0 * By1 * dx * dx * sigma_x ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 8 + r * eq_num ) = -( ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 * beta * newmark_B[5 + r * eq_num] * dt
				 + mesh[_x].Nk[5 + r * eq_num] ) ) / ( 72.0 * beta * dt * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 9 + i * eq_num ) = ( 1.0 / ( 3456.0 * beta * beta * dt * dt * dx * dx ) ) * ( h * ( eps_x_0 * h * h * ( 2.0 * ( 8.0 * mesh[_x].Nk[9 + i * eq_num] 
				 - mesh[_x].Nk[9 + r * eq_num] ) * ( 4.0 * mesh[_x].Nk[1 + i * eq_num] - mesh[_x].Nk[1 + r * eq_num] ) + 3.0 * By1 * ( -4.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + r * eq_num] ) ) * mesh[_x].Nk[5 + r * eq_num] + 2.0 * beta * dt * ( newmark_B[5 + r * eq_num] * eps_x_0 * h * h
				 * ( 2.0 * ( 8.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + r * eq_num] ) * ( 4.0 * mesh[_x].Nk[1 + i * eq_num] - mesh[_x].Nk[1 + r * eq_num] )
				 + 3.0 * By1 * ( -4.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] ) ) - 12.0 * newmark_B[4 + i * eq_num] * By1 * eps_x_0 * h * h
				 * mesh[_x].Nk[5 + r * eq_num]
				 + 3.0 * newmark_B[4 + r * eq_num] * By1 * eps_x_0 * h * h * mesh[_x].Nk[5 + r * eq_num] + 64.0 * newmark_B[1 + i * eq_num] * eps_x_0 * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] - 16.0 * newmark_B[1 + r * eq_num] * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[5 + r * eq_num]
				 - 8.0 * newmark_B[1 + i * eq_num] * eps_x_0 * h * h * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[5 + r * eq_num] + 2.0 * newmark_B[1 + r * eq_num]
				 * eps_x_0 * h * h * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[5 + r * eq_num] - 864.0 * dx * dx * eps_x_0 * mesh[_x].Nk[5 + i * eq_num]
				 * mesh[_x].Nk[8 + i * eq_num]
				 + 144.0 * eps_x_0 * h * h * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] - 24.0 * eps_x_0 * h * h * mesh[_x].Nk[5 + r * eq_num]
				 * mesh[_x].Nk[8 + i * eq_num] - 24.0 * eps_x_0 * h * h * mesh[_x].Nk[5 + r * eq_num] * mesh[_x].Nk[8 + r * eq_num]
				 - 432.0 * By1 * dx * dx * mesh[_x].Nk[1 + i * eq_num] * sigma_x + 48.0 * h * h * ( 6.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num]
				 - ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] ) * mesh[_x].Nk[4 + r * eq_num] ) * sigma_y )
				 - 4.0 * beta * beta * dt * dt * ( 2.0 * eps_x_0 * ( 432.0 * newmark_B[5 + i * eq_num] * dx * dx * mesh[_x].Nk[8 + i * eq_num] + h * h
				 * ( -72.0 * newmark_B[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
				 + newmark_B[5 + r * eq_num] * ( ( -( 4.0 * newmark_B[1 + i * eq_num] - newmark_B[1 + r * eq_num] ) ) * ( 8.0 * mesh[_x].Nk[9 + i * eq_num]
				 - mesh[_x].Nk[9 + r * eq_num] ) + 12.0 * ( mesh[_x].Nk[8 + i * eq_num] + mesh[_x].Nk[8 + r * eq_num] ) ) ) ) + 432.0 * newmark_B[1 + i * eq_num]
				 * By1 * dx * dx * sigma_x
				 + 12.0 * newmark_B[4 + i * eq_num] * h * h * ( newmark_B[5 + r * eq_num] * By1 * eps_x_0 - 24.0 * mesh[_x].Nk[9 + i * eq_num] * sigma_y )
				 + 3.0 * newmark_B[4 + r * eq_num] * h * h
				 * ( ( -newmark_B[5 + r * eq_num] ) * By1 * eps_x_0 + 16.0 * ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] ) * sigma_y ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 9 + r * eq_num ) = -( ( h * h * h * ( eps_x_0 * ( 2.0 * beta * newmark_B[5 + r * eq_num] * dt + mesh[_x].Nk[5 + r * eq_num] )
				* ( mesh[_x].Nk[9 + i * eq_num] * ( 8.0 * beta * newmark_B[1 + i * eq_num] * dt 
				- 2.0 * beta * newmark_B[1 + r * eq_num] * dt + 4.0 * mesh[_x].Nk[1 + i * eq_num] - mesh[_x].Nk[1 + r * eq_num] ) + 24.0 * beta * dt
				* mesh[_x].Nk[8 + i * eq_num] ) + 48.0 * beta * dt * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 * beta * newmark_B[4 + r * eq_num] * dt
				+ mesh[_x].Nk[4 + r * eq_num] ) * sigma_y ) ) / ( 1728.0 * beta * beta * dt * dt * dx * dx ) );

	matrA.coeffRef( 8 + i * eq_num, 0 + i * eq_num ) = ( -( ( 4.0 * mesh[_x].Nk[9 + i * eq_num] ) / 3.0 ) + ( 4.0 * mesh[_x].Nk[9 + r * eq_num] ) / 3.0 ) / ( 4.0 * beta * dt * dx );
	matrA.coeffRef( 8 + i * eq_num, 0 + r * eq_num ) = mesh[_x].Nk[9 + i * eq_num] / ( 4.0 * beta * dt * dx );
	matrA.coeffRef( 8 + i * eq_num, 9 + i * eq_num ) = 1.0 / ( 2.0 * beta * dt ) - ( 2.0 * ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2.0 * beta * dt ) ) )
				/ ( 3.0 * dx )
				+ ( newmark_B[0 + r * eq_num] + mesh[_x].Nk[0 + r * eq_num] / ( 2.0 * beta * dt ) ) / ( 2.0 * dx ) + 2.0 / ( 3.0 * dx * dx * mu * sigma_y );
	matrA.coeffRef( 8 + i * eq_num, 9 + r * eq_num ) = ( 2.0 * ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2.0 * beta * dt ) ) ) / ( 3.0 * dx ) - 2.0
				/ ( 3.0 * dx * dx * mu * sigma_y );

	matrA.coeffRef( 9 + i * eq_num, 1 + i * eq_num ) = ( mu * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 2.0 * beta * dt );
	matrA.coeffRef( 9 + i * eq_num, 4 + i * eq_num ) = -( ( By1 * mu * sigma_x ) / ( 4.0 * beta * dt ) );
	matrA.coeffRef( 9 + i * eq_num, 8 + i * eq_num ) = sigma_x_mu;
	matrA.coeffRef( 9 + i * eq_num, 9 + i * eq_num ) = mu * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2.0 * beta * dt ) ) * sigma_x;

	vectF( 2 + i * eq_num ) = ( h * ( 3.0 * eps_x_0 * mu * ( -4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + r * eq_num] + By1 * mesh[_x].Nk[4 + r * eq_num] )
			* mesh[_x].Nk[8 + i * eq_num] + 2.0 * beta * dt * ( 8.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			- 2.0 * mesh[_x].Nk[9 + i * eq_num] * ( 4.0 * mesh[_x].Nk[9 + r * eq_num] + 3.0 * newmark_B[1 + r * eq_num] * eps_x_0 * mu * mesh[_x].Nk[8 + i * eq_num] )
			+ 3.0 * dx * mu * ( 4.0 * newmark_A[0 + i * eq_num] * rho + newmark_B[0 + i * eq_num] * By1 * By1 * sigma_z ) ) ) )
			/ ( 24.0 * beta * dt * dx * mu );

	vectF( 3 + i * eq_num ) = ( 1.0 / ( 4.0 * beta * B22 * dt * dx ) ) * ( eps_x_0 * ( 2.0 * B12 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[0 + r * eq_num]
			- 4.0 * dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[3 + i * eq_num]
			+ B22 * By1 * dx * h * mesh[_x].Nk[5 + i * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] + B22 * dx * h * mesh[_x].Nk[9 + i * eq_num]
			* ( -4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] + By1 * mesh[_x].Nk[4 + i * eq_num] ) * sigma_x
			+ 2.0 * beta * dt * ( -2.0 * newmark_B[3 + i * eq_num] * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + h
			* ( B12 * newmark_B[0 + r * eq_num] * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 2.0 * newmark_A[1 + i * eq_num] * B22 * dx * rho
			- 2.0 * B22 * dx * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[8 + i * eq_num] ) * sigma_x ) ) );

	vectF( 6 + i * eq_num ) = ( 1.0 / ( 12.0 * beta * B22 * dt * dx * dx ) ) * ( mesh[_x].Nk[9 + i * eq_num] * ( B12 * eps_x_0 * h * h * h * ( 2.0
			* mesh[_x].Nk[4 + i * eq_num] - mesh[_x].Nk[4 + r * eq_num] ) * mesh[_x].Nk[8 + i * eq_num]
			+ dx * dx * ( -12.0 * eps_x_0 * mesh[_x].Nk[6 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + B22 * h * h * h * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[5 + i * eq_num] * sigma_x ) ) - beta * dt * ( 12.0 * newmark_B[6 + i * eq_num] * dx * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[8 + i * eq_num]
			+ h * h * h * ( B12 * ( -2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + r * eq_num] ) * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			+ B22 * dx * dx * ( newmark_A[5 + i * eq_num] * rho - newmark_B[5 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) ) ) );

	vectF( 7 + i * eq_num ) = ( 1.0 / ( 1728.0 * beta * beta * dt * dt * dx * dx ) ) * ( -3.0 * eps_x_0 * h * h * h
			* mesh[_x].Nk[9 + i * eq_num] * ( ( 4.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + r * eq_num] ) * ( 4.0 * mesh[_x].Nk[1 + i * eq_num]
			- mesh[_x].Nk[1 + r * eq_num]) + By1 * ( -4.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + r * eq_num] ) ) * mesh[_x].Nk[5 + r * eq_num]
			+ beta * dt * h * ( newmark_B[5 + r * eq_num] * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 16.0 * mesh[_x].Nk[9 + r * eq_num]
			* mesh[_x].Nk[1 + i * eq_num] - 4.0 * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[1 + r * eq_num]
			+ 16.0 * mesh[_x].Nk[9 + i * eq_num] * ( -4.0 * mesh[_x].Nk[1 + i * eq_num] + mesh[_x].Nk[1 + r * eq_num] ) + 12.0 * By1 * mesh[_x].Nk[4 + i * eq_num]
			- 3.0 * By1 * mesh[_x].Nk[4 + r * eq_num] ) + 12.0 * newmark_B[4 + i * eq_num] * By1 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num]
			- 3.0 * newmark_B[4 + r * eq_num] * By1 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] - 64.0 * newmark_B[1 + i * eq_num]
			* eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] + 16.0 * newmark_B[1 + r * eq_num] * eps_x_0 * h * h
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] + 16.0 * newmark_B[1 + i * eq_num] * eps_x_0 * h * h
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[5 + r * eq_num] - 4.0 * newmark_B[1 + r * eq_num] * eps_x_0 * h * h
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[5 + r * eq_num]
			+ 1728.0 * dx * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] - 288.0 * eps_x_0 * h * h
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 48.0 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[5 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			+ 48.0 * eps_x_0 * h * h * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[5 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 48.0 * eps_x_0 * h * h
			* mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + r * eq_num] * mesh[_x].Nk[8 + r * eq_num] + 432.0 * By1 * dx * dx * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[1 + i * eq_num] * sigma_x
			+ 48.0 * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -6.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[4 + r * eq_num] + 2.0 * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[4 + r * eq_num] ) * sigma_y ) + 4.0 * beta * beta * dt * dt
			* ( -216.0 * By1 * dx * dx * h * Jx + 4.0 * newmark_B[1 + r * eq_num] * newmark_B[5 + r * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[9 + i * eq_num]
			- newmark_B[1 + r * eq_num] * newmark_B[5 + r * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num] + 4.0
			* newmark_B[1 + i * eq_num] * newmark_B[5 + r * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -4.0 * mesh[_x].Nk[9 + i * eq_num]
			+ mesh[_x].Nk[9 + r * eq_num] )
			+ 432.0 * newmark_B[5 + i * eq_num] * dx * dx * eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] - 72.0 * newmark_B[5 + i * eq_num]
			* eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 12.0 * newmark_B[5 + r * eq_num] * eps_x_0
			* h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 12.0 * newmark_B[5 + r * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + r * eq_num]
			* mesh[_x].Nk[8 + i * eq_num] + 12.0 * newmark_B[5 + r * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + r * eq_num]
			+ 432.0 * newmark_A[4 + i * eq_num] * dx * dx * h * rho + 72.0 * newmark_A[4 + i * eq_num] * h * h * h * rho - 36.0 * newmark_A[4 + r * eq_num] * h * h * h * rho
			+ 432.0 * dx * dx * Pimp
			+ 12.0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -6.0 * newmark_B[4 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + newmark_B[4 + r * eq_num]
			* mesh[_x].Nk[9 + i * eq_num] + 2.0 * newmark_B[4 + r * eq_num] * mesh[_x].Nk[9 + r * eq_num] ) * sigma_y
			+ 9.0 * By1 * By1 * ( 12.0 * newmark_B[4 + i * eq_num] * dx * dx * h * sigma_x + ( 2.0 * newmark_B[4 + i * eq_num] - newmark_B[4 + r * eq_num] )
			* h * h * h * sigma_z ) ) );

	vectF( 8 + i * eq_num ) = ( 12.0 * beta * newmark_B[9 + i * eq_num] * dt * dx + 4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[0 + i * eq_num]
			- 4.0 * mesh[_x].Nk[9 + r * eq_num] * mesh[_x].Nk[0 + i * eq_num] - 3.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[0 + r * eq_num] )
			/ ( 12.0 * beta * dt * dx );

	vectF( 9 + i * eq_num ) = -( ( mu * ( beta * newmark_B[4 + i * eq_num] * By1 * dt + mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] ) * sigma_x )
			/ ( 2.0 * beta * dt ) );

	i = nx - 1;
	r = i + 1;
	rr = i + 2;
	j = i - 1;
	jj = i - 2;
	
	//for the right line:
	matrA.coeffRef( 0 + i * eq_num, 1 + j * eq_num ) = 1.0 / ( 2.0 * dx );
	matrA.coeffRef( 0 + i * eq_num, 2 + i * eq_num ) = 1.0 / ( h * B66 );

	matrA.coeffRef( 1 + i * eq_num, 0 + j * eq_num ) = B12 / ( B22 * 2.0 * dx );
	matrA.coeffRef( 1 + i * eq_num, 3 + i * eq_num ) = 1.0 / ( h * B22 );

	matrA.coeffRef( 2 + i * eq_num, 0 + j * eq_num ) = ( ( 3.0 * B12 * B12 - 4.0 * B11 * B22 ) * h ) / ( 4.0 * B22 * dx * dx );
	matrA.coeffRef( 2 + i * eq_num, 0 + i * eq_num ) = ( 2.0 * B11 * h ) / dx / dx - ( B12 * B12 * h ) / ( B22 * dx * dx ) + ( h * ( 8.0 * rho + By1 * By1 * dt * sigma_z ) )
				/ ( 8.0 * beta * dt * dt );
	matrA.coeffRef( 2 + i * eq_num, 1 + j * eq_num ) = -( ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt * dx ) );
	matrA.coeffRef( 2 + i * eq_num, 3 + j * eq_num ) = B12 / ( 2.0 * B22 * dx );
	matrA.coeffRef( 2 + i * eq_num, 4 + j * eq_num ) = ( By1 * eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] ) / ( 8.0 * beta * dt * dx );
	matrA.coeffRef( 2 + i * eq_num, 8 + i * eq_num ) = ( eps_x_0 * h * ( 2.0 * newmark_B[4 + j * eq_num] * beta * By1 * dt - 2.0 * mesh[_x].Nk[9 + i * eq_num]
				* ( 2.0 * newmark_B[1 + j * eq_num] * beta * dt + mesh[_x].Nk[1 + j * eq_num] )
				+ By1 * mesh[_x].Nk[4 + j * eq_num] ) ) / ( 8.0 * beta * dt * dx);
	matrA.coeffRef( 2 + i * eq_num, 9 + i * eq_num ) = ( h * ( 8.0 * beta * dt * ( 2.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) - 3.0 * eps_x_0 * mu
				* ( 2.0 * newmark_B[1 + j * eq_num] * beta * dt + mesh[_x].Nk[1 + j * eq_num] )
				* mesh[_x].Nk[8 + i * eq_num] ) ) / ( 12.0 * beta * dt * dx * mu );
	matrA.coeffRef( 2 + i * eq_num, 9 + j * eq_num ) = -( ( 2.0 * h * mesh[_x].Nk[9 + i * eq_num] ) / ( 3.0 * dx * mu ) );

	matrA.coeffRef( 3 + i * eq_num, 0 + j * eq_num ) = ( B12 * eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * B22 * beta * dt * dx );
	matrA.coeffRef( 3 + i * eq_num, 1 + j * eq_num ) = -( ( B66 * h ) / ( 4.0 * dx * dx ) );
	matrA.coeffRef( 3 + i * eq_num, 1 + i * eq_num ) = ( B66 * h ) / dx / dx + ( h * ( 2.0 * rho + dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) )
				 / ( 2.0 * beta * dt * dt );
	matrA.coeffRef( 3 + i * eq_num, 2 + j * eq_num ) = 1.0 / ( 2.0 * dx );
	matrA.coeffRef( 3 + i * eq_num, 3 + i * eq_num ) = ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * B22 * beta * dt );
	matrA.coeffRef( 3 + i * eq_num, 4 + i * eq_num ) = -( ( By1 * h * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 4.0 * beta * dt ) );
	matrA.coeffRef( 3 + i * eq_num, 5 + i * eq_num ) = -( ( By1 * eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt ) );
	matrA.coeffRef( 3 + i * eq_num, 8 + i * eq_num ) = ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( 4.0 * newmark_B[3 + i * eq_num] * beta * dt * dx + B12 * h
				* ( 2.0 * newmark_B[0 + j * eq_num] * beta * dt + mesh[_x].Nk[0 + j * eq_num] )
				+ 2.0 * dx * mesh[_x].Nk[3 + i * eq_num] ) - B22 * dx * h * ( By1 * eps_x_0 * ( 2.0 * newmark_B[5 + i * eq_num] * beta * dt + mesh[_x].Nk[5 + i * eq_num] )
				- 4.0 * beta * dt * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) )
				/ ( 4.0 * B22 * beta * dt * dx );
	matrA.coeffRef( 3 + i * eq_num, 9 + i * eq_num ) = ( 1.0 / ( 4.0 * B22 * beta * dt * dx ) ) * ( eps_x_0 * ( 4.0 * newmark_B[3 + i * eq_num] * beta * dt * dx
				+ B12 * h * ( 2.0 * newmark_B[0 + j * eq_num] * beta * dt + mesh[_x].Nk[0 + j * eq_num] ) + 2.0 * dx * mesh[_x].Nk[3 + i * eq_num] )
				* mesh[_x].Nk[8 + i * eq_num] + B22 * dx * h * ( ( 4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] - By1 * mesh[_x].Nk[4 + i * eq_num] )
				* sigma_x + 2.0 * beta * dt * ( 2.0 * Jx - newmark_B[4 + i * eq_num] * By1 * sigma_x + 4.0 * newmark_B[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x
				+ 2.0 * mesh[_x].Nk[8 + i * eq_num] * sigma_x ) ) );

	matrA.coeffRef( 4 + i * eq_num, 5 + i * eq_num ) = 1.0;

	matrA.coeffRef( 5 + i * eq_num, 4 + j * eq_num ) = -B12 / ( B22 * dx * dx );
	matrA.coeffRef( 5 + i * eq_num, 4 + i * eq_num ) = 2.0 * B12 / ( B22 * dx * dx );
	matrA.coeffRef( 5 + i * eq_num, 6 + i * eq_num ) = -12.0 / ( B22 * h * h * h );

	matrA.coeffRef( 6 + i * eq_num, 4 + j * eq_num ) = ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 24.0 * B22 * beta * dt * dx * dx );
	matrA.coeffRef( 6 + i * eq_num, 4 + i * eq_num ) =  -( ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 12.0 * B22 * beta * dt * dx * dx ) );
	matrA.coeffRef( 6 + i * eq_num, 5 + j * eq_num ) = h * h * h * B66 / ( 6.0 * dx * dx );
	matrA.coeffRef( 6 + i * eq_num, 5 + i * eq_num ) = -( ( B66 * h * h * h ) / ( 3.0 * dx * dx ) ) - ( h * h * h * rho ) / ( 12.0 * beta * dt * dt )
				- ( h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 24.0 * beta * dt );
	matrA.coeffRef( 6 + i * eq_num, 6 + i * eq_num ) = ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * B22 * beta * dt );
	matrA.coeffRef( 6 + i * eq_num, 7 + i * eq_num ) = 1.0;
	matrA.coeffRef( 6 + i * eq_num, 8 + i * eq_num ) = ( B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num]
				+ ( -2.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + j * eq_num] )
				/ ( 2.0 * beta * dt ) ) ) / ( 12.0 * B22 * dx * dx ) + (eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num]
				/ ( 2.0 * beta * dt ) ) ) / B22;
	matrA.coeffRef( 6 + i * eq_num, 9 + i * eq_num ) = ( B12 * eps_x_0 * h * h * h * ( -2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] + ( -2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) / ( 2.0 * beta * dt ) )
				 * mesh[_x].Nk[8 + i * eq_num] ) / ( 12.0 * B22 * dx * dx ) + ( eps_x_0 * ( newmark_B[6 + i * eq_num] + mesh[_x].Nk[6 + i * eq_num] / ( 2.0 * beta * dt ) )
				 * mesh[_x].Nk[8 + i * eq_num] ) / B22 - ( 1.0 / 6.0 ) * h * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2.0 * beta * dt ) ) * sigma_x;
	

	matrA.coeffRef( 7 + i * eq_num, 1 + j * eq_num ) = -( ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 4.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 * ( 2.0 * newmark_B[5 + j * eq_num] * beta * dt + mesh[_x].Nk[5 + j * eq_num] ) )
				 / ( 1728.0 * beta * beta * dt * dt * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 1 + i * eq_num ) = ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 4.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 * ( 2.0 * newmark_B[5 + j * eq_num] * beta * dt + mesh[_x].Nk[5 + j * eq_num])
				 - 108.0 * beta * By1 * dt * dx * dx * h * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 432.0 * beta * beta * dt * dt * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 4 + ( jj - 1 ) * eq_num ) = -( ( ( B12 * B12 - B11 * B22 + 2.0 * B12 * B66 ) * h * h * h ) / ( 12.0 * B22 * dx * dx * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 4 + jj * eq_num ) = ( ( 3.0 * B12 * B12 * B12 - 3.0 * B11 * B12 * B22 + 6.0 * B12 * B12 * B66 + 2.0 * B11 * B22 * B66 ) * h * h * h )
				 / ( 12.0 * B12 * B22 * dx * dx * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 4 + j * eq_num ) = -( ( 1.0 / ( 1152.0 * B12 * B22 * beta * beta * dt * dt * dx * dx * dx * dx ) )
				 * ( h * h * h * ( 192.0 * B12 * B12 * B12 * beta * beta * dt * dt + 384.0 * B12 * B12 * B66 * beta * beta * dt * dt
				 + 768.0 * B11 * B22 * B66 * beta * beta * dt * dt + B12 * B22 * ( -192.0 * B11 * beta * beta * dt * dt
				 + dx * dx * ( -2.0 * newmark_B[5 + j * eq_num] * beta * By1 * dt * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] - By1 * eps_x_0 * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[5 + j * eq_num] + 16.0 * beta * ( 6.0 * rho + dt
				 * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + i * eq_num] + 2.0 * mesh[_x].Nk[9 + j * eq_num] ) * sigma_y ) + 12.0 * beta * By1 * By1 * dt * sigma_z ) ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 4 + i * eq_num ) = ( 1.0 / ( 288.0 * B12 * B22 * beta * beta * dt * dt * dx * dx * dx * dx ) ) * ( -24.0 * B12 * B12 * B12
				 * beta* beta * dt * dt * h * h * h - 48.0 * B12 * B12 * B66 * beta * beta * dt * dt * h * h * h + 240.0 * B11 * B22 * B66 * beta * beta
				 * dt* dt * h * h * h + B12 * B22 * h * ( 24.0 * B11 * beta * beta * dt * dt * h * h + dx * dx * ( -2.0 * newmark_B[5 + j * eq_num] * beta * By1 * dt
				 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] - By1 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + j * eq_num] + 24.0 * beta
				 * ( 2.0 * ( 6.0 * dx * dx + h * h ) * rho
				 + dt * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * sigma_y ) + 6.0 * beta * By1 * By1 * dt * ( 6.0 * dx * dx * sigma_x + h * h * sigma_z ) ) ) );
	matrA.coeffRef( 7 + i * eq_num, 5 + j * eq_num ) = -( ( 1.0 / ( 3456.0 * beta * beta * dt * dt * dx * dx ) ) * ( eps_x_0 * h * h * h * ( 24.0
				 * newmark_B[4 + i * eq_num] * beta * By1 * dt * mesh[_x].Nk[9 + i * eq_num] - 6.0 * newmark_B[4 + j * eq_num] * beta * By1 * dt * mesh[_x].Nk[9 + i * eq_num]
				 - 64.0 * newmark_B[1 + i * eq_num] * beta * dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + 16.0 * newmark_B[1 + j * eq_num] * beta
				 * dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + 16.0 * newmark_B[1 + i * eq_num] * beta * dt * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[9 + j * eq_num] - 4.0 * newmark_B[1 + j * eq_num] * beta * dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num]
				 - 32.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num]
				 + 8.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[1 + i * eq_num] + 8.0 * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + j * eq_num] - 2.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[1 + j * eq_num]
				 + 12.0 * By1 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num]
				 - 3.0 * By1 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + j * eq_num] + 48.0 * beta * dt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
				 + 48.0 * beta * dt * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 48.0 * beta * dt * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[8 + j * eq_num] ) ) );
	matrA.coeffRef( 7 + i * eq_num, 5 + i * eq_num ) = ( eps_x_0 * ( -6.0 * dx * dx * h + h * h * h ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] )
				 / ( 12.0 * beta * dt * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 6 + j * eq_num ) = -( ( B12 + 2.0 * B66 ) / ( B22 * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 6 + i * eq_num ) = ( 2.0 * ( B12 + 2.0 * B66 ) ) / ( B22 * dx * dx );
	matrA.coeffRef( 7 + i * eq_num, 8 + j * eq_num ) =  -( ( eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 * newmark_B[5 + j * eq_num] * beta * dt
				 + mesh[_x].Nk[5 + j * eq_num] ) ) / ( 72.0 * beta * dt * dx * dx ) );
	matrA.coeffRef( 7 + i * eq_num, 8 + i * eq_num ) =  -( ( 1.0 / ( 72.0 * beta * dt * dx * dx ) ) * ( h * ( eps_x_0 * ( 12.0 * newmark_B[5 + i * eq_num] * beta * dt * ( 6.0
				 * dx * dx - h * h )
				 * mesh[_x].Nk[9 + i * eq_num] + 2.0 * newmark_B[5 + j * eq_num] * beta * dt * h * h * ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + j * eq_num] )
				 - 6.0 * ( -6.0 * dx * dx + h * h ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + i * eq_num] + h * h * ( mesh[_x].Nk[9 + i * eq_num]
				 + mesh[_x].Nk[9 + j * eq_num] )
				 * mesh[_x].Nk[5 + j * eq_num] ) + 36.0 * beta * By1 * dt * dx * dx * sigma_x ) ) );
	matrA.coeffRef( 7 + i * eq_num, 9 + j * eq_num ) = -( ( 1.0 / ( 1728.0 * beta * beta * dt * dt * dx * dx ) ) * ( h * h * h * ( eps_x_0 * ( 2.0 * newmark_B[5 + j * eq_num] * beta * dt
				 + mesh[_x].Nk[5 + j * eq_num] ) * ( 8.0 * newmark_B[1 + i * eq_num] * beta * dt * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + i * eq_num]
				 * ( 2.0 * newmark_B[1 + j * eq_num] * beta * dt - 4.0 * mesh[_x].Nk[1 + i * eq_num] + mesh[_x].Nk[1 + j * eq_num] ) + 24.0 * beta * dt
				 * mesh[_x].Nk[8 + i * eq_num] )
				 + 48.0 * beta * dt * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 * newmark_B[4 + j * eq_num] * beta * dt + mesh[_x].Nk[4 + j * eq_num] ) * sigma_y ) ) );
	matrA.coeffRef( 7 + i * eq_num, 9 + i * eq_num ) = ( 1.0 / ( 3456.0 * beta * beta * dt * dt ) ) * ( h * ( -1728.0 * beta * dt * eps_x_0
				 * ( 2.0 * newmark_B[5 + i * eq_num] * beta * dt + mesh[_x].Nk[5 + i * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] + ( 144.0 * beta * dt * eps_x_0 * h
				 * h * ( 4.0 * newmark_B[5 + i * eq_num] * beta * dt - 2.0 * newmark_B[5 + j * eq_num] * beta * dt + 2.0 * mesh[_x].Nk[5 + i * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / dx / dx + ( 1.0 / dx / dx) * ( eps_x_0 * h * h * ( 2.0 * newmark_B[5 + j * eq_num]
				 * beta * dt + mesh[_x].Nk[5 + j * eq_num] ) * ( -24.0 * newmark_B[4 + i * eq_num] * beta * By1 * dt
				 + 6.0 * newmark_B[4 + j * eq_num] * beta * By1 * dt + 3.0 * By1 * ( -4.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + j * eq_num] )
				 + 2.0 * ( 8.0 * newmark_B[1 + i * eq_num] * beta * dt * ( 8.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 + 2.0 * newmark_B[1 + j * eq_num] * beta * dt * ( -8.0 * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + j * eq_num] ) + ( 8.0 * mesh[_x].Nk[9 + i * eq_num]
				 - mesh[_x].Nk[9 + j * eq_num] ) * ( 4.0 * mesh[_x].Nk[1 + i * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) + 24.0 * beta * dt
				 * ( 2.0 * mesh[_x].Nk[8 + i * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) ) ) ) - 864.0 * beta * By1 * dt * ( 2.0 * newmark_B[1 + i * eq_num] * beta * dt
				 + mesh[_x].Nk[1 + i * eq_num] ) * sigma_x + ( 96.0 * beta * dt
				 * h * h * ( 12.0 * newmark_B[4 + i * eq_num] * beta * dt * mesh[_x].Nk[9 + i * eq_num] - 2.0 * newmark_B[4 + j * eq_num] * beta * dt
				 * ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + j * eq_num] ) + 6.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num]
				 - ( mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[4 + j * eq_num] ) * sigma_y ) / dx / dx ) );
	
	matrA.coeffRef( 8 + i * eq_num, 0 + j * eq_num ) = -( mesh[_x].Nk[9 + i * eq_num] / ( 4.0 * beta * dt * dx ) );
	matrA.coeffRef( 8 + i * eq_num, 0 + i * eq_num ) = ( ( 4.0 * mesh[_x].Nk[9 + i * eq_num] ) / 3.0 - ( 4.0 * mesh[_x].Nk[9 + j * eq_num] ) / 3.0 ) / ( 4.0 * beta * dt * dx );
	matrA.coeffRef( 8 + i * eq_num, 9 + j * eq_num ) = -( ( 2.0 * ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2.0 * beta * dt ) ) ) / ( 3.0 * dx ) ) - 2.0
				 / ( 3.0 * dx * dx * mu * sigma_y );
	matrA.coeffRef( 8 + i * eq_num, 9 + i * eq_num ) = 1.0 / ( 2.0 * beta * dt ) + ( 2.0 * ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2.0 * beta * dt ) ) )
				 / ( 3.0 * dx )
				 + ( -newmark_B[0 + j * eq_num] - mesh[_x].Nk[0 + j * eq_num] / ( 2.0 * beta * dt ) ) / ( 2.0 * dx ) + 2.0 / ( 3.0 * dx * dx * mu * sigma_y );

	matrA.coeffRef( 9 + i * eq_num, 1 + i * eq_num ) = ( mu * mesh[_x].Nk[9 + i * eq_num] * sigma_x ) / ( 2.0 * beta * dt );
	matrA.coeffRef( 9 + i * eq_num, 4 + i * eq_num ) = -( ( By1 * mu * sigma_x ) / ( 4.0 * beta * dt ) );
	matrA.coeffRef( 9 + i * eq_num, 8 + i * eq_num ) = sigma_x_mu;
	matrA.coeffRef( 9 + i * eq_num, 9 + i * eq_num ) = mu * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2.0 * beta * dt ) ) * sigma_x;

	vectF( 2 + i * eq_num ) = ( h * ( 3.0 * eps_x_0 * mu * ( 4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + j * eq_num] - By1 * mesh[_x].Nk[4 + j * eq_num] )
				 * mesh[_x].Nk[8 + i * eq_num] + 2.0 * beta * dt * ( -8.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
				 + mesh[_x].Nk[9 + i * eq_num] * ( 8.0 * mesh[_x].Nk[9 + j * eq_num] + 6.0 * newmark_B[1 + j * eq_num] * eps_x_0 * mu * mesh[_x].Nk[8 + i * eq_num] )
				 + 3.0 * dx * mu * ( 4.0 * newmark_A[0 + i * eq_num] * rho + newmark_B[0 + i * eq_num] * By1 * By1 * sigma_z ) ) ) )
				 / ( 24.0 * beta * dt * dx * mu );

	vectF( 3 + i * eq_num ) = ( 1.0 / ( 4.0 * B22 * beta * dt * dx ) ) * ( -4.0 * newmark_B[3 + i * eq_num] * beta * dt * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[8 + i * eq_num]
				 - 2.0 * B12 * eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[0 + j * eq_num] * beta * dt + mesh[_x].Nk[0 + j * eq_num] )
				 * mesh[_x].Nk[8 + i * eq_num] + dx * ( -4.0 * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[3 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
				 + B22 * By1 * eps_x_0 * h * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 4.0 * newmark_A[1 + i * eq_num] * B22 * beta * dt * h * rho
				 + B22 * h * mesh[_x].Nk[9 + i * eq_num] * ( By1 * mesh[_x].Nk[4 + i * eq_num]
				 - 4.0 * ( newmark_B[1 + i * eq_num] * beta * dt * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] + beta * dt
				 * mesh[_x].Nk[8 + i * eq_num] ) ) * sigma_x ) );

	vectF( 6 + i * eq_num ) = ( 1.0 / ( 12.0 * B22 * beta * dt * dx * dx ) ) * ( -12.0 * newmark_B[6 + i * eq_num] * beta * dt * dx * dx * eps_x_0
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
				 + B12 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 * newmark_B[4 + i * eq_num] * beta * dt - newmark_B[4 + j * eq_num] * beta * dt
				 + 2.0 * mesh[_x].Nk[4 + i * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num]
				 + dx * dx * ( -12.0 * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[6 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + B22 * h * h * h
				 * ( ( -newmark_A[5 + i * eq_num] ) * beta * dt * rho + mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
				 * ( newmark_B[5 + i * eq_num] * beta * dt + mesh[_x].Nk[5 + i * eq_num] ) * sigma_x ) ) );

	vectF( 7 + i * eq_num ) = ( 1.0 / ( 1728.0 * beta * beta * dt * dt * dx * dx ) ) * ( -3.0 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num]
				 * ( ( 4.0 * mesh[_x].Nk[9 + i * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 * ( 4.0 * mesh[_x].Nk[1 + i * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) + By1 * ( -4.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + j * eq_num] ) )
				 * mesh[_x].Nk[5 + j * eq_num] + beta * dt * h * ( newmark_B[5 + j * eq_num] * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * ( 16.0
				 * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[1 + i * eq_num] - 4.0 * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[1 + j * eq_num] + 16.0 * mesh[_x].Nk[9 + i * eq_num]
				 * ( -4.0 * mesh[_x].Nk[1 + i * eq_num] + mesh[_x].Nk[1 + j * eq_num] ) + 12.0 * By1 * mesh[_x].Nk[4 + i * eq_num] - 3.0 * By1 * mesh[_x].Nk[4 + j * eq_num] )
				 + 12.0 * newmark_B[4 + i * eq_num] * By1 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + j * eq_num] - 3.0 * newmark_B[4 + j * eq_num] * By1
				 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + j * eq_num] - 64.0 * newmark_B[1 + i * eq_num] * eps_x_0 * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + j * eq_num] + 16.0 * newmark_B[1 + j * eq_num] * eps_x_0 * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + j * eq_num] + 16.0 * newmark_B[1 + i * eq_num] * eps_x_0 * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[5 + j * eq_num]
				 - 4.0 * newmark_B[1 + j * eq_num] * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[5 + j * eq_num]
				 + 1728.0 * dx * dx * eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] - 288.0 * eps_x_0 * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 48.0 * eps_x_0 * h * h * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[5 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 48.0 * eps_x_0 * h * h * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[5 + j * eq_num]
				 * mesh[_x].Nk[8 + i * eq_num] + 48.0 * eps_x_0
				 * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[5 + j * eq_num] * mesh[_x].Nk[8 + j * eq_num] + 432.0 * By1 * dx * dx * mesh[_x].Nk[9 + i * eq_num]
				 * mesh[_x].Nk[1 + i * eq_num] * sigma_x + 48.0 * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -6.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + j * eq_num] + 2.0 * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[4 + j * eq_num] ) * sigma_y ) + 4.0 * beta
				 * beta * dt * dt * ( -216.0 * By1 * dx * dx * h * Jx + 4.0 * newmark_B[1 + j * eq_num]
				 * newmark_B[5 + j * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - newmark_B[1 + j * eq_num]
				 * newmark_B[5 + j * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] + 4.0 * newmark_B[1 + i * eq_num]
				 * newmark_B[5 + j * eq_num] * eps_x_0 * h * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * ( -4.0 * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + j * eq_num] ) + 432.0 * newmark_B[5 + i * eq_num] * dx * dx
				 * eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] - 72.0 * newmark_B[5 + i * eq_num] * eps_x_0 * h * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
				 + 12.0 * newmark_B[5 + j * eq_num] * eps_x_0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 12.0 * newmark_B[5 + j * eq_num]
				 * eps_x_0 * h * h * h * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 12.0 * newmark_B[5 + j * eq_num] * eps_x_0 * h * h * h
				 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + j * eq_num] + 432.0 * newmark_A[4 + i * eq_num] * dx * dx * h * rho + 72.0 * newmark_A[4 + i * eq_num]
				 * h * h * h * rho - 36.0 * newmark_A[4 + j * eq_num] * h * h * h * rho + 432.0 * dx * dx * Pimp
				 + 12.0 * h * h * h * mesh[_x].Nk[9 + i * eq_num] * ( -6.0 * newmark_B[4 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + newmark_B[4 + j * eq_num]
				 * mesh[_x].Nk[9 + i * eq_num] + 2.0 * newmark_B[4 + j * eq_num] * mesh[_x].Nk[9 + j * eq_num] ) * sigma_y + 9.0 * By1 * By1
				 * ( 12.0 * newmark_B[4 + i * eq_num]
				 * dx * dx * h * sigma_x + ( 2.0 * newmark_B[4 + i * eq_num] - newmark_B[4 + j * eq_num] ) * h * h * h * sigma_z ) ) );

	vectF( 8 + i * eq_num ) = ( 12.0 * newmark_B[9 + i * eq_num] * beta * dt * dx - 4.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[0 + i * eq_num]
				 + 4.0 * mesh[_x].Nk[9 + j * eq_num] * mesh[_x].Nk[0 + i * eq_num] + 3.0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[0 + j * eq_num] )
				 / ( 12.0 * beta * dt * dx );

	vectF( 9 + i * eq_num ) = -( ( mu * ( newmark_B[4 + i * eq_num] * beta * By1 * dt + mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] ) * sigma_x )
				 / ( 2.0 * beta * dt ) );

	for( int i = 1; i < nx - 1; ++i )
	{	
		if( stressType == stressCentered )
		{
			N_PRES rad2 = ( ( Km - 1 ) / 2 - _x ) * dy * ( ( Km - 1 ) / 2 - _x ) * dy + ( ( nx - 1 ) / 2 - i ) * dx * ( ( nx - 1 ) / 2 - i ) * dx;
			if( rad2 < impRadSq && cur_t < tauP )
			{
				Pimp = p0 * sqrt( 1 - rad2 / impRadSq ) * sin( M_PI * cur_t / tauP );
			}
			else if( _x == ( Km - 1 ) / 2 )
			{
				//cout << " == line " << i << " is out\n";
			}
		}

		r = i + 1;
		rr = i + 2;
		j = i - 1;
		jj = i - 2;

		matrA.coeffRef( 0 + i * eq_num, 1 + r * eq_num ) = -1.0 / ( 2.0 * dx );
		matrA.coeffRef( 0 + i * eq_num, 1 + j * eq_num ) = 1.0 / ( 2.0 * dx );
		matrA.coeffRef( 0 + i * eq_num, 2 + i * eq_num ) = 1.0 / ( h * B66 );

		matrA.coeffRef( 1 + i * eq_num, 0 + r * eq_num ) = - B12 / ( B22 * 2.0 * dx );
		matrA.coeffRef( 1 + i * eq_num, 0 + j * eq_num ) = B12 / ( B22 * 2.0 * dx );
		matrA.coeffRef( 1 + i * eq_num, 3 + i * eq_num ) = 1.0 / ( h * B22 );

		matrA.coeffRef( 2 + i * eq_num, 0 + r * eq_num ) = ( B12 * B12 / B22 - B11 ) * h / ( dx * dx );
		matrA.coeffRef( 2 + i * eq_num, 0 + j * eq_num ) = ( B12 * B12 / B22 - B11 ) * h / ( dx * dx );
		matrA.coeffRef( 2 + i * eq_num, 0 + i * eq_num ) = rho * h * 2.0 / ( Btdt * dt ) - ( B12 * B12 / B22 - B11 ) * 2.0 * h / ( dx * dx ) + h * sigma_z * By1 * By1 / ( 4.0 * Btdt );
		matrA.coeffRef( 2 + i * eq_num, 1 + r * eq_num ) = eps_x_0 * h / ( 2.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num];
		matrA.coeffRef( 2 + i * eq_num, 1 + j * eq_num ) = -eps_x_0 * h / ( 2.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num];
		matrA.coeffRef( 2 + i * eq_num, 3 + r * eq_num ) = -B12 / ( B22 * 2.0 * dx );
		matrA.coeffRef( 2 + i * eq_num, 3 + j * eq_num ) = B12 / ( B22 * 2.0 * dx );
		matrA.coeffRef( 2 + i * eq_num, 4 + r * eq_num ) = -eps_x_0 * h * By1 / ( 4.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num];
		matrA.coeffRef( 2 + i * eq_num, 4 + j * eq_num ) = eps_x_0 * h * By1 / ( 4.0 * Btdt * dx ) * mesh[_x].Nk[8 + i * eq_num];
		matrA.coeffRef( 2 + i * eq_num, 8 + i * eq_num ) = ( eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) 
				+newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] ) - eps_x_0 * h * By1 / ( 4.0 * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) 
				+newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] ) );
		matrA.coeffRef( 2 + i * eq_num, 9 + i * eq_num ) = ( h / ( mu * 2.0 * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) + eps_x_0 * h / ( 2.0 * dx ) * mesh[_x].Nk[8 + i * eq_num] 
				* ( 1.0 / Btdt * ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) 
				+ newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] ) );
		matrA.coeffRef( 2 + i * eq_num, 9 + r * eq_num ) = ( h / ( mu * 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 2 + i * eq_num, 9 + j * eq_num ) = ( -h / ( mu * 2.0 * dx ) * mesh[_x].Nk[9 + i * eq_num] );

		matrA.coeffRef( 3 + i * eq_num, 0 + r * eq_num ) = ( -B12 * eps_x_0 * h / ( 2.0 * Btdt * dx * B22 ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 3 + i * eq_num, 0 + j * eq_num ) = ( B12 * eps_x_0 * h / ( 2.0 * Btdt * dx * B22 ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 3 + i * eq_num, 1 + i * eq_num ) = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / Btdt * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 3 + i * eq_num, 2 + r * eq_num ) = -1.0 / ( 2.0 * dx );
		matrA.coeffRef( 3 + i * eq_num, 2 + j * eq_num ) = 1.0 / ( 2.0 * dx );
		matrA.coeffRef( 3 + i * eq_num, 3 + i * eq_num ) = ( eps_x_0 / ( B22 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 3 + i * eq_num, 4 + i * eq_num ) = ( -sigma_x * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 3 + i * eq_num, 5 + i * eq_num ) = ( -eps_x_0 * h * By1 / ( 2.0 * Btdt ) * mesh[_x].Nk[8 + i * eq_num] );
		matrA.coeffRef( 3 + i * eq_num, 8 + i * eq_num ) = ( mesh[_x].Nk[9 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) + newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] ) )
				- eps_x_0 * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) );
		matrA.coeffRef( 3 + i * eq_num, 9 + i * eq_num ) = ( mesh[_x].Nk[8 + i * eq_num] * ( sigma_x * h + eps_x_0 / B22 * ( 1.0 / Btdt * mesh[_x].Nk[3 + i * eq_num] + newmark_B[3 + i * eq_num] )
				- B12 * eps_x_0 * h / ( 2.0 * dx * B22 ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) + newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] ) )
				+ 2.0 * sigma_x * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				- sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[4 + i * eq_num] + newmark_B[4 + i * eq_num] ) + h * Jx );

		matrA.coeffRef( 4 + i * eq_num, 5 + i * eq_num ) = 1.0;

		matrA.coeffRef( 5 + i * eq_num, 4 + r * eq_num ) = -B12 / ( B22 * dx * dx );
		matrA.coeffRef( 5 + i * eq_num, 4 + i * eq_num ) = 2.0 * B12 / ( B22 * dx * dx );
		matrA.coeffRef( 5 + i * eq_num, 4 + j * eq_num ) = -B12 / ( B22 * dx * dx );
		matrA.coeffRef( 5 + i * eq_num, 6 + i * eq_num ) = -12.0 / ( B22 * h * h * h );

		matrA.coeffRef( 6 + i * eq_num, 4 + i * eq_num ) = ( -h * h * h * B12 / B22 * eps_x_0 / ( 6.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
		matrA.coeffRef( 6 + i * eq_num, 4 + r * eq_num ) = ( h * h * h * B12 / B22 * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
		matrA.coeffRef( 6 + i * eq_num, 4 + j * eq_num ) = ( h * h * h * B12 / B22 * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * 1.0 / Btdt );
		matrA.coeffRef( 6 + i * eq_num, 5 + i * eq_num ) = ( -rho * h * h * h / ( 6.0 * Btdt * dt ) - 2.0 * h * h * h * B66 / ( 6.0 * dx * dx )
				 - h * h * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( 12.0 * Btdt ) );
		matrA.coeffRef( 6 + i * eq_num, 5 + r * eq_num ) = h * h * h * B66 / ( 6.0 * dx * dx );
		matrA.coeffRef( 6 + i * eq_num, 5 + j * eq_num ) = h * h * h * B66 / ( 6.0 * dx * dx );
		matrA.coeffRef( 6 + i * eq_num, 6 + i * eq_num ) = ( eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( Btdt * B22 ) );
		matrA.coeffRef( 6 + i * eq_num, 7 + i * eq_num ) = 1.0;
		matrA.coeffRef( 6 + i * eq_num, 8 + i * eq_num ) = ( eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * B12 * eps_x_0 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) );
		matrA.coeffRef( 6 + i * eq_num, 9 + i * eq_num ) = ( -h * h * h / 6.0 * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 + eps_x_0 / B22 * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[6 + i * eq_num] + newmark_B[6 + i * eq_num] )
				 + h * h * h * eps_x_0 * B12 / ( 12.0 * B22 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) );

		matrA.coeffRef( 7 + i * eq_num, 1 + i * eq_num ) = ( -sigma_x * h * By1 / Btdt / 2.0 * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 7 + i * eq_num, 4 + i * eq_num ) = ( rho * h * 2.0 / ( Btdt * dt ) + sigma_x * h / ( 4.0 * Btdt ) * ( By1 * By1 + 1.0 / 3.0 * By2 * By2 ) + rho * h * h * h / ( 3.0 * dx * dx * Btdt * dt )
				 + h * h * h / 6.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) );
		matrA.coeffRef( 7 + i * eq_num, 4 + r * eq_num ) = ( -rho * h * h * h / ( 6.0 * dx * dx * Btdt * dt )
				 - h * h * h / ( 24.0 * dx * dx ) * sigma_y / Btdt * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 - h * h * h / 12.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) );
		matrA.coeffRef( 7 + i * eq_num, 4 + j * eq_num ) = ( -rho * h * h * h / ( 6.0 * dx * dx * Btdt * dt )
				 + h * h * h / ( 24.0 * dx * dx * Btdt ) * sigma_y * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
				 - h * h * h / 12.0 * ( sigma_y * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] + sigma_z * By1 * By1 / 4.0 ) / ( dx * dx * Btdt ) );              
		matrA.coeffRef( 7 + i * eq_num, 4 + i * eq_num ) = matrA.coeffRef( 7 + i * eq_num, 4 + i * eq_num ) - ( h * h * h / 2.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		matrA.coeffRef( 7 + i * eq_num, 4 + r * eq_num ) = matrA.coeffRef( 7 + i * eq_num, 4 + r * eq_num ) + ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		matrA.coeffRef( 7 + i * eq_num, 4 + j * eq_num ) = matrA.coeffRef( 7 + i * eq_num, 4 + j * eq_num ) + ( h * h * h / 3.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		if( i != nx - 2 )
		{
			matrA.coeffRef( 7 + i * eq_num, 4 + rr * eq_num ) = ( -h * h * h / 12.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );
		}
		if( i != 1 )
		{
			matrA.coeffRef( 7 + i * eq_num, 4 + jj * eq_num ) = ( -h * h * h / 12.0 * ( 2.0 * B66 * B12 / B22 - B11 + B12 * B12 / B22 ) / ( dx * dx * dx * dx ) );        
		}
		matrA.coeffRef( 7 + i * eq_num, 5 + i * eq_num ) = ( -eps_x_0 * h / Btdt * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
				 + h * h * h / 6.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt ) );
		matrA.coeffRef( 7 + i * eq_num, 5 + r * eq_num ) = ( -h * h * h / 12.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt )
				 - h * h * h / 12.0 * eps_x_0 * ( (mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] + ( mesh[_x].Nk[9 + r * eq_num]
				 - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( dx * dx * 4.0 * Btdt ) );
		matrA.coeffRef( 7 + i * eq_num, 5 + j * eq_num ) = ( -h * h * h / 12.0 * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] / ( dx * dx * Btdt )
				 + h * h * h / 12.0 * eps_x_0 * ( ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] + ( mesh[_x].Nk[9 + r * eq_num]
				 - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] ) / ( dx * dx * 4.0 * Btdt ) );
		matrA.coeffRef( 7 + i * eq_num, 6 + i * eq_num ) = ( 2.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
		matrA.coeffRef( 7 + i * eq_num, 6 + r * eq_num ) = ( -1.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
		matrA.coeffRef( 7 + i * eq_num, 6 + j * eq_num ) = ( -1.0 * ( B12 + 2.0 * B66 ) / ( dx * dx * B22 ) );
		matrA.coeffRef( 7 + i * eq_num, 8 + i * eq_num ) = ( -sigma_x * h * By1 / 2.0 - eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
				 - h * h * h * eps_x_0 / ( 12.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num]
				 + mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] + newmark_B[5 + j * eq_num] ) );
		matrA.coeffRef( 7 + i * eq_num, 8 + r * eq_num ) = ( -h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) );
		matrA.coeffRef( 7 + i * eq_num, 8 + j * eq_num ) = ( h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) );           
		matrA.coeffRef( 7 + i * eq_num, 9 + i * eq_num ) = ( -sigma_x * h * By1 / 2.0 * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] )
				 - eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
				 - h * h * h * sigma_y / ( 24.0 * dx * dx ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num]
				 - mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] )
				 - h * h * h * sigma_y * mesh[_x].Nk[9 + i * eq_num] / ( 6.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
				 + mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] )
				 - h * h * h * eps_x_0 * mesh[_x].Nk[8 + i * eq_num] / ( 12.0 * dx * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num]
				 + mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] + newmark_B[5 + j * eq_num] ) );
		matrA.coeffRef( 7 + i * eq_num, 9 + r * eq_num ) = ( -h * h * h * sigma_y / ( 24.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[4 + r * eq_num]
				 - mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] )
				 - h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * (mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) );
		matrA.coeffRef( 7 + i * eq_num, 9 + j * eq_num ) = ( h * h * h * sigma_y / ( 24.0 * dx * dx ) * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / Btdt * (mesh[_x].Nk[4 + r * eq_num]
				 - mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] ) 
				 + h * h * h * eps_x_0 / ( 48.0 * dx * dx ) * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / Btdt * ( mesh[_x].Nk[5 + r * eq_num]
				 - mesh[_x].Nk[5 + j * eq_num] ) + newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] ) ); 

		matrA.coeffRef( 8 + i * eq_num, 0 + i * eq_num ) = ( 1.0 / ( 2.0 * dx * Btdt ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) );
		matrA.coeffRef( 8 + i * eq_num, 0 + r * eq_num ) = ( 1.0 / ( 2.0 * dx * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 8 + i * eq_num, 0 + j * eq_num ) = ( -1.0 / ( 2.0 * dx * Btdt ) * mesh[_x].Nk[9 + i * eq_num] );
		matrA.coeffRef( 8 + i * eq_num, 9 + i * eq_num ) = ( 2.0 / ( sigma_y_mu * dx * dx ) + 1 / Btdt + 1.0 / ( 2.0 * dx ) * ( 1.0 / Btdt * ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] )
				 + newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] ) );
		matrA.coeffRef( 8 + i * eq_num, 9 + r * eq_num ) = ( -1.0 / ( sigma_y_mu * dx * dx ) + 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );
		matrA.coeffRef( 8 + i * eq_num, 9 + j * eq_num ) = ( -1.0 / ( sigma_y_mu * dx * dx ) - 1.0 / 2.0 / dx * ( mesh[_x].Nk[0 + i * eq_num] / Btdt + newmark_B[0 + i * eq_num] ) );

		matrA.coeffRef( 9 + i * eq_num, 1 + i * eq_num ) = sigma_x_mu / Btdt * mesh[_x].Nk[9 + i * eq_num];
		matrA.coeffRef( 9 + i * eq_num, 4 + i * eq_num ) = -sigma_x_mu * By1 / 4.0 / beta / dt;
		matrA.coeffRef( 9 + i * eq_num, 8 + i * eq_num ) = sigma_x_mu;
		matrA.coeffRef( 9 + i * eq_num, 9 + i * eq_num ) = sigma_x_mu * ( 1.0 / Btdt * mesh[_x].Nk[1 + i * eq_num] + newmark_B[1 + i * eq_num] );

		vectF( 2 + i * eq_num ) = ( - ( ( 2 * (B11 - B12 * B12 / B22) * h) / dx / dx + ( h * rho ) / ( beta * dt * dt ) + ( By1 * By1 * h * sigma_z ) / ( 8 * beta * dt ) ) ) 
			* mesh[_x].Nk[0 + i * eq_num] +	h * rho * ( newmark_A[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( beta * dt * dt) ) + ( 1.0 / 4.0 ) * By1 * By1 * h * sigma_z 
			* ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2 * beta * dt ) ) + ( ( B11 - B12 * B12 / B22 ) * h * mesh[_x].Nk[0 + r * eq_num] ) / dx / dx 
			+ ( ( B11 - B12 * B12 / B22 ) * h * mesh[_x].Nk[0 + j * eq_num] ) / dx / dx - ( ( B11 - B12 * B12 / B22 ) * h * ( -2.0 * mesh[_x].Nk[0 + i * eq_num] + mesh[_x].Nk[0 + r * eq_num] 
			+ mesh[_x].Nk[0 + j * eq_num] ) ) / dx / dx - ( h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + r * eq_num] ) / ( 2.0 * dx * mu ) 
			+ ( h * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 2.0 * dx * mu) 
			+ ( h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + j * eq_num] ) / ( 2 * dx * mu ) + ( B12 * mesh[_x].Nk[3 + r * eq_num] ) / ( 2.0 * B22 * dx )
			- ( B12 * ( mesh[_x].Nk[3 + r * eq_num] - mesh[_x].Nk[3 + j * eq_num] ) ) / ( 2 * B22 * dx ) - ( B12 * mesh[_x].Nk[3 + j * eq_num] ) / ( 2.0 * B22 * dx )
			- ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt * dx ) + ( eps_x_0 * h 
			* mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] + ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) / ( 2.0 * beta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num]) / ( 2 * dx ) + ( eps_x_0 * h * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt * dx) 
			+ ( By1 * eps_x_0 * h * mesh[_x].Nk[4 + r * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 8 * beta * dt * dx ) - ( ( eps_x_0 
			* h * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + r * eq_num] - newmark_B[1 + j * eq_num] + ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) 
			/ ( 2 * beta * dt ) ) ) / ( 2 * dx ) - ( By1 * eps_x_0 * h * ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] 
			+ ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2 * beta * dt ) ) ) / ( 4.0 * dx ) ) * mesh[_x].Nk[8 + i * eq_num] 
			- ( By1 * eps_x_0 * h * ( newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] + ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] ) / ( 2.0 * beta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 4 * dx ) - ( By1 * eps_x_0 * h * mesh[_x].Nk[4 + j * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 8 * beta * dt * dx )
			- mesh[_x].Nk[9 + i * eq_num] * ( ( h * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 2 * dx * mu ) + ( eps_x_0 * h * ( newmark_B[1 + r * eq_num] 
			- newmark_B[1 + j * eq_num] + ( mesh[_x].Nk[1 + r * eq_num] - mesh[_x].Nk[1 + j * eq_num] ) / ( 2.0 * beta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * dx ) );

		vectF( 3 + i * eq_num ) = h * Jx * mesh[_x].Nk[9 + i * eq_num] - ( ( h * rho ) / ( beta * dt * dt ) + ( h *sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] )
			/ ( 2 * beta * dt ) ) * mesh[_x].Nk[1 + i * eq_num] + h * rho * ( newmark_A[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( beta * dt * dt ) ) 
			+ h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2.0 * beta * dt ) ) 
			+ mesh[_x].Nk[2 + r * eq_num] / ( 2.0 * dx ) - ( mesh[_x].Nk[2 + r * eq_num] - mesh[_x].Nk[2 + j * eq_num] ) / ( 2.0 * dx ) - mesh[_x].Nk[2 + j * eq_num]
			/ ( 2.0 * dx ) + ( By1 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[4 + i * eq_num] ) / ( 4.0 * beta * dt ) - ( 1.0 / 2.0 ) * By1 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num]
			* ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2.0 * beta * dt ) ) + h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] + 
			( B12 * eps_x_0 * h * mesh[_x].Nk[0 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * B22 * dt * dx ) 
			- ( B12 * eps_x_0 * h * ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) / ( 2.0 * beta * dt ) ) * mesh[_x].Nk[9 + i * eq_num] 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * B22 * dx ) - ( B12 * eps_x_0 * h * mesh[_x].Nk[0 + j * eq_num] * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) 
			/ ( 4 * beta * B22 * dt * dx ) - ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[3 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * beta * B22 * dt ) 
			+ ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * (newmark_B[3 + i * eq_num] + mesh[_x].Nk[3 + i * eq_num] / ( 2.0 * beta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] ) / B22 
			+ ( By1 * eps_x_0 * h * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] ) / ( 4.0 * beta * dt ) - ( 1.0 / 2.0 ) * By1 * eps_x_0 * h 
			* ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2 * beta * dt ) ) * mesh[_x].Nk[8 + i * eq_num] - ( h * sigma_x 
			* mesh[_x].Nk[9 + i * eq_num] - ( B12 * eps_x_0 * h * ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) / ( 2 * beta * dt ) ) 
			* mesh[_x].Nk[9 + i * eq_num] ) / ( 2 * B22 * dx ) + ( eps_x_0 * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[3 + i * eq_num] + mesh[_x].Nk[3 + i * eq_num] / ( 2.0 * beta * dt ) ) )
			/ B22 - ( 1.0 / 2.0 ) * By1 * eps_x_0 * h * ( newmark_B[5 + i * eq_num] + mesh[_x].Nk[5 + i * eq_num] / ( 2.0 * beta * dt ) ) ) * mesh[_x].Nk[8 + i * eq_num] 
			- mesh[_x].Nk[9 + i * eq_num] * ( h * Jx + 2.0 * h * sigma_x * mesh[_x].Nk[9 + i * eq_num] * ( newmark_B[1 + i * eq_num] + mesh[_x].Nk[1 + i * eq_num] / ( 2.0 * beta * dt ) )
			- ( 1.0 / 2.0 ) * By1 * h * sigma_x * ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2.0 * beta * dt ) ) + h * sigma_x * mesh[_x].Nk[8 + i * eq_num] -
			( B12 * eps_x_0 * h * ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) / ( 2.0 * beta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / ( 2.0 * B22 * dx ) + ( eps_x_0 * ( newmark_B[3 + i * eq_num] + mesh[_x].Nk[3 + i * eq_num] / ( 2.0 * beta * dt ) ) 
			* mesh[_x].Nk[8 + i * eq_num] ) / B22 );

		vectF( 6 + i * eq_num ) = ( -rho * h * h * h / 12.0 * newmark_A[5 + i * eq_num] + sigma_x * h * h * h / 12.0 * ( 1.0 / beta / dt * mesh[_x].Nk[5 + i * eq_num]
			+ newmark_B[5 + i * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num]
			* mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[6 + i * eq_num] / 2.0 / beta / dt - eps_x_0 / B22 * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num]
			* ( mesh[_x].Nk[6 + i * eq_num] / 2.0 / beta / dt + newmark_B[6 + i * eq_num] ) - eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num]
			* mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx / 2.0 / beta / dt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
			+ mesh[_x].Nk[4 + j * eq_num] ) - eps_x_0 * h * h * h / 12.0 * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * B12 / B22 / dx / dx
			* ( 1.0 / 2.0 / beta / dt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num] + mesh[_x].Nk[4 + j * eq_num] )
			+ newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) ) / al;

		vectF( 7 + i * eq_num ) = ( rho * h * newmark_A[4 + i * eq_num] + Pimp + sigma_x * h / 4.0 * By1 * By1 * newmark_B[4 + i * eq_num]
			+ sigma_x * h * By1 / 4.0 / beta / dt * mesh[_x].Nk[1 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h / 2.0 / beta / dt * mesh[_x].Nk[5 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num]
			+ eps_x_0 * h * mesh[_x].Nk[8 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / beta / dt * mesh[_x].Nk[5 + i * eq_num] + newmark_B[5 + i * eq_num] ) - h / 2.0 * By1 * Jx
			- rho * h * h * h / 12.0 / dx / dx * ( newmark_A[4 + r * eq_num] - 2.0 * newmark_A[4 + i * eq_num] + newmark_A[4 + j * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / 2.0 / beta / dt * ( mesh[_x].Nk[4 + r * eq_num] - mesh[_x].Nk[4 + j * eq_num] )
			+ h * h * h / 12.0 / dx * sigma_y / 2.0 / dx * mesh[_x].Nk[9 + i * eq_num] * ( 1.0 / 2.0 / beta / dt * ( mesh[_x].Nk[4 + r * eq_num]
			- mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - newmark_B[4 + j * eq_num] ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] )
			+ h * h * h / 12.0 * sigma_y / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[9 + i * eq_num] * ( 2.0 / 2.0 / beta / dt * ( mesh[_x].Nk[4 + r * eq_num] - 2.0 * mesh[_x].Nk[4 + i * eq_num]
			+ mesh[_x].Nk[4 + j * eq_num] ) + newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * mesh[_x].Nk[8 + i * eq_num] / 2.0 / beta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * mesh[_x].Nk[9 + i * eq_num] / 2.0 / beta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * ( mesh[_x].Nk[8 + r * eq_num] - mesh[_x].Nk[8 + j * eq_num] ) * ( 1.0 / 2.0 / beta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / 4.0 / dx / dx * mesh[_x].Nk[8 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) * ( 1.0 / 2.0 / beta / dt * ( mesh[_x].Nk[5 + r * eq_num] - mesh[_x].Nk[5 + j * eq_num] )
			+ newmark_B[5 + r * eq_num] - newmark_B[5 + j * eq_num] )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / 2.0 / beta / dt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + j * eq_num] ) )
			+ h * h * h / 12.0 * eps_x_0 / dx / dx * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[8 + i * eq_num] * ( 1.0 / 2.0 / beta / dt * ( mesh[_x].Nk[5 + r * eq_num] - 2.0 * mesh[_x].Nk[5 + i * eq_num] + mesh[_x].Nk[5 + j * eq_num] )
			+ newmark_B[5 + r * eq_num] - 2.0 * newmark_B[5 + i * eq_num] + newmark_B[5 + j * eq_num] ) ) / al
			- ( By1 * By1 * h * h * h * sigma_z / 48.0 * ( newmark_B[4 + r * eq_num] - 2.0 * newmark_B[4 + i * eq_num] + newmark_B[4 + j * eq_num] ) ) / dx / dx / al;

		vectF( 8 + i * eq_num ) = newmark_B[9 + i * eq_num] + mesh[_x].Nk[9 + i * eq_num] / ( 2.0 * beta * dt ) - ( mesh[_x].Nk[0 + r * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) 
			/ ( 4.0 * beta * dt * dx ) - ( 1.0 / ( 2.0 * beta * dt ) + 2 / ( dx * dx * mu * sigma_y ) + ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] 
			+ ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] )	/ ( 2.0 * beta * dt ) ) / ( 2.0 * dx ) ) * mesh[_x].Nk[9 + i * eq_num] 
			+ ( ( newmark_B[0 + r * eq_num] - newmark_B[0 + j * eq_num] + ( mesh[_x].Nk[0 + r * eq_num] - mesh[_x].Nk[0 + j * eq_num] ) 
			/ ( 2.0 * beta * dt ) ) * mesh[_x].Nk[9 + i * eq_num] ) / ( 2 * dx ) + ( mesh[_x].Nk[0 + j * eq_num] * mesh[_x].Nk[9 + i * eq_num] ) / ( 4 * beta * dt * dx ) 
			- ( -( 1 / ( dx * dx * mu * sigma_y ) ) + ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2.0 * beta * dt ) ) / ( 2 * dx ) ) * mesh[_x].Nk[9 + r * eq_num] 
			- ( mesh[_x].Nk[0 + i * eq_num] * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 4.0 * beta * dt * dx ) + ( ( newmark_B[0 + i * eq_num] 
			+ mesh[_x].Nk[0 + i * eq_num] / ( 2 * beta * dt ) ) * ( mesh[_x].Nk[9 + r * eq_num] - mesh[_x].Nk[9 + j * eq_num] ) ) / ( 2 * dx ) - ( -( 1 / ( dx * dx * mu * sigma_y ) ) 
			- ( newmark_B[0 + i * eq_num] + mesh[_x].Nk[0 + i * eq_num] / ( 2 * beta * dt ) ) / ( 2 * dx ) ) * mesh[_x].Nk[9 + j * eq_num] 
			- ( -2.0 * mesh[_x].Nk[9 + i * eq_num] + mesh[_x].Nk[9 + r * eq_num] + mesh[_x].Nk[9 + j * eq_num] ) / ( dx * dx * mu * sigma_y );

		vectF( 9 + i * eq_num ) = -( ( mu * sigma_x * mesh[_x].Nk[9 + i * eq_num] * mesh[_x].Nk[1 + i * eq_num] ) / ( 2 * beta * dt ) ) + ( By1 * mu * sigma_x 
			* mesh[_x].Nk[4 + i * eq_num] ) / ( 4.0 * beta * dt ) - ( 1.0 / 2.0 ) * By1 * mu * sigma_x * ( newmark_B[4 + i * eq_num] + mesh[_x].Nk[4 + i * eq_num] / ( 2.0 * beta * dt ) );
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::rgkCalc( const Matrix<PL_NUM, Dynamic, Dynamic>& x, Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>* x1 )	
{
	//f1 = dy * Fhi( y, x )//see the theory, I'm not sure about these comments
	rgkF1 = dy * matrAprev * x;
	rgkF1.col( varNum / 2 ) += dy * vectFprev;	

	//f2 = dy * Fhi( y, x + d21 * f1 )
	rgkF2 = dy * matrAprev * ( x + rgkD21 *  rgkF1 );
	rgkF2.col( varNum / 2 ) += dy * vectFprev;	

	//f3 = dy * Fhi( y, x + d31 * f1 + d32 * f2 )
	rgkF3 = dy * matrAprev * ( x + rgkD31 * rgkF1 + rgkD32 * rgkF2 );
	rgkF3.col( varNum / 2 ) += dy * vectFprev;	

	//f4 = dy * Fhi( y + dy, x + d41 * f1 + d42 * f2 + d43 * f3 )
	rgkF4 = dy * matrA * ( x + rgkD41 * rgkF1 + rgkD42 * rgkF2 + rgkD43 * rgkF3 );
	rgkF4.col( varNum / 2 ) += dy * vectF;

	( *x1 ) = x + rgkC1 * rgkF1 + rgkC2 * rgkF2 + rgkC3 * rgkF3 + rgkC4 * rgkF4;
}

template<class PL_NUM>
void Solver<PL_NUM>::walkthrough( int mode )	//sequential version
{
	//Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2> orthoCheck1;

	time_t tBeg;
	time_t rgkT = 0;
	time_t orthoT = 0;

	calcNewmarkAB( 0, mode );
	calc_system( 0 );

	orthoBuilder->resetOrthoDoneInfo();

	int active = 1;		//how many nodes in a row we have, on which orthonormalization was not performed. 
						//We need this to know whem we can switch to ABM method

	//integrate and orthonorm
	int PhiInd = 0;
	int totalRgkSteps = 0;
	tBeg = 0;

	for( int _x = 0; _x < Km - 1; ++_x )
	{
		int rgkIsDone = 0;

		matrAprev = matrA;
		vectFprev = vectF;

		calcNewmarkAB( _x + 1, mode );
		calc_system( _x + 1 );
		tBeg = time( 0 );

		tempPhi = matrAprev * orthoBuilder->zi[_x];
		tempPhi.col( varNum / 2 ) += vectFprev;

		PhiInd = 0;
		if( active <= ABM_STAGE_NUM )
		{
			PhiInd = active - 1;
		}
		else
		{
			for( int stage = 0; stage < ABM_STAGE_NUM - 1; ++stage )
			{
				Phi[stage] = Phi[stage + 1];
			}
			PhiInd = ABM_STAGE_NUM - 1;
		}
		Phi[PhiInd] = tempPhi;

		if( active >= ABM_STAGE_NUM )
		{
			//use ABM method
			//predictor
			tempPhi = orthoBuilder->zi[_x] + dy / 24.0l * ( 55.0l * Phi[3] - 59.0l * Phi[2] + 37.0l * Phi[1] - 9.0l * Phi[0] );
			//corrector
			tempPhi2 = matrA * tempPhi;
			tempPhi2.col( varNum / 2 ) += vectF;
			//assignment
			decompVect = orthoBuilder->zi[_x] + dy / 24.0l * ( 9.0l * tempPhi2 + 19.0l * Phi[3] - 5.0l * Phi[2] + Phi[1] );
		}
		else
		{
			rgkCalc( orthoBuilder->zi[_x], &decompVect );
		}

		rgkT += time( 0 ) - tBeg;

		//tBeg = time( 0 );
		if( rgkIsDone )
		{
			++totalRgkSteps;
		}

		for( int i = 0; i < EQ_NUM * NUMBER_OF_LINES / 2 + 1; ++i )
		{
			for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
			{
				decompVectOrtho[i][j] = decompVect( j, i );
			}
		}

		//decompVectOrtho = decompVect;

		tBeg = time( 0 );
		orthoBuilder->orthonorm( _x, decompVectOrtho );
		orthoT += time( 0 ) - tBeg;

		if( orthoBuilder->checkOrtho( _x, decompVectOrtho, decompVect ) == 1 )			
		{
			active = 1;		//if orthonormalization has been performed, we have to restart the ABM method
			for( int i = 0; i < EQ_NUM * NUMBER_OF_LINES / 2 + 1; ++i )
			{
				for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
				{
					decompVect( j, i ) = decompVectOrtho[i][j];
				}
			}
			orthoBuilder->setOrthoDoneInfo( _x );
			orthoBuilder->setNextSolVects( _x, decompVectOrtho );
			//cout << " --- at x = " << _x << " ortho is needed\n";
		}
		else
		{
			++active;	//if no orthonormalization has been done, we have one more solution that can be used in ABM method
			orthoBuilder->zi[_x + 1] = decompVect;
		}

	}

	cout << " == rgkT \t" << rgkT << endl;
	cout << " == orthoT \t" << orthoT << endl;
	cout << " == rgk to total ratio: " << ( float )totalRgkSteps / ( float )( Km - 1 ) << endl;

	orthoBuilder->buildSolution( &mesh );

	orthoBuilder->setOmegasZero();
}

template<class PL_NUM>
void Solver<PL_NUM>::updateDerivs()
{
	for( int i = 0; i < Km; ++i )
	{
		for( int j = 0; j < varNum; ++j )
		{
			mesh[i].d2N[j] = ( mesh[i].Nk1[j] - mesh[i].Nk0[j] ) / beta / dt / dt - mesh[i].d1N0[j] / beta / dt - ( 0.5 - beta ) / beta * mesh[i].d2N0[j];
			mesh[i].d1N[j] = mesh[i].d1N0[j] + 0.5 * dt * ( mesh[i].d2N0[j] + mesh[i].d2N[j] );
			mesh[i].Nk0[j] = mesh[i].Nk1[j];
			mesh[i].d1N0[j] = mesh[i].d1N[j];
			mesh[i].d2N0[j] = mesh[i].d2N[j];
		}
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::pre_step()
{
	for( int i = 0; i < nx; ++i )			//TODO check indexes
	{
		orthoBuilder->zi[0]( i * eq_num + 2, i * eq_num / 2 + 0 ) = 1.0;
		orthoBuilder->zi[0]( i * eq_num + 3, i * eq_num / 2 + 1 ) = 1.0;
		orthoBuilder->zi[0]( i * eq_num + 5, i * eq_num / 2 + 2 ) = 1.0;
		orthoBuilder->zi[0]( i * eq_num + 7, i * eq_num / 2 + 3 ) = 1.0;
		orthoBuilder->zi[0]( ( i + 1 ) * eq_num - 1, ( i + 1 ) * eq_num / 2 - 1 ) = -1.0;
	}
	//z5s are already zeros

	calcNewmarkAB( 0, 0 );
	calc_system( 0 );

	walkthrough( 0 );
}

template<class PL_NUM>
PL_NUM Solver<PL_NUM>::do_step()
{	
	int cont = 1;
	while( cont == 1 )
	{
		cout << " = walk\n";
		calcNewmarkAB( 0, 1 );

		orthoBuilder->zi[0] = Matrix<PL_NUM, Dynamic, Dynamic>::Zero( EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1 );

		for( int i = 0; i < nx; ++i )			//TODO check indexes
		{
			orthoBuilder->zi[0]( i * eq_num + 2, i * eq_num / 2 + 0 ) = 1.0;
			orthoBuilder->zi[0]( i * eq_num + 3, i * eq_num / 2 + 1 ) = 1.0;
			orthoBuilder->zi[0]( i * eq_num + 5, i * eq_num / 2 + 2 ) = 1.0;
			orthoBuilder->zi[0]( i * eq_num + 7, i * eq_num / 2 + 3 ) = 1.0;
			orthoBuilder->zi[0]( ( i + 1 ) * eq_num - 2, ( i + 1 ) * eq_num / 2 - 1 ) = mesh[0].Nk1[i * eq_num + 1] / beta / 2.0 / dt + newmark_B[i * eq_num + 1];
			orthoBuilder->zi[0]( ( i + 1 ) * eq_num - 1, ( i + 1 ) * eq_num / 2 - 1 ) = -1.0;

			orthoBuilder->zi[0]( i * eq_num + 8, varNum / 2 ) = newmark_B[i * eq_num + 1] * mesh[0].Nk1[i * eq_num + 9] - newmark_B[i * eq_num + 4] * By0;
			orthoBuilder->zi[0]( i * eq_num + 9, varNum / 2 ) = -mesh[0].Nk1[i * eq_num + 9];
		}

		for( int i = 0; i < Km; ++i )
		{
			for( int j = 0; j < varNum; ++j )
			{
				mesh[i].Nk[j] = mesh[i].Nk1[j];
			}
		}
		walkthrough( 0 );
		cont = checkConv();
	}
	updateDerivs();

	return mesh[( Km - 1 ) / 2].Nk1[4 + (nx-1)/2 * eq_num];
}

template<class PL_NUM>
int Solver<PL_NUM>::checkConv()
{
	//if( newtonIt >= maxNewtonIt )		//stopping criterion with a fixed number of iterations
	//{
	//	newtonIt = 0;
	//	return 0;
	//}
	//else
	//{
	//	int i = 15;
	//	cout << " newton iteration " << newtonIt << endl;
	//	cout << " divergence in " << Km / 2 << " " << i << " " << fabsl( ( mesh[Km / 2].Nk1[i] - mesh[Km / 2].Nk[i] ) / mesh[Km / 2].Nk[i] ) << endl;
	//	++newtonIt;
	//	return 1;
	//}

	for( int x = 0; x < Km; ++x ) //old and weird stopping criterion. I think it works only because 1-2 iterations are almost always suffice. otherwise this doed not make sense to me
	{
		for( int i = 0; i < varNum; ++i )
		{
			if( mesh[x].Nk[i] != 0.0 )
			{
				if( fabs( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) < ALMOST_ZERO )
				{
					cout << " divergence " << x << " " << i << " " << fabs( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) << " delta is " << ALMOST_ZERO << endl;
					return 0;
				}
			}
			else
			{
				if( fabs( mesh[x].Nk1[i] ) < ALMOST_ZERO )
				{
					cout << " divergence " << x << " " << i << " " << fabs( mesh[x].Nk1[i] ) << " delta is " << ALMOST_ZERO << endl;
					return 0;
				}
			}
		}
	}
	return 1;

	//PL_NUM maxDiff = fabsl( mesh[0].Nk1[0] - mesh[0].Nk[0] );
	//for( int y = 1; y < Km; ++y )
	//{
	//	for( int i = 0; i < varNum; ++i )
	//	{
	//		if( fabsl( mesh[y].Nk1[i] - mesh[y].Nk[i] ) > maxDiff )
	//		{
	//			maxDiff = fabsl( mesh[y].Nk1[i] - mesh[y].Nk[i] );
	//		}
	//	}
	//}
	//if( prevVectDiff < 0.0 )
	//{
	//	prevVectDiff = maxDiff;
	//}
	//else
	//{
	//	PL_NUM ratio = maxDiff / prevVectDiff;
	//	if( ratio > 0.99 )
	//	{
	//		cout << " --- diverging: " << ratio << endl;
	//	}
	//	else
	//	{
	//		cout << " --- converging: " << ratio << endl;
	//	}
	//	prevVectDiff = maxDiff;
	//}

	//for( int y = 0; y < Km; ++y ) //new stopping criterion: the max of the absolute value of the relative difference + just the abs value check
	//{
	//	for( int i = 0; i < varNum; ++i )
	//	{
	//		if( fabsl( mesh[y].Nk[i] ) > ALMOST_ZERO * 1000.0 )
	//		{
	//			if( fabsl( ( mesh[y].Nk1[i] - mesh[y].Nk[i] ) / mesh[y].Nk[i] ) > QUASILIN_CHECK && fabsl( mesh[y].Nk1[i] - mesh[y].Nk[i] ) > QUASILIN_CHECK )
	//			{
	//				cout << " :: divergence " << y << " " << i << " " << fabsl( ( mesh[y].Nk1[i] - mesh[y].Nk[i] ) / mesh[y].Nk[i] ) << " delta is " << QUASILIN_CHECK << endl;
	//				cout << " :: " << mesh[y].Nk1[i] << " " << mesh[y].Nk[i] << endl;
	//				return 1;
	//			}
	//		}
	//		else
	//		{
	//			if( fabsl( mesh[y].Nk1[i] ) > ALMOST_ZERO * 1000.0 )
	//			{
	//				cout << " :: divergence -- 0 -- " << y << " " << i << " " << fabsl( mesh[y].Nk1[i] ) << " delta is " << QUASILIN_CHECK << endl;
	//				return 1;
	//			}
	//		}
	//	}
	//}
	//return 0;
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_sol()
{
	ofstream dumpSol;
	stringstream ss;
	ss << "sol_" << cur_t << ".txt";
	
	dumpSol.open ( ss.str() );
	
	for( int x = 0; x < Km; ++x )
	{
		for( int line = 0; line < nx; ++line )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				dumpSol << mesh[x].Nk1[line * eq_num + i] << " ";
			}
			dumpSol << "\n";
		}
		dumpSol << "\n\n";
	}

	dumpSol.close();
	return;
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_check_sol()	//dump numerical soln + the soln obtained analytically for 1d case
{
	N_PRES sum = 0.0;
	N_PRES h = hp;
	N_PRES a = bp;
	N_PRES t = cur_t - dt;

	int minusOne = -1;

	/*for( int i = 0; i <= 1000000; ++i )
	{
		PL_NUM omg = _MMM_PI * _MMM_PI * ( 2 * i + 1 ) * ( 2 * i + 1 ) * h / 2 / a / a * sqrtl( B22 / 3 / rho );

		minusOne = -minusOne;

		sum = sum + minusOne / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) * cosl( omg * t );
	}*/
	N_PRES wTheor;
	wTheor = - p0 * a * a * a * a / h / h / h / B22 * ( 5.0 / 32.0 - 48.0 / _MMM_PI / _MMM_PI / _MMM_PI / _MMM_PI / _MMM_PI * sum );
	wTheor = 1.0;

	ofstream of1( "test_sol.txt", ofstream::app );
	of1 << t << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] << " ; " << wTheor << " ; " << fabs( ( wTheor - mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] ) / wTheor ) << endl;
	of1.close();
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_check_sol2D()		//dump numerical soln + the soln obtained analytically for 2d case. WARNING: the plate should be square-shaped and isotropic !!!
{
	N_PRES sum = 0.0;
	N_PRES sum0 = 0.0;
	N_PRES h = hp;
	N_PRES a = bp;
	N_PRES t = cur_t - dt;

	N_PRES DD = E1 * h * h * h / 12 / ( 1 - nu21 * nu21 );
	N_PRES Om = 100.0 * _MMM_PI;

	for( int m = 1; m < 50; ++m )
	{
		for( int n = 1; n < 50; ++n )
		{
			N_PRES wmn = _MMM_PI * _MMM_PI * ( m * m / ap / ap + n * n / bp / bp ) * sqrt( DD / rho / h );
			N_PRES Wmn = -16.0 * p0 / _MMM_PI / _MMM_PI / m / n / rho / h / ( wmn * wmn - Om * Om );
			sum += ( sin( Om * t ) - Om / wmn * sin( wmn * t ) ) * Wmn * sin( m * _MMM_PI / 2.0 ) * sin( n * _MMM_PI / 2.0 );
			sum0 = sum;
		}
	}
	N_PRES wTheor = sum;

	ofstream of1( "test_sol.txt", ofstream::app );
	of1 << t << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] << " ; " << wTheor << " ; " << fabs( ( wTheor - mesh[ ( Km - 1 ) / 2 ].Nk1[4 + (nx-1)/2 * eq_num] ) / wTheor ) << endl;
	of1.close();
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_whole_sol( int var )
{
	stringstream ss;
	ss << "./res/" << var << "_sol_whole_" << curTimeStep << ".bin";
	ofstream of1( ss.str(), ofstream::out | ofstream::binary );
	for( int y = 0; y < Km; ++y )
	{
		for( int x = 0; x < nx; ++x )
		{
			of1.write( reinterpret_cast<char*>( &(mesh[y].Nk1[var + x * eq_num]) ), sizeof(PL_NUM) );
		}
	}
	of1.close();
}

#endif