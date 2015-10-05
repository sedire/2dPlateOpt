#include "Optimizer.h"

void optimizeASA( const vector<N_PRES>& params )
{
	time_t totOptStart = time( 0 );
	cout << "optimizeASA enter\n\n";

	const N_PRES threshold( 1.e-6 );

	N_PRES* x = new N_PRES[GRAD_SIZE];
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		x[i] = params[i];
	}

	N_PRES* lo = new N_PRES[GRAD_SIZE];
	N_PRES* hi = new N_PRES[GRAD_SIZE];
	lo[0] = -1.0;
	lo[1] = 0.002;
	hi[0] = 1.0;
	hi[1] = 1.e9;

	asacg_parm cgParm;
    asa_parm asaParm;
	asa_cg_default( &cgParm );
    asa_default( &asaParm );
    cgParm.PrintParms = TRUE;
    cgParm.PrintLevel = 3;
    asaParm.PrintParms = TRUE;
    asaParm.PrintLevel = 3;

	asa_cg( x, lo, hi, GRAD_SIZE, NULL, &cgParm, &asaParm, threshold, calcValASAcurSin, calcGradASAcurSin, calcValGradASAcurSin, 0, 0 );

	cout << "\n\n===============\nASA optimization complete. X is:\n";
	ofstream of( "solution.txt" );
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		of << x[i] << endl;
		cout << x[i] << endl;
	}
	of.close();

	cout << " total optimization time : " << time( 0 ) - totOptStart << endl;

	delete[] x;
	delete[] lo;
	delete[] hi;
}

N_PRES calcValTaus( N_PRES* x, long n )
{
	N_PRES ret = 0.0;
	time_t begin = time( 0 );

	cout << "try to calc at\n";
	vector<N_PRES> currentParms;
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		cout << "\t" << x[i] << endl;
		currentParms.push_back( x[i] );
	}
	cout << " ====\n";

	Solver<N_PRES>* solver = new Solver<N_PRES>( SOLVER_E1, SOLVER_E2, SOLVER_BY0, SOLVER_AP, SOLVER_BP );
	solver->setCurrentParms( currentSin, J0_SCALE, currentParms );

	vector<N_PRES> stressParms;
	stressParms.push_back( STRESS_P0 );
	stressParms.push_back( STRESS_TAU );
	stressParms.push_back( STRESS_RAD_FACTOR );
	solver->setStressParms( stressCentered, stressParms );

	cout << "\tcalculating func val\n";

	const N_PRES charTime( CHAR_TIME );

	N_PRES res = 0.0;
	N_PRES sum = 0.0;

	while( solver->getCurTime() <= charTime )
	{
		res = solver->do_step();
		sum += res * res; 

		solver->increaseTime();
	}

	cout << "\tfunc val done\n";
	cout << "\tval is " << sum << endl;
	cout << " -------------\n";

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return sum;
}

N_PRES calcValGradTaus( N_PRES* g, N_PRES* x, long n )
{
	N_PRES ret = 0.0;
	time_t begin = time( 0 );

	cout << "try to calc at\n";
	vector<HPD<N_PRES, GRAD_SIZE> > currentParms;
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		cout << "\t" << x[i] << endl;
		currentParms.push_back( x[i] );
		currentParms[i].elems[i + 1] = 1.0;
	}
	cout << " ====\n";

	Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >( SOLVER_E1, SOLVER_E2, SOLVER_BY0, SOLVER_AP, SOLVER_BP );
	solver->setCurrentParms( currentSin, J0_SCALE, currentParms );

	vector<N_PRES> stressParms;
	stressParms.push_back( STRESS_P0 );
	stressParms.push_back( STRESS_TAU );
	stressParms.push_back( STRESS_RAD_FACTOR );
	solver->setStressParms( stressCentered, stressParms );

	cout << "\tcalculating func val\n";

	const N_PRES charTime( CHAR_TIME );

	HPD<N_PRES, GRAD_SIZE> res = 0.0;
	HPD<N_PRES, GRAD_SIZE> sum = 0.0;

	while( solver->getCurTime() <= charTime )
	{
		res = solver->do_step();
		sum += res * res; 

		solver->increaseTime();
		cout << solver->getCurTime() << " -- step done\n";
	}

	cout << "\tfunc val done\n";
	cout << "\tval is " << sum << endl;
	cout << " -------------\n";

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	if( g != 0 )
	{
		for( int i = 0; i < GRAD_SIZE; ++i )
		{
			g[i] = sum.elems[i + 1];
		}
	}

	return sum.real();
}

N_PRES calcValASAcurSin( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Val\n";
	return calcValTaus( asa->x, asa->n );
}

void calcGradASAcurSin( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Grad\n";
	calcValGradTaus( asa->g, asa->x, asa->n );
}

N_PRES calcValGradASAcurSin( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Both\n";
	N_PRES ret = calcValGradTaus( asa->g, asa->x, asa->n );
	return ret;
}
