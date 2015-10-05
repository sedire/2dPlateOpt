#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include <time.h>
#include "Optimizer.h"
#include "hyperDual.h"
#include "plate_var_types.h"
#include <Eigen/Eigen>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";
	
	time_t beginTime;
	time_t endTime;

	N_PRES x[2];
	x[0] = 1.e5 / J0_SCALE;
	x[1] = 0.01;

	N_PRES val = calcValTaus( x, 0 );

	//Solver<HPD<N_PRES, 2> >* solver = new Solver<HPD<N_PRES, 2> >( 102970000000.0, 7550000000.0, 1.0, 0.1524 * 2.0, 0.1524 );

	//vector<HPD<N_PRES, 2> > currentParms;
	//currentParms.push_back( 1.e5 );
	//currentParms.push_back( 0.01 );
	//currentParms[0].elems[1] = 1.0;
	//currentParms[1].elems[2] = 1.0;
	//solver->setCurrentParms( currentSin, 1.0, currentParms );
	//
	//vector<N_PRES> stressParms;
	//stressParms.push_back( 30000.0 );
	////stressParms.push_back( 100.0 );
	//stressParms.push_back( 0.01 );
	//stressParms.push_back( 8.0 );
	//solver->setStressParms( stressCentered, stressParms );

	//cout << "\n doing pre_step...\n";
	//beginTime = time( 0 );
	//solver->pre_step();
	//endTime = time( 0 );
	//cout << "\n pre_step done\n";
	//cout << "pre_step done in " << float( endTime - beginTime ) << " ~~\n";

	//HPD<N_PRES, 2> sum = 0.0;
	//HPD<N_PRES, 2> res = 0.0;

	//while( solver->getCurTime() <= 0.03 )
	//{
	//	for( int i = 0; i < 1; ++i )
	//	{
	//		beginTime = time( 0 );
	//		res = solver->do_step();
	//		endTime = time( 0 );
	//		cout << "step done in " << float( endTime - beginTime ) << " ~~\n";
	//		cout << "solver returned " << res << "\n\n";

	//		sum += res * res;

	//		solver->increaseTime();
	//		cout << solver->getCurTime() << " -- step done\n";
	//	}
	//	//solver->dump_whole_sol( 0 );
	//	//solver->dump_whole_sol( 1 );
	//	//solver->dump_whole_sol( 2 );
	//	//solver->dump_whole_sol( 3 );
	//	//solver->dump_whole_sol( 4 );
	//	//solver->dump_whole_sol( 5 );
	//	//solver->dump_whole_sol( 6 );
	//	//solver->dump_whole_sol( 7 );
	//	//solver->dump_whole_sol( 8 );
	//	//solver->dump_whole_sol( 9 );
	//	//solver->dump_sol();
	//	//solver->dump_check_sol();
	//	solver->dump_check_sol2D();
	//}
	//
	//delete( solver );

	cout << "\nobj func is " << val << endl;
	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}

