#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include <omp.h>
#include <time.h>
#include "hyperDual.h"
#include <Eigen/Eigen>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";
	
	time_t beginTime;
	time_t endTime;

	//Matrix<HPD<double, 3>, 3, 1> v1;
	//Matrix<HPD<double, 3>, 3, 1> v2;

	//HPD<double, 3> a1, a2, a3, a4, a5, a6;
	//a1 = 1.0;
	//a1.elems[1] = 1.0;
	//a2 = 2.0;
	//a2.elems[2] = 1.0;
	//a3 = 3.0;
	//a3.elems[3] = 1.0;
	//a4 = 4.0;
	//a4.elems[1] = 0.5;
	//a5 = 5.0;
	//a5.elems[2] = 2.0;
	//a6 = 6.0;
	//a6.elems[3] = 3.0;

	//v1 << a1, a2, a3;
	//v2 << a4, a5, a6;

	Solver<HPD<N_PRES, 1> >* solver = new Solver<HPD<N_PRES, 1> >();

	cout << "\n doing pre_step...\n";
	beginTime = time( 0 );
	solver->pre_step();
	solver->dump_sol();
	endTime = time( 0 );

	cout << "\n pre_step done\n";
	cout << "pre_step done in " << float( endTime - beginTime ) << " ~~\n";

	while( solver->getCurTime() <= 0.05 )
	{
		for( int i = 0; i < 1; ++i )
		{
			beginTime = time( 0 );
			solver->do_step();
			endTime = time( 0 );
			cout << "step done in " << float( endTime - beginTime ) << " ~~\n";

			solver->increaseTime();
			cout << solver->getCurTime() << " -- step done\n";
		}
		//solver->dump_whole_sol( 0 );
		//solver->dump_whole_sol( 1 );
		//solver->dump_whole_sol( 2 );
		//solver->dump_whole_sol( 3 );
		//solver->dump_whole_sol( 4 );
		//solver->dump_whole_sol( 5 );
		//solver->dump_whole_sol( 6 );
		//solver->dump_whole_sol( 7 );
		//solver->dump_whole_sol( 8 );
		//solver->dump_whole_sol( 9 );
		//solver->dump_sol();
		//solver->dump_check_sol();
		solver->dump_check_sol2D();
	}
	
	delete( solver );

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}