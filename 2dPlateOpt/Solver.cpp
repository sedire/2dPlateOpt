#include "Solver.h"

template<>
void Solver<HPD<N_PRES, 2> >::calc_system( int _x )
{
	N_PRES h = hp;
	N_PRES Btdt = 2 * dt * beta;

	HPD<N_PRES, 2> Jx = 0.0;
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
	//HPD<N_PRES, 2> cur_X = _x * dy - bp / 2.0;
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
