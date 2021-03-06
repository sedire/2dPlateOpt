#include "OrthoBuilder.h"

//SolInfo::SolInfo() {}
//SolInfo::~SolInfo() {}

//void SolInfo::setup( int _varNum )
//{
//	o.resize( ( _varNum / 2 + 1 ) * ( _varNum / 2 + 2 ) / 2, 0.0 );		//the matrix is upper triangular. such size is to store less
//}
//
//void SolInfo::flushO()
//{
//	for( int i = 0; i < o.size(); ++i )
//	{
//		o[i] = 0.0;
//	}
//}

OrthoBuilder::~OrthoBuilder() 
{ 
	delete[] zi; 
}

OrthoBuilder::OrthoBuilder( const int _varNum, const int _Km ) :
	varNum( _varNum ),
	Km( _Km ),
	zi( 0 )
{
	zi = new Matrix<PL_NUM, Dynamic, Dynamic>[_Km];
	for( int i = 0; i < _Km; ++i )
	{
		zi[i] = Matrix<PL_NUM, Dynamic, Dynamic>::Zero( _varNum, _varNum / 2 + 1 );
	}
	omega.resize( _Km, vector<PL_NUM>( ( _varNum / 2 + 1 ) * ( _varNum / 2 + 2 ) / 2, 0.0 ) );
}

OrthoBuilderGSh::OrthoBuilderGSh( int _varNum, int _Km ) :
	OrthoBuilder( _varNum, _Km )
{

}

//void OrthoBuilder::flushO( int x )
//{
//	solInfoMap[x].flushO();
//}

void OrthoBuilder::setOmegasZero()
{
	for( int i = 0; i < Km; ++i )
	{
		for( int j = 0; j < ( varNum / 2 + 1 ) * ( varNum / 2 + 2 ) / 2; ++j )
		{
			omega[i][j] = 0.0;
		}
	}
}

void OrthoBuilder::setParams()
{
	//try 
	//{
	//	solInfoMap.resize( Km );
	//	for( int i = 0; i < solInfoMap.size(); ++i )
	//	{
	//		solInfoMap[i].setup( varNum );
	//	}
	//}
	//catch( bad_alloc &ba )
	//{
	//	cout << ba.what() << endl;
	//	std::abort();
	//}

	orthoDone.resize( Km - 1, 0 );
}

void OrthoBuilder::setInitVects( const vector<PL_NUM>& N1, const vector<PL_NUM>& N2, const vector<PL_NUM>& N3, const vector<PL_NUM>& N4, const vector<PL_NUM>& N5 )
{
	cout << "Warning! I am void!\n";
}

void OrthoBuilder::setOrthoDoneInfo( int y )
{
	orthoDone[y] = true;
}

void OrthoBuilder::resetOrthoDoneInfo()
{
	for( int y = 0; y < orthoDone.size(); ++y )
	{
		orthoDone[y] = false;
	}
}

inline void OrthoBuilderGSh::setNextSolVects( int n, const PL_NUM decompVect[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] )
{
	for( int vNum = 0; vNum < EQ_NUM * NUMBER_OF_LINES / 2 + 1; ++vNum )
	{
		for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
		{
			zi[n + 1]( j, vNum ) = decompVect[vNum][j];
		}
	}
}

inline PL_NUM OrthoBuilder::getInfNorm( PL_NUM* vect, int vectSize )
{
	if( vect != 0 )
	{
		PL_NUM ret = fabs( vect[0] );
		for( int i = 1; i < vectSize; ++i )
		{
			if( fabs( vect[i] ) > ret )
			{
				ret = fabs( vect[i] );
			}
		}
		return ret;
	}
	else
	{
		return -1;
	}
}

inline int OrthoBuilderGSh::checkOrtho( int n, 
									PL_NUM vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES], 
									//const Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>& vectSetOrtho,
									const Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>& vectSetOrig  )
{
	int ret = 0;
	PL_NUM eps = ORTHONORM_CHECK_EPS;

	//check for the main basis vectors
	for( int vNum = 1; vNum < EQ_NUM * NUMBER_OF_LINES / 2; ++vNum )
	{
		if( getInfNorm( vectSetOrtho[vNum], EQ_NUM * NUMBER_OF_LINES )
			/*vectSetOrtho.col( vNum ).lpNorm<Infinity>() * solInfoMap[n + 1].o[vNum * ( vNum + 1 ) / 2 + vNum]*/ <
			eps * vectSetOrig.col( vNum ).lpNorm<Infinity>() )
		{
			ret = 1;
			break;
		}
	}

	if( ret != 1 )	//check for the vector of particular solution
	{
		if( getInfNorm( vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2], EQ_NUM * NUMBER_OF_LINES )
			/*vectSetOrtho.col( varNum / 2 ).lpNorm<Infinity>()*/ <
			eps * vectSetOrig.col( varNum / 2 ).lpNorm<Infinity>() )
		{
			ret = 1;
		}
	}

	return ret;
}

void OrthoBuilderGSh::orthonorm( int y, PL_NUM NtoOrt[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] )
{
	PL_NUM k11 = sqrtl( 2.0 );
	PL_NUM norm( 0.0 );

	//theory is on p.45-46 of Scott, Watts article
	for( int baseVNum = 0; baseVNum < varNum / 2; ++baseVNum )
	{
		norm = 0.0;
		for( int i = 0; i < varNum; ++i )
		{
			norm += NtoOrt[baseVNum][i] * NtoOrt[baseVNum][i];
		}
		norm = sqrtl( norm );
		for( int i = 0; i < varNum / 2; ++i )
		{
			omega2[i] = 0.0;
		}

		for( int bvIt = 0; bvIt < baseVNum; ++bvIt )
		{
			for( int k = 0; k < varNum; ++k )
			{
				//solInfoMap[y + 1].o[baseVNum * ( baseVNum + 1 ) / 2 + bvIt] += NtoOrt[baseVNum][k] * NtoOrt[bvIt][k];			
				omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + bvIt] += NtoOrt[baseVNum][k] * NtoOrt[bvIt][k];			
			}
			for( int k = 0; k < varNum; ++k )
			{
				//NtoOrt[baseVNum][k] -= solInfoMap[y + 1].o[baseVNum * ( baseVNum + 1 ) / 2 + bvIt] * NtoOrt[bvIt][k];
				NtoOrt[baseVNum][k] -= omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + bvIt] * NtoOrt[bvIt][k];
			}
		}

		for( int k = 0; k < varNum; ++k )
		{
			omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + baseVNum] += NtoOrt[baseVNum][k] * NtoOrt[baseVNum][k];			
		}
		omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + baseVNum] = sqrtl( omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + baseVNum] );

		if( norm / omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + baseVNum] <= k11 )
		{
			for( int k = 0; k < varNum; ++k )
			{
				NtoOrt[baseVNum][k] /= omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + baseVNum];
			}
		}
		else
		{
			//cout << " === no ortho!!\n";
			for( int bvIt = 0; bvIt < baseVNum; ++bvIt )
			{
				for( int k = 0; k < varNum; ++k )
				{
					omega2[bvIt] += NtoOrt[baseVNum][k] * NtoOrt[bvIt][k];
				}
				for( int k = 0; k < varNum; ++k )
				{
					NtoOrt[baseVNum][k] -= omega2[bvIt] * NtoOrt[bvIt][k];
				}
			}
			for( int k = 0; k < varNum; ++k )
			{
				omega2[baseVNum] += NtoOrt[baseVNum][k] * NtoOrt[baseVNum][k];
			}
			omega2[baseVNum] = sqrtl( omega2[baseVNum] );
			for( int k = 0; k < varNum; ++k )
			{
				NtoOrt[baseVNum][k] /= omega2[baseVNum];
			}
			for( int bvIt = 0; bvIt < baseVNum; ++bvIt )
			{
				omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + bvIt] += omega2[bvIt];
			}
			omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + baseVNum] = omega2[baseVNum];
		}
	}
	//orthogonalization of the particular solution
	int partSolNum = varNum / 2;
	for( int bvIt = 0; bvIt < varNum / 2; ++bvIt )
	{
		for( int k = 0; k < varNum; ++k )
		{
			omega[y + 1][partSolNum * ( partSolNum + 1 ) / 2 + bvIt] += NtoOrt[partSolNum][k] * NtoOrt[bvIt][k];
		}
		for( int k = 0; k < varNum; ++k )
		{
			NtoOrt[partSolNum][k] -= omega[y + 1][partSolNum * ( partSolNum + 1 ) / 2 + bvIt] * NtoOrt[bvIt][k];
		}
	}
}

void OrthoBuilderGSh::orthonorm( int y, Matrix<PL_NUM, EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>* NtoOrt )
{
	qr.compute( *NtoOrt );

	qrQ = qr.householderQ();
	( *NtoOrt ) = qrQ.block<EQ_NUM * NUMBER_OF_LINES, EQ_NUM * NUMBER_OF_LINES / 2 + 1>( 0, 0 );
	
	qrR = qr.matrixQR().triangularView<Upper>();

	( *NtoOrt ).col( varNum / 2 ) *= qrR( varNum / 2, varNum / 2 );

	for( int baseVNum = 0; baseVNum < varNum / 2; ++baseVNum )
	{
		for( int bvIt = 0; bvIt < baseVNum + 1; ++bvIt )
		{
			omega[y + 1][baseVNum * ( baseVNum + 1 ) / 2 + bvIt] = qrR( bvIt, baseVNum );
		}
	}
	for( int bvIt = 0; bvIt < varNum / 2; ++bvIt )
	{
		omega[y + 1][varNum / 2 * ( varNum / 2 + 1 ) / 2 + bvIt] = qrR( bvIt, varNum / 2 );
	}
}

void OrthoBuilderGSh::buildSolution( vector<VarVect>* _mesh )
{
	static const int msize = NUMBER_OF_LINES * EQ_NUM / 2;
	Matrix<PL_NUM, msize, msize, RowMajor> M;
	Matrix<PL_NUM, msize, 1> f11;
	Matrix<PL_NUM, msize, 1> x1;
	Matrix<PL_NUM, msize, 1> res;
	Matrix<PL_NUM, msize, 1> res2;
	Matrix<PL_NUM, msize, 1> dx1;

	for( int i = 0; i < msize; ++i )
	{
		for( int j = 0; j < msize; ++j )
		{
			M( j, i ) = 0.0;
		}
		f11( i ) = 0.0;
		x1( i ) = 0.0;
		res( i ) = 0.0;
		res2( i ) = 0.0;
		dx1( i ) = 0.0;
	}
	
	//simply supported plate NO CURRENT PASSING THROUGH THE BOUNDARY

	int totLines = varNum / EQ_NUM;
	int _a = EQ_NUM / 2;
	for( int line = 0; line < totLines; ++line )
	{
		for( int vNum = 0; vNum < varNum / 2; ++vNum )
		{
			M( line * _a + 0, vNum ) = zi[Km - 1]( line * EQ_NUM + 0, vNum );		//TODO potential lags here!
			M( line * _a + 1, vNum ) = zi[Km - 1]( line * EQ_NUM + 1, vNum );
			M( line * _a + 2, vNum ) = zi[Km - 1]( line * EQ_NUM + 4, vNum );
			M( line * _a + 3, vNum ) = zi[Km - 1]( line * EQ_NUM + 6, vNum );
			M( line * _a + 4, vNum ) = zi[Km - 1]( line * EQ_NUM + 8, vNum );
		}
		f11( line * _a + 0 ) = -zi[Km - 1]( line * EQ_NUM + 0, varNum / 2 );
		f11( line * _a + 1 ) = -zi[Km - 1]( line * EQ_NUM + 1, varNum / 2 );
		f11( line * _a + 2 ) = -zi[Km - 1]( line * EQ_NUM + 4, varNum / 2 );
		f11( line * _a + 3 ) = -zi[Km - 1]( line * EQ_NUM + 6, varNum / 2 );
		f11( line * _a + 4 ) = -zi[Km - 1]( line * EQ_NUM + 8, varNum / 2 );
	}

	x1 = M.fullPivLu().solve( f11 );

	//refinement. I do not know the theoretical source of this procedure yet. just rewrote it
	//TODO test this
	res = f11 - M * x1;
	dx1 = M.fullPivLu().solve( res );
	x1 = x1 + dx1;

	PL_NUM nx = 0.0;
	PL_NUM ndx = 0.0;
	PL_NUM ndx2 = 0.0;
	PL_NUM temp;

	ndx = dx1.lpNorm<Infinity>();
	temp = ndx;							//FIXME may be we don't need temp
	int iterCount = 0;

	do
	{
		ndx = temp;
		++iterCount;

		res2 = res - M * dx1;
		res = res2;
		dx1 = M.fullPivLu().solve( res2 );
		x1 = x1 + dx1;


		res2 = res - M * dx1;
		res = res2;
		dx1 = M.fullPivLu().solve( res2 );
		x1 = x1 + dx1;

		nx = x1.lpNorm<Infinity>();
		ndx2 = dx1.lpNorm<Infinity>();
		temp = ndx2;
	} while( ndx2 < 0.9 * ndx && ndx2 / nx >= 2 * EPS_W );
	cout << " " << iterCount << " refnmt iterations\n";
	cout << " refnmt relative error is " << (PL_NUM)( ( M * x1 - f11 ).norm() / f11.norm() ) << endl;
	//refinement is over

	//now we determine coefficients for base solutions
	//the right-hand side:
	for( int i = 0; i < varNum / 2; ++i )
	{
		Cx1[i] = x1( i );
	}
	//calc the soln at the rhs
	for( int i = 0; i < varNum; ++i )
	{
		(*_mesh)[Km - 1].Nk1[i] = 0.0;
		for( int vNum = 0; vNum < varNum / 2; ++vNum )
		{
			(*_mesh)[Km - 1].Nk1[i] += Cx1[vNum] * zi[Km - 1]( i, vNum );			//FIXME lags may happen here
		}
		(*_mesh)[Km - 1].Nk1[i] += zi[Km - 1]( i, varNum / 2 );
	}

	//all the other points:
	for( int _x = Km - 2; _x >= 0; --_x )
	{
		if( orthoDone[_x] == true )
		{
			for( int i = varNum / 2 - 1; i >= 0; --i )
			{
				Cx[i] = Cx1[i] - omega[_x + 1][varNum / 2 * ( varNum / 2 + 1 ) / 2 + i];
				for( int j = varNum / 2 - 1; j > i; --j )
				{
					Cx[i] -= omega[_x + 1][j * ( j + 1 ) / 2 + i] * Cx[j];
				}
				Cx[i] /= omega[_x + 1][i * ( i + 1 ) / 2 + i];
			}
			for( int i = 0; i < varNum / 2; ++i )
			{
				Cx1[i] = Cx[i];
			}
		}
		else
		{
			for( int i = 0; i < varNum / 2; ++i )
			{
				Cx[i] = Cx1[i];
			}
			//... and there is no need to update the coefficients for the next interval
		}

		//now using the coefficients we write down the solution for current x
		for( int i = 0; i < varNum; ++i )
		{
			(*_mesh)[_x].Nk1[i] = 0.0;
			for( int vNum = 0; vNum < varNum / 2; ++vNum )
			{
				(*_mesh)[_x].Nk1[i] += Cx[vNum] * zi[_x]( i, vNum );			//FIXME lags may happen here
			}
			/*(*_mesh)[_x].Nk1[i] += solInfoMap[_x].z5[i];*/
			(*_mesh)[_x].Nk1[i] += zi[_x]( i, varNum / 2 );
		}
	}

	//force the BCs to be zero at y == a/2
	//TODO why do we need this??
	for( int line = 0; line < varNum / EQ_NUM; ++line )
	{
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 0] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 1] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 4] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 6] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 8] = 0.0;
	}
}