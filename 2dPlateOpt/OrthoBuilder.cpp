#include "OrthoBuilder.h"

template <>
N_PRES OrthoBuilder<HPD<N_PRES, 2> >::getInfNorm( HPD<N_PRES, 2>* vect, int vectSize )
{
	if( vect != 0 )
	{
		N_PRES ret = fabs( vect[0].real() );
		for( int i = 1; i < vectSize; ++i )
		{
			if( fabs( vect[i].real() ) > ret )
			{
				ret = fabs( vect[i].real() );
			}
		}
		return ret;
	}
	else
	{
		return -1;
	}
}