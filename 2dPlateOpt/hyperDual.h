#ifndef _PLATE_1D_HYPERDUAL_
#define _PLATE_1D_HYPERDUAL_ 1

#include <vector>
#include <iostream>

using std::vector;
using std::ostream;

template<class D_PRES, int NN>
class HPD;

//--------------------------------------------------------------------
// addition

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator+( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] + rhs.elems[i];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator+( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = rhs.elems[i];
	}
	ret.elems[0] += lhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator+( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i];
	}
	ret.elems[0] += rhs;
	return ret;
}

//----------------------------------------------------------------
// subtraction

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] - rhs.elems[i];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = -rhs.elems[i];
	}
	ret.elems[0] += lhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i];
	}
	ret.elems[0] -= rhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = -rhs.elems[i];
	}
	return ret;
}

//--------------------------------------------------------------------
// multiplication

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator*( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	ret.elems[0] = lhs.elems[0] * rhs.elems[0];
	for( int i = 1; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[0] * rhs.elems[i] + lhs.elems[i] * rhs.elems[0];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator*( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = rhs.elems[i] * lhs;
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator*( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] * rhs;
	}
	return ret;
}

//--------------------------------------------------------------------
// division

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator/( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	ret.elems[0] = lhs.elems[0] / rhs.elems[0];
	for( int i = 1; i < NN + 1; ++i )
	{
		ret.elems[i] = ( lhs.elems[i] * rhs.elems[0] - lhs.elems[0] * rhs.elems[i] ) / ( rhs.elems[0] * rhs.elems[0] ); 
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator/( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	ret.elems[0] = lhs / rhs.elems[0];
	for( int i = 1; i < NN + 1; ++i )
	{
		ret.elems[i] = ( -lhs * rhs.elems[i] ) / ( rhs.elems[0] * rhs.elems[0] ); 
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator/( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] / rhs;
	}
	return ret;
}

//----------------------------------------------------------------
// comparison

template<class D_PRES, int NN>
const bool operator<( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	return lhs.real() < rhs;
}

template<class D_PRES, int NN>
const bool operator<( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs < rhs.real();
}

template<class D_PRES, int NN>
const bool operator<( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs.real() < rhs.real();
}

template<class D_PRES, int NN>
const bool operator<=( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	return lhs.real() <= rhs;
}

template<class D_PRES, int NN>
const bool operator<=( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs <= rhs.real();
}

template<class D_PRES, int NN>
const bool operator<=( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs.real() <= rhs.real();
}

template<class D_PRES, int NN>
const bool operator>( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	return lhs.real() > rhs;
}

template<class D_PRES, int NN>
const bool operator>( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs > rhs.real();
}

template<class D_PRES, int NN>
const bool operator>( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs.real() > rhs.real();
}

template<class D_PRES, int NN>
const bool operator>=( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	return lhs.real() >= rhs;
}

template<class D_PRES, int NN>
const bool operator>=( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs >= rhs.real();
}

template<class D_PRES, int NN>
const bool operator>=( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	return lhs.real() >= rhs.real();
}

template<class D_PRES, int NN>
const bool operator==( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	if( lhs.real() != rhs )
	{
		return false;
	}
	return true;
}

template<class D_PRES, int NN>
const bool operator==( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	if( rhs.real() != lhs )
	{
		return false;
	}
	return true;
}

template<class D_PRES, int NN>
const bool operator==( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	for( int i = 0; i < NN + 1; ++i )
	{
		if( lhs.elems[i] != rhs.elems[i] )
		{
			return false;
		}
	}
	return true;
}

template<class D_PRES, int NN>
const bool operator!=( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	if( lhs.real() != rhs )
	{
		return true;
	}
	return false;
}

template<class D_PRES, int NN>
const bool operator!=( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )	//we are doing it this way probably on purpose (convergence checks, etc.)
{
	if( rhs.real() != lhs )
	{
		return true;
	}
	return false;
}

template<class D_PRES, int NN>
const bool operator!=( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	for( int i = 0; i < NN + 1; ++i )
	{
		if( lhs.elems[i] != rhs.elems[i] )
		{
			return true;
		}
	}
	return false;
}

//----------------------------------------------------------------
// cout

template<class D_PRES, int NN>
ostream& operator<<( ostream& os, const HPD<D_PRES, NN>& item )
{
	for( int i = 0; i < NN; ++i )
	{
		os << item.elems[i] << " ";
	}
	os << item.elems[NN];
	return os;
}

//----------------------------------------------------------------
// class itself

template<class D_PRES, int NN>
class HPD		//HPD = hyper dual
{
public:
	//vector<D_PRES> elems;
	D_PRES elems[NN + 1];

	HPD() 
	{	
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] = 0.0;
		}
	}
	//	elems.resize( NN + 1, 0.0 ); }
	HPD( D_PRES a )
	//{	elems.resize( NN + 1, 0.0 ); elems[0] = a; }
	{
		for( int i = 1; i < NN + 1; ++i )
		{
			elems[i] = 0.0;
		}
		elems[0] = a;
	}
	D_PRES real() const { return elems[0]; }

	template<class D_PRES_, int NN_>
	friend ostream& operator<<( ostream& os, const HPD& item );
	
	HPD& operator=( const HPD& rhs );
	HPD& operator=( const D_PRES& rhs );

	HPD& operator+=( const HPD& rhs );
	HPD& operator+=( const D_PRES& rhs );

	HPD& operator-=( const HPD& rhs );
	HPD& operator-=( const D_PRES& rhs );

	HPD& operator*=( const HPD& rhs );
	HPD& operator*=( const D_PRES& rhs );

	HPD& operator/=( const HPD& rhs );
	HPD& operator/=( const D_PRES& rhs );

	friend const HPD operator+<>( const HPD& lhs, const HPD& rhs );
	friend const HPD operator+<>( const D_PRES& lhs, const HPD& rhs );
	friend const HPD operator+<>( const HPD& lhs, const D_PRES& rhs );

	friend const HPD operator-<>( const HPD& lhs, const HPD& rhs );
	friend const HPD operator-<>( const D_PRES& lhs, const HPD& rhs );
	friend const HPD operator-<>( const HPD& lhs, const D_PRES& rhs );

	friend const HPD operator-<>( const HPD& rhs );

	friend const HPD operator*<>( const HPD& lhs, const HPD& rhs );
	friend const HPD operator*<>( const D_PRES& lhs, const HPD& rhs );
	friend const HPD operator*<>( const HPD& lhs, const D_PRES& rhs );

	friend const HPD operator/<>( const HPD& lhs, const HPD& rhs );
	friend const HPD operator/<>( const D_PRES& lhs, const HPD& rhs );
	friend const HPD operator/<>( const HPD& lhs, const D_PRES& rhs );

	friend const bool operator< <>( const HPD& lhs, const D_PRES& rhs );
	friend const bool operator< <>( const D_PRES& lhs, const HPD& rhs );
	friend const bool operator< <>( const HPD& lhs, const HPD& rhs );

	friend const bool operator<=<>( const HPD& lhs, const D_PRES& rhs );
	friend const bool operator<=<>( const D_PRES& lhs, const HPD& rhs );
	friend const bool operator<=<>( const HPD& lhs, const HPD& rhs );

	friend const bool operator> <>( const HPD& lhs, const D_PRES& rhs );
	friend const bool operator> <>( const D_PRES& lhs, const HPD& rhs );
	friend const bool operator> <>( const HPD& lhs, const HPD& rhs );

	friend const bool operator>=<>( const HPD& lhs, const D_PRES& rhs );
	friend const bool operator>=<>( const D_PRES& lhs, const HPD& rhs );
	friend const bool operator>=<>( const HPD& lhs, const HPD& rhs );

	friend const bool operator==<>( const HPD& lhs, const D_PRES& rhs );
	friend const bool operator==<>( const D_PRES& lhs, const HPD& rhs );
	friend const bool operator==<>( const HPD& lhs, const HPD& rhs );

	friend const bool operator!=<>( const HPD& lhs, const D_PRES& rhs );
	friend const bool operator!=<>( const D_PRES& lhs, const HPD& rhs );
	friend const bool operator!=<>( const HPD& lhs, const HPD& rhs );

//private:
	//HPD( const HPD& rhs );
};

//template<class D_PRES, int NN>
//HPD<D_PRES, NN>::HPD( const HPD& rhs )
//{
//	cout << " d-----\n";
//}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator=( const HPD& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] = rhs.elems[i];
		}
	}
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator=( const D_PRES& rhs )
{
	for( int i = 1; i < NN + 1; ++i )
	{
		elems[i] = 0.0;
	}
	elems[0] = rhs;
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator+=( const HPD& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] += rhs.elems[i];
		}
	}
	else
	{
		( *this ) *= 2.0;
	}
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator+=( const D_PRES& rhs )
{
	elems[0] += rhs;
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator-=( const HPD& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] -= rhs.elems[i];
		}
	}
	else
	{
		( *this ) = 0.0;
	}
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator-=( const D_PRES& rhs )
{
	elems[0] -= rhs;
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator*=( const HPD& rhs )
{
	HPD<D_PRES, NN> tmp;
	tmp = ( *this ) * rhs;
	//( *this ) = ( *this ) * rhs;
	( *this ) = tmp;

	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator*=( const D_PRES& rhs )
{
	for( int i = 0; i < NN + 1; ++i )
	{
		elems[i] *= rhs;
	}
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator/=( const HPD& rhs )
{
	HPD<D_PRES, NN> tmp;
	tmp = ( *this ) / rhs;
	//( *this ) = ( *this ) / rhs;
	( *this ) = tmp;

	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator/=( const D_PRES& rhs )
{
	for( int i = 0; i < NN + 1; ++i )
	{
		elems[i] /= rhs;
	}
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN> sin( const HPD<D_PRES, NN>& arg )
{
	HPD<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = sin( x0 ) + cos( x0 ) * ( arg - x0 );
	return ret;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN> cos( const HPD<D_PRES, NN>& arg )
{
	HPD<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = cos( x0 ) - sin( x0 ) * ( arg - x0 );
	return ret;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN> exp( const HPD<D_PRES, NN>& arg )
{
	HPD<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = exp( x0 ) + exp( x0 ) * ( arg - x0 );
	return ret;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN> sqrt( const HPD<D_PRES, NN>& arg )
{
	HPD<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = sqrt( x0 ) + 0.5 / sqrt( x0 ) * ( arg - x0 );
	return ret;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN> fabs( const HPD<D_PRES, NN>& arg )
{
	HPD<D_PRES, NN> ret = arg;
	if( arg.elems[0] < 0.0l )
	{
		ret = -ret;
	}
	return ret;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN> abs( const HPD<D_PRES, NN>& arg )
{
	return fabs( arg );
}

#endif