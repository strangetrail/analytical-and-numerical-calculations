#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <wchar.h>
#include <locale.h>
#include <vector>
//
#define N 255
#define M 255
#define L1_LENGTH 0.3e-6
#define L2_LENGTH L1_LENGTH
#define L_UNIT L1_LENGTH
//
#define K_SPLAY
#define K_TWIST
#define K_BEND
//
#define K_11 K_SPLAY
#define K_22 K_TWIST
#define K_33 K_BEND
//
#define SCALAR_ORDER_PARAMETER 1.0//????????????????????
#define CHIRALITY_LC_PARAMETER 1.0//????????????????????
#define ROTATIONAL_VISCOSITY 1.0//???????????????
#define eort 1.0//??????????
#define epar 1.0//??????????
//
enum axisIndex
{
	Z_IDX = 2,
	Y_IDX = 1,
	X_IDX = 0
}
enum edgeIndex
{
	M_IDX = 1,
	N_IDX = 0
};
enum sideMark
{
	SM_RIGHT = 6, SM_LEFT   = 5,
	SM_REAR  = 4, SM_FRONT  = 3,
	SM_TOP   = 2, SM_BOTTOM = 1,
	SM_INDIFFERENT = 0
};
typedef const char* EnumStrings_t [7];

static EnumStrings_t esSideMark =
                                  {
                                   "SM_INDIFFERENT",
                                   "SM_BOTTOM", "SM_TOP",
                                   "SM_FRONT",  "SM_REAR",
                                   "SM_LEFT",   "SM_RIGHT"
                                  };

enum sideRelation
{
	SR_PAR,
	SR_ORT_L1,
	SR_ORT_L2
};
int8_t eLeviCivita [3][3][3] =
{
 {
  { 0,  0,  0},
  { 0,  0,  1},
  { 0, -1,  0}
 },
 {
  { 0,  0, -1},
  { 0,  0,  0},
  { 1,  0,  0}
 },
 {
  { 0,  1,  0},
  {-1,  0,  0},
  { 0,  0,  0}
 }
}
int ***maxSideInd;//HERE
struct UnitSquareSimple
{
	double q;
	long long int mark;
	sideMark side;
};
//
struct UnitSquare
{
	double     z,
	          l1,// x
	          l2,// y
	           n,// l1
	           m,// l2
	          x0,
	          y0,
	          z0;
	uint8_t facet,
	         par;
//
	UnitSquare() :  z(0), l1(0), l2(0), n(0), m(0), x0(0), y0(0), z0(0),facet(0) {}// constructor
	UnitSquare                  \
	(                           \
	 double     cz,             \
	 double    cl1, double cl2, \
	 double     cn, double  cm, \
	 uint8_t cfacet             \
	)// Constructor.
	{
		z = cz;
		l1 = cl1; l2 = cl2;
		n = cn; m = cm;
		facet = cfacet;
	}
	uint8_t idx_facet()
	{
		return facet - 1;
	}
	UnitSquare                  \
	(                           \
	 double    cx0,             \
	 double    cy0,             \
	 double    cz0,             \
	 double     cz,             \
	 double    cl1, double cl2, \
	 double     cn, double  cm, \
	 uint8_t cfacet,             \
	 uint8_t cpar \
	)// constructor
	{
		z = cz;
		l1 = cl1; l2 = cl2;
		n = cn; m = cm;
		z0 = cz0; x0 = cx0; y0 = cy0;
		facet = cfacet;
		par = cpar;
	}
	double getLength(  int i )
	{
		if (i)
			return maxSideInd[par][facet][i] * l2;
		else
			return maxSideInd[par][facet][i] * l1;
	}
	double getL( int i )
	{
		if (i)
			return l2;
		else
			return l1;
	}
	double getNM(int i)
	{
		if (i)
			return m;
		else
			return n;
	}
	void setNM( int i, double val )
	{
		if (i)
			m = val;
		else
			n = val;
	}
};
enum idx_disambside
{
	DST = 1,
	SRC = 0
};
/*
struct DisambSides
{
	int srcSide,
	    dstSide;
};
*/
//union sideDisambiquation { int iSide; DisambSides uSide; };
//DisambSides equalEdges[6][6]/*indices are: source facet, destination facet*/ =
//{
// {{0,0}, {0,0}, {1,2}, {1,2}, {2,2}, {2,2}},/*SRC == BOTTOM*/
// {{0,0}, {0,0}, {1,2}, {1,2}, {2,2}, {2,2}},/*SRC == TOP   */
// {{2,1}, {2,1}, {0,0}, {0,0}, {1,1}, {1,1}},/*SRC == FRONT */
// {{2,1}, {2,1}, {0,0}, {0,0}, {1,1}, {1,1}},/*SRC == REAR  */
// {{2,2}, {2,2}, {1,1}, {1,1}, {0,0}, {0,0}},/*SRC == LEFT  */
// {{2,2}, {2,2}, {1,1}, {1,1}, {0,0}, {0,0}}/* SRC == RIGHT */
//};
//
//CUDA????? nope.., use braaaains!!!
//
inline uint32_t getFacetEdgeUnitCount( uint8_t i, uint8_t j, uint8_t k )
{
	return maxSideInd[i][j][k] + 1;
}
uint8_t mapFacet2AxisIdx( uint8_t idx_facet )
{
	return ( idx_facet + 1 ) / 3;
}
int8_t getAddOrSub( uint8_t idx_facet )
{
	return 1 - 2 * ( idx_facet % 2 );
}
//
	uint8_t getIsOrt( UnitSquare usSrc, UnitSquare usDst )
	{
		return (                                                               \
		     (( usSrc.facet == SM_BOTTOM ) || ( usSrc.facet == SM_TOP   )) && \
		     (( usDst.facet == SM_REAR   ) || ( usDst.facet == SM_FRONT ))    \
		    ) ||                                                            \
		    (                                                               \
		     (( usDst.facet == SM_BOTTOM ) || ( usDst.facet == SM_TOP   )) && \
		     (( usSrc.facet == SM_RIGHT  ) || ( usSrc.facet == SM_LEFT  ))    \
		    ) ||                                                               \
		    (  \
		     (( usDst.facet == SM_FRONT ) || ( usDst.facet == SM_REAR   )) && \
		     (( usSrc.facet == SM_RIGHT  ) || ( usSrc.facet == SM_LEFT  ))    \
		    );
	}
	double *getUnitSquareCoords( UnitSquare usOrigin )
	{
	double usX_src, usY_src, usZ_src, *usResults;
	switch ( usOrigin.facet )
	{
		case SM_BOTTOM:
		case SM_TOP:
				usX_src = usOrigin.x0  - L_UNIT * ( getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( usOrigin.facet ), N_IDX ) + usOrigin.n );
				usY_src =  usOrigin.y0 - L_UNIT * ( getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( usOrigin.facet ), M_IDX ) + usOrigin.m );
				usZ_src = usOrigin.z0 + getAddOrSub( usOrigin.facet ) * L_UNIT * getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( SM_FRONT ), M_IDX );
			break;
		case SM_LEFT:
		case SM_RIGHT:
				usY_src = usOrigin.y0 - L_UNIT * ( getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( usOrigin.facet ), M_IDX ) + usOrigin.m );
				usZ_src = usOrigin.z0 - L_UNIT * ( getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( usOrigin.facet ), N_IDX ) + usOrigin.n );
				usX_src = usOrigin.x0 + getAddOrSub( usOrigin.facet ) * L_UNIT * getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( SM_FRONT ), N_IDX );
			break;
		case SM_FRONT:
		case SM_REAR:
				usX_src = usOrigin.x0 - L_UNIT * ( getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( usOrigin.facet ), N_IDX ) + usOrigin.n );
				usZ_src = usOrigin.z0 - L_UNIT * ( getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( usOrigin.facet ), M_IDX ) + L_UNIT * usOrigin.m );
				usY_src = usOrigin.y0 + getAddOrSub( usOrigin.facet ) * L_UNIT * getFacetEdgeUnitCount( usOrigin.par, mapFacet2AxisIdx( SM_BOTTOM ), M_IDX );
			break;
	}
	usResults = new double[3];
	usResults[0] = usX_src;
	usResults[1] = usY_src;
	usResults[2] = usZ_src;
	return usResults ;
	}
//
UnitSquare getTargetUnitSquare ( UnitSquare usSrc, UnitSquare usDst )
{
	UnitSquare getResult;
	double *srcUnitPivot = getUnitSquareCoords( usSrc );
	double *dstUnitPivot = getUnitSquareCoords( usDst );
	switch ( usSrc.facet )
	{
		case SM_BOTTOM:
		case SM_TOP:
				getResult.z = fabs( srcUnitPivot[2]- dstUnitPivot[2]);
				getResult.z += getIsOrt( usSrc, usDst ) ? 0.5 * L_UNIT : 0.0;
				getResult.n = fabs( srcUnitPivot[0]- dstUnitPivot[0]) / L_UNIT;
				getResult.n += (( usDst.facet == SM_LEFT ) || ( usDst.facet == SM_RIGHT ) )? 0.5 * L_UNIT : 0.0;
				getResult.m = fabs( srcUnitPivot[1]- dstUnitPivot[1]) / L_UNIT;
				getResult.m += (( usDst.facet == SM_FRONT ) || ( usDst.facet == SM_REAR ) )? 0.5 * L_UNIT : 0.0;
			break;
		case SM_FRONT:
		case SM_REAR:
				getResult.z = fabs( srcUnitPivot[1]- dstUnitPivot[1]);
				getResult.z += getIsOrt( usSrc, usDst ) ? 0.5 * L_UNIT : 0.0;
				getResult.n = fabs( srcUnitPivot[0]- dstUnitPivot[0]) / L_UNIT;
				getResult.n += (( usDst.facet == SM_LEFT ) || ( usDst.facet == SM_RIGHT ) )? 0.5 * L_UNIT : 0.0;
				getResult.m = fabs( srcUnitPivot[2]- dstUnitPivot[2]) / L_UNIT;
				getResult.m += (( usDst.facet == SM_BOTTOM ) || ( usDst.facet == SM_TOP ) )? 0.5 * L_UNIT : 0.0;
			break;
		case SM_LEFT:
		case SM_RIGHT:
				getResult.z = fabs( srcUnitPivot[0]- dstUnitPivot[0]);
				getResult.z += getIsOrt( usSrc, usDst ) ? 0.5 * L_UNIT : 0.0;
				getResult.n = fabs( srcUnitPivot[2]- dstUnitPivot[2]) / L_UNIT;
				getResult.n += (( usDst.facet == SM_BOTTOM ) || ( usDst.facet == SM_TOP ) )? 0.5 * L_UNIT : 0.0;
				getResult.m = fabs( srcUnitPivot[1]- dstUnitPivot[1]) / L_UNIT;
				getResult.m += (( usDst.facet == SM_FRONT ) || ( usDst.facet == SM_REAR ) )? 0.5 * L_UNIT : 0.0;
			break;
	}
	getResult.l1 = L_UNIT;
	getResult.l2 = L_UNIT;
	return getResult;
}
//
//
double R0(double arg1, double arg2, double z)
{
	return sqrt( pow(arg1,2) + pow(arg2,2) + 4*pow(z,2) );
}
double Rarg( double arg0, double arg1, double arg2, double z )
{
	return arg0 + R0( arg1, arg2, z );
}
double asinh_arg(double A, double B, double C, double U)
{
	return asinh( ( 2*A*U + B ) / sqrt( 4*A*C - B*B ) );
}
//
double asin__arg(double A, double B, double C, double U)
{
	return asin( ( B*U + 2*C ) / ( fabs(U) * sqrt( B*B - 4*A*C ) ) );
}
//
#define RA(U) R0(Aa,Ba,Ca,U)
#define RB(U) R0(Ab,Bb,Cb,U)
//
double __getAsinh_LogExpansion_Complex( double lnX, double lmX, double z )
{
	return log(
	           ( lnX + sqrt( pow(lnX,2.0) + pow(lmX,2.0) + 4.0*pow(z,2.0) ) ) /
	           sqrt( pow(lmX,2.0) + 4.0*pow(z,2.0) )
	          );
}
//
double __getAsin( double lnX, double lmX, double z )
{
	return asin(
	            ( lmX * ( lmX + sqrt( pow(lnX,2.0) + pow(lmX,2.0) + 4.0*pow(z,2.0) ) ) + 4.0*pow(z,2.0) ) /
	            ( lmX + sqrt( pow(lnX,2.0) + pow(lmX,2.0) + 4*pow(z,2.0) ) ) *
	            sqrt( pow(lmX,2.0) + 4*pow(z,2.0) )
	           );
}
//
#define MAX_SERIES_ORDER 255
//
//
typedef std::vector<seriesMultipliedTerm> seriesMultipliedTermSum;
//
struct seriesMultipliedTerm
{
	private:
	uint16_t idxI;
	public:
	double                  a,
	                        b,
	                        z;
	uint64_t            pow_a,
	                    pow_b,
	                    pow_z,
	                 pow_sign;
	int64_t  pows_intProducts [MAX_SERIES_ORDER - 1];
	seriesMultipliedTerm() :
	                                a(0),     b(0),     z(0),
	                            pow_a(0), pow_b(0), pow_z(0),
	                         pow_sign(0)
	{
		memset( (void *)this->pows_intProducts, 0, (MAX_SERIES_ORDER + 1) * (sizeof ( uint64_t )) );
	}
	seriesMultipliedTerm( uint64_t pow_sign ) : seriesMultipliedTerm()
	{
		this->pow_sign = pow_sign;
	}
	seriesMultipliedTerm \
	(                    \
	 double a,           \
	 double b,           \
	 double z,           \
	 uint64_t pow_a,     \
	 uint64_t pow_b,     \
	 uint64_t pow_z      \
	) :
	    seriesMultipliedTerm()// Constructor
	{
		this->a = a;
		this->b = b;
		this->z = z;
		this->pow_a = pow_a;
		this->pow_b = pow_b;
		this->pow_z = pow_z;
	}
	seriesMultipliedTerm \
	(                    \
	 uint64_t num_bgn,   \
	 uint64_t num_end    \
	) :
	    seriesMultipliedTerm()// Constructor
	{
		memset                                                                       \
		(                                                                            \
		 ((void *)this->pows_intProducts) + (sizeof ( uint64_t )) * ( num_bgn - 2 ), \
		 1,                                                                          \
		 (num_end - num_bgn + 1) * (sizeof ( uint64_t ))                             \
		);
	}
	seriesMultipliedTerm                       \
	(                                          \
	 int16_t (*getNextNumIndex)( uint8_t idx, uint8_t maxIdx ), \
	 uint8_t max_local;
	) :
	    seriesMultipliedTerm()// Constructor
	{
		uint8_t idx_local = 0,
		         idx_term;
//
		while ( ( idx_term = getNextNumIndex( idx_local, max_local ) ) > -1 )
		{
			this->pows_intProducts[idx_term] = 1;
		}
	}
////
	friend seriesMultipliedTermSum operator +                                    \
	                                        (                                    \
	                                               seriesMultipliedTermSum  lhs, \
	                                         const seriesMultipliedTerm    &rhs  \
	                                        )
	{
		return lhs.push_back( rhs );
	}
////
	friend seriesMultipliedTermSum operator +                                    \
	                                        (                                    \
	                                         const seriesMultipliedTerm    &lhs, \
	                                               seriesMultipliedTermSum  rhs  \
	                                        )
	{
		return rhs + lhs;
	}
////
	friend seriesMultipliedTermSum operator *                                    \
	                                        (                                    \
	                                               seriesMultipliedTermSum  lhs, \
	                                         const seriesMultipliedTerm    &rhs  \
	                                        )
	{
		for ( idxI = 0; idxI < lhs.size(); idxI++ )
			lhs[idxI] = lhs[idxI] * rhs;
		return lhs;
	}
////
	friend seriesMultipliedTermSum operator *                                    \
	                                        (                                    \
	                                         const seriesMultipliedTerm    &rhs, \
	                                               seriesMultipliedTermSum  lhs  \
	                                        )
	{
		return rhs + lhs;
	}
}
//
int16_t getSqRtBinomCoeffNumerProdIdx( uint8_t idx, uint8_t max_L )
{
	uint8_t idx_L = idx + 2;
//
	if ( idx_L > max_L )
		return -1;
	else
		return 2 * idx_L - 3;
}
//
int16_t getSqRtBinomCoeffDenomProdIdx( uint8_t idx, uint8_t max_L )
{
	uint8_t idx_L = idx + 1;
//
	if ( idx_L > max_L )
		return -1;
	else
		return 2 * idx_L;
}
//
double __getIntLogSeries \
(                        \
 double  arg,            \
 uint8_t idx_last,       \
 double  value_a,        \
 double  value_b,        \
 double  value_z         \
)
{
	uint8_t idx_N,
	        idx_K,     idx_J, idx_Q, idx_R, idx_P, idx_H,
	        idx_M,     idx_S, idx_T, idx_U,
	        idx_F, idx_ALPHA,
	        idx_V,     idx_W;
////
	seriesMultipliedTerm minusOne_N;
////
	seriesMultipliedTermSum smts__Term_F;
////
	for ( idx_N = 1; idx_N < idx_last + 1; idx_N++ )
	{
		minusOne_N = seriesMultipliedTerm( 2 * idx_N );
		smt__Term_N = minusOne_N;
		for ( idx_K = 0; idx_K < idx_N; idx_K++ )
		{
//		//HERE 2 !
//		//
//		// Term with constant upper index of summation:
//		//
			minusOneStandalone = seriesMultipliedTerm( 1 );
//		//
			iSingleTermArg = 2 * (idx_N - idx_K) - 1;
			smt__SingleTerm_2tNm2tKm1 =                                 \
			                            seriesMultipliedTerm            \
			                            (                               \
			                             iSingleTermArg, iSingleTermArg \
			                            );
//		//
			smt__TermStandalone_AZ.z = value_z;
			smt__TermStandalone_AZ.pow_z = 1;
			smt__TermStandalone_AZ.a = value_a;
			smt__TermStandalone_AZ.pow_a = 1;
//		//
			smt__TermStandalone_B2.b = value_b;
			smt__TermStandalone_B2.pow_b = 2;
//		//
			smt__TermStandalone_Z2.z = value_z;
			smt__TermStandalone_Z2.pow_z = 2;
//		//
			smt__TermStandaloneMult = minusOneStandalone / smt__SingleTerm_2tNm2tKm1;
//		//
			smts__Term_NK =                           \
			                (                         \
			                 smts__Term_NK +          \
			                 smt__TermStandalone_AZ + \
			                 smt__TermStandalone_B2 + \
			                 smt__TermStandalone_Z2   \
			                ) *                       \
			                smt__TermStandaloneMult;
//		//
			smt__Prod_NK = seriesMultipliedTerm( 2 * (idx_N - idx_K) - 1, 2*idx_N - 1 );
			smt__Prod_K = seriesMultipliedTerm( idx_K + 1, 2 * idx_K );
//		//
			iSingleTermArg = 2 * idx_K + 1;
			smt__SingleTerm_2tKp1 = seriesMultipliedTerm( iSingleTermArg, iSingleTermArg );
//		//
			for ( idx_J = 0; idx_J < 2 * ( idx_N - idx_K ) - 1; idx_J++ )
			{
				for ( idx_Q = 0; idx_Q < 2 * ( idx_N - idx_K ) - idx_J - 1; idx_Q++ )
				{
					minusOne_Q = seriesMultipliedTerm( idx_Q );
					smt__Prod_NKJQ =                                           \
					                 seriesMultipliedTerm                      \
					                 (                                         \
					                  2 * (idx_N - idx_K) - idx_J - idx_Q - 1, \
					                  2 * (idx_N - idx_K - 1)                  \
					                 );
					smt__Prod_Q = seriesMultipliedTerm( 2, idx_Q );
					for ( idx_R = 0; idx_R < idx_J + 1; idx_R++ )
					{
						smt__Prod_R = seriesMultipliedTerm( 2, idx_R );
						smt__Prod_JR = seriesMultipliedTerm( 2, idx_J - idx_R );
						for ( idx_P = 0; idx_P < idx_K + 1; idx_P++ )
						{
							z_PJR.z = value_z;
							z_PJR.pow_z = 2 * idx_P + idx_J + idx_R;
							smt__Prod_P = seriesMultipliedTerm( 2, idx_P );
							for ( idx_H = 0; idx_H < idx_K - idx_P + 1; idx_H++ )
							{
								a_KJRPH.a = value_a;
								a_KJRPH.pow_a = 2 * ( idx_K - idx_P - idx_H ) + idx_J + idx_R;
								b_NKJQRPH.b = value_b;
								b_NKJQRPH.pow_b = 4 * ( idx_N - 1 ) + 2 * ( idx_H - idx_Q - idx_J - idx_K );
//							//HERE 6 !!! // Fix this error
								smt__Prod_H = seriesMultipliedTerm( 2, idx_H );
								smt__Prod_KPH = seriesMultipliedTerm( 2, idx_K - idx_P - idx_H );
								smt__Term_KJQRPH = minusOne_Q *                                              \
								                   z_PJR * a_NJQPHR * b_H *                                  \
								                   smt__Prod_NK * smt__Prod_NKJQ /                           \
								                   (                                                         \
								                    smt__Prod_Q * smt__Prod_R * smt__Prod_JR * smt__Prod_K * \
								                    smt__Prod_P * smt__Prod_H * smt__Prod_KPH                \
								                   );
								for ( idx_M = 0; idx_M < idx_last + 1; idx_M++ )
								{
									smt__Prod_Mnumer = seriesMultipliedTerm( getSqRtBinomCoeffNumerProdIdx, idx_M );
									smt__Prod_Mdenom = seriesMultipliedTerm( getSqRtBinomCoeffDenomProdIdx, idx_M );
									for ( idx_S = 0; idx_S < idx_M + 1; idx_S++ )
									{
										z_S.z = value_z;
										z_S.pow_z = 2 * idx_S;
										smt__Prod_S = seriesMultipliedTerm( 2, idx_S );
										for ( idx_T = 0; idx_T < idx_M - idx_S + 1; idx_T++ )
										{
											a_MST.a = value_a;
											a_MST.pow_a = 2 * ( idx_M - idx_S - idx_T ) + 1;
											smt__Prod_MST = seriesMultipliedTerm( idx_M - idx_S - idx_T, idx_M );
											for ( idx_U = 0; idx_U < 2 * idx_T + 1; idx_U++ )
											{
												minusOne_TU = seriesMultipliedTerm( 2 * idx_T - idx_U );
												b_U.b = value_b;
												b_U.pow_b = value_b;
												smt__Prod_TU = seriesMultipliedTerm( 2 * idx_T - idx_U + 1, 2 * idx_T );
												smt__Prod_U = seriesMultipliedTerm( 2, idx_U );
												smt__Term_MSTU = minusOne_TU *                                     \
												                 z_S * a_MST * b_U *                               \
												                 smt__Prod_TU * smt__Prod_MST * smt__Prod_Mnumer / \
												                 (                                                 \
												                  smt__Prod_S * smt__Prod_U * smt__Prod_Mdenom *   \
												                  smt__SingleTerm_2tKp1                            \
												                 );
												smts__Summand_NKJQRPHMSTUF = smt__Term_N * smt__Term_KJQRPH * \
												                             smt__Term_MSTU * smts__Term_F;
												sstVector__IntLogSeries = sstVector__IntLogSeries + smts__Summand_NKJQRPHMSTUF;
//											//HERE 3
												iSingleTermArg = 4 * idx_N - 2 * idx_K - 2;
												smt__SingleTerm_4tNm2tKm2 = seriesMultipliedTerm            \
												                            (                               \
												                             iSingleTermArg, iSingleTermArg \
												                            );
												iSingleTermArg = 2 * idx_N - 1;
												smt__SingleTerm_2tNm1 = seriesMultipliedTerm            \
												                        (                               \
												                         iSingleTermArg, iSingleTermArg \
												                        );
												smt__TermFree_0 = smt__SingleTerm_4tNm2tKm2 /                        \
												                  (                                                  \
												                   smt__SingleTerm_2tNm1 * smt__SingleTerm_2tNm2tKm1 \
												                  );
												smt__Summand_NKJQRPHMSTU = smt__Term_N * smt__Term_KJQRPH * \
												                            smt__Term_MSTU * smt_TermFree_0;
												ssts__IntLogSeries = ssts__IntLogSeries + smt__Summand_NKJQRPHMSTU;
//											// HERE 4 !!!
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return log_expansion;
}
//
double getDistancePotentialConstantAnisotropic( UnitCube ucArea )
{
	double Int_f1_dz, Int_f2_dz, Int_f3_dz, Int_f4_dz, Int_f5_dz;
//
	Int_f3_dz = (
	             2 * pow( z, 2 ) * atanh( (x * y) /
	             ( z * sqrt( pow( x, 2 ) + pow( y, 2 ) + pow( z, 2 ) ) ) ) -
	            
	            ) / 4
	Int_f4_dz = ( ( pow( y, 2 ) + pow( z, 2 ) ) * atan( y / z ) + y * z ) / 2
}
//
double getDistancePotentialConstant( UnitCube ucArea )
{
	double lmA, lmB, lnA, lnB,
	        Ba,  Bb,
	         V;
//
	lmA = usArea.l2 * (2.0*usArea.m - 1.0);
	lmB = usArea.l2 * (2.0*usArea.m + 1.0);
	lnA = usArea.l1 * (2.0*usArea.n - 1.0);
	lnB = usArea.l1 * (2.0*usArea.n + 1.0);
//
	Ba = (-1.0) * lmA;//Ba = (-2.0) * lmA;
	Bb = (-1.0) * lmB;//Bb = (-2.0) * lmB;
//
	V = (
	     (
	      lnA * log(
	                ( lmA + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnA,2.0) + pow(lmA,2.0) ) ) /
	                ( lmB + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnA,2.0) + pow(lmB,2.0) ) )
	               ) +
	      lnB * log(
	                ( lmB + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnB,2.0) + pow(lmB,2.0) ) ) /
	                ( lmA + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnB,2.0) + pow(lmA,2.0) ) )
	               )
	     ) / 2.0 +
	     /*asinh_log_expansion_complex_numbers*/
	     (
	      Ba * (
	            __getAsinh_LogExpansion_Complex( lnB, lmA, usArea.z ) -
	            __getAsinh_LogExpansion_Complex( lnA, lmA, usArea.z )
	           ) +
	      Bb * (
	            __getAsinh_LogExpansion_Complex( lnA, lmB, usArea.z ) -
	            __getAsinh_LogExpansion_Complex( lnB, lmB, usArea.z )
	           )
	     ) / 2.0 +
	     /*asin*/
	     (
	          __getAsin( lnA, lmB, usArea.z ) -
	          __getAsin( lnB, lmB, usArea.z ) +
	          __getAsin( lnB, lmA, usArea.z ) -
	          __getAsin( lnA, lmA, usArea.z )
	         )
	    );
//
	return V;
}
//
double getDistancePotentialConstant( UnitSquare usArea )
{
	double lmA, lmB, lnA, lnB,
	        Ba,  Bb,
	         V;
//
	lmA = usArea.l2 * (2.0*usArea.m - 1.0);
	lmB = usArea.l2 * (2.0*usArea.m + 1.0);
	lnA = usArea.l1 * (2.0*usArea.n - 1.0);
	lnB = usArea.l1 * (2.0*usArea.n + 1.0);
//
	Ba = (-2.0) * lmA;
	Bb = (-2.0) * lmB;
//
	V = (
	     (
	      lnA * log(
	                ( lmA + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnA,2.0) + pow(lmA,2.0) ) ) /
	                ( lmB + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnA,2.0) + pow(lmB,2.0) ) )
	               ) +
	      lnB * log(
	                ( lmB + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnB,2.0) + pow(lmB,2.0) ) ) /
	                ( lmA + sqrt( 4.0*pow(usArea.z,2.0) + pow(lnB,2.0) + pow(lmA,2.0) ) )
	               )
	     ) / 2.0 +
	     /*asinh_log_expansion_complex_numbers*/
	     (
	      Ba * (
	            __getAsinh_LogExpansion_Complex( lnB, lmA, usArea.z ) -
	            __getAsinh_LogExpansion_Complex( lnA, lmA, usArea.z )
	           ) +
	      Bb * (
	            __getAsinh_LogExpansion_Complex( lnA, lmB, usArea.z ) -
	            __getAsinh_LogExpansion_Complex( lnB, lmB, usArea.z )
	           )
	     ) / 2.0 +
	     /*asin*/
	     usArea.z * (
	          __getAsin( lnA, lmB, usArea.z ) -
	          __getAsin( lnB, lmB, usArea.z ) +
	          __getAsin( lnB, lmA, usArea.z ) -
	          __getAsin( lnA, lmA, usArea.z )
	         )
	    );
//
	return V;
}
//
UnitSquareSimple **getMarkedGrid \
(                                \
  int n,                         \
  int m,                         \
  int &k,                        \
  UnitSquare *(&distances),         \
  sideMark currentSide           \
)// HERE
{
	uint8_t i, j, l;
	UnitSquareSimple **_mgRes;
	int distCount;
	double distance, nHalf, mHalf, *dDistances;
	bool newDist;
//
	nHalf = n/2;
	mHalf = m/2;
	distCount = nHalf*mHalf;
	dDistances = new double[distCount];
	distances = new UnitSquare[distCount];
	memset( (void *)dDistances, 0, distCount * (sizeof ( double )) );
	_mgRes = new UnitSquareSimple* [n];
	for ( i = 0; i < n; i++ )
		_mgRes[i] = new UnitSquareSimple[m];
//
	k = 0;
	for ( i = 0; i < nHalf; i++ )
		for ( j = 0; j < mHalf; j++ )
		{
			distance = pow(i, 2) + pow(j, 2);
			for ( l = 0, newDist = true; ( l < k ) && newDist; l++ )
				newDist = distance != dDistances[l];
			if ( newDist )
			{
				_mgRes[i][j].mark = l;
				dDistances[l] = distance;
				distances[l] = UnitSquare(0, 0, 0, i, j, currentSide);// CHECK syntax.// HERE
				k++;
			}
			else
			{
				_mgRes[i][j].mark = l-1;
				distances[l-1] = UnitSquare( 0, 0, 0, i, j, currentSide );// CHECK syntax.// HERE
			}
			_mgRes[n-i-1][j].mark = _mgRes[i][j].mark;
			_mgRes[i][m-j-1].mark = _mgRes[i][j].mark;
			_mgRes[n-i-1][m-j-1].mark = _mgRes[i][j].mark;
		}
//
	delete [] dDistances;
//
	return _mgRes;
}
//
//Check stack push/pull methods online.
double **getEquationMatrix                                                            \
(                                                                                     \
  sideMark *datMarks                                                                  \
  /*Sequence of marks respects sequense of matrices ((UnitSquareSimple**)(data*)).*/, \
  UnitSquareSimple ***data,                                                           \
  int dataItems,                                                                      \
  /*Number of sides*/                                                                 \
  int count,                                                                          \
  /*Number of marks*/                                                                 \
  uint8_t n,                                                                          \
  uint8_t m,                                                                          \
  UnitSquare **distances                                                              \
)//Required rigorous verification.
{
	uint8_t i, j, k, l;
	double **matrRresult;
	UnitSquare usSrc, usDst;
//
	matrRresult = new double*[count];
	for ( i = 0; i < count; i++ )
	{
		// One additional column for RHS/LHS in Gaussian elimination:
		matrRresult[i] = new double[count+1];
		memset( matrRresult[i], 0, count+1 );
	}
//
	for ( i = 0; i < dataItems; i++ )
	{
		for ( l = 0; l < count; l++ )
		{
			usDst.facet = (distances[i][l]).facet;
			usDst.l1 = L1_LENGTH;
			usDst.l2 = L2_LENGTH;
			usDst.n = (distances[i][l]).n;
			usDst.m = (distances[i][l]).m;
			for ( j = 0; j < n; j++ )
				for ( k = 0; k < m; k++ )
				{
					usSrc.facet = datMarks[i];
					usSrc.l1 = L1_LENGTH;
					usSrc.l2 = L2_LENGTH;
					usSrc.n = j;
					usSrc.m = k;
					// L indexing Destination; data[i][j][k].mark indexing Source.
					matrRresult[l][data[i][j][k].mark] += getDistancePotentialConstant( getTargetUnitSquare( usSrc, usDst ) );// HERE
				}
		}
	}
	return matrRresult;
}
//
void swapRows( double **A, short k, short i, short m )
{
	double tmp;
	for ( short j = 0; j < m + 1; j++ )
	{
		tmp = A[k][j];
		A[k][j] = A[i][j];
		A[i][j] = tmp;
	}
}
//
void getGaussianSolution( double **A, short m )
{
//	short ind_max;
	double /*Ai_max = 0.0,*/
	       tmp;
//
	for ( short k = 0; k < m; k++ )
	{
		// Find pivot for column k:
/*
		for ( short i = k; i < m; i++ )
		{
			if ( ( tmp = fabs( A[i][k] ) ) > Ai_max )
			{
				Ai_max = tmp;
				ind_max = i;
			}
		}
		if ( A[ind_max][k] == 0 )
		{
			printf( "Error: Matrix is singular!" );
		}
		swapRows( A, k, ind_max, m );
*/
		// Do for all rows below pivot:
		for ( short i = k + 1; i < m; i++ )
		{
			// Do for all remaining elements in current row:
				for ( short j = k+1; j < m + 1; j++ )
				{
					A[i][j] = A[i][j] - A[k][j] * ( A[i][k] / A[k][k] );
				}
				// Fill lower triangular matrix with zeros:
			 A[i][k] = 0;
		}
	}
	for ( short k = m - 1; k > -1; k-- )
	{
		tmp = 0;
		for ( short j = k + 1; j < m; j++ )
		{
			A[k][j] = A[k][j] * A[j][j];
			tmp += A[k][j];
			A[k][j] = 0;
		}
		A[k][k] = (A[k][m] - tmp) / A[k][k];
		A[k][m] = k+1;
	}
}
//
struct LCMolDirector
{
	uint64_t xidx, yidx, zidx;
	void x(cx) { values_[0] = cx; }
	void y(cy) { values_[1] = cy; }
	void z(cz) { values_[2] = cz; }
	float& x() { return values_[0]; }
	float& y() { return values_[1]; }
	float& z() { return values_[2]; }
	LCMolDirector( double cx, double cy, double cz ): z(cz), y(cy), x(cx) {}
	LCMolDirector(): z(0), y(0), x(0) {}
	float  operator [] (unsigned i) const { return this->values_[i]; }
	float& operator [] (unsigned i)       { return this->values_[i]; }
	operator float*() const { return this->values_; }
	private:
	float[3] values_;
	static LCMolDirector ***DirectorLCField;
	static double ***PotentialLCField;
	double diff_n( LCMolDirector ***n_vicinity, uint8_t idx1, uint8_t idx2 );
}
//
getIsCconverges( int t )
{
	bool bConverges = false;
//
//
	return bConverges;
}
//
double f_g( double n )
{
	double __result;
	return __result;
}
//
typedef double *( *fg_t )();
//
double *df_g__dn_x()
{
	return 3 * this->x() * SCALAR_ORDER_PARAMETER * df_g__dQ( X_IDX, X_IDX ) + \
	       3 / 2 * this->y() * SCALAR_ORDER_PARAMETER *                        \
	       ( df_g__dQ( X_IDX, Y_IDX ) + df_g__dQ( Y_IDX, X_IDX ) ) +           \
	       3 / 2 * this->z() * SCALAR_ORDER_PARAMETER *                        \
	       ( df_g__dQ( X_IDX, Z_IDX ) + df_g__dQ( Z_IDX, X_IDX ) );
}
//
double *df_g__dn_y()
{
	return 3 * this->y() * SCALAR_ORDER_PARAMETER * df_g__dQ( Y_IDX, Y_IDX ) + \
	       3 / 2 * this->x() * SCALAR_ORDER_PARAMETER *                        \
	       ( df_g__dQ( X_IDX, Y_IDX ) + df_g__dQ( Y_IDX, X_IDX ) ) +           \
	       3 / 2 * this->z() * SCALAR_ORDER_PARAMETER *                        \
	       ( df_g__dQ( Y_IDX, Z_IDX ) + df_g__dQ( Z_IDX, Y_IDX ) );
}
//
double *df_g__dn_z()
{
	return 3 * this->z() * SCALAR_ORDER_PARAMETER * df_g__dQ( Z_IDX, Z_IDX ) + \
	       3 / 2 * this->x() * SCALAR_ORDER_PARAMETER *                        \
	       ( df_g__dQ( X_IDX, Z_IDX ) + df_g__dQ( Z_IDX, X_IDX ) ) +           \
	       3 / 2 * this->y() * SCALAR_ORDER_PARAMETER *                        \
	       ( df_g__dQ( Y_IDX, Z_IDX ) + df_g__dQ( Z_IDX, Y_IDX ) );
}
//
// TODO -> BELOW  put everything in order !!!! complete D(f_g)/D(n[i]) -currently  you have only first summand !!! 
//
#define DX 0.0001
#define DY 0.0001
#define DZ 0.0001
//
double Dfreevar [3] = {DX, DY, DZ};
//
double diff_n( LCMolDirector ***n_vicinity, uint8_t aidx, uint8_t idx1 )
{
	return diff_n( n_vicinity, aidx, idx1, 0 );// This is a crutch.
}
double diff_n( LCMolDirector ***n_vicinity, uint8_t aidx, uint8_t idx1, uint8_t idx2 )
{
	if ( idx2 == 0 )
	{
		if ( idx1 == idx2 )
		{
			nx = !(1^idx1) + !(1^idx2);
			ny = !(2^idx1) + !(2^idx2);
			nz = !(3^idx1) + !(3^idx2);
			return                                                                 \
			(                                                                      \
			 n_vicinity[this->xidx + nx][this->yidx + ny][this->zidx + nz][aidx] + \
			 n_vicinity[this->xidx - nx][this->yidx - ny][this->zidx - nz][aidx] - \
			 2 * n_vicinity[this->xidx][this->yidx][this->zidx][aidx]              \
			) / ( 2 * Dfreevar[idx1 - 1] );
		}
		else
		{
			nx = !(1^idx1) + !(1^idx2);
			ny = !(2^idx1) + !(2^idx2);
			nz = !(3^idx1) + !(3^idx2);
			Dn = n_vicinity[this->xidx + nx][this->yidx + ny][this->zidx + nz][aidx] + \
			     n_vicinity[this->xidx - nx][this->yidx - ny][this->zidx - nz][aidx];
			nx = !(1^idx1) - !(1^idx2);
			ny = !(2^idx1) - !(2^idx2);
			nz = !(3^idx1) - !(3^idx2);
			Dn -= n_vicinity[this->xidx + nx][this->yidx + ny][this->zidx + nz][aidx];
			nx = !(1^idx2) - !(1^idx1);
			ny = !(2^idx2) - !(2^idx1);
			nz = !(3^idx2) - !(3^idx1);
			Dn -= n_vicinity[this->xidx + nx][this->yidx + ny][this->zidx + nz][aidx];
			return Dn / ( 4 * Dfreevar[idx1 - 1] * Dfreevar[idx2 - 1] );
		}
	}
	else
	{
		nx = !(1^idx1);
		ny = !(2^idx1);
		nz = !(3^idx1);
		return                                                           \
		      (                                                          \
		       n_vicinity[this->xidx+nx][this->yidx+ny][this->zidx+nz] - \
		       n_vicinity[this->xidx-nx][this->yidx-ny][this->zidx-nz]   \
		      ) / ( 2 * Dfreevar[idx1 - 1] );
	}
}
// Change this->/*...*/ to /*...*/_vicinity[center][center][center]
double diff_V( double ***V_vicinity, uint8_t idx1 )
{
	nx = !(1^idx1);
	ny = !(2^idx1);
	nz = !(3^idx1);
	return                                \
	      (                               \
	       V_vicinity[1+nx][1+ny][1+nz] - \
	       V_vicinity[1-nx][1-ny][1-nz]   \
	      ) / ( 2 * Dfreevar[idx1 - 1] );
}
//
double d_Q_common      \
       (               \
        uint8_t aidx1, \
        uint8_t aidx2  \
       )
{
	return d_Q_common( NULL, Q, 0, 0, aidx1, aidx2 ) + KronekerD[aidx1][aidx2];
}
//
double d_Q_common      \
       (               \
        n_vicinity,    \
        d_Q Func2Call, \
        uint8_t didx1, \
        uint8_t aidx1, \
        uint8_t aidx2  \
       )
{
	return d_Q_common( n_vicinity, Func2Call, didx1, 0, aidx1, aidx2 );// This is a crutch.
}
//
double d_Q_common      \
       (               \
        n_vicinity,    \
        d_Q Func2Call, \
        uint8_t didx1, \
        uint8_t didx2, \
        uint8_t aidx1, \
        uint8_t aidx2  \
       )
{
	return SCALAR_ORDER_PARAMETER / 2 *                                 \
	       ( 3 * Func2Call( n_vicinity, didx1, didx2, aidx1, aidx2 ) );
}
//
double d2Q                           \
       (                             \
        LCMolDirector ***n_vicinity, \
        uint8_t didx1,               \
        uint8_t didx2,               \
        uint8_t aidx1,               \
        uint8_t aidx2                \
       )
{
	return diff_n( n_vicinity, aidx1, didx1 + 1, didx2 + 1 ) * \
	       this->values_[aidx2] +                              \
	       diff_n( n_vicinity, aidx2, didx1 + 1, didx2 + 1 ) * \
	       this->values_[aidx1] +                              \
	       diff_n( n_vicinity, aidx1, didx1 + 1 ) *            \
	       diff_n( n_vicinity, aidx2, didx2 + 1 ) +            \
	       diff_n( n_vicinity, aidx1, didx2 + 1) *             \
	       diff_n( n_vicinity, aidx2, didx1 + 1);
}
//
double dQ                            \
       (                             \
        LCMolDirector ***n_vicinity, \
        uint8_t didx1,               \
        uint8_t didx2,               \
        uint8_t aidx1,               \
        uint8_t aidx2                \
       )
{
	(void)didx2;// This is a crutch.
	//
	return diff_n( n_vicinity, aidx1, didx1 + 1 ) * \
	       this->values_[aidx2] +                   \
	       diff_n( n_vicinity, aidx2, didx1 + 1 ) * \
	       this->values_[aidx1];
}
//
double Q                             \
       (                             \
        LCMolDirector ***n_vicinity, \
        uint8_t didx1,               \
        uint8_t didx1,               \
        uint8_t aidx1,               \
        uint8_t aidx2                \
       )
{
	return this->values_[aidx1] * this->values_[aidx2];
}
//
typedef double ( *d_Q )( LCMolDirector ***n, uint8_t didx1, uint8_t didx2, uint8_t aidx, uint8_t aidx );
//typedef double ( *d2Q )( LCMolDirector ***n );
//
//dQ diff2Q[][]
//
//
//TODO
double dG1_dQ( LCMolDirector ***n_vicinity, uint8_t aidx1, uint8_t aidx2 )
{
	uint8_t idx_l;
	double __result = 0;
	//
	for ( idx_l = 0; idx_l < 2; idx_l++ )
		__result += d_Q_common( n_vicinity, d2Q, idx_l, idx_l, aidx1, aidx2 );
	//
	return (-2) * __result;
}
//
double dG2_dQ( LCMolDirector ***n_vicinity, uint8_t aidx1, uint8_t aidx2 )
{
	uint8_t idx_l;
	double __result = 0;
	//
	for ( idx_l = 0; idx_l < 2; idx_l++ )
		__result += d_Q_common( n_vicinity, d2Q, idx_l, aidx2, aidx1, idx_l );
	//
	return (-2) * __result;
}
//
double dG6_dQ( LCMolDirector ***n_vicinity, uint8_t aidx1, uint8_t aidx2 )
{
	uint8_t idx_l,
	        idx_m;
	double  __result = 0;
	//
	for ( idx_l = 0; idx_l < 2; idx_l++ )
		for ( idx_m = 0; idx_m < 2; idx_m++ )
			__result += d_Q_common( n_vicinity, dQ, aidx1, idx_l, idx_m ) *         \
			            d_Q_common( n_vicinity, dQ, aidx2, idx_l, idx_m ) -         \
			            d_Q_common( n_vicinity, dQ, idx_l, idx_l, idx_m ) *         \
			            d_Q_common( n_vicinity, dQ, idx_m, aidx1, aidx2 ) -         \
			            d_Q_common( idx_l, idx_m ) *                                \
			            d_Q_common( n_vicinity, d2Q, idx_m, idx_l, aidx1, aidx2 ) - \
			            d_Q_common( n_vicinity, dQ, idx_m, idx_l, idx_m ) *         \
			            d_Q_common( n_vicinity, dQ, idx_l, aidx1, aidx2 ) -         \
			            d_Q_common( idx_l, idx_m ) *                                \
			            d_Q_common( n_vicinity, d2Q, idx_l, idx_m, aidx1, aidx2 );
	//
	return __result;
}
//
double dG4_dQ( LCMolDirector ***n_vicinity, uint8_t aidx1, uint8_t aidx2 )
{
	uint8_t idx_l,
	        idx_m;
	double __result = 0;
	//
	for ( idx_l = 0; idx_l < 2; idx_l++ )
		for ( idx_m = 0; idx_m < 2; idx_m++ )
			__result += eLeviCivita[aidx1][idx_l][idx_m] * \
			            d_Q_common( n_vicinity, dQ, idx_l, idx_m, aidx2 );
	//
	return (-2) * __result;
}
//
template <typename T>
T [3][3][3]getVicinity( T ***field, uint64_t x, uint64_t y, uint64_t z, bool hasBoundary )
{
	T __result [3][3][3];
	//
	for ( int8_t i = -1; i < 2; i++ )
		for ( int8_t j = -1; j < 2; j++ )
			for ( int8_t k = -1; k < 2; k++ )
				__result[i + 1][j + 1][k + 1] = field[x + i + hasBoundary][y + j + hasBoundary][z + k + hasBoundary];
	//
	return __result;
}
//
//
double df_g__dQ( axisIndex aidx1, axisIndex aidx2 )
{
	double Spow2 = pow( SCALAR_ORDER_PARAMETER, 2 );
	double ***V_vicinity = getVicinity<double>( this->PotentialLCField, this->x(), this->y(), this->z(), false );
	LCMolDirector ***n_vicinity = getVicinity<LCMolDirector>( this->DirectorLCField, this->x(), this->y(), this->z(), true );
	//
	return 1 / 27 * ( K_33 - K_11 + 3 * K_22 ) *                         \
	       dG1_dQ( n_vicinity, aidx1, aidx2 ) / Spow2 +                  \
	       2 /  9 * ( K_11 - K_22 ) *                                    \
	       dG2_dQ( n_vicinity, aidx1, aidx2 ) / Spow2 +                  \
	       2 / 27 * ( K_33 - K_11 ) *                                    \
	       dG6_dQ( n_vicinity, aidx1, aidx2 ) /                          \
	       ( SCALAR_ORDER_PARAMETER * Spow2 ) +                          \
	       4 /  9 * K_22 * ( 2*M_PI / CHIRALITY_LC_PARAMETER ) *         \
	       dG4_dQ( n_vicinity, aidx1, aidx2 ) / Spow2 -                  \
	       1 /  3 * e0 * ( epar - eort ) * diff_V( V_vicinity, aidx1 ) * \
	                                       diff_V( V_vicinity, aidx2 );
}
// FDTD.

LCMolDirector ***getDirectorField( double ***ElPotentialField, uint64_t X, uint64_t Y, uint64_t Z )
{
	LCMolDirector ***__result,
	                ***n;
	double dt,
	        n_tilde,
	         gamma1;
	fg_t f_g[3] = {df_g__dn_x, df_g__dn_y, df_g__dn_z};
// Add 3 `for' loops above.
	__result = new LCMolDirector** [X];
	do
	{
		for ( i = 0; i < X; i++ )
		{
			for ( j = 0; j < Y; j++ )
			{
				for ( k = 0; k < Z; k++ )
				{
					for ( idx = 0; idx < 2; idx++ )
					{
						n_tilde = this->values_[idx] - deltaTau / ROTATIONAL_VISCOSITY * ( f_g[idx] )( n );
						this->values_[idx] = n_tilde / fabs( n_tilde );
					}
				}
			}
		}
		updatePotentials( this->PotentialLCField );
	}
	while ( getIsCconverge( t ) );
	return __result;
}
//
int main( int argc, char *argv[] )
{
	int i, j, n, m, k, kl, l;
	double **A, tmp;
	UnitSquareSimple ***ussDataTest;
	UnitSquare usItem, usData1, usData2, **usDistTest;
	sideMark *smDataTest;
	double **B;
//
/**/
	(void)argc; (void)argv;
	(void)A; (void)tmp; (void)B;
	(void)k; (void)i; (void)j; (void)m; (void)n; (void)l; (void)kl;
	(void)usItem; (void)usData1; (void)usData2; (void)usDistTest;
	(void)ussDataTest;
	(void)smDataTest;
/**/
//
	srand( time( 0 ) );
//
//	setlocale( LC_CTYPE , "");
/*
	wchar_t wcTest = 0x254B;
	setlocale( LC_CTYPE , "");
	wprintf( L"%lc", wcTest );
*/
/*
	printf("\nDistance constant test:\t\t%E\n",                                  \
	 getDistancePotentialConstant( UnitSquare( 0.0, 1.0, 1.0, 0.0, 0.0, 0 ) ) );
*/
//
/*
	kl = 2;
	n = ceil( (double)rand() / RAND_MAX * 10 );
	n = 10;
	m = ceil( (double)rand() / RAND_MAX * 10 );
	m = 12;
	smDataTest = new sideMark [kl];
	usDistTest = new UnitSquare* [kl];
	ussDataTest = new UnitSquareSimple** [kl];
	for ( l = 0; l < kl; l++ )
	{
		smDataTest[l] = (sideMark)ceil( (double)rand() / RAND_MAX * 5 );
		ussDataTest[l] = getMarkedGrid( n, m, k, usDistTest[l], smDataTest[l] );
	}
	B = getEquationMatrix( smDataTest, ussDataTest, l, k, n, m, usDistTest );
*/
/*
	wprintf( L"\nk = %i;\n", k );
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < m; j++ )
		{
			if ( j % 4 == 0 )
				wprintf( L" \\ \n" );
			if ( i + 1 == n )
				wprintf( L"  ussTest[%i][%i] = %llu", i, j, ussTest[i][j].mark );
			else
				wprintf( L"  ussTest[%i][%i] = %llu,", i, j, ussTest[i][j].mark );
		}
		wprintf( L"\n" );
	}
	wprintf( L"\n" );
	uint8_t tableLength = 5*m + 2;
	wchar_t *hLine = new wchar_t [tableLength];
	hLine[0] = 0x2554;
	for ( i = 1; i < tableLength - 2; i++ )
		if ( i % 5 )
			hLine[i] = 0x2550;
		else
			hLine[i] = 0x2566;
	hLine[i] = 0x2557;
	hLine[++i] = 0;
	wprintf( L"%ls", hLine );
	wprintf( L"\n" );
	hLine[0] = 0x2560;
	for ( i = 1; i < tableLength - 2; i++ )
		if ( i % 5 )
			hLine[i] = 0x2550;
		else
			hLine[i] = 0x256C;
	hLine[i] = 0x2563;
	for ( i = 0; i < n; i++ )
	{
		wprintf( L"\u2551" );
		for ( j = 0; j < m; j++ )
			wprintf( L" %2i \u2551", ussTest[i][j].mark );
		wprintf( L"\n" );
		if ( i + 1 < n )
			wprintf( L"%ls\n", hLine );
	}
	hLine[0] = 0x255A;
	for ( i = 5; i < tableLength - 2; i+=5 )
		hLine[i] = 0x2569;
	hLine[i] = 0x255D;
	wprintf( L"%ls", hLine );
	wprintf( L"\n" );
	wprintf( L"\n" );
	for ( i = 0; i < k; i++ )
	{
		if ( i % 4 == 0 )
			wprintf( L" \\ \n" );
		if ( i + 1 == k )
			wprintf( L"  [%i]: n = %f, m = %f", i, usDistTest[i].n, usDistTest[i].m );
		else
			wprintf( L"  [%i]: n = %f, m = %f;", i, usDistTest[i].n, usDistTest[i].m );
	}
	wprintf( L"\n" );
*/
/*
	n = 4;
	printf( "\nn = %i\n", n );
	A = new double*[n];
	B = new double*[n];
	for ( i = 0; i < n; i++ )
	{
		A[i] = new double[n+1];
		B[i] = new double[n+1];
		for ( j = 0; j < n+1; j++ )
			A[i][j] = B[i][j] = (int8_t)ceil( (double)rand() / RAND_MAX * 256 );
	}
	printf( "\n\n" );
	for ( i = 0; i < n; i++ )
	{
		printf( " Eq(%f,", A[i][n] );
		for ( j = 0; j < n; j++ )
			if ( j + 1 == n )
				printf( " %f*x%i ),", A[i][j], j );
			else
				printf( " %f*x%i +", A[i][j], j );
		printf( "\n" );
	}
*/
/*
	for ( l = 0; l < kl; l++ )
	{
		// B = ussDataTest[l];
		printf( "\n\n" );
		for ( i = 0; i < n; i++ )
		{
			printf( " Eq(%E,", 1.0 * B[i][n] );
			for ( j = 0; j < m; j++ )
				if ( j + 1 == m )
					printf( " %E*x%i ),", B[i][j], j );
				else
					printf( " %E*x%i +", B[i][j], j );
			printf( "\n" );
		}
	}
*/
/*
	getGaussianSolution( A, n );
	printf( "\n\n" );
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n + 1; j++ )
			if ( j + 1 == n )
				printf( " %f*x%i", A[i][j], j );
			else
				printf( " %f*x%i +", A[i][j], j );
		printf( "\n" );
	}
	printf( "\n" );
	for ( i = 0; i < n; i++ )
	{
		tmp = 0;
		for ( j = 0; j < n; j++ )
		{
			tmp += B[i][j] * A[j][j];
		}
		printf( "\ntmp%i = %f; A[%i][n] = %f\n", i, tmp, i, B[i][n] );
	}
*/
//
	maxSideInd = new int** [2];
	for (i = 0; i < 2; i++)
	{
		maxSideInd[i] = new int*[6];
		for (j = 0; j < 6; j++)
		{
			maxSideInd[i][j] = new int[2];
			maxSideInd[i][j][0] = N;
			maxSideInd[i][j][1] = M;
		}
	}
	for ( i = 1; i < 7; i++ )
		for ( j = 1; j < 7; j++ )
		{
			usData1 = UnitSquare( 1.0, 2.0, 3.0, 0.0, 1.0, 1.0, ceil( (double)rand() / RAND_MAX * 256 ), ceil( (double)rand() / RAND_MAX * 256 ), i, 0 );
			usData2 = UnitSquare( -0.3, 0.1, -2.0, 0.0, 1.0, 1.0, ceil( (double)rand() / RAND_MAX * 256 ), ceil( (double)rand() / RAND_MAX * 256 ), j, 1 );
			printf \
			 ( \
			  "\nUnit square 1: n = %f, m = %f, side = %s\nUnit square 2: n = %f, m = %f, side = %s\n", \
			  usData1.n, usData1.m, esSideMark[usData1.facet], \
			  usData2.n, usData2.m, esSideMark[usData2.facet]\
			 );
			usItem = getTargetUnitSquare( usData1, usData2 );
			printf                                                                    \
			 (                                                                        \
			   "\n\nUnit square orientation test %i.%i:\t\tz = %E, m = %f, n = %f\n", \
			   i, j, usItem.z, usItem.m, usItem.n                                     \
			 );
		}
//
	return 0;
}
//
//
