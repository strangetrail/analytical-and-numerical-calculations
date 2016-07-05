//
/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                                                                          */
/*                               libUtils.h                                 */
/*                                                                          */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
//
/*
void getEmptyArrays( unsigned int typeSize, ... ); // [unsigned char typeSize; unsigned short [d1, d2, .., dN], 0; long address; ...]; unsigned char 0
*/
#define U_CHAR unsigned char
#define U_SHORT unsigned short
#define U_INT unsigned int
#define U_LONG unsigned long
#define L_DOUBLE double
/*
#define TOKEN_TO_STRING(TOK) # TOK
#define STRINGIZE_TOKEN(TOK,...) TOKEN_TO_STRING(TOK)
//
#define PP_NARG(...) PP_NARG_(__VA_ARGS__,PP_RSEQ_N())
#define PP_NARG_(...) PP_ARG_N(__VA_ARGS__)
#define PP_ARG_N(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,_12,_13,_14,_15,_16,_17,_18,_19,_20,_21,_22,_23,_24,_25,N,...) N
#define PP_RSEQ_N() 25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
//
#define EMPTY(...)
#define DEFER(...) __VA_ARGS__ EMPTY()
#define OBSTRUCT(...) __VA_ARGS__ DEFER(EMPTY)()
#define EXPAND(...) __VA_ARGS__
//
#define EVAL(...)  EVAL1(EVAL1(EVAL1(__VA_ARGS__)))
#define EVAL1(...) EVAL2(EVAL2(EVAL2(__VA_ARGS__)))
#define EVAL2(...) EVAL3(EVAL3(EVAL3(__VA_ARGS__)))
#define EVAL3(...) EVAL4(EVAL4(EVAL4(__VA_ARGS__)))
#define EVAL4(...) EVAL5(EVAL5(EVAL5(__VA_ARGS__)))
#define EVAL5(...) __VA_ARGS__
//
#define CAT(a, ...) PRIMITIVE_CAT(a, __VA_ARGS__)
#define PRIMITIVE_CAT(a, ...) a ## __VA_ARGS__
//
#define INC(x) PRIMITIVE_CAT(INC_, x)
#define INC_0 1
#define INC_1 2
#define INC_2 3
#define INC_3 4
#define INC_4 5
#define INC_5 6
#define INC_6 7
#define INC_7 8
#define INC_8 9
#define INC_9 9
//
#define DEC(x) PRIMITIVE_CAT(DEC_, x)
#define DEC_0 0
#define DEC_1 0
#define DEC_2 1
#define DEC_3 2
#define DEC_4 3
#define DEC_5 4
#define DEC_6 5
#define DEC_7 6
#define DEC_8 7
#define DEC_9 8
//
#define CHECK_N(x, n, ...) n
#define CHECK(...) CHECK_N(__VA_ARGS__, 0,)
//
#define NOT(x) CHECK(PRIMITIVE_CAT(NOT_, x))
#define NOT_0 ~, 1,
//
#define COMPL(b) PRIMITIVE_CAT(COMPL_, b)
#define COMPL_0 1
#define COMPL_1 0
//
#define BOOL(x) COMPL(NOT(x))
//
#define IIF(c) PRIMITIVE_CAT(IIF_, c)
#define IIF_0(t, ...) __VA_ARGS__
#define IIF_1(t, ...) t
//
#define IF(c) IIF(BOOL(c))
//
#define EAT(...)
#define EXPAND(...) __VA_ARGS__
#define WHEN(c) IF(c)(EXPAND, EAT)
//
#define REPEATX(X, count, macro, ...) \
	WHEN(count) \
	( \
		OBSTRUCT(PRIMITIVE_CAT(REPEAT_INDIRECT,X)) () \
		( \
			DEC(count), macro, __VA_ARGS__ \
		) \
		OBSTRUCT(macro) \
		( \
			DEC(count), __VA_ARGS__ \
		) \
	)
//
#define REPEAT1(count, macro, ...) \
	REPEATX(1,count,macro,__VA_ARGS__)
#define REPEAT2(count, macro, stub, ...) \
	REPEATX(2,count,macro,__VA_ARGS__)
#define REPEAT3(count, macro, c, stub, ...) \
	REPEATX(3,count,macro,c,__VA_ARGS__)
//
#define REPEAT_INDIRECT1() REPEAT1
#define REPEAT_INDIRECT2() REPEAT2
#define REPEAT_INDIRECT3() REPEAT3
//
// #define M(i, _) i
// EVAL(REPEAT(8, M, ~)) // 0 1 2 3 4 5 6 7
//
#define MAP(pref,data,suff) pref data suff
#define SUFFZZZ )
//
#define PREFTYPE (U_CHAR)(sizeof(
#define SUFFTYPE ))
#define PPT(a) MAP(PREFTYPE,a,SUFFTYPE),
//
#define PREFDIM (U_INT)(
#define PPD(a) MAP(PREFDIM,a,SUFFZZZ),
//
#define BGNINIT getEmptyArrays(
#define ENDINIT MAP(PREFDIM,0,SUFFZZZ) );
*/
/*
#define PPDI1(i, a) PPD(a)
#define PPDR(n,a) EVAL(REPEAT1(n, PPDI1, a))
#define PPDI2(i, a,...) PPD(a)
#define PPDL(...) EVAL(REPEAT2(PP_NARG(__VA_ARGS__), PPDI2, 0, __VA_ARGS__, 0))
*/
/*
#define PREFADDR PPD(0) (U_LONG)&(
#define PPA(a) MAP(PREFADDR,a,SUFFZZZ),
//
#define DD(i,PPL,a,...) PPL () PPA(a)
#define DDL(PPL,...) EVAL(REPEAT3(PP_NARG(__VA_ARGS__), DD, PPL, 0, __VA_ARGS__, 0))
//
// #define MACROSTOPREPEND() PPT(typeName), PPD(N),
// DDL(MACROSTOPREPEND,ptr1Name,ptr2Name,ptr3Name)
*/
