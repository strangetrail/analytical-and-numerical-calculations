#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "./libUtils.h"
#include "./libFFTcstm.h"

//using System;
//using System.Collections.Generic;
//using System.Text;
//using System.IO;

#define _USE_MATH_DEFINES
//#define RAND_MAX 100

//namespace _2dFFT
//{
//  class Program
//  {

const char *const_outFormatPrecision = "  % 1.4e %c %1.4ei",
           *mFileNameSuffix = ".m",
           *fftFileNameSuffix = ".txt",
           *mFileNamePrefix = "m_",
           *fftFileNamePrefix = "fft_";

//Матрица из "Сойфер. Методы компьютерной обработки изображений", стр. 317.
static double metric_tensor[4][4] = {{ 1,  1,  1,  1 },
                                     { 1,  1, -1, -1 },
                                     { 1, -1,  1, -1 },
                                     { 1, -1, -1,  1 }};
Complex operator* (const iComplex& x, const Complex& y)
{
	Complex c_res;
	if (( x.Re == 0.0 ) && ( x.Im == 1.0 ))
	{
		c_res.Re = -y.Im;
		c_res.Im = y.Re;
		return c_res;
	}
	if (( x.Re == 0.0 ) && ( x.Im == -1.0 ))
	{
		c_res.Re = y.Im;
		c_res.Im = -y.Re;
		return c_res;
	}
	if (( x.Im == 0.0 ) && ( x.Re == 1.0 ))
	{
		c_res.Re = y.Re;
		c_res.Im = y.Im;
		return c_res;
	}
	if (( x.Im == 0.0 ) && ( x.Re == -1.0 ))
	{
		c_res.Re = -y.Re;
		c_res.Im = -y.Im;
		return c_res;
	}
	c_res.Re = x.Re * y.Re - x.Im * y.Im;
	c_res.Im = x.Re * y.Im + x.Im * y.Re;
	return c_res;
}
Complex operator* (const Complex& x, const iComplex& y)
{
	return y * x;
}
Complex operator+ (const Complex& x, const Complex& y)
{
	Complex c_res;

	c_res.Re = x.Re + y.Re;
	c_res.Im = x.Im + y.Im;
	return c_res;
}
Complex operator+= (const Complex& x, const Complex& y)
{
	Complex c_res;

	c_res.Re = x.Re + y.Re;
	c_res.Im = x.Im + y.Im;
	return c_res;
}
Complex operator* (const Complex& x, const Complex& y)
{
	Complex c_res;
	c_res.Re = x.Re * y.Re - x.Im * y.Im;
	c_res.Im = x.Re * y.Im + x.Im * y.Re;
	return c_res;
}

Complex operator* (const Complex& x, const double& y)
{
	Complex c_res;
	c_res.Re = x.Re * y;
	c_res.Im = x.Im * y;
	return c_res;
}
Complex operator* (const double& x, const Complex& y)
{
	return y * x;
}
Complex operator/ (const Complex& x, const double& y)
{
	Complex c_res;
	c_res.Re = x.Re / y;
	c_res.Im = x.Im / y;
	return c_res;
}


        /*private static Complex[][] T_(short N, Complex[][] f_cplx)
        {
            Complex[][] local_temp = new Complex[N][];
            for (short i = 0; i < N; i++)
            {
                local_temp[i] = new Complex[N];
                for (short j = 0; j < N; j++)
                    local_temp[i][j] = new Complex();
            }
            for (short i = 0; i < N; i++)
                for (short j = 0; j < N - i; j++)
                {
                    local_temp[i][j].Re = f_cplx[N - 1 - j][N - 1 - i].Re;
                    local_temp[N - 1 - j][N - 1 - i].Re = f_cplx[i][j].Re;
                    local_temp[i][j].Im = f_cplx[i][j].Im;
                    local_temp[N - 1 - j][N - 1 - i].Im = f_cplx[N - 1 - j][N - 1 - i].Im;
                }
            return local_temp;
        }*/

        static void T_(short n, Complex **f_cplx)
        {
            double local_temp;
            for (short i = 0; i < n; i++)
                for (short j = i; j < n - 1; j++)
                {
                    local_temp = f_cplx[i][j + 1].Im;
                    f_cplx[i][j + 1].Im = f_cplx[j + 1][i].Im;
                    f_cplx[j + 1][i].Im = local_temp;
                }
        }

        static void Evaluation(short n, double w1, double w2, Complex **f_cplx)
        {
            short  a, b, n1, n2, length, i, j;
            double pow1, pow2, w1_temp, w2_temp, tmpln;
            Complex **quad_f_cplx,
                    m;
            if (n > 2)//Ограничение рек. вызовов
            {
                length = sizeof( Complex* );
                Complex **temp_f_cplx = (Complex **)calloc( n, length );//Для хранения оригинального параметра f_cplx.
                length = sizeof( Complex );
                for (i = 0; i < n; i++)
                    temp_f_cplx[i] = (Complex *)calloc( n, length );
                for (n1 = 0; n1 < n; n1++)
                    for (n2 = 0; n2 < n; n2++)
                    {
                        temp_f_cplx[n1][n2].Re = f_cplx[n1][n2].Re;
                        temp_f_cplx[n1][n2].Im = f_cplx[n1][n2].Im;
                    }
                length = sizeof( Complex * );
                quad_f_cplx = (Complex **)calloc( n / 2, length );
                for (i = 0; i < n / 2; i++)
                    quad_f_cplx[i] = (Complex *)calloc( n / 2, sizeof( Complex ) );
                U_CHAR count = 0;
                for(a = 0; a < 2; a++)
                    for (b = 0; b < 2; b++)
                    {
                        //Все замечания я делал для себя.
                        //Замечание 0
                        //Возможно есть способ обойтись без создания матрицы для дополнительных вычислений,
                        //но не факт, т.к. другой способ передать рекурсивно подмножество индексов,
                        //не имеющих строго прямоугольного геометрического отображения в исходной
                        //матрице мне не представляется.
                        for (n1 = 0; n1 < n / 2; n1++)
                            for (n2 = 0; n2 < n / 2; n2++)
                            {
                                quad_f_cplx[n1][n2].Re = temp_f_cplx[2 * n1 + a][2 * n2 + b].Re;
                                quad_f_cplx[n1][n2].Im = temp_f_cplx[2 * n1 + a][2 * n2 + b].Im;
                            }
                        //Замечание 0 end
                        w1_temp = 1;
                        Evaluation((short)trunc( n / 2 ), w1 * w1, w2 * w2, quad_f_cplx);
                        for (i = 0; i < n / 2; i++)
                        {
                            w2_temp = 1;
                            for (j = 0; j < n / 2; j++)
                            {
                                //w1_temp = Math.Pow(W1, i);
                                //W2_temp = Math.Pow(W2, j);
                                if (a == 0)
                                    pow1 = 1;
                                else
                                    pow1 = w1_temp;
                                if (b == 0)
                                    pow2 = 1;
                                else
                                    pow2 = w2_temp;
                                //Замечание 1
                                //Не забыть обнулить как-нибудь по-другому, если получится, f[i][j]
                                //и её зависимые f[i+...][j+...], если a=b=0.
                                if ((a == 0) && (b == 0))
                                {
                                    f_cplx[i][j].Re = 0;
                                    f_cplx[i + n / 2][j].Re = 0;
                                    f_cplx[i][j + n / 2].Re = 0;
                                    f_cplx[i + n / 2][j + n / 2].Re = 0;
                                    f_cplx[i][j].Im = 0;
                                    f_cplx[i + n / 2][j].Im = 0;
                                    f_cplx[i][j + n / 2].Im = 0;
                                    f_cplx[i + n / 2][j + n / 2].Im = 0;
                                }
                                //Замечание 1 end

                                //Далее - рассчёт по формуле 5.46 из
                                //"Сойфер. Методы компьютерной обработки изображений", стр. 317;
                                //при условии, что отдельные компоненты левого множителя вычисляются
                                //после каждого вызова процедуры Evaluation(...).

                                tmpln = log( pow1 * pow2 );
                                m.Re = cos( tmpln );
                                m.Im = sin( tmpln );
                                f_cplx[i][j].Re += quad_f_cplx[i][j].Re * m.Re - quad_f_cplx[i][j].Im * m.Im;
                                f_cplx[i][j].Im += quad_f_cplx[i][j].Re * m.Im + quad_f_cplx[i][j].Im * m.Re;
                                f_cplx[i + n / 2][j].Re += metric_tensor[1][count] * ( quad_f_cplx[i][j].Re * m.Re - quad_f_cplx[i][j].Im * m.Im );
                                f_cplx[i + n / 2][j].Im += metric_tensor[1][count] * (quad_f_cplx[i][j].Re * m.Im + quad_f_cplx[i][j].Im * m.Re);
                                f_cplx[i][j + n / 2].Re += metric_tensor[2][count] * (quad_f_cplx[i][j].Re * m.Re - quad_f_cplx[i][j].Im * m.Im);
                                f_cplx[i][j + n / 2].Im += metric_tensor[2][count] * (quad_f_cplx[i][j].Re * m.Im + quad_f_cplx[i][j].Im * m.Re);
                                f_cplx[i + n / 2][j + n / 2].Re += metric_tensor[3][count] * (quad_f_cplx[i][j].Re * m.Re - quad_f_cplx[i][j].Im * m.Im);
                                f_cplx[i + n / 2][j + n / 2].Im += metric_tensor[3][count] * (quad_f_cplx[i][j].Re * m.Im + quad_f_cplx[i][j].Im * m.Re);
                                w2_temp = w2_temp * w2;
                            }
                            w1_temp = w1_temp * w1;
                        }
                        count++;
                    }
            }
            else
            {
                quad_f_cplx = (Complex **)calloc( 2, sizeof (Complex *) );
                for (i = 0; i < 2; i++)
                    quad_f_cplx[i] = (Complex *)calloc( 2, sizeof (Complex) );
                for (n1 = 0; n1 < 2; n1++)
                    for (n2 = 0; n2 < 2; n2++)
                    {
                        quad_f_cplx[n1][n2].Re = f_cplx[n1][n2].Re;
                        quad_f_cplx[n1][n2].Im = f_cplx[n1][n2].Im;
                    }
                //Замечание 2
                //Не знаю стоит ли делать завершение рекурсивных вызовов как-нибудь иначе,
                //просто если останавливаться на матрице 1х1, можно получить лишние возвраты из
                //рекурсивных процедур, а также лишние вычисления в основной части алгоритма
                //при работе с матрицой размера 2х2.
                f_cplx[0][0].Re = quad_f_cplx[0][0].Re + quad_f_cplx[0][1].Re + quad_f_cplx[1][0].Re + quad_f_cplx[1][1].Re;
                f_cplx[0][0].Im = quad_f_cplx[0][0].Im + quad_f_cplx[0][1].Im + quad_f_cplx[1][0].Im + quad_f_cplx[1][1].Im;

                m.Re = -1;
                m.Im = 0;

                f_cplx[0][1].Re = quad_f_cplx[0][0].Re + (quad_f_cplx[0][1].Re * m.Re - quad_f_cplx[0][1].Im * m.Im);
                f_cplx[0][1].Re += quad_f_cplx[1][0].Re + (quad_f_cplx[1][1].Re * m.Re - quad_f_cplx[1][1].Im * m.Im);
                f_cplx[0][1].Im = quad_f_cplx[0][0].Im + (quad_f_cplx[0][1].Im * m.Re + quad_f_cplx[0][1].Re * m.Im);
                f_cplx[0][1].Im += quad_f_cplx[1][0].Im + (quad_f_cplx[1][1].Im * m.Re + quad_f_cplx[1][1].Re * m.Im);

                f_cplx[1][0].Re = quad_f_cplx[0][0].Re + (quad_f_cplx[1][0].Re * m.Re - quad_f_cplx[1][0].Im * m.Im);
                f_cplx[1][0].Re += quad_f_cplx[0][1].Re + (quad_f_cplx[1][1].Re * m.Re - quad_f_cplx[1][1].Im * m.Im);
                f_cplx[1][0].Im = quad_f_cplx[0][0].Im + (quad_f_cplx[1][0].Im * m.Re + quad_f_cplx[1][0].Re * m.Im);
                f_cplx[1][0].Im += quad_f_cplx[0][1].Im + (quad_f_cplx[1][1].Im * m.Re + quad_f_cplx[1][1].Re * m.Im);

                f_cplx[1][1].Re = quad_f_cplx[0][0].Re + (quad_f_cplx[1][0].Re * m.Re - quad_f_cplx[1][0].Im * m.Im);
                f_cplx[1][1].Im = quad_f_cplx[0][0].Im + (quad_f_cplx[1][0].Im * m.Re + quad_f_cplx[1][0].Re * m.Im);
                f_cplx[1][1].Re += (quad_f_cplx[0][1].Re * m.Re - quad_f_cplx[0][1].Im * m.Im);
                f_cplx[1][1].Im += (quad_f_cplx[0][1].Im * m.Re + quad_f_cplx[0][1].Re * m.Im);

                m.Re = 1;
                m.Im = 0;

                f_cplx[1][1].Re += (quad_f_cplx[1][1].Re * m.Re - quad_f_cplx[1][1].Im * m.Im);
                f_cplx[1][1].Im += (quad_f_cplx[1][1].Im * m.Re + quad_f_cplx[1][1].Re * m.Im);
                //Замечание 2 end
            }
        }

				void printResults( char *chrs_filename, char *outFormatPrecision, short ncols, int n, Complex **f_cplx )
				{
					char *fileName,
					      complexImSign;
					FILE *file_handler;
//
// printf( 
// STRINGIZE_TOKEN( 
// 					BGNINIT 
// 						PPT(char) PPD(9+strlen(chrs_filename)) PPA(fileName)
// 					ENDINIT 
// ) 
// );
//
//
fileName = ( char * )malloc( 9 + strlen( chrs_filename ) );
//
//
					// getEmptyVectors< char >( 9 + strlen( chrs_filename ), &fileName, NULL );
					sprintf( fileName, "%s%s%s", fftFileNamePrefix, chrs_filename, fftFileNameSuffix );
					file_handler = fopen( fileName, "w" );
					fprintf( file_handler, "\n" );
					for ( int k = 0; k < n / ncols + (( n % ncols ) > 0 ); k++ )
					{
						for ( int i = 0; i < n; i++ )
						{
							for ( int j = ncols * k; ( j < ncols * ( k + 1 )) && ( j < n ); j++ )
							{
								if ( f_cplx[i][j].Im >= 0 )
									complexImSign = '+';
								else
									complexImSign = '-';
								fprintf( file_handler, ( const char * )outFormatPrecision, f_cplx[j][i].Re, complexImSign, fabs( f_cplx[i][j].Im ));
							}
							fprintf( file_handler, "\n" );
						}
						fprintf( file_handler, "\n\n\n" );
					}
					fprintf( file_handler, "\n" );
					fclose( file_handler );
				}

        Complex **doFFT( short inv, int n, Complex **f_cplx, char *chrs_filename, char *expFormat )
        {
            try
            {
                char  arrayParenthesis,
                     *outFormatPrecision,
                     *fileName;
                short /*k,*/
                      length/*,
                      size = 36*/;
//              clock_t begin_time, end_time;
                int i,
                    j;
                FILE *file_handler;
                //Два параметра W1 и W2 нужны, если матрица не квадратная, в квадратной они равны.
                double w1 = exp ((-2 + 4 * inv) * M_PI / n),
                       w2 = exp ((-2 + 4 * inv) * M_PI / n);
                length = sizeof (Complex *);
                Complex **f_copy = (Complex **)calloc( n, length );
                length = sizeof (Complex);
                for (i = 0; i < n; i++)
                {
                    f_copy[i] = (Complex *)calloc( n, length );
                }
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                    {
                        f_copy[i][j].Re = f_cplx[i][j].Re;
                        f_copy[i][j].Im = f_cplx[i][j].Im;
                    }
                Evaluation( n, w1, w2, f_copy );
                T_( n, f_copy );
                printf( "Creating Log file ... \n" );
// 								printf( 
// STRINGIZE_TOKEN( 
// 								BGNINIT 
// 									PPA(char) PPD(1+strlen(chrs_filename)) PPA(fileName)
// 								ENDINIT 
// ) 
// );
//
//
fileName = ( char * )malloc( 1 + strlen( chrs_filename ) );
//
//
                // getEmptyVectors< char >( 1 + strlen( chrs_filename ), &fileName, NULL );
                sprintf( fileName, "%s%s%s", mFileNamePrefix, chrs_filename, mFileNameSuffix );
                file_handler = fopen( fileName, "w" );
                fprintf( file_handler, "format SHORTE\nI = [" );
                arrayParenthesis = ' ';
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                        fprintf( file_handler, "  % 1.4e%+1.4ei", f_cplx[j][i].Re, f_cplx[j][i].Im );
                    if ( i == n - 1 )
                      arrayParenthesis = ']';
                    fprintf( file_handler, "%c;\n", arrayParenthesis );
                }
                if ( inv )
                  fprintf( file_handler, "J = ifft2( I );\ndisp( J );\n" );
                else
                  fprintf( file_handler, "J = fft2( I );\ndisp( J );\n" );
                fclose( file_handler );
                if (expFormat != NULL)
                {
                  outFormatPrecision = ( char * )calloc( 2 * strlen( expFormat ) + 8, sizeof ( char ) );
                  sprintf( outFormatPrecision, "  %s %s %si", expFormat, "%c", expFormat );
                }
                else
                  outFormatPrecision = ( char * )const_outFormatPrecision;
                printResults( chrs_filename, outFormatPrecision, 6, n, f_copy );
                return f_copy;
            }
            catch( ... )
            {
                printf( "\nUnhandled exception.\n" );
                return NULL;
            }
        }
//  }
//}
