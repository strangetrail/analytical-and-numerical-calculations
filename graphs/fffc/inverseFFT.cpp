//using System;
//using System.Collections.Generic;
//using System.Text;
//using System.IO;

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define _USE_MATH_DEFINES
#define U_CHAR unsigned char
//#define RAND_MAX 100

//namespace _2dFFT
//{
//  class Program
//  {
        struct Complex
        {
            double Re;
            double Im;
        };

        //Матрица из "Сойфер. Методы компьютерной обработки изображений", стр. 317.
        static double metric_tensor[4][4] = {{ 1, 1, 1, 1 },
                                             { 1, 1, -1, -1 },
                                             { 1, -1, 1, -1 },
                                             { 1, -1, -1, 1 } };

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

        int main(int argc, char **argv)
        {
            try
            {
                const char *const_outFormatPrecision = "  % 1.4e %c %1.4ei";
                char  arrayParenthesis,
                      complexImSign,
                     *outFormatPrecision;
                short n, k,
                      length,
                      size = 36;
                clock_t begin_time, end_time;
                int i,
                    j;
                FILE *file_handler;
                //Чтобы при тестировании не вводить N в консоли.
                if (argc > 1)
                    n = atoi(argv[1]);
                else
                    n = 16;//Здесь задавать значение N (N - степень двойки).
                //Два параметра W1 и W2 нужны, если матрица не квадратная, в квадратной они равны.
                double w1 = exp (2 * M_PI / n), // test here
                       w2 = exp (2 * M_PI / n), // test here
                       elapsed_time;
                length = sizeof (Complex *);
                Complex **f_cplx = (Complex **)calloc( n, length ),
                        **f_orig = (Complex **)calloc( n, length );
                length = sizeof (Complex);
                for (i = 0; i < n; i++)
                {
                    f_cplx[i] = (Complex *)calloc( n, length );
                    f_orig[i] = (Complex *)calloc( n, length );
                }
                srand (time (NULL));
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                    {
                        f_cplx[i][j].Re = ((double) rand()) / RAND_MAX * 2000;
                        f_cplx[i][j].Im = ((double) rand()) / RAND_MAX * 2000;
                        f_orig[i][j].Re = f_cplx[i][j].Re;
                        f_orig[i][j].Im = f_cplx[i][j].Im;
                    }
                begin_time = clock();
                Evaluation( n, w1, w2, f_cplx );
                int n2 = n * n;
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                    {
                        f_cplx[i][j].Re = f_cplx[i][j].Re / n2;
                        f_cplx[i][j].Im = f_cplx[i][j].Im / n2;
                    }
                T_( n, f_cplx );
                end_time = clock() - begin_time;
                elapsed_time = double( end_time ) /  CLOCKS_PER_SEC;
                printf( "\n" );
                printf( "Elapsed time: %f\n", elapsed_time );
                printf( "Creating Log file ... \n" );
                file_handler = fopen("output.dat", "w");
                fprintf( file_handler, "Start here:\nformat SHORTE\nI = [" );
                arrayParenthesis = ' ';
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                        fprintf( file_handler, "  % 1.4e%+1.4ei", f_orig[j][i].Re, f_orig[j][i].Im );
                    if ( i == n - 1 )
                      arrayParenthesis = ']';
                    fprintf( file_handler, "%c;\n", arrayParenthesis );
                }
                fprintf( file_handler, "J = ifft2( I );\ndisp( J );\n\n" );
                fprintf( file_handler, "fft:\n" );
                if (argc > 2)
                {
                  outFormatPrecision = ( char * )calloc( 2 * strlen( argv[2] ) + 8, sizeof ( char ) );
                  sprintf( outFormatPrecision, "  %s %s %si", argv[2], "%c", argv[2] );
                }
                else
                  outFormatPrecision = ( char * )const_outFormatPrecision;
                short ncols = 6;
                for ( k = 0; k < n / ncols + ((n % ncols) > 0); k++ )
                {
                  for (i = 0; i < n; i++)
                  {
                      for ( j = ncols * k; (j < ncols * (k + 1)) && (j < n); j++ )
                      {
                        if (f_cplx[i][j].Im >= 0)
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
                printf( "done.\n" );
            }
            catch( ... )
            {
                printf( "\nUnhandled exception.\n" );
            }
        }
//  }
//}
