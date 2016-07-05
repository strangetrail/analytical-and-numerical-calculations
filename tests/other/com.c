#include <stdio.h>

template<int A,int B>
struct S
{
  float f1[A];
  long f0;
  float f2[B];
};

int main(int argc, char *argv[])
{
  S<3,5> test1;
  S<5,6> test2;
  test1.f0 = 3;
  test2.f0 = 4;
  printf ( "%d\n", test1.f0 );
  printf ( "%d\n", ((S<5,6>)test1).f0 );
  return 0;
}
