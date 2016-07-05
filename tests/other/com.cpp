#include <stdio.h>

struct SBaseNonTemplatePolymorphic
{
private:
  virtual void makeMePolymorphic () {};
};

struct SBaseNonTemplate : public SBaseNonTemplatePolymorphic
{
  double r;
};

template<int A>
struct SBase : public SBaseNonTemplate
{
  float f1[A];
};

template<int A,int B>
struct S : public SBase<A>
{
  long f0;
  float f2[B];
};

int main(int argc, char *argv[])
{
  (void)argc;
  (void)argv;

  S<3,5> test11;
  S<5,6> test22;

  test11.f0 = 3;
  test22.f0 = 4;

  SBaseNonTemplate& test1 = test11;
  SBaseNonTemplate &test2 = test22;

  S<3,5> *cast1 = dynamic_cast<S<3,5>*>(&test1);
  if ( cast1 != NULL )
    printf ( "%ld\n", cast1->f0 );
  S<5,6> *cast2 = dynamic_cast<S<5,6>*>(&test1);
  if ( cast2 != NULL )
    printf ( "%ld\n", cast2->f0 );

  printf ( "%ld\n", ((S<5,6> &)test2).f0 );

  return 0;
}
