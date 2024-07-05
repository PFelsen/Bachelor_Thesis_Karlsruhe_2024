/***************************************************************************/
/* Adaption of http://sprng.org/Version5.0/examples5/sprng-simple_mpi.cpp  */
/*           Demonstrates sprng use with one stream per process            */
/* A distinct stream is created on each process, then prints a few         */
/* random numbers. MPI is already initialized by TestEnvironment.          */
/***************************************************************************/

#include "TestEnvironment.hpp"

#include <cstdio>
#include <iostream>


#define SIMPLE_SPRNG /* simple interface                        */
#define USE_MPI      /* use MPI to find number of processes     */
#include "sprng.h"


#define SEED 985456376

using namespace std;

TEST(TestSPRNG5, SimpleInterface) {
  double rn;
  int i;
  int myid = PPM->Proc(0);
  int gtype = 0; // lfg: 0  lcg: 1  lcg64: 2  cmrg: 3  mlfg: 4  pmlcg: 5

  PPM->BCastOnCommSplit(gtype, 0);

  init_sprng(SEED, SPRNG_DEFAULT, gtype); /* initialize stream        */
  printf("\n\nProcess %d, print information about stream:\n", myid);
  print_sprng();

  /************************ print random numbers ***************************/

  for (i = 0; i < 3; i++) {
    rn = sprng(); /* generate double precision random number */
    printf("Process %d, random number %d: %.14f\n", myid, i + 1, rn);
  }
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv).WithPPM()).RUN_ALL_MPP_TESTS();
}
