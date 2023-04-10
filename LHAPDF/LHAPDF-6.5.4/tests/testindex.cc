// Example program to test LHAPDF index lookup and loading

#include "LHAPDF/PDFIndex.h"
#include <iostream>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
using namespace LHAPDF;
using namespace std;

void lookup(int id) {
  pair<string, int> set_id = lookupPDF(id);
  cout << "ID=" << id << " -> set=" << set_id.first << ", mem=" << set_id.second << endl;
}


int main(int argc, char* argv[]) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  lookup(10800);
  lookup(10801);
  lookup(10042);
  lookup(10041);
  lookup(10799);
  lookup(12346);
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
