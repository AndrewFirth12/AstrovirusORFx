#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
using namespace std;

int main(int argc, char* argv[])
{

  if (4 != argc) {
    cerr << "Usage '" << argv[0] << " seed nbootstrap ntrials'.\n";
    exit(EXIT_FAILURE);
  }

  int i, j, k, nlines, codonstats[3][10000], sum0, sum1, sum2, r1;
  int seed = atoi(argv[1]);
  int nbootstrap = atoi(argv[2]);
  int ntrials = atoi(argv[3]);

  for (k = 0; k < 10000; ++k) {
    codonstats[0][k] = codonstats[1][k] = codonstats[2][k] = 0;
  }
  
  ifstream statsfile("inputstats.txt");
  if (!statsfile) {
    cerr << "Aborting: Can't find file 'inputstats.txt'.\n";
    exit(EXIT_FAILURE);
  }
  nlines = -1;
  while (statsfile.ignore(1000, '\n')) {
    ++nlines;
  }
  statsfile.clear();
  statsfile.seekg(0);
  //  cout << "Found " << nlines << " lines in codon stats file.\n";

  for (k = 0; k < nlines; ++k) {
    statsfile >> codonstats[0][k];
    statsfile >> codonstats[1][k];
    statsfile >> codonstats[2][k];
    statsfile.ignore(1000, '\n');
    //    cout << codonstats[0][k] << " " << codonstats[1][k] << " "
    //	 << codonstats[2][k] << "\n" ;
  }
		     
  srand(seed);

  for (j = 0; j < ntrials; ++j) {
    sum0 = sum1 = sum2 = 0;
    // The random numbers range from 0 to nlines - 1 so can be used as indices
    // in the codonstats matrix.
    for (i = 0; i < nbootstrap; ++i) {
      r1 = int(nlines * float(rand())/float(RAND_MAX));
      sum0 += codonstats[0][r1];
      sum1 += codonstats[1][r1];
      sum2 += codonstats[2][r1];
    }
    cout << sum0 << " " << sum1 << " " << sum2 << "\n";
  }
  
}
