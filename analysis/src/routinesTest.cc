//
//  routinesTest.cpp
//  main program to test main functions for
//  high level usage of LatticeTester
//

// Include Header
#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

// Include LatticeTester Header
#include "latticetester/LatticeTesterRoutines.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;


//==================================================================================

#if 0
int main (int argc, char *argv[])
{
   // all the parameters
   BMat matrix;
   NormType norm = L2NORM;
   NormaType normalizerType = BESTLAT;
   long maxNodesBB = 5;
   PreReductionType preRed = BKZ; //LenstraLL;
   PrecisionType doublePrecision = QUADRUPLE; 
   double fact = 0.9;
   int blocksize = 10;

   // initialization of matrix B

   /*
   matrix.resize(2,2);
   matrix[0][0]=1; matrix[0][1]=12;
   matrix[1][0]=0; matrix[1][1]=101;
   */
   
   matrix.resize(8,8);
   //dual
   matrix[0][0]=9223372036854773561;    
   matrix[1][1]=9223372036854773561;   
   matrix[2][2]=9223372036854773561;
   matrix[3][3]=1;
   matrix[4][4]=1;
   matrix[5][5]=1;
   matrix[6][6]=1;
   matrix[7][7]=1;
   matrix[3][0]=-9222187883300163885; matrix[3][1]=0;                    matrix[3][2]=-1145902849652723;
   matrix[4][0]=-6629950407641799081; matrix[4][1]=-9222187883300163885; matrix[4][2]= -6821832294614679088;
   matrix[5][0]=-7721957468801369126; matrix[5][1]=-6629950407641799081; matrix[5][2]=-2092544303070727312;
   matrix[6][0]=-8814617570338055415; matrix[6][1]=-7721957468801369126; matrix[6][2]=-1127730827481453799;
   matrix[7][0]=-3020738662083495906; matrix[7][1]=-8814617570338055415; matrix[7][2]=-1975258172869346907;
   //primal
   /*
   matrix[0][0]=1;
   matrix[1][1]=1;
   matrix[2][2]=1;
   matrix[3][3]=9223372036854773561;
   matrix[4][4]=9223372036854773561;
   matrix[5][5]=9223372036854773561;
   matrix[6][6]=9223372036854773561;
   matrix[7][7]=9223372036854773561;
   matrix[0][3]=9222187883300163885; matrix[1][3]=0;                   matrix[2][3]=1145902849652723;
   matrix[0][4]=6629950407641799081; matrix[1][4]=9222187883300163885; matrix[2][4]=6821832294614679088;
   matrix[0][5]=7721957468801369126; matrix[1][5]=6629950407641799081; matrix[2][5]=2092544303070727312;
   matrix[0][6]=8814617570338055415; matrix[1][6]=7721957468801369126; matrix[2][6]=1127730827481453799;
   matrix[0][7]=3020738662083495906; matrix[1][7]=8814617570338055415; matrix[2][7]=1975258172869346907;
   */
   
   
   cout << "matrix = \n" << matrix << endl;

   // printing the length of shortest vector
   cout << "Length = " << ShortestVector(matrix, norm) << endl; //, preRed, doublePrecision, fact, blocksize) << endl;

   // printing the FoM
   cout << "FoM = " << FigureOfMerit(matrix, normalizerType) << endl; //, preRed, doublePrecision, blocksize) << endl;
   
   // printing the FoM Beyer
   cout << "FoM Beyer = " << FigureOfMeritBeyer(matrix) << endl; //, preRed, doublePrecision, blocksize) << endl;
   
   return 0;
}
#endif

