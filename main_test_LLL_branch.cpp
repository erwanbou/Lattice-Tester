//
//  main.cpp
//  Lattice Tester
// 

#include <iostream>
#include <iomanip>
#include <string>
#include <iterator>
#include <set>
#include <time.h>

#include "latticetester/Normalizer.h"
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Basis.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Reducer.h"
#include "latticetester/Coordinates.h"
#include "latticetester/ntlwrap.h"

#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include "NTL/vec_ZZ.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include <boost/progress.hpp>

using namespace std;
using namespace LatticeTester;
using namespace NTL;


void RandomMatrix (mat_ZZ& A, ZZ& det, int min, int max, int seed){
    
    int dim = (int) A.NumRows() ;
    srand(seed);
    
    do{
        for (int i = 1; i < dim; i++){
            for (int j = 1; j < dim; j++)
                A[i][j] = min + (rand() % (int)(max - min + 1));
        }
        // Richard implementation of basis include an empty first column
        // and an empty first lign, so we had a 1 on top left position
        // to preserve determinant calculation
        A[0][0] = 1;
        
        det = determinant(A);
        
    } while ( det == 0 );
    
    // we set the top left coefficient back to 0
    A[0][0] = 0;
}

void RealLatticePrint (IntLattice& lattice)
{
    int dim = lattice.getPrimalBasis().getDim();
    
    for (int i = 0; i < dim+1; i++){
        cout << "[";
        for (int j = 0; j < dim+1; j++)
            cout << lattice.getPrimalBasis()[i][j] << " ";
        cout << "]" << endl;
    }
    cout << endl;
}

template<typename Type, long Size>
void print(string name, Type const(& array)[Size], bool isIntegerOutput) {
    cout << name << " = ";
    for(int i=0; i<Size; i++){
        if (isIntegerOutput)
            cout << conv<ZZ>(array[i]) << " ";
        else
            cout << array[i] << " ";
    }
    //cout << endl;
}

template<typename Type, long Size>
Type Average(Type const(& array)[Size]) {
    Type sum (0);
    for(int i=0; i<Size; i++)
        sum += array[i];
    return sum / Size;
}




//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

int main()
{
    bool printMatricesDetails = false;
    bool printLists = false;
    
    // main parameters for the test 
    int dimension = 30;
    int min = 1;
    int max = 100;
    
    long a = 999999;
    long b = 1000000;
    double delta = (double) a/b;
    double epsilon = 1.0 - delta;
    
    int maxcpt = 1000000;
    int d = 0;
    long blocksize = 30; // for BKZ insertions
    
    // iteration loop over matrices of same dimension
    const int maxIteration = 10;
    
    // print important information
    cout << "epsilon = " << epsilon << endl;
    cout << "dimension = " << dimension << endl;
    cout << "nombre de matrices testées = " << maxIteration << endl;
    cout << endl;
    
    // to display progress bar
    boost::progress_display show_progress(maxIteration);
    
    // arrays to store values

    double timing_PairRedPrimal [maxIteration];
    double timing_PairRedPrimalRandomized [maxIteration];

    double timing_LLL [maxIteration];
    double timing_LLL_PairRedPrimal [maxIteration];
    double timing_LLL_PostPairRedPrimal [maxIteration];
    double timing_LLL_PairRedPrimalRandomized [maxIteration];
    double timing_LLL_PostPairRedPrimalRandomized [maxIteration];

    double timing_LLLNTL [maxIteration];
    double timing_LLLNTL_PairRedPrimal [maxIteration];
    double timing_LLLNTL_PostPairRedPrimal [maxIteration];
    double timing_LLLNTL_PairRedPrimalRandomized [maxIteration];
    double timing_LLLNTL_PostPairRedPrimalRandomized [maxIteration];

    //double timing_LLL_NTL_Exact [maxIteration];
    
    double timing_BKZNTL [maxIteration];
    double timing_BKZNTL_PairRedPrimal [maxIteration];
    double timing_BKZNTL_PostPairRedPrimal [maxIteration];
    double timing_BKZNTL_PairRedPrimalRandomized [maxIteration];
    double timing_BKZNTL_PostPairRedPrimalRandomized [maxIteration];


    NScal length_Initial [maxIteration];
    NScal length_PairRedPrimal [maxIteration];
    NScal length_PairRedPrimalRandomized [maxIteration];
    
    NScal length_LLL [maxIteration];
    NScal length_LLL_PostPairRedPrimal [maxIteration];
    NScal length_LLL_PostPairRedPrimalRandomized [maxIteration];
    
    NScal length_LLLNTL [maxIteration];
    NScal length_LLLNTL_PostPairRedPrimal [maxIteration];
    NScal length_LLLNTL_PostPairRedPrimalRandomized [maxIteration];

    //NScal length_LLL_NTL_Exact [maxIteration];

    NScal length_BKZNTL [maxIteration];
    NScal length_BKZNTL_PostPairRedPrimal [maxIteration];
    NScal length_BKZNTL_PostPairRedPrimalRandomized [maxIteration];


    for (int iteration = 0; iteration < maxIteration; iteration++){
        ++show_progress;
    
        int seed = (iteration+1) * (iteration+1) * 123456789;
        //int seed = (int) (iteration+1) * 12345 * time(NULL);
        
        // We create copies of the same basis for: LLL, 
        // pairwisePrimalRed + LLL, LLL NTL floating point, 
        // LLL NTL exact version, pairwiseRed + LLL, BKZ NTL, 
        Basis basis_PairRedPrimal (dimension, L2NORM);
        ZZ det;
        RandomMatrix(basis_PairRedPrimal, det, min, max,seed);
        

        Basis basis_PairRedPrimalRandomized (basis_PairRedPrimal);
        Basis basis_LLL (basis_PairRedPrimal);
        Basis basis_PairRedPrimal_LLL (basis_PairRedPrimal);
        Basis basis_PairRedPrimalRandomized_LLL (basis_PairRedPrimal);
        Basis basis_LLLNTL (basis_PairRedPrimal);
        Basis basis_PairRedPrimal_LLLNTL (basis_PairRedPrimal);
        Basis basis_PairRedPrimalRandomized_LLLNTL (basis_PairRedPrimal);
        //Basis basis_LLLNTL_Exact (basis_PairRedPrimal);
        Basis basis_BKZNTL (basis_PairRedPrimal);
        Basis basis_PairRedPrimal_BKZNTL (basis_PairRedPrimal);
        Basis basis_PairRedPrimalRandomized_BKZNTL (basis_PairRedPrimal);
        
        IntLattice lattice_PairRedPrimal (basis_PairRedPrimal,1);
        IntLattice lattice_PairRedPrimalRandomized (basis_PairRedPrimalRandomized,1);
        IntLattice lattice_LLL (basis_LLL, 1);
        IntLattice lattice_PairRedPrimal_LLL (basis_PairRedPrimal_LLL, 1);
        IntLattice lattice_PairRedPrimalRandomized_LLL (basis_PairRedPrimalRandomized_LLL, 1);
        IntLattice lattice_LLLNTL (basis_LLLNTL,1);
        IntLattice lattice_PairRedPrimal_LLLNTL (basis_PairRedPrimal_LLLNTL,1);
        IntLattice lattice_PairRedPrimalRandomized_LLLNTL (basis_PairRedPrimalRandomized_LLLNTL,1);
        //IntLattice lattice_LLLNTL_Exact (basis_LLLNTL_Exact,1);
        IntLattice lattice_BKZNTL (basis_BKZNTL,1);
        IntLattice lattice_PairRedPrimal_BKZNTL (basis_PairRedPrimal_BKZNTL,1);
        IntLattice lattice_PairRedPrimalRandomized_BKZNTL (basis_PairRedPrimalRandomized_BKZNTL,1);

        lattice_PairRedPrimal.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal.sort(0);
        NScal initialShortestVectorLength = lattice_PairRedPrimal.getPrimalBasis().getVecNorm(1);
        length_Initial [iteration] = initialShortestVectorLength;
        
        if (printMatricesDetails) {
            cout << "\n*** Initial basis ***" << endl;
            cout << "det = " << det << endl;
            cout << "Shortest vector = " << initialShortestVectorLength << endl;
            lattice_PairRedPrimal.getPrimalBasis().write();
        }

        Reducer reducer_PairRedPrimal (lattice_PairRedPrimal);
        Reducer reducer_PairRedPrimalRandomized (lattice_PairRedPrimalRandomized);
        Reducer reducer_LLL (lattice_LLL);
        Reducer reducer_PairRedPrimal_LLL (lattice_PairRedPrimal_LLL);
        Reducer reducer_PairRedPrimalRandomized_LLL (lattice_PairRedPrimalRandomized_LLL);
        Reducer reducer_LLLNTL (lattice_LLLNTL);
        Reducer reducer_PairRedPrimal_LLLNTL (lattice_PairRedPrimal_LLLNTL);
        Reducer reducer_PairRedPrimalRandomized_LLLNTL (lattice_PairRedPrimalRandomized_LLLNTL);
        //Reducer reducer_LLLNTL_Exact (lattice_LLLNTL_Exact);
        Reducer reducer_BKZNTL (lattice_BKZNTL);
        Reducer reducer_PairRedPrimal_BKZNTL (lattice_PairRedPrimal_BKZNTL);
        Reducer reducer_PairRedPrimalRandomized_BKZNTL (lattice_PairRedPrimalRandomized_BKZNTL);
        

        //------------------------------------------------------------------------------------
        // Pairwise reduction in primal basis only
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimal = clock();
        reducer_PairRedPrimal.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimal = clock();

        lattice_PairRedPrimal.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal.sort(0);

        if (printMatricesDetails) {
            cout << "*** Pairwise reduction in primal basis only ***" << endl;
            cout << "Shortest vector = ";
            cout << lattice_PairRedPrimal.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimal.getPrimalBasis().write();
        }


        //------------------------------------------------------------------------------------
        // Randomized pairwise reduction in primal basis only
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimalRandomized = clock();
        reducer_PairRedPrimalRandomized.preRedDieterPrimalOnlyRandomized(d);
        clock_t end_PairRedPrimalRandomized = clock();

        lattice_PairRedPrimalRandomized.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized.sort(0);

        if (printMatricesDetails) {
            cout << "*** Randomized pairwise reduction in primal basis only ***" << endl;
            cout << "Shortest vector = ";
            cout << lattice_PairRedPrimalRandomized.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimalRandomized.getPrimalBasis().write();
        }


        //------------------------------------------------------------------------------------
        // LLL Richard
        //------------------------------------------------------------------------------------

        clock_t begin_LLL = clock();
        reducer_LLL.redLLL(delta, maxcpt, dimension);
        clock_t end_LLL = clock();
        
        lattice_LLL.getPrimalBasis().setNegativeNorm(true);
        lattice_LLL.getPrimalBasis().updateVecNorm();
        lattice_LLL.sort(0);
        
        if (printMatricesDetails) {
            cout << "*** LLL only ***" << endl;
            cout << "Shortest vector = " << lattice_LLL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_LLL.getPrimalBasis().write();
        }
        
        
        //------------------------------------------------------------------------------------
        // Pairwise reduction (in primal basis only) and then LLL Richard
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimal_LLL1 = clock();
        reducer_PairRedPrimal_LLL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimal_LLL1 = clock();
        
        lattice_PairRedPrimal_LLL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal_LLL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal_LLL.sort(0);
        
        NScal intermediateLength = lattice_PairRedPrimal_LLL.getPrimalBasis().getVecNorm(1);
        
        clock_t begin_PairRedPrimal_LLL2 = clock();
        reducer_PairRedPrimal_LLL.redLLL(delta, maxcpt, dimension);
        clock_t end_PairRedPrimal_LLL2 = clock();
        
        lattice_PairRedPrimal_LLL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal_LLL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal_LLL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Pairwise reduction in primal and LLL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimal_LLL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimal_LLL.getPrimalBasis().write();
        }


        //------------------------------------------------------------------------------------
        // Randomized pairwise reduction (in primal basis only) and then LLL Richard
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimalRandomized_LLL1 = clock();
        reducer_PairRedPrimalRandomized_LLL.preRedDieterPrimalOnlyRandomized(d);
        clock_t end_PairRedPrimalRandomized_LLL1 = clock();
        
        lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized_LLL.sort(0);
        
        NScal intermediateLengthRandomized = lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().getVecNorm(1);
        
        clock_t begin_PairRedPrimalRandomized_LLL2 = clock();
        reducer_PairRedPrimalRandomized_LLL.redLLL(delta, maxcpt, dimension);
        clock_t end_PairRedPrimalRandomized_LLL2 = clock();
        
        lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized_LLL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Randomized pairwise reduction in primal and LLL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().write();
        }
        
        
        //------------------------------------------------------------------------------------
        // LLL NTL reduction (floating point version = proxy)
        //------------------------------------------------------------------------------------

        clock_t begin_LLLNTL = clock();
        reducer_LLLNTL.redLLLNTLProxy(delta);
        clock_t end_LLLNTL = clock();
        
        lattice_LLLNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_LLLNTL.getPrimalBasis().updateVecNorm();
        lattice_LLLNTL.sort(0);
        
        if (printMatricesDetails) {
            cout << "*** LLL NTL Proxy only ***" << endl;
            cout << "Shortest vector = " << lattice_LLLNTL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_LLLNTL.getPrimalBasis().write();
        }
        

        //------------------------------------------------------------------------------------
        // Pairwise reduction (in primal basis only) and then LLL NTL proxy
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimal_LLLNTL1 = clock();
        reducer_PairRedPrimal_LLLNTL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimal_LLLNTL1 = clock();
        
        lattice_PairRedPrimal_LLLNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal_LLLNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal_LLLNTL.sort(0);
        
        // useful ?
        //NScal intermediateLengthBis = lattice_PairRedPrimal_LLL_NTL.getPrimalBasis().getVecNorm(1);
        
        clock_t begin_PairRedPrimal_LLLNTL2 = clock();
        reducer_PairRedPrimal_LLLNTL.redLLLNTLProxy(delta);
        clock_t end_PairRedPrimal_LLLNTL2 = clock();
        
        lattice_PairRedPrimal_LLLNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal_LLLNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal_LLLNTL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Pairwise reduction in primal and LLL NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimal_LLLNTL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimal_LLLNTL.getPrimalBasis().write();
        }


        //------------------------------------------------------------------------------------
        // Randomized pairwise reduction (in primal basis only) and then LLL NTL proxy
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimalRandomized_LLLNTL1 = clock();
        reducer_PairRedPrimalRandomized_LLLNTL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimalRandomized_LLLNTL1 = clock();
        
        lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized_LLLNTL.sort(0);
        
        // usefull ?
        //NScal intermediateShortestVectorLengthBis = lattice_PairRedPrimal_LLL_NTL.getPrimalBasis().getVecNorm(1);
        
        clock_t begin_PairRedPrimalRandomized_LLLNTL2 = clock();
        reducer_PairRedPrimalRandomized_LLLNTL.redLLLNTLProxy(delta);
        clock_t end_PairRedPrimalRandomized_LLLNTL2 = clock();
        
        lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized_LLLNTL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Randomized pairwise reduction in primal and LLL NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().write();
        }
        

        //------------------------------------------------------------------------------------
        // LLL NTL Exact reduction only
        //------------------------------------------------------------------------------------

        /*ZZ det2;
        clock_t begin_LLLNTL_Exact = clock();
        reducer_LLLNTL_Exact.redLLLNTLExact(det2, a, b);
        clock_t end_LLLNTL_Exact = clock();
        
        lattice_LLLNTL_Exact.getPrimalBasis().setNegativeNorm(true);
        lattice_LLLNTL_Exact.getPrimalBasis().updateVecNorm();
        lattice_LLLNTL_Exact.sort(0);
        
        if (printMatricesDetails) {
            cout << "*** LLL NTL Exact only ***" << endl;
            cout << "Shortest vector = " << lattice_LLLNTL_Exact.getPrimalBasis().getVecNorm(1) << endl;
            lattice_LLLNTL_Exact.getPrimalBasis().write();
        }*/


        //------------------------------------------------------------------------------------
        // BKZ NTL reduction
        //------------------------------------------------------------------------------------
        
        clock_t begin_BKZNTL = clock();
        reducer_BKZNTL.redBKZ(delta, blocksize);
        clock_t end_BKZNTL = clock();

        lattice_BKZNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_BKZNTL.getPrimalBasis().updateVecNorm();
        lattice_BKZNTL.sort(0);

        if (printMatricesDetails) {
            cout << "*** BKZ NTL only ***" << endl;
            cout << "Shortest vector = " << lattice_BKZNTL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_BKZNTL.getPrimalBasis().write();
        }


        //------------------------------------------------------------------------------------
        // Pairwise reduction (in primal basis only) and then BKZ NTL proxy
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimal_BKZNTL1 = clock();
        reducer_PairRedPrimal_BKZNTL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimal_BKZNTL1 = clock();
        
        lattice_PairRedPrimal_BKZNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal_BKZNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal_BKZNTL.sort(0);
        
        // useful ?
        //NScal intermediateLengthBis = lattice_PairRedPrimal_LLL_NTL.getPrimalBasis().getVecNorm(1);
        
        clock_t begin_PairRedPrimal_BKZNTL2 = clock();
        reducer_PairRedPrimal_BKZNTL.redLLLNTLProxy(delta);
        clock_t end_PairRedPrimal_BKZNTL2 = clock();
        
        lattice_PairRedPrimal_BKZNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimal_BKZNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimal_BKZNTL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Pairwise reduction in primal and BKZ NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimal_BKZNTL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimal_BKZNTL.getPrimalBasis().write();
        }


        //------------------------------------------------------------------------------------
        // Randomized pairwise reduction (in primal basis only) and then BKZ NTL proxy
        //------------------------------------------------------------------------------------

        clock_t begin_PairRedPrimalRandomized_BKZNTL1 = clock();
        reducer_PairRedPrimalRandomized_BKZNTL.preRedDieterPrimalOnly(d);
        clock_t end_PairRedPrimalRandomized_BKZNTL1 = clock();
        
        lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized_BKZNTL.sort(0);
        
        // usefull ?
        //NScal intermediateShortestVectorLengthBis = lattice_PairRedPrimal_LLL_NTL.getPrimalBasis().getVecNorm(1);
        
        clock_t begin_PairRedPrimalRandomized_BKZNTL2 = clock();
        reducer_PairRedPrimalRandomized_BKZNTL.redLLLNTLProxy(delta);
        clock_t end_PairRedPrimalRandomized_BKZNTL2 = clock();
        
        lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().setNegativeNorm(true);
        lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().updateVecNorm();
        lattice_PairRedPrimalRandomized_BKZNTL.sort(0);
        
        if (printMatricesDetails){
            cout << "*** Randomized pairwise reduction in primal and BKZ NTL ***" << endl;
            cout << "Shortest vector = " << lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().getVecNorm(1) << endl;
            lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().write();
        }

        //------------------------------------------------------------------------------------
        // timing updating
        //------------------------------------------------------------------------------------
        

        double runningTime_PairRedPrimal = double (end_PairRedPrimal - begin_PairRedPrimal) / CLOCKS_PER_SEC;
        double runningTime_PairRedPrimalRandomized = double (end_PairRedPrimalRandomized - begin_PairRedPrimalRandomized) / CLOCKS_PER_SEC;

        double runningTime_LLL = double(end_LLL - begin_LLL) / CLOCKS_PER_SEC;
        double runningTime_LLL_PairRedPrimal = double(end_PairRedPrimal_LLL1 - begin_PairRedPrimal_LLL1) / CLOCKS_PER_SEC;
        double runningTime_LLL_PostPairRedPrimal = double(end_PairRedPrimal_LLL2 - begin_PairRedPrimal_LLL2) / CLOCKS_PER_SEC;
        double runningTime_LLL_PairRedPrimalRandomized = double(end_PairRedPrimalRandomized_LLL1 - begin_PairRedPrimalRandomized_LLL1) / CLOCKS_PER_SEC;
        double runningTime_LLL_PostPairRedPrimalRandomized = double (end_PairRedPrimalRandomized_LLL2 - begin_PairRedPrimalRandomized_LLL2) / CLOCKS_PER_SEC;

        double runningTime_LLLNTL = double(end_LLLNTL - begin_LLLNTL) / CLOCKS_PER_SEC;
        double runningTime_LLLNTL_PairRedPrimal = double(end_PairRedPrimal_LLLNTL1 - begin_PairRedPrimal_LLLNTL1) / CLOCKS_PER_SEC;
        double runningTime_LLLNTL_PostPairRedPrimal = double (end_PairRedPrimal_LLLNTL2 - begin_PairRedPrimal_LLLNTL2) / CLOCKS_PER_SEC;
        double runningTime_LLLNTL_PairRedPrimalRandomized = double(end_PairRedPrimalRandomized_LLLNTL1 - begin_PairRedPrimalRandomized_LLLNTL1) / CLOCKS_PER_SEC;
        double runningTime_LLLNTL_PostPairRedPrimalRandomized = double (end_PairRedPrimalRandomized_LLLNTL2 - begin_PairRedPrimalRandomized_LLLNTL2) / CLOCKS_PER_SEC;
        //double runningTime_LLLNTL_Exact = double(end_LLLNTL_Exact - begin_LLLNTL_Exact) / CLOCKS_PER_SEC;
        double runningTime_BKZNTL = double (end_BKZNTL - begin_BKZNTL) / CLOCKS_PER_SEC;
        double runningTime_BKZNTL_PairRedPrimal = double(end_PairRedPrimal_BKZNTL1 - begin_PairRedPrimal_BKZNTL1) / CLOCKS_PER_SEC;
        double runningTime_BKZNTL_PostPairRedPrimal = double (end_PairRedPrimal_BKZNTL2 - begin_PairRedPrimal_BKZNTL2) / CLOCKS_PER_SEC;
        double runningTime_BKZNTL_PairRedPrimalRandomized = double(end_PairRedPrimalRandomized_BKZNTL1 - begin_PairRedPrimalRandomized_BKZNTL1) / CLOCKS_PER_SEC;
        double runningTime_BKZNTL_PostPairRedPrimalRandomized = double (end_PairRedPrimalRandomized_BKZNTL2 - begin_PairRedPrimalRandomized_BKZNTL2) / CLOCKS_PER_SEC;

        //------------------------------------------------------------------------------------
        // timing and length arrays updating
        //------------------------------------------------------------------------------------

        timing_PairRedPrimal [iteration] = runningTime_PairRedPrimal;
        timing_PairRedPrimalRandomized [iteration] = runningTime_PairRedPrimalRandomized;

        timing_LLL [iteration] = runningTime_LLL;
        timing_LLL_PairRedPrimal [iteration] = runningTime_LLL_PairRedPrimal;
        timing_LLL_PostPairRedPrimal [iteration] = runningTime_LLL_PostPairRedPrimal;
        timing_LLL_PairRedPrimalRandomized [iteration] = runningTime_LLL_PairRedPrimalRandomized;
        timing_LLL_PostPairRedPrimalRandomized [iteration] = runningTime_LLL_PostPairRedPrimalRandomized;

        timing_LLLNTL [iteration] = runningTime_LLLNTL;
        timing_LLLNTL_PairRedPrimal [iteration] = runningTime_LLLNTL_PairRedPrimal;
        timing_LLLNTL_PostPairRedPrimal [iteration] = runningTime_LLLNTL_PostPairRedPrimal;
        timing_LLLNTL_PairRedPrimalRandomized [iteration] = runningTime_LLLNTL_PairRedPrimalRandomized;
        timing_LLLNTL_PostPairRedPrimalRandomized [iteration] = runningTime_LLLNTL_PostPairRedPrimalRandomized;

        //timing_LLLNTL_Exact [iteration] = runningTime_LLLNTL_Exact;

        timing_BKZNTL [iteration] = runningTime_BKZNTL;
        timing_BKZNTL_PairRedPrimal [iteration] = runningTime_BKZNTL_PairRedPrimal;
        timing_BKZNTL_PostPairRedPrimal [iteration] = runningTime_BKZNTL_PostPairRedPrimal;
        timing_BKZNTL_PairRedPrimalRandomized [iteration] = runningTime_BKZNTL_PairRedPrimalRandomized;
        timing_BKZNTL_PostPairRedPrimalRandomized [iteration] = runningTime_BKZNTL_PostPairRedPrimalRandomized;


        length_PairRedPrimal [iteration] = lattice_PairRedPrimal.getPrimalBasis().getVecNorm(1);
        length_PairRedPrimalRandomized [iteration] = lattice_PairRedPrimalRandomized.getPrimalBasis().getVecNorm(1);

        length_LLL [iteration] = lattice_LLL.getPrimalBasis().getVecNorm(1);
        length_LLL_PostPairRedPrimal [iteration] = lattice_PairRedPrimal_LLL.getPrimalBasis().getVecNorm(1);
        length_LLL_PostPairRedPrimalRandomized [iteration] = lattice_PairRedPrimalRandomized_LLL.getPrimalBasis().getVecNorm(1);

        length_LLLNTL [iteration] = lattice_LLLNTL.getPrimalBasis().getVecNorm(1);
        length_LLLNTL_PostPairRedPrimal [iteration] = lattice_PairRedPrimal_LLLNTL.getPrimalBasis().getVecNorm(1);
        length_LLLNTL_PostPairRedPrimalRandomized [iteration] = lattice_PairRedPrimalRandomized_LLLNTL.getPrimalBasis().getVecNorm(1);

        //length_LLLNTL_Exact [iteration] = lattice_LLLNTL_Exact.getPrimalBasis().getVecNorm(1);

        length_BKZNTL [iteration] = lattice_BKZNTL.getPrimalBasis().getVecNorm(1);
        length_BKZNTL_PostPairRedPrimal [iteration] = lattice_PairRedPrimal_BKZNTL.getPrimalBasis().getVecNorm(1);
        length_BKZNTL_PostPairRedPrimalRandomized [iteration] = lattice_PairRedPrimalRandomized_BKZNTL.getPrimalBasis().getVecNorm(1);

    } // end iteration loop over matrices of same dimension
    


    // print arrays
    
    // mettre à jour 

    /*
    if (printLists) {

        cout << "\nTIMING LIST ---------" << endl;
        print("          PairRedPrimal", timing_PairRedPrimal, false);
        print("PairRedPrimalRandomized", timing_PairRedPrimalRandomized, false);
        print("                    LLL", timing_LLL, false);
        print("          PairRedPrimal", timing_PairRedPrimal, false);
        print("      PostPairRedPrimal", timing_LLL_PostPairRedPrimal, false);
        print("          LLL_NTL_Proxy", timing_LLL_NTL_Proxy, false);
        print("      PairRedPrimal_NTL", timing_PairRedPrimal_NTL, false);
        print("  PostPairRedPrimal_NTL", timing_LLL_NTL_PostPairRedPrimal, false);
        //print("          LLL_NTL_Exact", timing_LLL_NTL_Exact, false);
        print("                BKZ_NTL", timing_BKZ_NTL, false);

        cout << "\nLENGTH LIST ---------" << endl;
        print("Initial",length_Initial,true);
        print("PairRedPrimalRandomized", length_PairRedPrimalRandomized, true);
        print("                    LLL", length_LLL, true);
        print("          PairRedPrimal", length_PairRedPrimal, true);
        print("      PostPairRedPrimal", length_LLL_PostPairRedPrimal, true);
        print("          LLL_NTL_Proxy", length_LLL_NTL_Proxy, true);
        print("      PairRedPrimal_NTL", length_PairRedPrimal_NTL, true);
        print("  PostPairRedPrimal_NTL", length_LLL_NTL_PostPairRedPrimal, true);
        //print("          LLL_NTL_Exact", length_LLL_NTL_Exact, true);
        print("                BKZ_NTL", length_BKZ_NTL, true);
    }
    */

    
    //------------------------------------------------------------------------------------
    // Results printing in console
    //------------------------------------------------------------------------------------

    // print parameters used
    cout << "\n" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "dimension = " << dimension << endl;
    cout << "nombre de matrices testées = " << maxIteration << endl;
    
    // print statistics
    

    cout << "\n---------------- TIMING AVG ----------------\n" << endl;

    cout << "       PairRedPrimal = " << Average(timing_PairRedPrimal) << endl;
    cout << " PairRedPrimalRandom = " << Average(timing_PairRedPrimalRandomized) << endl;
    cout << endl;

    cout << "                 LLL = " << Average(timing_LLL) << endl;
    cout << "         PairRed+LLL = " << Average(timing_LLL_PairRedPrimal) + Average(timing_LLL_PostPairRedPrimal);
    cout << " (" << Average(timing_LLL_PostPairRedPrimal) << ")" << endl;
    cout << "   PairRedRandom+LLL = " << Average(timing_LLL_PairRedPrimalRandomized) + Average(timing_LLL_PostPairRedPrimalRandomized);
    cout << " (" << Average(timing_LLL_PostPairRedPrimalRandomized) << ")" << endl;
    cout << endl;
    
    cout << "              LLLNTL = " << Average(timing_LLLNTL) << endl;
    cout << "      PairRed+LLLNTL = " << Average(timing_LLLNTL_PairRedPrimal) + Average(timing_LLLNTL_PostPairRedPrimal);
    cout << " (" << Average(timing_LLLNTL_PostPairRedPrimal) << ")" << endl;
    cout << "PairRedRandom+LLLNTL = " << Average(timing_LLLNTL_PairRedPrimalRandomized) + Average(timing_LLLNTL_PostPairRedPrimalRandomized);
    cout << " (" << Average(timing_LLLNTL_PostPairRedPrimalRandomized) << ")" << endl;
    cout << endl;

    //cout << "        LLLNTL_Exact = " << Average(timing_LLL_NTL_Exact) << endl;
    cout << endl;

    cout << "              BKZNTL = " << Average(timing_BKZNTL) << endl;
    cout << "      PairRed+BKZNTL = " << Average(timing_BKZNTL_PairRedPrimal) + Average(timing_BKZNTL_PostPairRedPrimal);
    cout << " (" << Average(timing_BKZNTL_PostPairRedPrimal) << ")" << endl;
    cout << "PairRedRandom+BKZNTL = " << Average(timing_BKZNTL_PairRedPrimalRandomized) + Average(timing_BKZNTL_PostPairRedPrimalRandomized);
    cout << " (" << Average(timing_BKZNTL_PostPairRedPrimalRandomized) << ")" << endl;
    cout << endl;

    cout << "\n--------------------------------------------" << endl;
    


    cout << "\n---------------- LENGTH AVG ----------------\n" << endl;
    
    cout << "             Initial = " << conv<ZZ>(Average(length_Initial)) << endl;
    cout << "       PairRedPrimal = " << conv<ZZ>(Average(length_PairRedPrimal)) << endl;
    cout << " PairRedPrimalRandom = " << conv<ZZ>(Average(length_PairRedPrimalRandomized)) << endl;    
    cout << endl;
    
    cout << "                 LLL = " << conv<ZZ>(Average(length_LLL)) << endl;
    cout << "         PairRed+LLL = " << conv<ZZ>(Average(length_LLL_PostPairRedPrimal)) << endl;
    cout << "   PairRedRandom+LLL = " << conv<ZZ>(Average(length_LLL_PostPairRedPrimalRandomized)) << endl;
    cout << endl;
    
    cout << "              LLLNTL = " << conv<ZZ>(Average(length_LLLNTL)) << endl;
    cout << "      PairRed+LLLNTL = " << conv<ZZ>(Average(length_LLLNTL_PostPairRedPrimal)) << endl;
    cout << "PairRedRandom+LLLNTL = " << conv<ZZ>(Average(length_LLLNTL_PostPairRedPrimalRandomized)) << endl;
    cout << endl;

    //cout << "       LLL_NTL_Exact = " << conv<ZZ>(Average(length_LLL_NTL_Exact)) << endl;
    cout << endl;

    cout << "              BKZNTL = " << conv<ZZ>(Average(length_BKZNTL)) << endl;
    cout << "      PairRed+BKZNTL = " << conv<ZZ>(Average(length_BKZNTL_PostPairRedPrimal)) << endl;
    cout << "PairRedRandom+BKZNTL = " << conv<ZZ>(Average(length_BKZNTL_PostPairRedPrimalRandomized)) << endl;
    cout << endl;

    cout << "\n--------------------------------------------" << endl;
    
    return 0;
}





