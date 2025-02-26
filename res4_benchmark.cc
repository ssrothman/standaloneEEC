#include <stdio.h>
#include <chrono>
#include <random>
#include <boost/histogram.hpp>

#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/EECs/src/Res4Calculator.h"
#include "SRothman/EECs/src/Res4Result.h"
#include "SRothman/EECs/src/Res4TransferResult.h"

#include "testCalculator.h"

int main(){
    JetFactory jetFactory;

    std::vector<double> bins(NBINS-1);
    for (int i=0; i<NBINS-1; ++i){
        bins[i] = double(i)/(NBINS-2);
    }

    EEC::Res4Calculator calculator(
        bins,
        bins, bins,
        bins, bins,
        bins, bins,
        0.05, 0.05,
        EEC::normType::RAWPT
    );

    printf("Benchmarking Res4Calculator...\n\n");
    printf("NPART = %d\n", NPART);
    printf("NITER = %d\n", NITER);
    printf("NBINS = %d\n", NBINS);
    printf("\n");

    printf("PLAIN\n");
    unsigned long t1 = benchmarkCalculator<EEC::Res4Result_Vector, EEC::Res4Calculator, false>(
            calculator, jetFactory);
    printf("\tJIT vector: %lu ms\n", t1);

    unsigned long t2 = benchmarkCalculator<EEC::Res4Result_MultiArray, EEC::Res4Calculator, false>(
            calculator, jetFactory);
    printf("\tJIT multiarray: %lu ms\n", t2);

    unsigned long t3 = benchmarkCalculator<EEC::Res4Result_Vector, EEC::Res4Calculator, true>(
            calculator, jetFactory);
    printf("\tprecomputed vector: %lu ms\n", t3);

    unsigned long t4 = benchmarkCalculator<EEC::Res4Result_MultiArray, EEC::Res4Calculator, true>(
            calculator, jetFactory);
    printf("\tprecomputed multiarray: %lu ms\n", t4);

    printf("WITH UNMATCHED\n");
    unsigned long t5 = benchmarkCalculator_unmatched<EEC::Res4Result_Vector, EEC::Res4Calculator, false>(
            calculator, jetFactory);
    printf("\tJIT vector: %lu ms\n", t5);

    unsigned long t6 = benchmarkCalculator_unmatched<EEC::Res4Result_MultiArray, EEC::Res4Calculator, false>(
            calculator, jetFactory);
    printf("\tJIT multiarray: %lu ms\n", t6);

    unsigned long t7 = benchmarkCalculator_unmatched<EEC::Res4Result_Vector, EEC::Res4Calculator, true>(
            calculator, jetFactory);
    printf("\tprecomputed vector: %lu ms\n", t7);

    unsigned long t8 = benchmarkCalculator_unmatched<EEC::Res4Result_MultiArray, EEC::Res4Calculator, true>(
            calculator, jetFactory);
    printf("\tprecomputed multiarray: %lu ms\n", t8);

    EEC::Res4TransferCalculator tcalculator(
        bins, 
        bins, bins,
        bins, bins,
        bins, bins,
        bins, 
        bins, bins,
        bins, bins,
        bins, bins,
        0.05, 0.05,
        EEC::normType::RAWPT
    );
    
    printf("TRANSFER\n");
    unsigned long t9 = benchmarkCalculator_transfer<EEC::Res4Result_Vector, 
                                                    EEC::Res4TransferResult_Vector,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    false>(
            tcalculator, jetFactory
    );
    printf("\tJIT vector - vector: %lu ms\n", t9);

    unsigned long t10 = benchmarkCalculator_transfer<EEC::Res4Result_Vector, 
                                                    EEC::Res4TransferResult_MultiArray,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    false>(
            tcalculator, jetFactory
    );
    printf("\tJIT vector - multiarray: %lu ms\n", t10);

    unsigned long t11 = benchmarkCalculator_transfer<EEC::Res4Result_MultiArray, 
                                                    EEC::Res4TransferResult_Vector,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    false>(
            tcalculator, jetFactory
    );
    printf("\tJIT multiarray - vector: %lu ms\n", t11);

    unsigned long t12 = benchmarkCalculator_transfer<EEC::Res4Result_MultiArray, 
                                                    EEC::Res4TransferResult_MultiArray,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    false>(
            tcalculator, jetFactory
    );
    printf("\tJIT multiarray - multiarray: %lu ms\n", t12);

    unsigned long t13 = benchmarkCalculator_transfer<EEC::Res4Result_Vector, 
                                                    EEC::Res4TransferResult_Vector,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    true>(
            tcalculator, jetFactory
    );
    printf("\tprecomputed vector - vector: %lu ms\n", t13);

    unsigned long t14 = benchmarkCalculator_transfer<EEC::Res4Result_Vector, 
                                                    EEC::Res4TransferResult_MultiArray,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    true>(
            tcalculator, jetFactory
    );
    printf("\tprecomputed vector - multiarray: %lu ms\n", t14);

    unsigned long t15 = benchmarkCalculator_transfer<EEC::Res4Result_MultiArray, 
                                                    EEC::Res4TransferResult_Vector,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    true>(
            tcalculator, jetFactory
    );
    printf("\tprecomputed multiarray - vector: %lu ms\n", t15);

    unsigned long t16 = benchmarkCalculator_transfer<EEC::Res4Result_MultiArray, 
                                                    EEC::Res4TransferResult_MultiArray,
                                                    EEC::Res4TransferCalculator, 
                                                    true,
                                                    true>(
            tcalculator, jetFactory
    );
    printf("\tprecomputed multiarray - multiarray: %lu ms\n", t16);

    return 0;
}
