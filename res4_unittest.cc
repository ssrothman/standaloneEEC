#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/EECs/src/usings.h"
#include "SRothman/EECs/src/Res4Calculator.cc"

#include <stdio.h>

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

    printf("UNIT TEST Res4Calculator...\n\n");
    printf("NPART = %d\n", NPART);
    printf("NBINS = %d\n", NBINS);
    printf("\n");

    printf("PLAIN\n");
    checkAllSame<EEC::Res4Result>(calculator, jetFactory);
    checkSum<EEC::Res4Result_MultiArray>(calculator, jetFactory);

    printf("WITH UNMATCHED\n");
    checkAllSame_matched<EEC::Res4Result>(calculator, jetFactory);
    checkSum_matched<EEC::Res4Result_MultiArray>(calculator, jetFactory);

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

    printf("WITH TRANSFER\n");
    checkAllSame_transfer<EEC::Res4Result, 
                          EEC::Res4TransferResult, 
                          true>(
            tcalculator, jetFactory
    );
    check_sum_transfer<EEC::Res4Result_MultiArray, 
                       EEC::Res4TransferResult_Vector, 
                       true>(
            calculator, tcalculator, jetFactory
    );

    return 0;
}
