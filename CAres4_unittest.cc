#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/EECs/src/usings.h"
#include "SRothman/EECs/src/CARes4Calculator.cc"

#include <stdio.h>

#include "testCalculator.h"

int main(){
    JetFactory jetFactory;

    std::vector<double> bins(NBINS-1);
    for (int i=0; i<NBINS-1; ++i){
        bins[i] = double(i)/(NBINS-2);
    }

    EEC::CARes4Calculator calculator(
        bins,
        bins, bins,
        bins, bins,
        EEC::normType::RAWPT
    );

    printf("UNIT TEST CARes4Calculator...\n\n");
    printf("NPART = %d\n", NPART);
    printf("NBINS = %d\n", NBINS);
    printf("\n");

    printf("PLAIN\n");
    checkAllSame<EEC::CARes4Result>(calculator, jetFactory);
    checkSum<EEC::CARes4Result_MultiArray>(calculator, jetFactory);

    printf("WITH UNMATCHED\n");
    checkAllSame_matched<EEC::CARes4Result>(calculator, jetFactory);
    checkSum_matched<EEC::CARes4Result_MultiArray>(calculator, jetFactory);

    EEC::CARes4TransferCalculator tcalculator(
        bins, 
        bins, bins,
        bins, bins,
        bins, 
        bins, bins,
        bins, bins,
        EEC::normType::RAWPT
    );

    printf("WITH TRANSFER\n");
    checkAllSame_transfer<EEC::CARes4Result, 
                          EEC::CARes4TransferResult, 
                          false>(
            tcalculator, jetFactory
    );
    check_sum_transfer<EEC::CARes4Result_MultiArray, 
                       EEC::CARes4TransferResult_Vector, 
                       false>(
            calculator, tcalculator, jetFactory
    );

    return 0;
}
