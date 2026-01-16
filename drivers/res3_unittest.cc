#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/EECs/src/usings.h"
#include "SRothman/EECs/src/Res3Calculator.h"

#include <stdio.h>

#include "testResCalculator.h"

int main(){
    JetFactory jetFactory;

    std::vector<double> bins(NBINS-1);
    for (int i=0; i<NBINS-1; ++i){
        bins[i] = double(i)/(NBINS-2);
    }

    EEC::Res3Calculator calculator(
        bins,
        bins, bins,
        EEC::normType::RAWPT
    );

    printf("UNIT TEST Res3Calculator...\n\n");
    printf("NPART = %d\n", NPART);
    printf("NBINS = %d\n", NBINS);
    printf("\n");

    printf("PLAIN\n");
    checkAllSame<EEC::Res3Result>(calculator, jetFactory);
    checkSum<EEC::Res3Result_MultiArray>(calculator, jetFactory);

    printf("WITH UNMATCHED\n");
    checkAllSame_matched<EEC::Res3Result>(calculator, jetFactory);
    checkSum_matched<EEC::Res3Result_MultiArray>(calculator, jetFactory);

    EEC::Res3TransferCalculator tcalculator(
        bins, bins, bins,
        bins, bins, bins,
        EEC::normType::RAWPT
    );

    printf("WITH TRANSFER\n");
    checkAllSame_transfer<EEC::Res3Result, 
                          EEC::Res3TransferResult, 
                          false>(
            tcalculator, jetFactory
    );
    check_sum_transfer<EEC::Res3Result_MultiArray, 
                       EEC::Res3TransferResult_Vector, 
                       false>(
            calculator, tcalculator, jetFactory
    );

    return 0;
}
