#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/EECs/src/usings.h"

#include "SRothman/EECs/src/ProjCalculator.h"

#include <stdio.h>

#include "testProjCalculator.h"

int main(){
    JetFactory jetFactory;

    std::vector<double> bins(NBINS-1);
    for (int i=0; i<NBINS-1; ++i){
        bins[i] = double(i)/(NBINS-2);
    }

    EEC::ProjCalculator calculator(
        bins,
        EEC::normType::RAWPT
    );

    printf("UNIT TEST ProjCalculator...\n\n");
    printf("NPART = %d\n", NPART);
    printf("NBINS = %d\n", NBINS);
    printf("\n");

    printf("PLAIN\n");
    checkSum(calculator, jetFactory);
    checkAllSame(calculator, jetFactory);

    printf("WITH UNMATCHED\n");
    checkSum_matched(calculator, jetFactory);
    checkAllSame_matched(calculator, jetFactory);

    EEC::ProjTransferCalculator tcalculator(
        bins, bins, 
        EEC::normType::RAWPT
    );

    printf("WITH TRANSFER\n");
    checkSum_transfer(calculator, tcalculator, jetFactory);
    checkAllSame_transfer(tcalculator, jetFactory);

    return 0;
}
