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

    EEC::CARes4TransferCalculator tcalculator(
        bins, 
        bins, bins,
        bins, bins,
        bins, 
        bins, bins,
        bins, bins,
        EEC::normType::RAWPT
    );

    printf("Check by hand.....\n");
    printf("NPART = %d\n", NPART);
    printf("NBINS = %d\n", NBINS);
    printf("\n");

    printf("bins: ");
    for (int i=0; i<NBINS-1; ++i){
        printf("%0.2f, ", bins[i]);
    }
    printf("\n");
    printf("\n");

    simon::jet J_gen;
    jetFactory.makeJet(J_gen, NPART);

    simon::jet J_reco;
    Eigen::MatrixXd tmat;
    jetFactory.makeTransferJet(J_gen, J_reco, tmat);

    EEC::CARes4Result_MultiArray gen(
        calculator
    );
    EEC::CARes4Result_MultiArray unmatched_gen(
        calculator
    );
    EEC::CARes4TransferResult_Vector transfer(
        tcalculator
    );

    tcalculator.compute_precomputed(
            J_gen, J_reco, tmat,
            gen, unmatched_gen, transfer);

}
