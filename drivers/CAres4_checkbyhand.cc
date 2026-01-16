#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"

#include "SRothman/EECs/src/usings.h"
#include "SRothman/EECs/src/CARes4Calculator.h"

#include <stdio.h>

#include "testResCalculator.h"

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

    jetFactory.initialize();
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

    printf("total weight gen:\n");
    printf("\tchain: %g\n", gen.total_chain_weight());
    printf("\tsymmetric_wrtR : %g\n", gen.total_symmetric_wrtR_weight());
    printf("\tsymmetric_wrtr : %g\n", gen.total_symmetric_wrtr_weight());
    printf("total weight gen_unmatched\n");
    printf("\tchain: %g\n", unmatched_gen.total_chain_weight());
    printf("\tsymmetric_wrtR : %g\n", unmatched_gen.total_symmetric_wrtR_weight());
    printf("\tsymmetric_wrtr : %g\n", unmatched_gen.total_symmetric_wrtr_weight());
    printf("total weight transfer\n");
    printf("\tchain -> chain: %g\n", transfer.total_chain_weight_reco());
    printf("\tchain -> symmetric_wrtR: %g\n", transfer.total_chain_to_symmetric_wrtR_weight_reco());
    printf("\tchain -> symmetric_wrtr: %g\n", transfer.total_chain_to_symmetric_wrtr_weight_reco());
    printf("\tsymmetric_wrtR -> chain: %g\n", transfer.total_symmetric_wrtR_to_chain_weight_reco());
    printf("\tsymmetric_wrtR -> symmetric_wrtR: %g\n", transfer.total_symmetric_wrtR_weight_reco());
    printf("\tsymmetric_wrtr -> symmetric_wrtr: %g\n", transfer.total_symmetric_wrtr_weight_reco());
}
