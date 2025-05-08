#include "SRothman/SimonTools/src/ToyShowerer.h"
#include "SRothman/EECs/src/Res4Calculator.h"

#include <iostream>

int main(){
    simon::ToyShowerer showerer("COS2PHI", "GLUON", "LNX", 0.01, 0.01, 0.4, true);

    simon::jet thejet;

    for (unsigned REP=0; REP<10; ++REP){
        showerer.shower(200, 2, 0, 40, thejet);

        std::vector<double> bins = {{0.0, 1.0}};

        EEC::Res4Calculator calculator(
            bins, 
            bins, bins,
            bins, bins,
            bins, bins,
            0.05, 0.05,
            EEC::normType::RAWPT
        );

        EEC::Res4Result_Unbinned result;
        calculator.compute_precomputed(thejet, result);

        for (const auto& entry : result.get_tee().get_data()){
            printf("%f, %f, %f, %g\n", entry.iR, entry.ir, entry.ic, entry.wt);
        }
    }
}
