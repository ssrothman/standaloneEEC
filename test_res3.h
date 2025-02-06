#ifndef TEST_RES3_H
#define TEST_RES3_H

#include "example_jet.h"
#include <boost/histogram.hpp>
#include "SRothman/EECs/src/theOnlyHeader.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"
#include "SRothman/SimonTools/src/ToyShowerer.h"
#include "axes.h"
#include <memory>

bool test_res3(unsigned Ntest, unsigned Npart) noexcept {
    auto RLax = get_RLax();
    auto RLax_coarse = get_RLax_coarse();
    auto xi_ax = get_xi_res3_ax();
    auto phi_ax = get_phi_res3_ax();

    fastEEC::normType norm = fastEEC::normType::SUMPT;

    fastEEC::result_t<double> accu;

    bool passed = true;

    ToyShowerer showerer("UNIFORM", "GLUON", "LNX", "OFF",
                         0.1, 0.1, 0.5, "");

    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        showerer.shower(100, 0, 0, 0, Npart, genJet);

        double acc_res3;

        if (REPEAT==0){
            runFastEEC(
                accu,
                genJet, RLax, norm,
                3, fastEEC::DORES3,
                RLax_coarse, xi_ax, phi_ax
            );

            acc_res3 = recursive_reduce(
                *accu.resolved3, 0.0
            );
        } else {
            fastEEC::result_t<double> EEC_gen;

            runFastEEC(
                EEC_gen,
                genJet, RLax, norm,
                3, fastEEC::DORES3,
                RLax_coarse, xi_ax, phi_ax
            );

            acc_res3 = recursive_reduce(
                *EEC_gen.resolved3, 0.0
            );
        }
        passed = passed && std::abs(acc_res3 - 1) < 1e-10;
    }
    return passed;
}

bool test_res3_PU(unsigned Ntest, unsigned Npart) noexcept {
    auto RLax = get_RLax();
    auto RLax_coarse = get_RLax_coarse();
    auto xi_ax = get_xi_res3_ax();
    auto phi_ax = get_phi_res3_ax();

    fastEEC::normType norm = fastEEC::normType::SUMPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        std::vector<bool> isPU(Npart, false);
        isPU[0] = true;

        fastEEC::result_t<double> EEC_gen;

        runFastEEC(
            EEC_gen,
            genJet, RLax, norm,
            3, fastEEC::DORES3 | fastEEC::DOPU,
            RLax_coarse, xi_ax, phi_ax,
            nullptr, nullptr, 
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            0.0,
            &isPU
        );

        double acc_res3 = recursive_reduce(
            *EEC_gen.resolved3, 0.0
        );
        passed = passed && std::abs(acc_res3 - 1) < 1e-10;

        double acc_res3_PU = recursive_reduce(
            *EEC_gen.resolved3_PU, 0.0
        );
        passed = passed && acc_res3_PU > 0.0 && acc_res3_PU < 1.0;
    }

    return passed;
}

bool test_res3_transfer(unsigned Ntest, unsigned Npart) noexcept {
    auto RLax = get_RLax();
    auto RLax_coarse = get_RLax_coarse();
    auto xi_ax = get_xi_res3_ax();
    auto phi_ax = get_phi_res3_ax();
    
    fastEEC::normType norm = fastEEC::normType::RAWPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        jet recoJet;
        make_random_jet(recoJet, Npart);

        Eigen::MatrixXd ptrans = Eigen::MatrixXd::Identity(Npart,Npart);

        /*for(unsigned iPart=0; iPart < Npart; ++iPart){
            ptrans(iPart, iPart) = recoJet.particles[iPart].pt / genJet.particles[iPart].pt;
        }*/
        //ptrans *= genJet.pt / recoJet.pt;

        std::vector<bool> PU(Npart, false);
        std::vector<bool> UM(Npart, false);

        PU[0] = true;
        UM[0] = true;
        ptrans(0,0) = 0.0;

        fastEEC::result_t<double> EEC_gen;

        runFastEEC(
            EEC_gen,
            genJet, RLax, norm,
            3, fastEEC::DORES3 | fastEEC::DOPU | fastEEC::DOTRANSFER,
            RLax_coarse, xi_ax, phi_ax,
            nullptr, nullptr, 
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            0.0,
            &UM, &recoJet, &ptrans
        );

        fastEEC::result_t<double> EEC_reco;

        runFastEEC(
            EEC_reco,
            recoJet, RLax, norm,
            3, fastEEC::DOPU | fastEEC::DORES3,
            RLax_coarse, xi_ax, phi_ax,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            0.0,
            &PU
        );

        double acc_res3 = recursive_reduce(
            *EEC_reco.resolved3, 0.0
        );

        double acc_res3_PU = recursive_reduce(
            *EEC_reco.resolved3_PU, 0.0
        );

        double acc_res3_transfer = recursive_reduce(
            *EEC_gen.transfer_res3, 0.0
        );

        passed = passed && std::abs(acc_res3-acc_res3_PU-acc_res3_transfer) < 1e-8;
    }

    return passed;
}

#endif
