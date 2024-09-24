#ifndef TEST_RES3_H
#define TEST_RES3_H

#include "example_jet.h"
#include <boost/histogram.hpp>
#include "SRothman/EECs/src/theOnlyHeader.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"
#include <memory>

bool test_res3(unsigned Ntest, unsigned Npart) noexcept {
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);

    std::vector<double> xiedges({1e-3, 0.5, 1.01});
    auto xi_ax = std::make_shared<boost::histogram::axis::variable<double>>(xiedges);

    std::vector<double> phi_edges({1e-3, 0.5, 1.0});
    auto phi_ax = std::make_shared<boost::histogram::axis::variable<double>>(phi_edges);
    
    fastEEC::normType norm = fastEEC::normType::SUMPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        fastEEC::result_t<double> EEC_gen;

        runFastEEC(
            EEC_gen,
            genJet, RLax, norm,
            3, fastEEC::DORES3,
            RLax_coarse, xi_ax, phi_ax
        );

        double acc_res3 = recursive_reduce(
            *EEC_gen.resolved3, 0.0
        );
        passed = passed && std::abs(acc_res3 - 1) < 1e-10;
    }
    return passed;
}

bool test_res3_PU(unsigned Ntest, unsigned Npart) noexcept {
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);

    std::vector<double> xiedges({1e-3, 0.5, 1.01});
    auto xi_ax = std::make_shared<boost::histogram::axis::variable<double>>(xiedges);

    std::vector<double> phi_edges({1e-3, 0.5, 1.0});
    auto phi_ax = std::make_shared<boost::histogram::axis::variable<double>>(phi_edges);
    
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
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);

    std::vector<double> xiedges({1e-3, 0.5, 1.01});
    auto xi_ax = std::make_shared<boost::histogram::axis::variable<double>>(xiedges);

    std::vector<double> phi_edges({1e-3, 0.5, 1.0});
    auto phi_ax = std::make_shared<boost::histogram::axis::variable<double>>(phi_edges);
    
    fastEEC::normType norm = fastEEC::normType::RAWPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        jet recoJet;
        make_random_jet(recoJet, Npart);

        arma::mat ptrans = arma::eye<arma::mat>(Npart, Npart);

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
