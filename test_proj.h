#ifndef TEST_PROJ_H
#define TEST_PROJ_H

#include "example_jet.h"
#include <boost/histogram.hpp>
#include "SRothman/EECs/src/theOnlyHeader.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"
#include <memory>

bool test_proj(unsigned Ntest,
              unsigned Npart,
              unsigned order) noexcept {

    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    
    fastEEC::normType norm = fastEEC::normType::SUMPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        fastEEC::result_t<double> EEC_gen;

        runFastEEC(
            EEC_gen,
            genJet, RLax, norm,
            order, 0
        );

        for(unsigned testorder=2; testorder<=order; ++testorder){
            double accumulation = recursive_reduce(
                *EEC_gen.wts[testorder-2], 0.0
            );
            passed = passed && std::abs(accumulation - 1.0) < 1e-8;
        }

   }

    return passed;
}

bool test_proj_PU(unsigned Ntest,
              unsigned Npart,
              unsigned order) noexcept {

    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    

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
            order, fastEEC::DOPU,
            nullptr, nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            0.0,
            &isPU
        );

        for(unsigned testorder=2; testorder<=order; ++testorder){
            double accumulation = recursive_reduce(
                *EEC_gen.wts[testorder-2], 0.0
            );
            passed = passed && std::abs(accumulation - 1.0) < 1e-8;

            double accumulation_PU = recursive_reduce(
                *EEC_gen.wts_PU[testorder-2], 0.0
            );
            passed = passed && accumulation_PU > 0 && accumulation_PU < 1.0;
        }

   }

    return passed;
}

bool test_proj_transfer(
        unsigned Ntest,
        unsigned Npart,
        unsigned order) noexcept {

    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    

    fastEEC::normType norm = fastEEC::normType::RAWPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        jet recoJet;
        make_random_jet(recoJet, Npart);

        Eigen::MatrixXd ptrans = Eigen::MatrixXd::Identity(Npart, Npart);

        /*for(unsigned iPart=0; iPart < Npart; ++iPart){
            ptrans(iPart, iPart) = recoJet.particles[iPart].pt / genJet.particles[iPart].pt;
        }*/
        //ptrans *= genJet.rawpt / recoJet.rawpt;

        std::vector<bool> PU(Npart, false);
        std::vector<bool> UM(Npart, false);

        PU[0] = true;
        UM[0] = true;
        ptrans(0,0) = 0.0;

        fastEEC::result_t<double> EEC_gen;

        runFastEEC(
            EEC_gen,
            genJet, RLax, norm,
            order, fastEEC::DOPU | fastEEC::DOTRANSFER,
            nullptr, nullptr, nullptr,
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
            order, fastEEC::DOPU,
            nullptr, nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            nullptr, nullptr,
            0.0,
            &PU
        );

        for(unsigned testorder=2; testorder<=order; ++testorder){

            std::vector<double> reco_PU = *EEC_reco.wts_PU[testorder-2];
            std::vector<double> reco = *EEC_reco.wts[testorder-2];
            std::vector<double> recopure(reco.size(), 0.0);

            for(unsigned iReco=0; iReco<reco.size(); ++iReco){
                recopure[iReco] = reco[iReco] - reco_PU[iReco];
            }

            std::vector<double> wts_transfer(reco.size(), 0.0);
            for (unsigned iGen=0; iGen<reco.size(); ++iGen){
                for(unsigned iReco=0; iReco<reco.size(); ++iReco){
                    wts_transfer[iReco] += (*EEC_gen.transfer_wts[testorder-2])[iGen][iReco];
                }
            }

            double accumulation = recursive_reduce(
                *EEC_reco.wts[testorder-2], 0.0
            );

            double accumulation_PU = recursive_reduce(
                *EEC_reco.wts_PU[testorder-2], 0.0
            );

            double accumulation_transfer = recursive_reduce(
                *EEC_gen.transfer_wts[testorder-2], 0.0
            );

            passed = passed && std::abs(accumulation - accumulation_PU - accumulation_transfer) < 1e-8;
        }
   }

    return passed;
}

#endif
