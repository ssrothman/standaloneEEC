#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <unordered_map>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/printPart.h"

#include "SRothman/Matching/src/matcher.h"

#include "SRothman/EECs/src/run.h"

#include "SRothman/SimonTools/src/recursive_reduce.h"
#include "testcoords.h"

#include "SRothman/armadillo-12.2.0/include/armadillo"

void example_recojet(const jet& genJet, jet& recoJet, arma::mat& ptrans,
                    std::vector<bool>& PU, std::vector<bool>& UM){
    static constexpr float pmiss = 0.1;
    static constexpr float psplit = 0.1;
    static std::default_random_engine gen(12);
    static std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double totalPt=0.0;
    double weightedEta=0.0;
    double weightedPhi=0.0;

    PU.clear();
    UM.clear();

    std::unordered_map<unsigned, std::vector<unsigned>> genToReco;

    for (unsigned iGen =0; iGen < genJet.nPart; ++iGen){
        const auto& pGen = genJet.particles[iGen];

        if (uniform(gen) < pmiss){
            PU.emplace_back(1);
            UM.emplace_back(1);
            continue;
        } else {
            PU.emplace_back(0);
            UM.emplace_back(0);

            if(uniform(gen) < psplit){
                recoJet.particles.emplace_back(pGen.pt/2, 
                        pGen.eta+0.02, pGen.phi-0.02);

                genToReco[iGen].emplace_back(recoJet.nPart);
                ++recoJet.nPart;

                totalPt += pGen.pt/2;
                weightedEta += pGen.pt*(pGen.eta+0.02)/2;
                weightedPhi += pGen.pt*(pGen.phi-0.02)/2;

                recoJet.particles.emplace_back(pGen.pt/2, 
                        pGen.eta-0.02, pGen.phi+0.02);

                genToReco[iGen].emplace_back(recoJet.nPart);
                ++recoJet.nPart;

                totalPt += pGen.pt/2;
                weightedEta += pGen.pt*(pGen.eta-0.02)/2;
                weightedPhi += pGen.pt*(pGen.phi+0.02)/2;
            } else {
                recoJet.particles.emplace_back(pGen.pt, pGen.eta, pGen.phi);

                genToReco[iGen].emplace_back(recoJet.nPart);
                ++recoJet.nPart;

                totalPt += pGen.pt;
                weightedEta += pGen.pt*pGen.eta;
                weightedPhi += pGen.pt*pGen.phi;
            }
        }
    }

    weightedEta /= totalPt;
    weightedPhi /= totalPt;

    recoJet.pt = totalPt;
    recoJet.sumpt = totalPt;
    recoJet.rawpt = totalPt;

    recoJet.eta = weightedEta;
    recoJet.phi = weightedPhi;

    ptrans = arma::mat(recoJet.nPart, genJet.nPart, arma::fill::zeros);
    for (const auto& [iGen, iRecos] : genToReco){
        if (iRecos.empty()){
            continue;
        } else {
            const unsigned long len = iRecos.size();
            for (const auto& iReco : iRecos){
                ptrans(iReco, iGen) = 1.0/len;
            }
        }
    }

    std::cout << ptrans << std::endl;
}

void make_random_jet(jet& genJet, unsigned nPart){
    static std::default_random_engine gen(12);
    static std::normal_distribution<double> normal(0.0, 0.4);
    //static std::uniform_real_distribution<double> normal(-0.4, 0.4);

    double totalPt=0;
    double weightedEta=0;
    double weightedPhi=0;
    for(unsigned i=0; i<nPart; ++i){
        double pt = 1.0;
        double eta = normal(gen);
        double phi = normal(gen);

        genJet.particles.emplace_back(
            pt,
            eta,
            phi
        );
        ++genJet.nPart;

        totalPt += pt;
        weightedEta += pt*eta;
        weightedPhi += pt*phi;
    }

    weightedEta /= totalPt;
    weightedPhi /= totalPt;

    genJet.pt = totalPt;
    genJet.sumpt = totalPt;
    genJet.rawpt = totalPt;

    genJet.eta = weightedEta;
    genJet.phi = weightedPhi;
}

int main(){
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);

    std::vector<double> xiedges({1e-3, 0.5, 1.01});
    auto xi_ax = std::make_shared<boost::histogram::axis::variable<double>>(xiedges);

    std::vector<double> phi_edges({1e-3, 0.5, 1.0});
    auto phi_ax = std::make_shared<boost::histogram::axis::variable<double>>(phi_edges);

    std::vector<double> r_dipole_edges({0.0, 0.05, 0.1, 0.15, 0.2,
                                        0.25, 0.3, 0.35, 0.4, 0.45,
                                        0.5, 0.55, 0.6, 0.65, 0.7,
                                        0.75, 0.8, 0.85, 0.9, 0.95,
                                        1.0});
    auto r_dipole_ax = std::make_shared<boost::histogram::axis::variable<double>>(r_dipole_edges);

    std::vector<double> ct_dipole_edges({0.0, M_PI/32, 2*M_PI/32, 3*M_PI/32,
                                         4*M_PI/32, 5*M_PI/32, 6*M_PI/32,
                                         7*M_PI/32, 8*M_PI/32, 9*M_PI/32,
                                         10*M_PI/32, 11*M_PI/32, 12*M_PI/32,
                                         13*M_PI/32, 14*M_PI/32, 15*M_PI/32,
                                         16*M_PI/32});
                                
    auto ct_dipole_ax = std::make_shared<boost::histogram::axis::variable<double>>(ct_dipole_edges);

    std::vector<double> r_tee_edges = r_dipole_edges;
    auto r_tee_ax = std::make_shared<boost::histogram::axis::variable<double>>(r_tee_edges);

    std::vector<double> ct_tee_edges = ct_dipole_edges;
    auto ct_tee_ax = std::make_shared<boost::histogram::axis::variable<double>>(ct_tee_edges);

    std::vector<double> r_tri_edges({0.0, 0.5, 1.0});
    auto r_tri_ax = std::make_shared<boost::histogram::axis::variable<double>>(r_tri_edges);

    std::vector<double> ct_tri_edges({0.0, 0.5, 1.0});
    auto ct_tri_ax = std::make_shared<boost::histogram::axis::variable<double>>(ct_tri_edges);

    fastEEC::normType norm = fastEEC::normType::RAWPT;


    std::shared_ptr<fastEEC::result_t<double>> EEC_accu = nullptr;
    for(int REP=0; REP<2; ++REP){
        jet genJet;
        make_random_jet(genJet, 10);

        jet recoJet;
        arma::mat ptrans;
        std::vector<bool> PU, UM;
        example_recojet(genJet, recoJet, ptrans, PU, UM);

        auto EEC_reco = std::make_shared<fastEEC::result_t<double>>();

        fastEEC::runSuperSpecific<double>(
            *EEC_reco,
            genJet, RLax, norm,
            6, fastEEC::DORES4 | fastEEC::DORES3 | fastEEC::DOPU | fastEEC::DOTRANSFER,
            RLax_coarse, xi_ax, phi_ax,
            r_dipole_ax, ct_dipole_ax,
            r_tee_ax, ct_tee_ax,
            r_tri_ax, ct_tri_ax,
            0.05,
            &UM, &recoJet, &ptrans
        );

        if (REP==0){
            EEC_accu = EEC_reco;
        } else {
            *EEC_accu += *EEC_reco;
        }

        printf("Ran %u:\n", REP);
        //printf("\tdipole: %g\n", recursive_reduce(*(EEC_reco->resolved4_shapes->dipole), 0.0));
        //printf("\ttee   : %g\n", recursive_reduce(*(EEC_reco->resolved4_shapes->tee), 0.0));
        EEC_reco->summarize();
        fflush(stdout);
    }

    printf("\n\n");
    printf("Final:\n");
    EEC_accu->summarize();
    fflush(stdout);

    fastEEC::dumpToFile(*(EEC_accu->resolved4_shapes->dipole), "dipole.dat");
    fastEEC::dumpToFile(*(EEC_accu->resolved4_shapes->tee), "tee.dat");
    return 0;
}
