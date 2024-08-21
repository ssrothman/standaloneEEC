#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <unordered_map>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/printPart.h"

#include "SRothman/Matching/src/matcher.h"

#include "SRothman/EECs/src/theOnlyHeader.h"

#include "SRothman/SimonTools/src/recursive_reduce.h"
#include "testcoords.h"

#include "cov.h"

#include "SRothman/armadillo-12.2.0/include/armadillo"

void example_recojet(const jet& genJet, jet& recoJet, arma::mat& ptrans,
                    std::vector<bool>& PU, std::vector<bool>& UM) noexcept {
    static constexpr float pmiss = 0.00;
    static constexpr float psplit = 0.00;
    static std::default_random_engine gen(12);
    static std::uniform_real_distribution<double> uniform(0.0, 1.0);
    static std::normal_distribution<double> smearpt(1.0, 0.03);
    static std::normal_distribution<double> smearep(0.0, 0.03);

    double totalPt=0.0;
    double weightedEta=0.0;
    double weightedPhi=0.0;

    PU.clear();
    UM.clear();

    std::unordered_map<unsigned, std::vector<unsigned>> genToReco;

    for (unsigned iGen =0; iGen < genJet.nPart; ++iGen){
        const auto& pGen = genJet.particles[iGen];

        if (uniform(gen) < pmiss){
            UM.emplace_back(1);
            continue;
        } else {
            UM.emplace_back(0);

            if(uniform(gen) < psplit){
                double newpt = pGen.pt/2 * smearpt(gen);
                double neweta = pGen.eta + smearep(gen);
                double newphi = pGen.phi + smearep(gen);

                recoJet.particles.emplace_back(
                        newpt, neweta, newphi
                );

                genToReco[iGen].emplace_back(recoJet.nPart);
                ++recoJet.nPart;

                totalPt += newpt;
                weightedEta += newpt * neweta;
                weightedPhi += newpt * newphi;

                double newpt2 = pGen.pt/2 * smearpt(gen);
                double neweta2 = pGen.eta + smearep(gen);
                double newphi2 = pGen.phi + smearep(gen);

                recoJet.particles.emplace_back(
                        newpt2, neweta2, newphi2
                );

                genToReco[iGen].emplace_back(recoJet.nPart);
                ++recoJet.nPart;

                totalPt += newpt2;
                weightedEta += newpt2 * neweta2;
                weightedPhi += newpt2 * newphi2;
                PU.emplace_back(0);
                PU.emplace_back(0);
            } else {
                double newpt = pGen.pt * smearpt(gen);
                double neweta = pGen.eta + smearep(gen);
                double newphi = pGen.phi + smearep(gen);

                recoJet.particles.emplace_back(
                        newpt, neweta, newphi
                );

                genToReco[iGen].emplace_back(recoJet.nPart);
                ++recoJet.nPart;

                totalPt += newpt;
                weightedEta += newpt * neweta;
                weightedPhi += newpt * newphi;
                PU.emplace_back(0);
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
            for (const auto& iReco : iRecos){
                ptrans(iReco, iGen) = recoJet.particles[iReco].pt / genJet.particles[iGen].pt * genJet.pt / recoJet.pt;
            }
        }
    }

    //std::cout << ptrans << std::endl;
}

void make_random_jet(jet& genJet, unsigned nPart) noexcept {
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

int main() noexcept {
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

    std::shared_ptr<fastEEC::result_t<double>> EEC_reco_accu = nullptr;
    std::shared_ptr<fastEEC::result_t<double>> EEC_gen_accu = nullptr;
    /*auto cov = std::make_shared<EECcov<double>>();
    cov->setup(fastEEC::AXextent(*RLax), 5,
               fastEEC::AXextent(*RLax_coarse),
               fastEEC::AXextent(*xi_ax),
               fastEEC::AXextent(*phi_ax),
               fastEEC::AXextent(*RLax_coarse),
               fastEEC::AXextent(*r_dipole_ax),
               fastEEC::AXextent(*ct_dipole_ax),
               fastEEC::AXextent(*RLax_coarse),
               fastEEC::AXextent(*r_tee_ax),
               fastEEC::AXextent(*ct_tee_ax));*/

    for(int REP=0; REP<1; ++REP){
        jet genJet;
        make_random_jet(genJet, 50);

        jet recoJet;
        arma::mat ptrans;
        std::vector<bool> PU, UM;
        example_recojet(genJet, recoJet, ptrans, PU, UM);

        auto EEC_reco = std::make_shared<fastEEC::result_t<double>>();
        auto EEC_gen = std::make_shared<fastEEC::result_t<double>>();

        printf("RUNNING GEN\n");
        fflush(stdout);
        runFastEEC(
            *EEC_gen,
            genJet, RLax, norm,
            4, fastEEC::DORES4 | fastEEC::DORES3 | fastEEC::DOPU | fastEEC::DOTRANSFER,
            RLax_coarse, xi_ax, phi_ax,
            r_dipole_ax, ct_dipole_ax,
            r_tee_ax, ct_tee_ax,
            r_tri_ax, ct_tri_ax,
            0.05,
            &UM, &recoJet, &ptrans
        );

        printf("RUNNING RECO\n");
        fflush(stdout);
        runFastEEC(
                *EEC_reco,
                recoJet, RLax, norm,
                4, fastEEC::DORES4 | fastEEC::DORES3 | fastEEC::DOPU,
                RLax_coarse, xi_ax, phi_ax,
                r_dipole_ax, ct_dipole_ax,
                r_tee_ax, ct_tee_ax,
                r_tri_ax, ct_tri_ax,
                0.05,
                &PU
        );

        /*printf("FILLING COV\n");
        fflush(stdout);
        cov->fill(*EEC_reco);*/

        printf("ACCUMULATING\n");
        fflush(stdout);
        if (REP==0){
            EEC_reco_accu = EEC_reco;
            EEC_gen_accu = EEC_gen;
        } else {
            *EEC_reco_accu += *EEC_reco;
            *EEC_gen_accu += *EEC_gen;
        }


        printf("Ran %u:\n", REP);
        //printf("\tdipole: %g\n", recursive_reduce(*(EEC_reco->resolved4_shapes->dipole), 0.0));
        //printf("\ttee   : %g\n", recursive_reduce(*(EEC_reco->resolved4_shapes->tee), 0.0));
        printf("RECO:\n");
        EEC_reco->summarize();
        printf("GEN:\n");
        EEC_gen->summarize();
        fflush(stdout);
    }

    printf("\n\n");
    printf("Final:\n");
    printf("RECO:\n");
    EEC_reco_accu->summarize();
    printf("GEN:\n");
    EEC_gen_accu->summarize();
    fflush(stdout);

    //fastEEC::dumpToFile(*(EEC_accu->resolved4_shapes->dipole), "dipole.dat");
    //fastEEC::dumpToFile(*(EEC_accu->resolved4_shapes->tee), "tee.dat");
    return 0;
}
