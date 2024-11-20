#ifndef TEST_JETS_H
#define TEST_JETS_H

#include <unordered_map>
#include <random>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/printPart.h"

void example_recojet(const jet& genJet, jet& recoJet, Eigen::MatrixXd& ptrans,
                    std::vector<bool>& PU, std::vector<bool>& UM) noexcept {
    static constexpr float pmiss = 0.20;
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

    ptrans = Eigen::MatrixXd(recoJet.nPart, genJet.nPart);
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
    static std::uniform_real_distribution<double> uniform(0.0, 1.0);
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

    double extrapt = uniform(gen);
    double JEC = 1 + 0.2*uniform(gen)-0.05;
    genJet.pt = (totalPt + extrapt) * JEC;
    genJet.sumpt = totalPt;
    genJet.rawpt = (totalPt + extrapt);

    genJet.eta = weightedEta;
    genJet.phi = weightedPhi;
}



#endif
