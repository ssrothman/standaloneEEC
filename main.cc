#include <stdio.h>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/printPart.h"

#include "SRothman/Matching/src/matcher.h"

#include "SRothman/EECs/src/fast.h"

void setup_example_genjet(jet& result, const jet& recojet){
    static std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
    static std::uniform_real_distribution<double> rand(0,1);
    static std::normal_distribution<double> normal(0, 1);
    static std::exponential_distribution<double> exp(0.2);
    static std::poisson_distribution<int> poisson(20);

    result.eta = recojet.eta;
    result.phi = recojet.phi;

    for (unsigned iPart=0; iPart < recojet.nPart; ++iPart){
        const particle& recoPart = recojet.particles[iPart];
        double pt = recoPart.pt * (1 + 0.001*normal(rng));
        double eta = recoPart.eta + 0.001*normal(rng);
        double phi = recoPart.phi + 0.001*normal(rng);

        //in principle theres also some chance
        //that it isn't perfect one-to-one
        //but whatever
        result.particles.emplace_back(
                pt, eta, phi,
                0,  0,   0,
                recoPart.pdgid, recoPart.charge);
        ++result.nPart;
        result.pt += pt;
        result.sumpt += pt;
        result.rawpt += pt;
    }
}

void setup_example_recojet(jet& result){
    static std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
    static std::uniform_real_distribution<double> rand(0,1);
    static std::normal_distribution<double> normal(0, 1);
    static std::exponential_distribution<double> exp(0.2);
    static std::poisson_distribution<int> poisson(20);
    
    int nPart = poisson(rng);

    for (unsigned iPart=0; iPart<nPart; ++iPart){
        double pt = exp(rng);
        double eta = normal(rng);
        double phi = normal(rng);

        double pdgid_f = rand(rng);
        int pdgid;
        int charge;
        if(pdgid_f < 0.6){
            pdgid = 211;
            if(pdgid_f < 0.3){
                charge = 1;
            } else {
                charge = -1;
            }
        } else if(pdgid_f < 0.65){
            pdgid = 11;
            if(pdgid_f < 0.625){
                charge = 1;
            } else {
                charge = -1;
            }
        } else if(pdgid_f < 0.7){
            pdgid = 13;
            if(pdgid_f < 0.675){
                charge = 1;
            } else {
                charge = -1;
            }
        } else if(pdgid_f < 0.9){
            pdgid = 22;
            charge = 0;
        } else {
            pdgid = 130;
            charge = 0;
        }

        result.particles.emplace_back(
                pt, eta, phi,
                0,  0,   0,
                pdgid, charge);
        ++result.nPart;

        result.pt += pt;
        result.sumpt += pt;
        result.rawpt += pt;
    }
    result.pt *= 1.1; //fake JEC
    result.eta = 0;
    result.phi = 0;

    std::sort(result.particles.begin(), result.particles.end(),
            [](const particle& a, const particle& b){
                return a.pt > b.pt;
            });
}

int main(){
    //setup example reco jet
    jet recoJet;
    setup_example_recojet(recoJet);
    printf("Reco Jet: (%0.3f, %0.3f, %0.3f)\n", 
            recoJet.pt, recoJet.eta, recoJet.phi);
    printf("PARTICLES:\n");
    for (const auto& part : recoJet.particles){
        printPart(part);
    }

    //setup example gen jet
    jet genJet;
    setup_example_genjet(genJet, recoJet);
    printf("Gen Jet: (%0.3f, %0.3f, %0.3f)\n", 
            genJet.pt, genJet.eta, genJet.phi);
    printf("PARTICLES:\n");
    for (const auto& part : genJet.particles){
        printPart(part);
    }

    //setup matcher
    //there's a ton of configuration values here
    matcher match(recoJet, genJet,

                  0.05,

                  spatialLoss::TYPE1,

                  {1.0, 1.0, 1.0, 1.0, 1.0},
                  {4.0, 4.0, 4.0, 4.0, 4.0},
                  {2e12, 2e12, 2e12, 2e12, 2e12},

                  "Standard",

                  {"AnyPhoton", "AnyNeutralHadron", 
                  "AnyCharged", "AnyCharged", "AnyCharged"},
                  {"AnyPhoton", "AnyNeutralHadron", 
                  "AnyCharged", "AnyCharged", "AnyCharged"},
                  {0.0, 0.0, 0.0, 0.0, 0.0},

                  {"ChargeSign", "ChargeSign",
                  "ChargeSign", "ChargeSign", "ChargeSign"},
                  {"Fixed", "Fixed", 
                   "Tracking", "Tracking", "Tracking"},
                  {"Best", "Best", "Best", "Best", "Best"},

                  "OneGenOneRecoPerType",

                  "NONE", 
                  "NONE",
                  
                  false, true,
                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0},
                  4.2,

                  false,
                  {0.0, 0.0, 0.0},

                  {0.17, 0.18, 0.80},
                  {0.0074, 0.0253, 0.0253},
                  {0.05, 0.05, 0.07},
                  {0.05, 0.05, 0.07},
                  {0.0, 0.9, 1.4, 3.0},

                  {1.63, 3.90, 6.44},
                  {0.21, 0.14, 0.10},
                  {0.10, 0.10, 0.15},
                  {0.10, 0.10, 0.15},
                  {0.0, 0.9, 0.14, 3.0},

                  {0.000069, 0.000072, 0.000072},
                  {0.0076, 0.014, 0.018},
                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0},
                  {0.001, 0.001, 0.001},
                  {0.001, 0.001, 0.001},
                  {0.0, 0.9, 1.4, 3.0},

                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0},
                  {0.0, 0.0, 0.0},

                  {0.050, 0.050, 0.050},
                  {0.000, 0.000, 0.000},
                  {0.000, 0.000, 0.000},

                  {0.100, 0.100, 0.150},
                  {0.000, 0.000, 0.000},
                  {0.000, 0.000, 0.000},

                  {0.002, 0.002, 0.003},
                  {0.010, 0.010, 0.015},
                  {0.050, 0.050, 0.070},

                  {0.002, 0.002, 0.003},
                  {0.010, 0.010, 0.015},
                  {0.050, 0.050, 0.070},

                  {0.002, 0.002, 0.003},
                  {0.010, 0.010, 0.015},
                  {0.050, 0.050, 0.070},

                  50,

                  100);
    match.minimize();

    //matching matrix
    //can be either ptrans() -> includes pT smearing
    //           or rawmat() -> ignores pT smearing
    arma::mat matchMatrix = match.ptrans();

    std::vector<bool> isPU(recoJet.nPart, true);
    std::vector<bool> isUmmatched(genJet.nPart, true);
    for(unsigned iPReco=0; iPReco<recoJet.nPart; ++iPReco){
        for(unsigned iPGen=0; iPGen<genJet.nPart; ++iPGen){
            if(matchMatrix(iPReco, iPGen) > 0){
                isPU[iPReco] = false;
                isUmmatched[iPGen] = false;
            }
        }
    }

    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    int maxOrder = 6;
    fastEEC::normType norm = fastEEC::normType::RAWPT;

    fastEEC::result<double> EEC_reco;
    fastEEC::fastEEC<double, true, false, false, false, false>(
        EEC_reco,
        recoJet, RLax, maxOrder, norm,
        RLax, RLax, RLax,//these are a bunch of other axes
                         //for the resolved correlators
                         //which we aren't calculating
                         //but you still have to pass something
                         //so whatever
        RLax, RLax, 
        RLax, RLax, 
        RLax, RLax, 
        0.0,
        &isPU
    );

    fastEEC::result<double> EEC_gen;
    fastEEC::fastEEC<double, true, true, false, false, false, false>(
        EEC_gen,
        genJet, RLax, maxOrder, norm,
        RLax, RLax, RLax,
        RLax, RLax,
        RLax, RLax,
        RLax, RLax,
        0.0,
        &isUmmatched, &recoJet, &matchMatrix
    );

    //binned EEC results are in EEC_reco and EEC_gen
    //and also the transfer matrix entries are in EEC_gen

    return 0;
}
