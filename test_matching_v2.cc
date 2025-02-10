#include "SRothman/Matching/src/v2/TrackMatcher.h"
#include "SRothman/SimonTools/src/jet.h"
#include <stdio.h>
#include <random>
#include <vector>
#include <iostream>

constexpr int NPART = 10;

class JetFactory {
public:
    JetFactory(){
        initialize();
    }

    void initialize(){
        rng.seed(0);
        norm = std::normal_distribution<double>(0, 0.4);
        gamma = std::gamma_distribution<double>(2, 2);
    }

    void makeJet(simon::jet& J, const int nPart){
        J.nPart = nPart;

        J.sumpt = 0;
        J.pt = 0;
        J.rawpt = 0;

        J.eta = 0;
        J.phi = 0;

        J.particles.clear();
        for(int i=0; i<nPart; ++i){
            double pt = gamma(rng);
            double phi = norm(rng);
            double eta = norm(rng);
            J.particles.emplace_back(pt, phi, eta);
            J.sumpt += pt;
            J.pt += pt;
            J.rawpt += pt;
        }
    }

    void makeTransferJet(const simon::jet& J1, simon::jet& J2, Eigen::MatrixXd& tmat){
        J2.nPart = 0;

        J2.sumpt = 0;
        J2.pt = 0;
        J2.rawpt = 0;

        J2.eta = J1.eta;
        J2.phi = J1.phi;

        for (unsigned i=0; i<J1.nPart; ++i){
            double pt = J1.particles[i].pt;
            double phi = J1.particles[i].phi;
            double eta = J1.particles[i].eta;

            ++J2.nPart;
            J2.particles.emplace_back(pt, eta, phi);
            J2.sumpt += pt;
            J2.pt += pt;
            J2.rawpt += pt;

            tmat(i, i) = 1;
        }
    }
private:
    std::default_random_engine rng;
    std::normal_distribution<double> norm;
    std::gamma_distribution<double> gamma;
};


int main(){
    matching::TrackMatcher matcher(
        "Const",   //jet_dr_mode
        0.4,       //jet_dr_param1
        0.0,       //jet_dr_param2
        0.0,       //jet_dr_param3
        "Const",   //jet_ptres_mode
        30.0,      //jet_ptres_param1
        0.0,       //jet_ptres_param2
        "Const",   //jet_angres_mode
        0.4,       //jet_angres_param1
        0.0,       //jet_angres_param2
        "TrackPt", //particle_dr_mode
        0.1,       //particle_dr_param1
        0.1,       //particle_dr_param2
        0.4,       //particle_dr_param3
        "TrackPt", //particle_ptres_mode
        0.1,       //particle_ptres_param1
        0.1,       //particle_ptres_param2
        "TrackAng",//particle_angres_mode
        0.1,       //particle_angres_param1
        0.1,       //particle_angres_param2
        10.0,      //opp_charge_penalty
        10.0       //no_charge_penalty
    );

    simon::jet genjet, recojet;
    Eigen::MatrixXd tmat(NPART, NPART);
    JetFactory jetFactory;
    jetFactory.initialize();
    jetFactory.makeJet(genjet, NPART);
    jetFactory.makeTransferJet(genjet, recojet, tmat);

    Eigen::MatrixXd tmat2(NPART, NPART);
    matcher.matchParticles(genjet, recojet, tmat2);

    printf("DONE\n");

    std::cout << tmat2 << std::endl;

    return 0;
}
