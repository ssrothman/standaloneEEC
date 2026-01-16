#include "SRothman/Matching/src/TrackMatcher.h"
#include "SRothman/SimonTools/src/JetFactory.h"
#include "SRothman/SimonTools/src/jet.h"

#include <stdio.h>
#include <random>
#include <vector>
#include <iostream>

constexpr int NPART = 10;

int main(){
    matching::TrackMatcher matcher(
        0.4,
        999,
        "TrackPt", //electron_dr_mode
        0.1,       //electron_dr_param1
        0.1,       //electron_dr_param2
        0.4,       //electron_dr_param3
        "TrackPt", //electron_ptres_mode
        0.1,       //electron_ptres_param1
        0.1,       //electron_ptres_param2
        "TrackAng",//electron_angres_mode
        0.1,       //electron_angres_param1
        0.1,       //electron_angres_param2
        10.0,      //electron_opp_charge_penalty
        10.0,      //electron_no_charge_penalty
        "Any",
        "Any",
        "TrackPt", //muon_dr_mode
        0.1,       //muon_dr_param1
        0.1,       //muon_dr_param2
        0.4,       //muon_dr_param3
        "TrackPt", //muon_ptres_mode
        0.1,       //muon_ptres_param1
        0.1,       //muon_ptres_param2
        "TrackAng",//muon_angres_mode
        0.1,       //muon_angres_param1
        0.1,       //muon_angres_param2
        10.0,      //muon_opp_charge_penalty
        10.0,      //muon_no_charge_penalty
        "Any",
        "Any",
        "TrackPt", //hadch_dr_mode
        0.1,       //hadch_dr_param1
        0.1,       //hadch_dr_param2
        0.4,       //hadch_dr_param3
        "TrackPt", //hadch_ptres_mode
        0.1,       //hadch_ptres_param1
        0.1,       //hadch_ptres_param2
        "TrackAng",//hadch_angres_mode
        0.1,       //hadch_angres_param1
        0.1,       //hadch_angres_param2
        10.0,      //hadch_opp_charge_penalty
        10.0,      //hadch_no_charge_penalty
        "Any",
        "Any" 
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
