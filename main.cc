#include <stdio.h>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/printPart.h"

#include "SRothman/Matching/src/matcher.h"

#include "SRothman/EECs/src/run.h"

#include "testcoords.h"

void setup_example_recojet(jet& recoJet){
    double totalPt=0.0;
    double weightedEta=0.0;
    double weightedPhi=0.0;

    for(unsigned i=0; i<etas.size(); ++i){
        recoJet.particles.emplace_back(
            pts.at(i),
            etas.at(i),
            phis.at(i)
        );
        ++recoJet.nPart;

        totalPt += pts.at(i);
        weightedEta += pts.at(i)*etas.at(i);
        weightedPhi += pts.at(i)*phis.at(i);
    }
    weightedEta /= totalPt;
    weightedPhi /= totalPt;

    recoJet.pt = totalPt;
    recoJet.sumpt = totalPt;
    recoJet.rawpt = totalPt;

    recoJet.eta = weightedEta;
    recoJet.phi = weightedPhi;
}

void make_random_jet(jet& recoJet, unsigned nPart){
    static std::default_random_engine gen(12);
    static std::normal_distribution<double> normal(0.0, 0.4);

    double totalPt=0;
    double weightedEta=0;
    double weightedPhi=0;
    for(unsigned i=0; i<nPart; ++i){
        double pt = 1.0;
        double eta = normal(gen);
        double phi = normal(gen);

        recoJet.particles.emplace_back(
            pt,
            eta,
            phi
        );
        ++recoJet.nPart;

        totalPt += pt;
        weightedEta += pt*eta;
        weightedPhi += pt*phi;
    }

    weightedEta /= totalPt;
    weightedPhi /= totalPt;

    recoJet.pt = totalPt;
    recoJet.sumpt = totalPt;
    recoJet.rawpt = totalPt;

    recoJet.eta = weightedEta;
    recoJet.phi = weightedPhi;
}

int main(){
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.001, 0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);

    std::vector<double> xiedges({1e-3, 0.5, 1.01});
    auto xi_ax = std::make_shared<boost::histogram::axis::variable<double>>(xiedges);

    std::vector<double> phi_edges({1e-3, 0.5, 1.0});
    auto phi_ax = std::make_shared<boost::histogram::axis::variable<double>>(phi_edges);

    std::vector<double> r_dipole_edges({0.0, 0.5, 1.0});
    auto r_dipole_ax = std::make_shared<boost::histogram::axis::variable<double>>(r_dipole_edges);

    std::vector<double> ct_dist_edges({0.0, 0.5, 1.0});
    auto ct_dist_ax = std::make_shared<boost::histogram::axis::variable<double>>(ct_dist_edges);

    std::vector<double> r_tee_edges({0.0, 0.5, 1.0});
    auto r_tee_ax = std::make_shared<boost::histogram::axis::variable<double>>(r_tee_edges);

    std::vector<double> ct_tee_edges({0.0, 0.5, 1.0});
    auto ct_tee_ax = std::make_shared<boost::histogram::axis::variable<double>>(ct_tee_edges);

    std::vector<double> r_tri_edges({0.0, 0.5, 1.0});
    auto r_tri_ax = std::make_shared<boost::histogram::axis::variable<double>>(r_tri_edges);

    std::vector<double> ct_tri_edges({0.0, 0.5, 1.0});
    auto ct_tri_ax = std::make_shared<boost::histogram::axis::variable<double>>(ct_tri_edges);

    fastEEC::normType norm = fastEEC::normType::RAWPT;


    for(int REP=0; REP<3; ++REP){
        //printf("\n\n\n\n\n\n");
        //setup example reco jet
        jet recoJet;
        //setup_example_recojet(recoJet);
        make_random_jet(recoJet, 50);
        //for(const auto& part : recoJet.particles){
        //    printf("%0.3f, %0.3f, %0.3f\n", part.pt, part.eta, part.phi);
        //}
        //printf("Reco Jet: (%0.3f, %0.3f, %0.3f)\n", 
        //        recoJet.pt, recoJet.eta, recoJet.phi);
        //printf("PARTICLES:\n");
        //for (const auto& part : recoJet.particles){
        //    printPart(part);
        //}

        fastEEC::result_t<double> EEC_reco;
        fastEEC::run<double, false, false, true, true, false, 4>(
            EEC_reco,
            recoJet, RLax, norm,
            RLax_coarse, xi_ax, phi_ax,
            r_dipole_ax, ct_dist_ax,
            r_tee_ax, ct_tee_ax,
            r_tri_ax, ct_tri_ax,
            0.1
        );
    }

    return 0;
}
