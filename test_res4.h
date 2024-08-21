#ifndef TEST_RES4_H
#define TEST_RES4_H

#include "example_jet.h"
#include <boost/histogram.hpp>
#include "SRothman/EECs/src/theOnlyHeader.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"
#include <memory>

bool test_res4(unsigned Ntest, unsigned Npart) noexcept {
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);
    
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


    fastEEC::normType norm = fastEEC::normType::SUMPT;

    bool passed = true;
    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        make_random_jet(genJet, Npart);

        fastEEC::result_t<double> EEC_gen;

        runFastEEC(
            EEC_gen,
            genJet, RLax, norm,
            4, fastEEC::DORES4,
            RLax_coarse, nullptr, nullptr,
            r_dipole_ax, ct_dipole_ax,
            r_tee_ax, ct_tee_ax,
            r_tri_ax, ct_tri_ax,
            0.05
        );

        double acc_res4_dipole = recursive_reduce(
            *EEC_gen.resolved4_shapes->dipole, 0.0
        );
        double acc_res4_tee = recursive_reduce(
            *EEC_gen.resolved4_shapes->tee, 0.0
        );
        passed = passed && acc_res4_dipole > 0 && acc_res4_tee > 0;
    }

    return passed;

}

bool test_res4_PU(unsigned Ntest, unsigned Npart) noexcept {
    std::vector<double> RLedges({1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0});
    auto RLax = std::make_shared<boost::histogram::axis::variable<double>>(RLedges);
    std::vector<double> RLedges_coarse({0.10, 0.20});
    auto RLax_coarse = std::make_shared<boost::histogram::axis::variable<double>>(RLedges_coarse);
    
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
            4, fastEEC::DORES4 | fastEEC::DOPU,
            RLax_coarse, nullptr, nullptr,
            r_dipole_ax, ct_dipole_ax,
            r_tee_ax, ct_tee_ax,
            r_tri_ax, ct_tri_ax,
            0.05,
            &isPU
        );

        double acc_res4_dipole = recursive_reduce(
            *EEC_gen.resolved4_shapes->dipole, 0.0
        );
        double acc_res4_tee = recursive_reduce(
            *EEC_gen.resolved4_shapes->tee, 0.0
        );
        passed = passed && acc_res4_dipole > 0 && acc_res4_tee > 0;

        double acc_res4_dipole_PU = recursive_reduce(
            *EEC_gen.resolved4_shapes_PU->dipole, 0.0
        );
        double acc_res4_tee_PU = recursive_reduce(
            *EEC_gen.resolved4_shapes_PU->tee, 0.0
        );
        passed = passed && acc_res4_dipole_PU > 0 && acc_res4_tee_PU > 0;
        passed = passed && acc_res4_dipole_PU < acc_res4_dipole;
        passed = passed && acc_res4_tee_PU < acc_res4_tee;
    }

    return passed;
}

#endif
