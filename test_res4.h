#ifndef TEST_RES4_H
#define TEST_RES4_H

#include "example_jet.h"
#include "axes.h"
#include <boost/histogram.hpp>
#include "SRothman/EECs/src/theOnlyHeader.h"
#include "SRothman/SimonTools/src/recursive_reduce.h"
#include <memory>
#include "libnpy/include/npy.hpp"
#include "SRothman/SimonTools/src/ToyShowerer.h"

bool test_res4(unsigned Ntest, unsigned Npart) noexcept {
    auto RLax = get_RLax();
    auto RLax_coarse = get_RLax_coarse();
    auto r_dipole_ax = get_r_dipole_ax();
    auto ct_dipole_ax = get_phi_dipole_ax();
    auto r_tee_ax = get_r_tee_ax();
    auto ct_tee_ax = get_phi_tee_ax();
    auto r_tri_ax = get_r_tri_ax();
    auto ct_tri_ax = get_phi_tri_ax();
    auto r_minR_ax = get_r_minR_ax();
    auto phi_minR_ax = get_phi_minR_ax();
    auto phidiff_minR_ax = get_phidiff_minR_ax();
    auto theta_minR_ax = get_theta_minR_ax();

    fastEEC::normType norm = fastEEC::normType::SUMPT;

    fastEEC::result_t<double> accu;

    bool passed = true;

    ToyShowerer showerer("UNIFORM", "GLUON", "LNX", "OFF",
                         0.1, 0.1, 0.5, "");

    for (unsigned REPEAT=0; (REPEAT<Ntest) && passed; ++REPEAT){
        jet genJet;
        showerer.shower(100, 0, 0, 0, Npart, genJet);

        //make_random_jet(genJet, Npart);

        double acc_res4_dipole, acc_res4_tee, acc_res4_tri;// acc_res4_minR;

        if(REPEAT==0){
            runFastEEC(
                accu,
                genJet, RLax, norm,
                4, fastEEC::DORES4,
                RLax_coarse, nullptr, nullptr,
                r_dipole_ax, ct_dipole_ax,
                r_tee_ax, ct_tee_ax,
                r_tri_ax, ct_tri_ax,
                r_minR_ax, phi_minR_ax,
                phidiff_minR_ax, theta_minR_ax,
                0.05
            );

            acc_res4_dipole = recursive_reduce(
                *accu.resolved4_shapes->dipole, 0.0
            );
            acc_res4_tee = recursive_reduce(
                *accu.resolved4_shapes->tee, 0.0
            );
            acc_res4_tri = recursive_reduce(
                *accu.resolved4_shapes->triangle, 0.0
            );
            /*acc_res4_minR = recursive_reduce(
                *accu.resolved4_shapes->minR_phi1, 0.0
            );*/
        
        } else {
            fastEEC::result_t<double> EEC_gen;

            runFastEEC(
                EEC_gen,
                genJet, RLax, norm,
                4, fastEEC::DORES4,
                RLax_coarse, nullptr, nullptr,
                r_dipole_ax, ct_dipole_ax,
                r_tee_ax, ct_tee_ax,
                r_tri_ax, ct_tri_ax,
                r_minR_ax, phi_minR_ax,
                phidiff_minR_ax, theta_minR_ax,
                0.05
            );

            accu += EEC_gen;

            acc_res4_dipole = recursive_reduce(
                *EEC_gen.resolved4_shapes->dipole, 0.0
            );
            acc_res4_tee = recursive_reduce(
                *EEC_gen.resolved4_shapes->tee, 0.0
            );
            acc_res4_tri = recursive_reduce(
                *EEC_gen.resolved4_shapes->triangle, 0.0
            );
            /*acc_res4_minR = recursive_reduce(
                *EEC_gen.resolved4_shapes->minR_phi1, 0.0
            );*/
        }
        passed = passed && acc_res4_dipole > 0 && acc_res4_tee > 0 && acc_res4_tri > 0; //&& acc_res4_minR > 0;
    }

    std::vector<double> proj;
    std::vector<unsigned long> proj_shape;

    std::vector<double> res4_dipole;
    std::vector<unsigned long> res4_dipole_shape;

    std::vector<double> res4_tee;
    std::vector<unsigned long> res4_tee_shape;

    std::vector<double> res4_triangle;
    std::vector<unsigned long> res4_triangle_shape;

    proj.reserve(5*accu.wts[0]->size());
    proj.insert(proj.end(), accu.wts[0]->begin(), accu.wts[0]->end());
    proj.insert(proj.end(), accu.wts[1]->begin(), accu.wts[1]->end());
    proj.insert(proj.end(), accu.wts[2]->begin(), accu.wts[2]->end());
    proj.insert(proj.end(), accu.wts[3]->begin(), accu.wts[3]->end());
    proj.insert(proj.end(), accu.wts[4]->begin(), accu.wts[4]->end());

    proj_shape.push_back(5);
    proj_shape.push_back(accu.wts[0]->size());

    unsigned long Ndipole = accu.resolved4_shapes->dipole->num_elements();
    res4_dipole.reserve(Ndipole);
    res4_dipole.insert(res4_dipole.end(), accu.resolved4_shapes->dipole->data(), accu.resolved4_shapes->dipole->data()+Ndipole);
    
    res4_dipole_shape.push_back(accu.resolved4_shapes->dipole->shape()[0]);
    res4_dipole_shape.push_back(accu.resolved4_shapes->dipole->shape()[1]);
    res4_dipole_shape.push_back(accu.resolved4_shapes->dipole->shape()[2]);

    unsigned long Ntee = accu.resolved4_shapes->tee->num_elements();
    res4_tee.reserve(Ntee);
    res4_tee.insert(res4_tee.end(), accu.resolved4_shapes->tee->data(), accu.resolved4_shapes->tee->data()+Ntee);

    res4_tee_shape.push_back(accu.resolved4_shapes->tee->shape()[0]);
    res4_tee_shape.push_back(accu.resolved4_shapes->tee->shape()[1]);
    res4_tee_shape.push_back(accu.resolved4_shapes->tee->shape()[2]);

    unsigned long Ntriangle = accu.resolved4_shapes->triangle->num_elements();
    res4_triangle.reserve(Ntriangle);
    res4_triangle.insert(res4_triangle.end(), accu.resolved4_shapes->triangle->data(), accu.resolved4_shapes->triangle->data()+Ntriangle);
    
    res4_triangle_shape.push_back(accu.resolved4_shapes->triangle->shape()[0]);
    res4_triangle_shape.push_back(accu.resolved4_shapes->triangle->shape()[1]);
    res4_triangle_shape.push_back(accu.resolved4_shapes->triangle->shape()[2]);

    npy::npy_data<double> d_proj;
    d_proj.data = proj;
    d_proj.shape = proj_shape;
    npy::write_npy("test_res4_proj.npy", d_proj);

    npy::npy_data<double> d_res4_dipole;
    d_res4_dipole.data = res4_dipole;
    d_res4_dipole.shape = res4_dipole_shape;
    npy::write_npy("test_res4_dipole.npy", d_res4_dipole);

    npy::npy_data<double> d_res4_tee;
    d_res4_tee.data = res4_tee;
    d_res4_tee.shape = res4_tee_shape;
    npy::write_npy("test_res4_tee.npy", d_res4_tee);

    npy::npy_data<double> d_res4_triangle;
    d_res4_triangle.data = res4_triangle;
    d_res4_triangle.shape = res4_triangle_shape;
    npy::write_npy("test_res4_triangle.npy", d_res4_triangle);

    return passed;

}

bool test_res4_PU(unsigned Ntest, unsigned Npart) noexcept {
    auto RLax = get_RLax();
    auto RLax_coarse = get_RLax_coarse();
    auto r_dipole_ax = get_r_dipole_ax();
    auto ct_dipole_ax = get_phi_dipole_ax();
    auto r_tee_ax = get_r_tee_ax();
    auto ct_tee_ax = get_phi_tee_ax();
    auto r_tri_ax = get_r_tri_ax();
    auto ct_tri_ax = get_phi_tri_ax();
    auto r_minR_ax = get_r_minR_ax();
    auto phi_minR_ax = get_phi_minR_ax();
    auto phidiff_minR_ax = get_phidiff_minR_ax();
    auto theta_minR_ax = get_theta_minR_ax();

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
            r_minR_ax, phi_minR_ax,
            phidiff_minR_ax, theta_minR_ax,
            0.05,
            &isPU
        );

        double acc_res4_dipole = recursive_reduce(
            *EEC_gen.resolved4_shapes->dipole, 0.0
        );
        double acc_res4_tee = recursive_reduce(
            *EEC_gen.resolved4_shapes->tee, 0.0
        );
        double acc_res4_tri = recursive_reduce(
            *EEC_gen.resolved4_shapes->triangle, 0.0
        );
        double acc_res4_minR = recursive_reduce(
            *EEC_gen.resolved4_shapes->minR_phi1, 0.0
        );
        passed = passed && acc_res4_dipole > 0 && acc_res4_tee > 0 && acc_res4_tri > 0 && acc_res4_minR > 0;

        double acc_res4_dipole_PU = recursive_reduce(
            *EEC_gen.resolved4_shapes_PU->dipole, 0.0
        );
        double acc_res4_tee_PU = recursive_reduce(
            *EEC_gen.resolved4_shapes_PU->tee, 0.0
        );
        double acc_res4_tri_PU = recursive_reduce(
            *EEC_gen.resolved4_shapes_PU->triangle, 0.0
        );
        double acc_res4_minR_PU = recursive_reduce(
            *EEC_gen.resolved4_shapes_PU->minR_phi1, 0.0
        );

        passed = passed && acc_res4_dipole_PU > 0 && acc_res4_tee_PU > 0 && acc_res4_tri_PU > 0 && acc_res4_minR_PU > 0;
        passed = passed && acc_res4_dipole_PU < acc_res4_dipole;
        passed = passed && acc_res4_tee_PU < acc_res4_tee;
        passed = passed && acc_res4_tri_PU < acc_res4_tri;
        //passed = passed && acc_res4_minR_PU < acc_res4_minR;
    }

    return passed;
}

#endif
