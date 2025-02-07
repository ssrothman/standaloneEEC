#include <stdio.h>
#include <chrono>
#include <random>
#include <boost/histogram.hpp>

#include "SRothman/EECs/src/standalones/res4_standalone.h"
#include "SRothman/SimonTools/src/simon_jet.h"
#include "SRothman/EECs/src/theOnlyHeader.h"

int main(){
    printf("Hello, world!\n");

    simon_jet J;
    
    J.nPart = 7;
    J.sumpt = 7;
    J.pt = 7;
    J.rawpt = 7;

    J.particles.emplace_back(1, 0, 0);
    J.particles.emplace_back(1, 0, 1);
    J.particles.emplace_back(1, 0.75, 0);
    J.particles.emplace_back(1, 1.1, 0.8);
    J.particles.emplace_back(1, -0.1, 0.7);
    J.particles.emplace_back(1, -0.2, -0.1);
    J.particles.emplace_back(1, +0.1, -0.6);

    standaloneEEC::axis ax_R(std::vector<double>({{
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
    }}));
    standaloneEEC::axis ax_r(std::vector<double>({{
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
        0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
        0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00
    }}));
    standaloneEEC::axis ax_c(std::vector<double>({{
        0.00*M_PI, 0.05*M_PI, 0.10*M_PI,
        0.15*M_PI, 0.20*M_PI, 0.25*M_PI, 0.30*M_PI,
        0.35*M_PI, 0.40*M_PI, 0.45*M_PI,
        0.50*M_PI, 0.55*M_PI, 0.60*M_PI, 0.65*M_PI,
        0.70*M_PI, 0.75*M_PI, 0.80*M_PI, 
        0.85*M_PI, 0.90*M_PI, 0.95*M_PI, 1.00*M_PI

    }}));

    standaloneEEC::res4_result_multi_array result(8, 
                                                  8, 8,
                                                  8, 8,
                                                  8, 8);

    standaloneEEC::res4_result_vector result_vec(8, 
                                                  8, 8,
                                                  8, 8,
                                                  8, 8);

    standaloneEEC::res4_standalone_multi_array(
            result, 
            J, standaloneEEC::normType::RAWPT,
            ax_R, 
            ax_r, ax_c,
            ax_r, ax_c,
            ax_r, ax_c,
            0.1, 0.1);
    printf("ran multi_array\n");

    standaloneEEC::res4_standalone_vector(
            result_vec, 
            J, standaloneEEC::normType::RAWPT,
            ax_R, 
            ax_r, ax_c,
            ax_r, ax_c,
            ax_r, ax_c,
            0.1, 0.1);
    printf("ran vector\n");

    standaloneEEC::res4_standalone_multi_array_precomputed(
            result, 
            J, standaloneEEC::normType::RAWPT,
            ax_R, 
            ax_r, ax_c,
            ax_r, ax_c,
            ax_r, ax_c,
            0.1, 0.1);
    printf("ran multi_arry precomputed\n");

    standaloneEEC::res4_standalone_vector_precomputed(
            result_vec, 
            J, standaloneEEC::normType::RAWPT,
            ax_R, 
            ax_r, ax_c,
            ax_r, ax_c,
            ax_r, ax_c,
            0.1, 0.1);
    printf("ran vector precomputed\n");

    printf("START BENCHMARK\n");

    std::default_random_engine rng(0);
    std::normal_distribution<double> norm(0.0, 0.4);
    std::gamma_distribution<double> gam(2.0, 2.0);

    constexpr size_t NUM_TRIALS = 100;
    constexpr size_t NPART = 80;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t trial=0; trial<NUM_TRIALS; ++trial){
        simon_jet J;
        J.nPart = NPART;
        J.eta = 0;
        J.phi = 0;
        for (unsigned i=0; i<J.nPart; ++i){
            double eta = norm(rng);
            double phi = norm(rng);
            double pt = gam(rng);
            J.particles.emplace_back(pt, eta, phi);
            J.sumpt += pt;
            J.pt += pt;
            J.rawpt += pt;
        }

        standaloneEEC::res4_result_multi_array result(
            12, 22, 22, 22, 22, 22, 22
        );

        //standaloneEEC::res4_standalone_multi_array(
        //        result,
        //        J, standaloneEEC::normType::RAWPT,
        //        ax_R,
        //        ax_r, ax_c,
        //        ax_r, ax_c,
        //        ax_r, ax_c,
        //        0.1, 0.1);
       
        standaloneEEC::res4_standalone_multi_array_precomputed(
                result,
                J, standaloneEEC::normType::RAWPT,
                ax_R,
                ax_r, ax_c,
                ax_r, ax_c,
                ax_r, ax_c,
                0.1, 0.1);
       
        //standaloneEEC::res4_result_vector result(
        //    12, 22, 22, 22, 22, 22, 22
        //);

        //standaloneEEC::res4_standalone_vector(
        //        result,
        //        J, standaloneEEC::normType::RAWPT,
        //        ax_R,
        //        ax_r, ax_c,
        //        ax_r, ax_c,
        //        ax_r, ax_c,
        //        0.1, 0.1);
       
        //standaloneEEC::res4_standalone_vector_precomputed(
        //        result,
        //        J, standaloneEEC::normType::RAWPT,
        //        ax_R,
        //        ax_r, ax_c,
        //        ax_r, ax_c,
        //        ax_r, ax_c,
        //        0.1, 0.1);
       
        //fastEEC::axisptr ax_R_ptr = std::make_shared<fastEEC::axis_t>(ax_R);
        //fastEEC::axisptr ax_r_ptr = std::make_shared<fastEEC::axis_t>(ax_r);
        //fastEEC::axisptr ax_c_ptr = std::make_shared<fastEEC::axis_t>(ax_c);

        //fastEEC::result_t<double> result;
        //runFastEEC(
        //    result,
        //    J, 
        //    ax_R_ptr, 
        //    fastEEC::normType::RAWPT,
        //    4,
        //    fastEEC::DORES4,

        //    ax_R_ptr, 
        //    ax_R_ptr, ax_R_ptr,
        //    ax_r_ptr, ax_c_ptr,
        //    ax_r_ptr, ax_c_ptr,
        //    ax_r_ptr, ax_c_ptr,
        //    ax_R_ptr, ax_R_ptr,
        //    ax_R_ptr, ax_R_ptr,

        //    0.1
        //); 

    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    printf("multi_array precomputed: %lu ms\n", duration);

    return 0;
}
