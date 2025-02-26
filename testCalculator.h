#ifndef TESTCALCULATOR_H
#define TESTCALCULATOR_H

#include "SRothman/SimonTools/src/JetFactory.h"
#include <chrono>
#include <random>
#include <Eigen/Dense>

constexpr int NPART = 60;
constexpr int NITER = 10;
constexpr int NBINS = 20;

template <typename RESULT, typename CALCULATOR, bool precomputed>
unsigned long benchmarkCalculator(const CALCULATOR& calculator,
                   JetFactory& jetFactory){

    jetFactory.initialize();

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i=0; i<NITER; ++i){
        simon::jet J;
        jetFactory.makeJet(J, NPART);

        RESULT result(
            calculator
        );
        if constexpr (precomputed){
            calculator.compute_precomputed(J, result);
        } else {
            calculator.compute_JIT(J, result);
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
}

template <typename RESULT, typename CALCULATOR, bool precomputed>
unsigned long benchmarkCalculator_unmatched(const CALCULATOR& calculator,
                                     JetFactory& jetFactory){

    jetFactory.initialize();


    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i=0; i<NITER; ++i){
        simon::jet J;
        jetFactory.makeJet(J, NPART);
        std::vector<bool> matched(NPART, false);

        RESULT result(
            calculator
        );
        RESULT result_unmatched(
            calculator
        );

        if constexpr (precomputed){
            calculator.compute_precomputed_matched(
                    J, matched,
                    result, 
                    result_unmatched);
        } else {
            calculator.compute_JIT_matched(
                    J, matched,
                    result,
                    result_unmatched);
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
}

template <typename RESULT,
          typename RESULT_TRANSFER, 
          typename CALCULATOR, 
          bool has_untransferred,
          bool precomputed>
unsigned long benchmarkCalculator_transfer(const CALCULATOR& calculator,
                   JetFactory& jetFactory){

    jetFactory.initialize();

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i=0; i<NITER; ++i){
        simon::jet J_gen;
        jetFactory.makeJet(J_gen, NPART);

        simon::jet J_reco;
        Eigen::MatrixXd tmat;
        jetFactory.makeTransferJet(J_gen, J_reco, tmat); 

        RESULT result(
            calculator.get_axes_gen()
        );
        RESULT unmatched_gen(
            calculator.get_axes_gen()
        );
        RESULT_TRANSFER tresult(
            calculator
        );
        if constexpr (has_untransferred){
            RESULT untransferred_gen(
                calculator.get_axes_gen()
            );
            RESULT untransferred_reco(
                    calculator.get_axes_reco()
            );
            if constexpr (precomputed){
                calculator.compute_precomputed(
                        J_reco, J_gen,
                        tmat,
                        result,
                        unmatched_gen,
                        tresult,
                        untransferred_reco,
                        untransferred_gen);
            } else {
                calculator.compute_JIT(
                        J_reco, J_gen,
                        tmat,
                        result,
                        unmatched_gen,
                        tresult,
                        untransferred_reco,
                        untransferred_gen);
            }
        } else {
            if constexpr (precomputed){
                calculator.compute_precomputed(
                        J_reco, J_gen,
                        tmat,
                        result,
                        unmatched_gen,
                        tresult);
            } else {
                calculator.compute_JIT(
                        J_reco, J_gen,
                        tmat,
                        result,
                        unmatched_gen,
                        tresult);
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
}


#endif
