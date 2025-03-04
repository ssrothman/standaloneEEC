#ifndef TESTCALCULATOR_H
#define TESTCALCULATOR_H

#include "SRothman/SimonTools/src/JetFactory.h"
#include <chrono>
#include <random>
#include <Eigen/Dense>

#include "SRothman/EECs/src/Res3Result.h"
#include "SRothman/EECs/src/Res4Result.h"

#include "SRothman/EECs/src/ResResultContainers.h"
#include "SRothman/EECs/src/ResTransferResultContainers.h"

#include "SRothman/EECs/src/Res3Calculator.h"
#include "SRothman/EECs/src/Res4Calculator.h"
#include "SRothman/EECs/src/CARes3Calculator.h"
#include "SRothman/EECs/src/CARes4Calculator.h"

#ifndef CHECK_BY_HAND
constexpr int NPART = 10;
#else
constexpr int NPART = 50;
#endif

constexpr int NITER = 10;
constexpr int NBINS = 20;

template <typename RESULT, typename CALCULATOR>
void checkSum(const CALCULATOR& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    RESULT result(
        calculator
    );
    calculator.compute_JIT(J, result);

    if constexpr(std::is_same_v<CALCULATOR, EEC::Res3Calculator>){
        printf("Total weight = %g\n", result.total_weight());
    } else if constexpr(std::is_same_v<CALCULATOR, EEC::Res4Calculator>){
        printf("Total tee weight = %g\n", result.total_tee_weight());
        printf("Total dipole weight = %g\n", result.total_dipole_weight());
        printf("Total triangle weight = %g\n", result.total_triangle_weight());
    } else if constexpr(std::is_same_v<CALCULATOR, EEC::CARes3Calculator>){
        printf("Total weight = %g\n", result.total_weight());
    } else if constexpr(std::is_same_v<CALCULATOR, EEC::CARes4Calculator>){
        printf("Total chain weight = %g\n", result.total_chain_weight());
        printf("Total symmetric_wrtR weight = %g\n", result.total_symmetric_wrtR_weight());
        printf("Total symmetric_wrtr weight = %g\n", result.total_symmetric_wrtr_weight());
    }
}

template <typename RESULT, typename CALCULATOR>
void checkSum_matched(const CALCULATOR& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);
    std::vector<bool> matched(NPART);
    jetFactory.makeMatchedVec(J, matched, 0.5);

    RESULT result(
        calculator
    );
    RESULT result_unmatched(
        calculator
    );
    calculator.compute_JIT_matched(J, matched, result, result_unmatched);

    if constexpr(std::is_same_v<CALCULATOR, EEC::Res3Calculator>){
        printf("Total weight = %g\n", result.total_weight());
        printf("Total weight unmatched = %g\n", result_unmatched.total_weight());
    } else if constexpr(std::is_same_v<CALCULATOR, EEC::Res4Calculator>){
        printf("Total tee weight = %g\n", result.total_tee_weight());
        printf("\tunmatched = %g\n", result_unmatched.total_tee_weight());
        printf("Total dipole weight = %g\n", result.total_dipole_weight());
        printf("\tunmatched = %g\n", result_unmatched.total_dipole_weight());
        printf("Total triangle weight = %g\n", result.total_triangle_weight());
        printf("\tunmatched = %g\n", result_unmatched.total_triangle_weight());
    } else if constexpr(std::is_same_v<CALCULATOR, EEC::CARes3Calculator>){
        printf("Total weight = %g\n", result.total_weight());
        printf("Total weight unmatched = %g\n", result_unmatched.total_weight());
    } else if constexpr(std::is_same_v<CALCULATOR, EEC::CARes4Calculator>){
        printf("Total chain weight = %g\n", result.total_chain_weight());
        printf("\tunmatched = %g\n", result_unmatched.total_chain_weight());
        printf("Total symmetric_wrtR weight = %g\n", result.total_symmetric_wrtR_weight());
        printf("\tunmatched = %g\n", result_unmatched.total_symmetric_wrtR_weight());
        printf("Total symmetric_wrtr weight = %g\n", result.total_symmetric_wrtr_weight());
        printf("\tunmatched = %g\n", result_unmatched.total_symmetric_wrtr_weight());
    }
}

template <typename RESULT, typename TRESULT, bool has_untransferred, typename CALCULATOR, typename TCALCULATOR>
void check_sum_transfer(const CALCULATOR& calculator, 
                        const TCALCULATOR& tcalculator,
                        JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J_gen;
    jetFactory.makeJet(J_gen, NPART);

    simon::jet J_reco;
    Eigen::MatrixXd tmat;
    jetFactory.makeTransferJet(J_gen, J_reco, tmat);

    RESULT gen(
        tcalculator.get_axes_gen()
    );
    RESULT unmatched_gen(
        tcalculator.get_axes_gen()
    );
    TRESULT tresult(
        tcalculator
    );
    RESULT untransferred_gen(
        tcalculator.get_axes_gen()
    );
    RESULT untransferred_reco(
        tcalculator.get_axes_reco()
    );

    RESULT reco(
        tcalculator.get_axes_reco()
    );
    RESULT reco_unmatched(
        tcalculator.get_axes_reco()
    );
    std::vector<bool> matched_reco(NPART, false);
    for (unsigned iReco=0; iReco<J_reco.nPart; ++iReco){
        for (unsigned iGen=0; iGen<J_gen.nPart; ++iGen){
            if (tmat(iReco, iGen) > 0){
                matched_reco[iReco] = true;
                break;
            }
        }
    }
    calculator.compute_precomputed_matched(
        J_reco, matched_reco,
        reco, reco_unmatched
    );

    if constexpr(has_untransferred){
        tcalculator.compute_precomputed(
            J_reco, J_gen, tmat,
            gen, unmatched_gen, tresult,
            untransferred_reco, untransferred_gen
        );
    } else {
        tcalculator.compute_precomputed(
            J_reco, J_gen, tmat,
            gen, unmatched_gen, tresult
        );
    }

    printf("total transfer chain: %g\n", tresult.total_chain_weight_reco());

    auto sum_over_gen = tresult.get_sum_over_gen();
    bool pass = sum_over_gen + untransferred_reco == reco - reco_unmatched;
    printf("sum_over_gen + untransferred_reco == reco - reco_unmatched? %s\n", pass ? "true" : "false");
    if constexpr(std::is_same_v<CALCULATOR, EEC::CARes4Calculator>){
        bool pass_chain = sum_over_gen.get_chain() + untransferred_reco.get_chain() == reco.get_chain() - reco_unmatched.get_chain();
        bool pass_symmetric_wrtR = sum_over_gen.get_symmetric_wrtR() + untransferred_reco.get_symmetric_wrtR() == reco.get_symmetric_wrtR() - reco_unmatched.get_symmetric_wrtR();
        bool pass_symmetric_wrtr = sum_over_gen.get_symmetric_wrtr() + untransferred_reco.get_symmetric_wrtr() == reco.get_symmetric_wrtr() - reco_unmatched.get_symmetric_wrtr();
        printf("\tchain: %s\n", pass_chain ? "true" : "false");
        printf("\t\ttotal_sum_over_gen = %g\n", sum_over_gen.total_chain_weight());
        printf("\t\ttotal_chain_to_chain = %g\n", tresult.total_chain_weight_reco());
        printf("\t\ttotal_symmetric_wrtR_to_chain = %g\n", tresult.total_symmetric_wrtR_to_chain_weight_reco());
        printf("\t\tuntransferred_reco = %g\n", untransferred_reco.total_chain_weight());
        printf("\t\treco = %g\n", reco.total_chain_weight());
        printf("\t\treco_unmatched = %g\n", reco_unmatched.total_chain_weight());
        printf("\tsymmetric_wrtR: %s\n", pass_symmetric_wrtR ? "true" : "false");
        printf("\t\ttotal_transfered = %g\n", sum_over_gen.total_symmetric_wrtR_weight());
        printf("\t\tuntransferred_reco = %g\n", untransferred_reco.total_symmetric_wrtR_weight());
        printf("\t\treco = %g\n", reco.total_symmetric_wrtR_weight());
        printf("\t\treco_unmatched = %g\n", reco_unmatched.total_symmetric_wrtR_weight());
        printf("\tsymmetric_wrtr: %s\n", pass_symmetric_wrtr ? "true" : "false");
        printf("\t\ttotal_transfered = %g\n", sum_over_gen.total_symmetric_wrtr_weight());
        printf("\t\tuntransferred_reco = %g\n", untransferred_reco.total_symmetric_wrtr_weight());
        printf("\t\treco = %g\n", reco.total_symmetric_wrtr_weight());
        printf("\t\treco_unmatched = %g\n", reco_unmatched.total_symmetric_wrtr_weight());
    }
}

template <template <typename CONTAINER> typename RESULT_T, typename CALCULATOR>
void checkAllSame(const CALCULATOR& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    RESULT_T<EEC::ResVectorContainer> result_vec_JIT(
        calculator
    );

    RESULT_T<EEC::ResVectorContainer> result_vec_precomputed(
        calculator
    );

    RESULT_T<EEC::ResMultiArrayContainer> result_arr_JIT(
        calculator
    );
    RESULT_T<EEC::ResMultiArrayContainer> result_arr_precomputed(
        calculator
    );

    calculator.compute_JIT(J, result_vec_JIT);
    calculator.compute_precomputed(J, result_vec_precomputed);
    calculator.compute_JIT(J, result_arr_JIT);
    calculator.compute_precomputed(J, result_arr_precomputed);

    printf("vec == vec_precomputd? %s\n", result_vec_JIT == result_vec_precomputed ? "true" : "false");
    printf("arr == arr_precomputd? %s\n", result_arr_JIT == result_arr_precomputed ? "true" : "false");
    printf("vec == arr? %s\n", result_vec_JIT == result_arr_JIT ? "true" : "false");
}

template <template <typename CONTAINER> typename RESULT_T, typename CALCULATOR>
void checkAllSame_matched(const CALCULATOR& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;   
    jetFactory.makeJet(J, NPART);
    std::vector<bool> matched(NPART);
    jetFactory.makeMatchedVec(J, matched, 0.5);


    RESULT_T<EEC::ResVectorContainer> result_vec_JIT(
        calculator
    );
    RESULT_T<EEC::ResVectorContainer> unmatched_vec_JIT(
        calculator
    );
    RESULT_T<EEC::ResVectorContainer> result_vec_precomputed(
        calculator
    );
    RESULT_T<EEC::ResVectorContainer> unmatched_vec_precomputed(
        calculator
    );
    RESULT_T<EEC::ResMultiArrayContainer> result_arr_JIT(
        calculator
    );
    RESULT_T<EEC::ResMultiArrayContainer> unmatched_arr_JIT(
        calculator
    );
    RESULT_T<EEC::ResMultiArrayContainer> result_arr_precomputed(
        calculator
    );
    RESULT_T<EEC::ResMultiArrayContainer> unmatched_arr_precomputed(
        calculator
    );

    calculator.compute_JIT_matched(J, matched, result_vec_JIT, unmatched_vec_JIT);
    calculator.compute_precomputed_matched(J, matched, result_vec_precomputed, unmatched_vec_precomputed);
    calculator.compute_JIT_matched(J, matched, result_arr_JIT, unmatched_arr_JIT);
    calculator.compute_precomputed_matched(J, matched, result_arr_precomputed, unmatched_arr_precomputed);

    bool same_vec = result_vec_JIT == result_vec_precomputed;
    same_vec &= unmatched_vec_JIT == unmatched_vec_precomputed;

    bool same_arr = result_arr_JIT == result_arr_precomputed;
    same_arr &= unmatched_arr_JIT == unmatched_arr_precomputed;

    bool same = result_vec_JIT == result_arr_JIT;
    same &= unmatched_vec_JIT == unmatched_arr_JIT;

    printf("vec == vec_precomputd? %s\n", same_vec ? "true" : "false");
    printf("arr == arr_precomputd? %s\n", same_arr ? "true" : "false");
    printf("vec == arr? %s\n", same ? "true" : "false");
}

template <template <typename CONTAINER> typename RESULT_T,
          template <typename CONATINER> typename TRESULT_T,
          bool has_untransferred,
          typename CALCULATOR>
void checkAllSame_transfer(const CALCULATOR& calculator,
                           JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J_gen;
    jetFactory.makeJet(J_gen, NPART);

    simon::jet J_reco;
    Eigen::MatrixXd tmat;
    jetFactory.makeTransferJet(J_gen, J_reco, tmat);

    RESULT_T<EEC::ResVectorContainer> result_vec_JIT(
        calculator.get_axes_gen()
    );
    RESULT_T<EEC::ResVectorContainer> unmatched_vec_JIT(
        calculator.get_axes_gen()
    );
    TRESULT_T<EEC::ResTransferVectorContainer> tresult_vec_JIT(
        calculator
    );

    RESULT_T<EEC::ResVectorContainer> result_vec_precomputed(
        calculator.get_axes_gen()
    );
    RESULT_T<EEC::ResVectorContainer> unmatched_vec_precomputed(
        calculator.get_axes_gen()
    );
    TRESULT_T<EEC::ResTransferVectorContainer> tresult_vec_precomputed(
        calculator
    );

    RESULT_T<EEC::ResMultiArrayContainer> result_arr_JIT(
        calculator.get_axes_gen()
    );
    RESULT_T<EEC::ResMultiArrayContainer> unmatched_arr_JIT(
        calculator.get_axes_gen()
    );
    TRESULT_T<EEC::ResTransferMultiArrayContainer> tresult_arr_JIT(
        calculator
    );

    RESULT_T<EEC::ResMultiArrayContainer> result_arr_precomputed(
        calculator.get_axes_gen()
    );
    RESULT_T<EEC::ResMultiArrayContainer> unmatched_arr_precomputed(
        calculator.get_axes_gen()
    );
    TRESULT_T<EEC::ResTransferMultiArrayContainer> tresult_arr_precomputed(
        calculator
    );

    if constexpr (has_untransferred){
        RESULT_T<EEC::ResVectorContainer> ut_reco_vec_JIT(
            calculator.get_axes_reco()
        );
        RESULT_T<EEC::ResVectorContainer> ut_gen_vec_JIT(
            calculator.get_axes_gen()
        );

        RESULT_T<EEC::ResVectorContainer> ut_reco_vec_precomputed(
            calculator.get_axes_reco()
        );
        RESULT_T<EEC::ResVectorContainer> ut_gen_vec_precomputed(
            calculator.get_axes_gen()
        );

        RESULT_T<EEC::ResMultiArrayContainer> ut_reco_arr_JIT(
            calculator.get_axes_reco()
        );
        RESULT_T<EEC::ResMultiArrayContainer> ut_gen_arr_JIT(
            calculator.get_axes_gen()
        );

        RESULT_T<EEC::ResMultiArrayContainer> ut_reco_arr_precomputed(
            calculator.get_axes_reco()
        );
        RESULT_T<EEC::ResMultiArrayContainer> ut_gen_arr_precomputed(
            calculator.get_axes_gen()
        );

        calculator.compute_JIT(
            J_reco, J_gen, tmat,
            result_vec_JIT,
            unmatched_vec_JIT,
            tresult_vec_JIT,
            ut_reco_vec_JIT,
            ut_gen_vec_JIT
        );

        calculator.compute_precomputed(
            J_reco, J_gen, tmat,
            result_vec_precomputed,
            unmatched_vec_precomputed,
            tresult_vec_precomputed,
            ut_reco_vec_precomputed,
            ut_gen_vec_precomputed
        );

        calculator.compute_JIT(
            J_reco, J_gen, tmat,
            result_arr_JIT,
            unmatched_arr_JIT,
            tresult_arr_JIT,
            ut_reco_arr_JIT,
            ut_gen_arr_JIT
        );

        calculator.compute_precomputed(
            J_reco, J_gen, tmat,
            result_arr_precomputed,
            unmatched_arr_precomputed,
            tresult_arr_precomputed,
            ut_reco_arr_precomputed,
            ut_gen_arr_precomputed
        );

        bool same_vec = result_vec_JIT == result_vec_precomputed;
        same_vec &= unmatched_vec_JIT == unmatched_vec_precomputed;
        same_vec &= tresult_vec_JIT == tresult_vec_precomputed;
        same_vec &= ut_reco_vec_JIT == ut_reco_vec_precomputed;
        same_vec &= ut_gen_vec_JIT == ut_gen_vec_precomputed;

        bool same_arr = result_arr_JIT == result_arr_precomputed;
        same_arr &= unmatched_arr_JIT == unmatched_arr_precomputed;
        same_arr &= tresult_arr_JIT == tresult_arr_precomputed;
        same_arr &= ut_reco_arr_JIT == ut_reco_arr_precomputed;
        same_arr &= ut_gen_arr_JIT == ut_gen_arr_precomputed;

        bool same = result_vec_JIT == result_arr_JIT;
        same &= unmatched_vec_JIT == unmatched_arr_JIT;
        same &= tresult_vec_JIT == tresult_arr_JIT;
        same &= ut_reco_vec_JIT == ut_reco_arr_JIT;
        same &= ut_gen_vec_JIT == ut_gen_arr_JIT;

        printf("vec == vec_precomputd? %s\n", same_vec ? "true" : "false");
        printf("arr == arr_precomputd? %s\n", same_arr ? "true" : "false");
        printf("vec == arr? %s\n", same ? "true" : "false");
    } else {
        calculator.compute_JIT(
            J_reco, J_gen, tmat,
            result_vec_JIT,
            unmatched_vec_JIT,
            tresult_vec_JIT
        );

        calculator.compute_precomputed(
            J_reco, J_gen, tmat,
            result_vec_precomputed,
            unmatched_vec_precomputed,
            tresult_vec_precomputed
        );

        calculator.compute_JIT(
            J_reco, J_gen, tmat,
            result_arr_JIT,
            unmatched_arr_JIT,
            tresult_arr_JIT
        );

        calculator.compute_precomputed(
            J_reco, J_gen, tmat,
            result_arr_precomputed,
            unmatched_arr_precomputed,
            tresult_arr_precomputed
        );

        bool same_vec = result_vec_JIT == result_vec_precomputed;
        same_vec &= unmatched_vec_JIT == unmatched_vec_precomputed;
        same_vec &= tresult_vec_JIT == tresult_vec_precomputed;

        bool same_arr = result_arr_JIT == result_arr_precomputed;
        same_arr &= unmatched_arr_JIT == unmatched_arr_precomputed;
        same_arr &= tresult_arr_JIT == tresult_arr_precomputed;

        bool same = result_vec_JIT == result_arr_JIT;
        same &= unmatched_vec_JIT == unmatched_arr_JIT;
        same &= tresult_vec_JIT == tresult_arr_JIT;

        printf("vec == vec_precomputd? %s\n", same_vec ? "true" : "false");
        printf("arr == arr_precomputd? %s\n", same_arr ? "true" : "false");
        printf("vec == arr? %s\n", same ? "true" : "false");
    }

}

template <typename RESULT, typename CALCULTOR, bool precomputed>
void runCalculator(RESULT& result, const CALCULTOR& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    if constexpr (precomputed){
        calculator.compute_precomputed(J, result);
    } else {
        calculator.compute_JIT(J, result);
    }
    return result;
}

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
void runCalculator_matched(RESULT& result, RESULT& result_unmatched, const CALCULATOR& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);
    std::vector<bool> matched(NPART);
    jetFactory.makeMatchedVec(J, matched, 0.5);

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
    return result;
}

template <typename RESULT, typename CALCULATOR, bool precomputed>
unsigned long benchmarkCalculator_unmatched(const CALCULATOR& calculator,
                                     JetFactory& jetFactory){

    jetFactory.initialize();


    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i=0; i<NITER; ++i){
        simon::jet J;
        jetFactory.makeJet(J, NPART);
        std::vector<bool> matched(NPART);
        jetFactory.makeMatchedVec(J, matched, 0.5);

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
         typename RESULT_TRANFSER,
         typename CALCULATOR,
         bool has_untransferred,
         bool precomputed>
void runCalculator_transfer(RESULT& result,
                            RESULT& unmatched_gen,
                            RESULT_TRANFSER& tresult,
                            RESULT* untransferred_gen,
                            RESULT* untransferred_reco,
                            const CALCULATOR& calculator,
                            JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J_gen;
    jetFactory.makeJet(J_gen, NPART);

    simon::jet J_reco;
    Eigen::MatrixXd tmat;
    jetFactory.makeTransferJet(J_gen, J_reco, tmat); 

    if constexpr (has_untransferred){
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
