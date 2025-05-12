#ifndef TESTPROJCALCULATOR_H
#define TESTPROJCALCULATOR_H

#include "SRothman/SimonTools/src/JetFactory.h"
#include <chrono>
#include <random>
#include <Eigen/Dense>

#include "SRothman/EECs/src/ProjResult.h"
#include "SRothman/EECs/src/ProjTransferResult.h"
#include "SRothman/EECs/src/ProjCalculator.h"

#ifdef CHECK_BY_HAND
constexpr int NPART = 10;
#else
constexpr int NPART = 70;
#endif

constexpr int NITER = 10;
constexpr int NBINS = 20;

inline void checkSum(const EEC::ProjCalculator& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    EEC::ProjResult_Array result(calculator);
    calculator.compute(J, result);

    for (unsigned order = 2; order <= 6; ++order){
        printf("Sum(EEC) for order %u = %g\n", order, result.total_weight(order));
    }
}

inline void checkSum_matched(const EEC::ProjCalculator& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    EEC::ProjResult_Array result(calculator);
    EEC::ProjResult_Array unmatched(calculator);

    std::vector<bool> matched(NPART, true);
    for ( unsigned i = 0; i < NPART; ++i){
        matched[i] = (i % 5 != 0);
    }

    calculator.compute_matched(J, matched, result, unmatched);
    for (unsigned order = 2; order <= 6; ++order){
        printf("Order %u\n", order);
        printf("\tTotal = %g\n", result.total_weight(order));
        printf("\tUnmatched = %g\n", unmatched.total_weight(order));
    }
}

inline void checkSum_transfer(const EEC::ProjCalculator& calculator,
                              const EEC::ProjTransferCalculator& tcalculator,
                              JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J_gen;
    jetFactory.makeJet(J_gen, NPART);

    simon::jet J_reco;
    Eigen::MatrixXd tmat;
    jetFactory.makeTransferJet(J_gen, J_reco, tmat);

    EEC::ProjResult_Vector gen(tcalculator.get_axes_gen());
    EEC::ProjResult_Vector unmatched_gen(tcalculator.get_axes_gen());
    EEC::ProjResult_Vector reco(tcalculator.get_axes_reco());
    EEC::ProjResult_Vector unmatched_reco(tcalculator.get_axes_reco());
    EEC::ProjTransferResult_Vector tresult(tcalculator);

    std::vector<bool> matched_reco(NPART, false);
    for (unsigned iReco=0; iReco<J_reco.nPart; ++iReco){
        for (unsigned iGen=0; iGen<J_gen.nPart; ++iGen){
            if (tmat(iReco, iGen) > 0){
                matched_reco[iReco] = true;
                break;
            }
        }
    }

    calculator.compute_matched(J_reco, matched_reco, reco, unmatched_reco);
    tcalculator.compute(J_reco, J_gen, tmat, gen, unmatched_gen, tresult);
    
    auto sum_over_gen = tresult.get_sum_over_gen();
    auto sum_over_reco = tresult.get_sum_over_reco();

    auto RHS = reco - unmatched_reco;
    printf("sum_over_gen == reco - unmatched_reco? %s\n",
           (sum_over_gen == RHS) ? "yes" : "no");
    for (unsigned order = 2; order <= 6; ++order){
        printf("Order %u\n", order);
        printf("\ttotal reco = %g\n", reco.total_weight(order));
        printf("\ttotal unmatched reco = %g\n", unmatched_reco.total_weight(order));
        printf("\tsum over gen = %g\n", sum_over_gen.total_weight(order));
    }
    auto RHS2 = gen - unmatched_gen;
    printf("sum_over_reco == gen - unmatched_gen? %s\n",
           (sum_over_reco == RHS2) ? "yes" : "no");
    for (unsigned order = 2; order <= 6; ++order){
        printf("Order %u\n", order);
        printf("\ttotal gen = %g\n", gen.total_weight(order));
        printf("\ttotal unmatched gen = %g\n", unmatched_gen.total_weight(order));
        printf("\tsum over reco = %g\n", sum_over_reco.total_weight(order));
    }
}

inline void checkAllSame(const EEC::ProjCalculator& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    EEC::ProjResult_Array result(calculator);
    EEC::ProjResult_Vector result2(calculator);

    calculator.compute(J, result);
    calculator.compute(J, result2);

    bool allSame = result == result2;
    printf("All same? %s\n", allSame ? "yes" : "no");
}

inline void checkAllSame_matched(const EEC::ProjCalculator& calculator, JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J;
    jetFactory.makeJet(J, NPART);

    EEC::ProjResult_Array result(calculator);
    EEC::ProjResult_Vector result2(calculator);

    std::vector<bool> matched(NPART, true);
    for ( unsigned i = 0; i < NPART; ++i){
        matched[i] = (i % 5 != 0);
    }

    EEC::ProjResult_Array unmatched(calculator);
    EEC::ProjResult_Vector unmatched2(calculator);

    calculator.compute_matched(J, matched, result, unmatched);
    calculator.compute_matched(J, matched, result2, unmatched2);

    bool allSame = result == result2 &&
                     unmatched == unmatched2;
    printf("All same? %s\n", allSame ? "yes" : "no");
}

inline void checkAllSame_transfer(const EEC::ProjTransferCalculator& tcalculator,
                                  JetFactory& jetFactory){
    jetFactory.initialize();
    simon::jet J_gen;
    jetFactory.makeJet(J_gen, NPART);

    simon::jet J_reco;
    Eigen::MatrixXd tmat;
    jetFactory.makeTransferJet(J_gen, J_reco, tmat);

    EEC::ProjResult_Array gen_array(tcalculator.get_axes_gen());
    EEC::ProjResult_Array unmatched_array(tcalculator.get_axes_gen());
    EEC::ProjTransferResult_Array tresult_array(tcalculator);

    EEC::ProjResult_Vector gen_vector(tcalculator.get_axes_gen());
    EEC::ProjResult_Vector unmatched_vector(tcalculator.get_axes_gen());
    EEC::ProjTransferResult_Vector tresult_vector(tcalculator);

        tcalculator.compute(J_reco, J_gen, tmat, 
                gen_array, unmatched_array, tresult_array);
    tcalculator.compute(J_reco, J_gen, tmat,
            gen_vector, unmatched_vector, tresult_vector);
    bool allSame = gen_array == gen_vector &&
                   unmatched_array == unmatched_vector &&
                   tresult_array == tresult_vector;
    printf("All same? %s\n", allSame ? "yes" : "no");
}

template <class Container>
inline unsigned long benchmarkCalculator(const EEC::ProjCalculator& calculator,
                                  JetFactory& jetFactory){
    jetFactory.initialize();

    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < NITER; ++i){
        simon::jet J;
        jetFactory.makeJet(J, NPART);

        EEC::ProjResult<Container> result(calculator);
        calculator.compute(J, result);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
}

inline unsigned long benchmarkCalculator_matched(
        const EEC::ProjCalculator& calculator,
        JetFactory& jetFactory){

    jetFactory.initialize();

    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < NITER; ++i){
        simon::jet J;
        jetFactory.makeJet(J, NPART);

        EEC::ProjResult_Array result(calculator);
        EEC::ProjResult_Array unmatched(calculator);

        std::vector<bool> matched(NPART, true);
        for ( unsigned i = 0; i < NPART; ++i){
            matched[i] = (i % 5 != 0);
        }

        calculator.compute_matched(J, matched, result, unmatched);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
}

inline unsigned long benchmarkCalculator_transfer(){
    //TODO
    return 0;
}

#endif
