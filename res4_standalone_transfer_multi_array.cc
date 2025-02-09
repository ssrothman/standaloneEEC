#include <stdio.h>
#include <chrono>
#include <random>
#include <boost/histogram.hpp>
#include <Eigen/Dense>

#include "SRothman/EECs/src/standalones/res4_standalone.h"
#include "SRothman/SimonTools/src/simon_jet.h"
#include "SRothman/EECs/src/theOnlyHeader.h"

constexpr int NPART = 50;
constexpr int NITER = 10;
constexpr int NBINS = 20;

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

    void makeJet(simon_jet& J, const int nPart){
        J.nPart = nPart;

        J.sumpt = 0;
        J.pt = 0;
        J.rawpt = 0;

        J.eta = 0;
        J.phi = 0;

        J.particles.clear();
        for(int i=0; i<NPART; ++i){
            double pt = gamma(rng);
            double phi = norm(rng);
            double eta = norm(rng);
            J.particles.emplace_back(pt, phi, eta);
            J.sumpt += pt;
            J.pt += pt;
            J.rawpt += pt;
        }
    }

    void makeTransferJet(const simon_jet& J1, simon_jet& J2, Eigen::MatrixXd& tmat){
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

    JetFactory jetFactory;
    jetFactory.initialize();

    std::vector<double> bins(NBINS-1);
    for (int i=0; i<NBINS-1; ++i){
        bins[i] = double(i)/(NBINS-2);
    }

    standaloneEEC::axis ax(bins);
    auto ax_ptr = std::make_shared<standaloneEEC::axis>(ax);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i=0; i<NITER; ++i){
        simon_jet J1, J2;
        Eigen::MatrixXd tmat(NPART, NPART);
        jetFactory.makeJet(J1, NPART);
        jetFactory.makeTransferJet(J1, J2, tmat);

        standaloneEEC::res4_result_multi_array result(
            NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS
        );
        standaloneEEC::res4_transfer_result_multi_array transfer_result(
            NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS,
            NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS
        );

        standaloneEEC::res4_standalone_transfer_multi_array(
                result,
                transfer_result,
                J1, J2,
                standaloneEEC::normType::RAWPT,
                tmat,
                ax,
                ax, ax,
                ax, ax,
                ax, ax,
                0.1, 0.1);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    auto deltat = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    printf("standalone multi array: %lu ms\n", deltat);
}
