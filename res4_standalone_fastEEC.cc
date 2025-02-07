#include <stdio.h>
#include <chrono>
#include <random>
#include <boost/histogram.hpp>

#include "SRothman/EECs/src/standalones/res4_standalone.h"
#include "SRothman/SimonTools/src/simon_jet.h"
#include "SRothman/EECs/src/theOnlyHeader.h"

constexpr int NPART = 70;
constexpr int NITER = 50;
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

    void makeJet(simon_jet& J, int nPart){
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
        simon_jet J;
        jetFactory.makeJet(J, NPART);

        fastEEC::result_t<double> result;
        runFastEEC(
            result,
            J,
            ax_ptr,
            fastEEC::normType::RAWPT,
            4,
            fastEEC::DORES4,

            ax_ptr,
            ax_ptr, ax_ptr,
            ax_ptr, ax_ptr,
            ax_ptr, ax_ptr,
            ax_ptr, ax_ptr,
            ax_ptr, ax_ptr,
            ax_ptr, ax_ptr,

            0.1
            );
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    auto deltat = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    printf("fastEEC: %lu ms\n", deltat);

}
