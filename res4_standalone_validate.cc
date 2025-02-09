#include <stdio.h>
#include <chrono>
#include <random>
#include <boost/histogram.hpp>

#include "SRothman/EECs/src/standalones/res4_standalone.h"
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/EECs/src/theOnlyHeader.h"

constexpr int NPART = 9;
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

    void makeJet(simon::jet& J, int nPart){
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

    simon::jet J;
    jetFactory.makeJet(J, NPART);

    standaloneEEC::res4_result_multi_array result_multi_array(
        NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS
    );
    standaloneEEC::res4_standalone_multi_array(
            result_multi_array,
            J, standaloneEEC::normType::RAWPT,
            ax,
            ax, ax,
            ax, ax,
            ax, ax,
            0.1, 0.1);

    printf("ran standalone multi array\n");
    /*printf("CONTENTS:\n");
    printf("dipole:\n");
    standaloneEEC::print_nonzero(result_multi_array.get_dipole());
    printf("tee\n");
    standaloneEEC::print_nonzero(result_multi_array.get_tee());
    printf("triangle\n");
    standaloneEEC::print_nonzero(result_multi_array.get_triangle());
    printf("\n");*/

    /*standaloneEEC::res4_result_multi_array result_multi_array_precomputed(
        NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS
    );
    standaloneEEC::res4_standalone_multi_array_precomputed(
            result_multi_array_precomputed,
            J, standaloneEEC::normType::RAWPT,
            ax,
            ax, ax,
            ax, ax,
            ax, ax,
            0.1, 0.1);
    printf("ran standalone multi array precomputed\n");*/
    /*printf("CONTENTS:\n");
    printf("dipole:\n");
    standaloneEEC::print_nonzero(result_multi_array_precomputed.get_dipole());
    printf("tee\n");
    standaloneEEC::print_nonzero(result_multi_array_precomputed.get_tee());
    printf("triangle\n");
    standaloneEEC::print_nonzero(result_multi_array_precomputed.get_triangle());
    printf("\n");*/

    /*printf("COMPARING multi array VS multi array precomputed\n");
    bool agrees = result_multi_array == result_multi_array_precomputed;
    printf("%s\n", agrees? "YES" : "NO");
    printf("\n");

    standaloneEEC::res4_result_vector result_vector(
        NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS
    );
    standaloneEEC::res4_standalone_vector(
            result_vector,
            J, standaloneEEC::normType::RAWPT,
            ax,
            ax, ax,
            ax, ax,
            ax, ax,
            0.1, 0.1);
    printf("ran standalone vector\n");*/
    /*printf("CONTENTS:\n");
    printf("dipole:\n");
    standaloneEEC::print_nonzero(result_vector.get_dipole());
    printf("tee\n");
    standaloneEEC::print_nonzero(result_vector.get_tee());
    printf("triangle\n");
    standaloneEEC::print_nonzero(result_vector.get_triangle());
    printf("\n");*/

    /*printf("COMPARING vector VS multi array\n");
    agrees = result_vector == result_multi_array;
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");
    printf("COMPARING vector VS multi array precomputed\n");
    agrees = result_vector == result_multi_array_precomputed;
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");

    standaloneEEC::res4_result_vector result_vector_precomputed(
        NBINS, NBINS, NBINS, NBINS, NBINS, NBINS, NBINS
    );
    standaloneEEC::res4_standalone_vector_precomputed(
            result_vector_precomputed,
            J, standaloneEEC::normType::RAWPT,
            ax,
            ax, ax,
            ax, ax,
            ax, ax,
            0.1, 0.1);
    printf("ran standalone vector precomputed\n");*/
    /*printf("CONTENTS:\n");
    printf("dipole:\n");
    standaloneEEC::print_nonzero(result_vector_precomputed.get_dipole());
    printf("tee\n");
    standaloneEEC::print_nonzero(result_vector_precomputed.get_tee());
    printf("triangle\n");
    standaloneEEC::print_nonzero(result_vector_precomputed.get_triangle());
    printf("\n");*/

    /*printf("COMPARING vector precomputed VS multi array\n");
    agrees = result_vector_precomputed == result_multi_array;
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");
    printf("COMPARING vector precomputed VS multi array precomputed\n");
    agrees = result_vector_precomputed == result_multi_array_precomputed;
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");
    printf("COMPARING vector precomputed VS vector\n");
    agrees = result_vector_precomputed == result_vector;
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");*/

    fastEEC::result_t<double> result_fastEEC;
    runFastEEC(
        result_fastEEC,
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
    printf("ran fastEEC\n");

    const auto& dipole = *(result_fastEEC.resolved4_shapes->dipole); 
    const auto& tee = *(result_fastEEC.resolved4_shapes->tee);
    const auto& triangle = *(result_fastEEC.resolved4_shapes->triangle);
    /*printf("CONTENTS:\n");
    printf("dipole:\n");
    for(unsigned iR=0; iR<dipole.shape()[0]; ++iR){
        for(unsigned ir=0; ir<dipole.shape()[1]; ++ir){
            for(unsigned ic=0; ic<dipole.shape()[2]; ++ic){
                if (dipole[iR][ir][ic] > 0){
                    printf("(%u, %u, %u): %g\n", iR, ir, ic, dipole[iR][ir][ic]/24.0);
                }
            }
        }
    }
    printf("tee\n");
    for(unsigned iR=0; iR<tee.shape()[0]; ++iR){
        for(unsigned ir=0; ir<tee.shape()[1]; ++ir){
            for(unsigned ic=0; ic<tee.shape()[2]; ++ic){
                if (tee[iR][ir][ic] > 0){
                    printf("(%u, %u, %u): %g\n", iR, ir, ic, tee[iR][ir][ic]/24.0);
                }
            }
        }
    }
    printf("triangle\n");
    for(unsigned iR=0; iR<triangle.shape()[0]; ++iR){
        for(unsigned ir=0; ir<triangle.shape()[1]; ++ir){
            for(unsigned ic=0; ic<triangle.shape()[2]; ++ic){
                if (triangle[iR][ir][ic] > 0){
                    printf("(%u, %u, %u): %g\n", iR, ir, ic, triangle[iR][ir][ic]/24.0);
                }
            }
        }
    }
    printf("\n");*/

    printf("COMPARING fastEEC VS multi array\n");
    
    printf("DIPOLE\n");
    for(size_t iR=0; iR<dipole.shape()[0]; ++iR){
        for(size_t ir=0; ir<dipole.shape()[1]; ++ir){
            for(size_t ic=0; ic<dipole.shape()[2]; ++ic){
                double val1 = dipole[iR][ir][ic]/24.0;
                double val2 = result_multi_array.get_dipole().get_data()[iR][ir][ic];
                if (std::abs(val1 - val2) > 1e-6){
                    printf("(%lu, %lu, %lu): %g != %g\n", 
                            iR, ir, ic, 
                            val1, val2);
                }
            }
        }
    }
    printf("TEE\n");
    for(size_t iR=0; iR<tee.shape()[0]; ++iR){
        for(size_t ir=0; ir<tee.shape()[1]; ++ir){
            for(size_t ic=0; ic<tee.shape()[2]; ++ic){
                double val1 = tee[iR][ir][ic]/24.0;
                double val2 = result_multi_array.get_tee().get_data()[iR][ir][ic];
                if (std::abs(val1 - val2) > 1e-6){
                    printf("(%lu, %lu, %lu): %g != %g\n", 
                            iR, ir, ic, 
                            val1, val2);
                }
            }
        }
    }
    printf("TRIANGLE\n");
    for(size_t iR=0; iR<triangle.shape()[0]; ++iR){
        for(size_t ir=0; ir<triangle.shape()[1]; ++ir){
            for(size_t ic=0; ic<triangle.shape()[2]; ++ic){
                double val1 = triangle[iR][ir][ic]/24.0;
                double val2 = result_multi_array.get_triangle().get_data()[iR][ir][ic];
                if (std::abs(val1 - val2) > 1e-6){
                    printf("(%lu, %lu, %lu): %g != %g\n", 
                            iR, ir, ic, 
                            val1, val2);
                }
            }
        }
    }


    /*agrees = std::equal(
            dipole.data(), 
            dipole.data() + dipole.num_elements(),
            result_multi_array.get_dipole().get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            tee.data(), 
            tee.data() + tee.num_elements(),
            result_multi_array.get_tee().get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            triangle.data(), 
            triangle.data() + triangle.num_elements(),
            result_multi_array.get_triangle().get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");

    printf("COMPARING fastEEC VS multi array precomputed\n");
    agrees = std::equal(
            dipole.data(), 
            dipole.data() + dipole.num_elements(),
            result_multi_array_precomputed.get_dipole().get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            tee.data(), 
            tee.data() + tee.num_elements(),
            result_multi_array_precomputed.get_tee().get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            triangle.data(), 
            triangle.data() + triangle.num_elements(),
            result_multi_array_precomputed.get_triangle().get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");

    printf("COMPARING fastEEC VS vector\n");
    agrees = std::equal(
            dipole.data(), 
            dipole.data() + dipole.num_elements(),
            standaloneEEC::res4_multi_array_container(result_vector.get_dipole()).get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            tee.data(), 
            tee.data() + tee.num_elements(),
            standaloneEEC::res4_multi_array_container(result_vector.get_tee()).get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            triangle.data(), 
            triangle.data() + triangle.num_elements(),
            standaloneEEC::res4_multi_array_container(result_vector.get_triangle()).get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");

    printf("COMPARING fastEEC VS vector precomputed\n");
    agrees = std::equal(
            dipole.data(), 
            dipole.data() + dipole.num_elements(),
            standaloneEEC::res4_multi_array_container(result_vector_precomputed.get_dipole()).get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            tee.data(), 
            tee.data() + tee.num_elements(),
            standaloneEEC::res4_multi_array_container(result_vector_precomputed.get_tee()).get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    agrees &= std::equal(
            triangle.data(), 
            triangle.data() + triangle.num_elements(),
            standaloneEEC::res4_multi_array_container(result_vector_precomputed.get_triangle()).get_data().data(),
            [](double a, double b){
                return std::abs(a/24.0 - b) < 1e-6;
            });
    printf("%s\n", agrees ? "YES" : "NO");
    printf("\n");*/

    return 0;
}
