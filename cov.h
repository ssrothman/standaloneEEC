#ifndef EEC_COV_H
#define EEC_COV_H

#include "SRothman/EECs/src/fastStructs.h"
#include <boost/multi_array.hpp>

template <typename T>
struct EECcov {
    std::shared_ptr<boost::multi_array<T, 4>> covProj;

    std::shared_ptr<boost::multi_array<T, 6>> covRes3;

    std::shared_ptr<boost::multi_array<T, 6>> covRes4dipole;
    std::shared_ptr<boost::multi_array<T, 6>> covRes4tee;

    void setup(const unsigned nProj, const unsigned nOrder,
               const unsigned nRL_res3, 
               const unsigned nxi_res3,
               const unsigned nphi_res3,
               const unsigned nRL_res4dipole,
               const unsigned nr_res4dipole,
               const unsigned nphi_res4dipole,
               const unsigned nRL_res4tee,
               const unsigned nr_res4tee,
               const unsigned nphi_res4tee){
        covProj = std::make_shared<boost::multi_array<T, 4>>(
                boost::extents[nOrder][nProj][nOrder][nProj]
        );

        covRes3 = std::make_shared<boost::multi_array<T, 6>>(
                boost::extents[nRL_res3][nxi_res3][nphi_res3][nRL_res3][nxi_res3][nphi_res3]
        );

        covRes4dipole = std::make_shared<boost::multi_array<T, 6>>(
                boost::extents[nRL_res4dipole][nr_res4dipole][nphi_res4dipole][nRL_res4dipole][nr_res4dipole][nphi_res4dipole]
        );

        covRes4tee = std::make_shared<boost::multi_array<T, 6>>(
                boost::extents[nRL_res4tee][nr_res4tee][nphi_res4tee][nRL_res4tee][nr_res4tee][nphi_res4tee]
        );
    }

    void fill(const fastEEC::result_t<T>& ans){
        //cov proj
        for (unsigned order1=0; order1 < 5; ++order1){
            for (unsigned RL1=0; RL1 < ans.wts[order1]->size(); ++RL1){
                for (unsigned order2=0; order2 < 5; ++order2){
                    for(unsigned RL2=0; RL2 < ans.wts[order2]->size(); ++RL2){
                        (*covProj)[order1][RL1][order2][RL2] +=
                            ans.wts[order1]->at(RL1) * 
                            ans.wts[order2]->at(RL2);
                    }
                }
            }
        }

        if (ans.resolved3){
            //cov res3
            for (unsigned RL1=0; RL1 < ans.resolved3->shape()[0]; ++RL1){
                for(unsigned xi1=0; xi1 < ans.resolved3->shape()[1]; ++xi1){
                    for(unsigned phi1=0; phi1 < ans.resolved3->shape()[2]; ++ phi1){
                        for (unsigned RL2=0; RL2 < ans.resolved3->shape()[0]; ++ RL2){
                            for(unsigned xi2=0; xi2 < ans.resolved3->shape()[1]; ++ xi2){
                                for(unsigned phi2=0; phi2 < ans.resolved3->shape()[2]; ++ phi2){
                                    (*covRes3)[RL1][xi1][phi1][RL2][xi2][phi2] +=
                                        (*ans.resolved3)[RL1][xi1][phi1] *
                                        (*ans.resolved3)[RL2][xi2][phi2];
                                }
                            }
                        }
                    }
                }
            }
        }

        if(ans.resolved4_shapes){
            //cov res4 dipole
            for(unsigned R1=0; R1<ans.resolved4_shapes->dipole->shape()[0]; ++R1){
                for(unsigned r1=0; r1<ans.resolved4_shapes->dipole->shape()[1]; ++ r1){
                    for(unsigned ct1=0; ct1<ans.resolved4_shapes->dipole->shape()[2]; ++ ct1){
                        for(unsigned R2=0; R2<ans.resolved4_shapes->dipole->shape()[0]; ++ R2){
                            for(unsigned r2=0; r2<ans.resolved4_shapes->dipole->shape()[1]; ++ r2){
                                for(unsigned ct2=0; ct2<ans.resolved4_shapes->dipole->shape()[2]; ++ ct2){
                                    (*covRes4dipole)[R1][r1][ct1][R2][r2][ct2] +=
                                        (*ans.resolved4_shapes->dipole)[R1][r1][ct1] *
                                        (*ans.resolved4_shapes->dipole)[R2][r2][ct2];
                                }
                            }
                        } 
                    }
                }
            }

            //cov res4 tee
            for (unsigned R1=0; R1<ans.resolved4_shapes->tee->shape()[0]; ++ R1){
                for (unsigned r1=0; r1<ans.resolved4_shapes->tee->shape()[1]; ++ r1){
                    for (unsigned ct1=0; ct1<ans.resolved4_shapes->tee->shape()[2]; ++ ct1){
                        for (unsigned R2=0; R2<ans.resolved4_shapes->tee->shape()[0]; ++ R2){
                            for (unsigned r2=0; r2<ans.resolved4_shapes->tee->shape()[1]; ++ r2){
                                for (unsigned ct2=0; ct2<ans.resolved4_shapes->tee->shape()[2]; ++ ct2){
                                    (*covRes4tee)[R1][r1][ct1][R2][r2][ct2] +=
                                        (*ans.resolved4_shapes->tee)[R1][r1][ct1] *
                                        (*ans.resolved4_shapes->tee)[R2][r2][ct2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif
