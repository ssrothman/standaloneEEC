#ifndef TEST_AXES_H
#define TEST_AXES_H

#include <boost/histogram.hpp>

using axis_t = boost::histogram::axis::variable<double>;
using axisptr_t = std::shared_ptr<axis_t>;

inline axisptr_t get_xi_res3_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0, 0.2, 0.4, 0.6, 0.8, 1.0
    }}));
}

inline axisptr_t get_phi_res3_ax(){
#ifdef ALT_RES3
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32
    }}));
#else
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0, 0.2, 0.4, 0.6, 0.8, 1.0
    }}));
#endif
}

inline axisptr_t get_RLax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0
    }}));
}

inline axisptr_t get_RLax_coarse(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.2, 0.3, 0.4, 0.5
    }}));
}

inline axisptr_t get_r_dipole_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
    }}));
}

inline axisptr_t get_phi_dipole_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32
    }}));
}

inline axisptr_t get_r_tee_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
    }}));
}

inline axisptr_t get_phi_tee_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32
    }}));
}

inline axisptr_t get_r_tri_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0
    }}));
}

inline axisptr_t get_phi_tri_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32,
             17*M_PI/32, 18*M_PI/32, 19*M_PI/32, 20*M_PI/32,
             21*M_PI/32, 22*M_PI/32, 23*M_PI/32, 24*M_PI/32,
             25*M_PI/32, 26*M_PI/32, 27*M_PI/32, 28*M_PI/32,
             29*M_PI/32, 30*M_PI/32, 31*M_PI/32, 32*M_PI/32,
             33*M_PI/32, 34*M_PI/32, 35*M_PI/32, 36*M_PI/32,
             37*M_PI/32, 38*M_PI/32, 39*M_PI/32, 40*M_PI/32,
             41*M_PI/32, 42*M_PI/32, 43*M_PI/32, 44*M_PI/32,
             45*M_PI/32, 46*M_PI/32, 47*M_PI/32, 48*M_PI/32,
             49*M_PI/32, 50*M_PI/32, 51*M_PI/32, 52*M_PI/32,
             53*M_PI/32, 54*M_PI/32, 55*M_PI/32, 56*M_PI/32,
             57*M_PI/32, 58*M_PI/32, 59*M_PI/32, 60*M_PI/32,
             61*M_PI/32, 62*M_PI/32, 63*M_PI/32, 64*M_PI/32
    }}));
}

inline axisptr_t get_r_minR_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
    }}));
}

inline axisptr_t get_phi_minR_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32
    }}));
}

inline axisptr_t get_phidiff_minR_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32
    }}));
}

inline axisptr_t get_theta_minR_ax(){
    return std::make_shared<axis_t>(std::vector<double>({{
        0.0,    M_PI/32,  2*M_PI/32,  3*M_PI/32,  4*M_PI/32, 
              5*M_PI/32,  6*M_PI/32,  7*M_PI/32,  8*M_PI/32, 
              9*M_PI/32, 10*M_PI/32, 11*M_PI/32, 12*M_PI/32, 
             13*M_PI/32, 14*M_PI/32, 15*M_PI/32, 16*M_PI/32
    }}));
}
#endif
