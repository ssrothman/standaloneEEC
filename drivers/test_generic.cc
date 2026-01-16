#include "SRothman/EECs/src/generic/backend.h"
#include "SRothman/EECs/src/Adjacency.h"
#include "SRothman/EECs/src/EECjet.h"
#include "SRothman/SimonTools/src/jet.h"

int main(){
    simon::jet simonjet;
    simonjet.particles.emplace_back(10, 0, 0);
    simonjet.particles.emplace_back(20, 0, 1);
    simonjet.particles.emplace_back(5, -1, 0);
    simonjet.particles.emplace_back(15, 1, 0.5);
    simonjet.particles.emplace_back(12, 1, 0.5);
    simonjet.particles.emplace_back(11, 1, 0.5);
    simonjet.particles.emplace_back(4, 1, 0.5);
    simonjet.nPart = 7;
    simonjet.pt = simonjet.sumpt = simonjet.rawpt = 50;

    auto thejet = std::make_shared<EEC::EECjet_JIT>(simonjet, EEC::normType::RAWPT);

    struct nothing{
        int pass=0;
    };

    auto result = std::make_shared<nothing>();
    auto axis_gen = std::make_shared<nothing>();
    nothing config;

    mainloop<nothing, nothing, EEC::EECjet_JIT, nothing, nothing, nothing, 0b111>(
        result,
        thejet,
        axis_gen,
        config,

        nullptr, nullptr, nullptr, nullptr,
        nullptr
    );
}