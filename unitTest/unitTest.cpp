#include <iostream>

#include "../interface.hpp"

int main() {
    auto face = interface::Interface();

    // calculate error
    auto histFile = hist_parser::Hist<double>("./test.hist.pmf");
    face.CalcError(histFile, {0.2, 90}, "./test.error");
    std::cout << "error calculation finished!" << std::endl;

    // calculate RMSD
    auto pmfFile = pmf_parser::Pmf<double>("./test.pmf");
    face.CalcRMSD(histFile, pmfFile, "./test.rmsd.pmfref");
    std::cout << "RMSD calculation 1 finished!" << std::endl;
    face.CalcRMSD(histFile, "./test.rmsd.selfref");
    std::cout << "RMSD calculation 2 finished!" << std::endl;

    // calculate pathway
    auto pmfFile_nanma = pmf_parser::Pmf<double>("./nanma_test.pmf");
    face.FindPathway(pmfFile_nanma, {-156, 160}, {78, -58}, {true, true},
                     "./nanma_test", {{0, 100}}, {{0.1, 0.1}}, true);
    std::cout << "pathway calculation finished!" << std::endl;

    // calculate free-energy difference
    auto diff = face.CalcFreeEnergyDiff(pmfFile_nanma, {-156, 160}, 3,
                                        {78, -58}, 5, 300);
    std::cout << "free-energy difference is " << diff << std::endl;
    std::cout << "free-energy difference calculation finished!" << std::endl;

    // calculate pait-interaction energies in CV space
    face.CalcEnergyInCVSpace({-11.0}, {0.2}, {11.0}, "./test_namdlog.log",
                             "./test_cvtrj.colvars.traj", {1},
                             "./test_pairEnergy", false, 4);
    std::cout << "pair interaction in CV space calculation finished"
              << std::endl;
    // face.CalcEnergyInCVSpace({-11.0}, {0.2}, {11.0},
    //                          "D:\\Project\\ringThroughRing\\008_onlysmr_Water\\pair_solvent.log",
    //                          "D:\\Project\\ringThroughRing\\008_onlysmr_Water\\output\\abf.colvars.traj",
    //                          {1},
    //                          "D:\\Project\\ringThroughRing\\008_onlysmr_Water\\pair_solvent",
    //                          true,
    //                          4);
    // std::cout << "self pair interaction in CV space calculation finished" <<
    // std::endl;

    // merge windows
    face.MergeWindows(
        {-18, -15}, {18, 15}, {0.1, 0.2},
        {"./grads\\win1-2.abf1.grad", "./grads\\win2-1.abf1.grad",
         "./grads\\win3-1.abf1.grad", "./grads\\win4-1.abf1.grad",
         "./grads\\win5-1.abf1.grad", "./grads\\win6.abf1.grad"},
        {"./grads\\win1-2.abf1.count", "./grads\\win2-1.abf1.count",
         "./grads\\win3-1.abf1.count", "./grads\\win4-1.abf1.count",
         "./grads\\win5-1.abf1.count", "./grads\\win6.abf1.count"},
        "./grads\\merged");
    std::cout << "merging windows finished" << std::endl;

    // reweight
    face.ReweightPmf(
        "E:/codes/MDTools/rewrite/unitTest/chig_ml_2o.abf1.czar.pmf",
        {"E:/codes/MDTools/rewrite/unitTest/chig_ml_1o.colvars.traj",
         "E:/codes/MDTools/rewrite/unitTest/chig_ml_2o.colvars.traj",
         "E:/codes/MDTools/rewrite/unitTest/chig_mlo.colvars.traj"},
        {1, 3}, {0, 0}, {0.2, 0.2}, {20, 20}, {5, 6}, 300,
        "E:/codes/MDTools/rewrite/unitTest/reweight.pmf");
    std::cout << "reweighting finished" << std::endl;

    return 0;
}
