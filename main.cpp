#include <string>
#include "interface_config_reader.hpp"

int main(int argc, char* argv[]) {
    std::cout << "MDTools v0.1 beta\n" << std::endl;

    if (argc < 3) {
        std::cerr << "Error, a config file and a job type must be provided!"
                  << std::endl;
        exit(1);
    }
    std::string ini_file(argv[1]);
    auto job_type = interface::JobType(std::stoi(std::string(argv[2])));

    auto reader = interface::ConfigReader(ini_file, job_type);

    if (job_type == interface::kCalcFreeEnergyDiff)
        std::cout << reader.GetFreeEnergyDifference() << std::endl;

    return 0;
}
