#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <cassert>
#include <string>

#include "common_pmf_parser.hpp"
#include "tool_grad_merger.hpp"
#include "tool_hist_parser.hpp"
#include "tool_pair_interaction.hpp"
#include "tool_path_finder.hpp"
#include "tool_projector.hpp"
#include "tool_reweightor.hpp"

// a general interface for functions of parseABF
//
// Usage:
//    auto f = interface::interface();
//
//    // calculate the error bar of a pmf
//    // note that this is only an approximation, it supposes that the sampling
//    // along the CV-space are evenly distributed
//    // it guarantees Err(align_position) = 0
//    f.CalcError(const hist_parser::Hist<double>& pmf_hist,
//                const std::vector<double>& align_position,
//                const std::string output_file))
//
//    // find the optimal pathway connecting two points
//    // one can set targetedPoints and forceConstants to set additional
//    // restraints (Manhatton Potential)
//    f.FindPathway(
//                  const pmf_parser::Pmf<double>& pmf_info,
//                  const std::vector<double>& initial_point,
//                  const std::vector<double>& end_point,
//                  const std::vector<bool>& pbc,
//                  const std::string& output_prefix,
//                  const std::vector<std::vector<double> >& targeted_points =
//                  {}, const std::vector<std::vector<double> >& force_constants
//                  = {}, bool write_explored_points = false
//                 )
//
//    // calculate time evolution of RMSD over a reference PMF
//    // if the reference PMF is not provided, then the the PMF of last frame in
//    // the hist file will be set as the reference
//    f.CalcRMSD(const std::string hist_file,
//               const std::string pmf_file, const std::string output_file)
//    f.CalcRMSD(const std::string hist_file, const std::string output_file)
//
//    // calculate free-energy difference between two basins
//    // one needs to provide two points in each basin and the corresponding
//    // cutoff energy
//    f.CalcFreeEnergyDiff(
//                         const pmf_parser::Pmf<double>& pmf,
//                         const std::vector<double>& point1,
//                         double cutoff_energy1,
//                         const std::vector<double>& point2,
//                         double cutoff_energy2,
//                         double temperature
//                        )
//
//    // calculate the pair-interaction energies in a CV space
//    // an namd log file and a cvtrj file need to be provided
//    // one needs to indicate which cols (0-based) correspond to the CV values
//    // the infomation (frames) in the cvtrj file must be more than or equal to
//    // that in namd log file
//    // one may set stride to "jump over" frames in the cvtrj file
//    // for instance, if the output freq of cvtrj file is 1000, for namd log
//    // file, 5000, then stride = 4
//    f.CalcEnergyInCVSpace(
//                          const std::vector<double>& lowerboundary,
//                          const std::vector<double>& width,
//                          const std::vector<double>& upperboundary,
//                          const std::string& namd_log_file,
//                          const std::string& cvtrj_file,
//                          const std::vector<int>& cols,
//                          const std::string& output_prefix,
//                          bool self_interaction = false,
//                          int stride = 0
//                          )
//
//    // merge NAMD grad-count file pairs
//    // one need to provide lowerboundary, upperboundary and width for the
//    // output results
//    f.MergeWindows(const std::vector<double>& lowerboundary,
//                   const std::vector<double>& upperboundary,
//                   const std::vector<double>& width,
//                   const std::vector<std::string>& grad_files,
//                   const std::vector<std::string>& count_files,
//                   const std::string& output_prefix)
//
//    // reweight PMF, such that on can get PMF along an alternative RC
//    // one requires the original PMF file and the corresponding cvtrjfile
//    // the columns (0-based) for the original CVs and the new CVs
//    // need to be properly set
//    f.ReweightPmf(const std::string& old_pmf,
//                  const std::vector<std::string>& cvtrj_files,
//                  const std::vector<int>& old_columns,
//                  const std::vector<double>& lowerboundary,
//                  const std::vector<double>& width,
//                  const std::vector<double>& upperboundary,
//                  const std::vector<int>& new_columns, double temperature,
//                  const std::string& output_file)

namespace interface {

    class Interface {
       public:
        // calculate the error bar of a pmf
        void CalcError(const hist_parser::Hist<double>& pmf_hist,
                       const std::vector<double>& align_position,
                       const std::string output_file) const {
            int max_frame = pmf_hist.GetNumberOfFrames() - 1;
            int half_frame = max_frame / 2;

            auto pmf_data = pmf_hist.GetData();

            auto end_pmf = pmf_data[max_frame]->GetPmfData();
            auto first_half_pmf = pmf_data[half_frame]->GetPmfData();
            auto real_pos = pmf_data[max_frame]->RCToInternal(align_position);

            auto last_half_pmf = (end_pmf - first_half_pmf * 0.5) / 0.5;
            last_half_pmf -= last_half_pmf[real_pos] - first_half_pmf[real_pos];

            auto error = pmf_data[0];
            error->SetPMFData((first_half_pmf - last_half_pmf) / 1.41421356);

            error->WritePmfFile(output_file);
        }

        // find optimized pathway
        // return the total number of points explored
        int FindPathway(
            const pmf_parser::Pmf<double>& pmf_info,
            const std::vector<double>& initial_point,
            const std::vector<double>& end_point, const std::vector<bool>& pbc,
            const std::string& output_prefix,
            const std::vector<std::vector<double> >& targeted_points = {},
            const std::vector<std::vector<double> >& force_constants = {},
            bool writeExploredPoints = false) const {
            auto path_find = path_finder::PathFinder(pmf_info, initial_point,
                                                     end_point, pbc);
            std::vector<std::vector<double> > results;
            std::vector<double> energy_results;

            if (targeted_points.size() != 0 && force_constants.size() != 0) {
                path_find.SetTargetedPoints(targeted_points, force_constants);
                path_find.Dijkstra(
                    &path_finder::PathFinder::ManhattonPotential);
            } else {
                path_find.Dijkstra();
            }

            path_find.GetResults(results, energy_results);

            std::string traj_file = output_prefix + ".traj";
            std::string energy_file = output_prefix + ".energy";

            // if one wants to write explored points
            std::vector<std::vector<double> > explored_points;
            std::string explored_points_file = output_prefix + ".explored";
            if (writeExploredPoints) {
                path_find.GetExploredPoints(explored_points);
            }

            // write traj and energy
            WriteData(traj_file, results);
            WriteData(energy_file, energy_results);

            // write explored points
            if (writeExploredPoints) {
                WriteData(explored_points_file, explored_points);
            }

            return path_find.GetExploredPointNum();
        }

        // calculate RMSD over a pmf
        void CalcRMSD(const hist_parser::Hist<double>& hist,
                      const pmf_parser::Pmf<double>& pmf,
                      const std::string& output_file) const {
            auto result = hist.CalcRMSD(pmf);
            WriteData(output_file, result);
        }

        void CalcRMSD(const hist_parser::Hist<double>& hist,
                      const std::string& output_file) const {
            auto result = hist.CalcRMSD();
            WriteData(output_file, result);
        }

        // calculate the free-energy difference between two regions
        double CalcFreeEnergyDiff(const pmf_parser::Pmf<double>& pmf,
                                  const std::vector<double>& point1,
                                  double cutoff_energy1,
                                  const std::vector<double>& point2,
                                  double cutoff_energy2,
                                  double temperature) const {
            auto r1 = pmf.DetermineRegion(point1, cutoff_energy1);
            auto r2 = pmf.DetermineRegion(point2, cutoff_energy2);
            return pmf.CalcFreeEnergyDiff(r1, r2, temperature);
        }

        // do a pair interaction in CV space
        void CalcEnergyInCVSpace(const std::vector<double>& lowerboundary,
                                 const std::vector<double>& width,
                                 const std::vector<double>& upperboundary,
                                 const std::string& namd_log_file,
                                 const std::string& cvtrj_file,
                                 const std::vector<int>& cols,
                                 const std::string& output_prefix,
                                 bool self_interaction = false,
                                 int stride = 0) {
            auto pair_int = pair_interaction::PairInteraction(
                lowerboundary, width, upperboundary, namd_log_file, cvtrj_file,
                cols, self_interaction, stride);
            pair_int.WriteFiles(output_prefix);
        }

        // merge a set of grad-count file pairs
        void MergeWindows(const std::vector<double>& lowerboundary,
                          const std::vector<double>& upperboundary,
                          const std::vector<double>& width,
                          const std::vector<std::string>& grad_files,
                          const std::vector<std::string>& count_files,
                          const std::string& output_prefix) {
            auto MergeFiles =
                grad_merger::GradMerge(lowerboundary, upperboundary, width,
                                       grad_files, count_files, output_prefix);
        }

        void ReweightPmf(const std::string& old_pmf,
                         const std::vector<std::string>& cvtrj_files,
                         const std::vector<int>& old_columns,
                         const std::vector<double>& lowerboundary,
                         const std::vector<double>& width,
                         const std::vector<double>& upperboundary,
                         const std::vector<int>& new_columns,
                         double temperature, bool is_gawtm_eabf,
                         const std::vector<std::string>& namd_log_files,
                         int ignore_samples, const std::string& output_file) {
            auto reweight = reweightor::Reweightor(
                old_pmf, cvtrj_files, old_columns, lowerboundary, width,
                upperboundary, new_columns, temperature, is_gawtm_eabf,
                namd_log_files, ignore_samples);
            reweight.WriteFiles(output_file);
        }

        void ProjectPmf(const pmf_parser::Pmf<double>& pmf_data,
                        const std::vector<int>& final_axis, double temperature,
                        const std::string& output_file) {
            auto project =
                projector::Projector(pmf_data, final_axis, temperature);
            project.WriteFile(output_file);
        }

       private:
        // write a serial of vectors to a file
        // used for writing a trajectory
        void WriteData(const std::string& file,
                       const std::vector<std::vector<double> >& points) const {
            std::ofstream write_file;
            write_file.open(file, std::ios::out);
            if (!write_file.is_open()) {
                std::cerr << "Cannot open " << file << std::endl;
                exit(1);
            }
            for (const auto& result : points) {
                for (const auto& item : result) {
                    write_file << item << " ";
                }
                write_file << std::endl;
            }
            write_file.close();
        }

        // write a serial of numbers to a file
        // used for writing energies
        void WriteData(const std::string& file,
                       const std::vector<double>& data) const {
            std::ofstream write_file;
            write_file.open(file, std::ios::out);
            if (!write_file.is_open()) {
                std::cerr << "Cannot open " << file << std::endl;
                exit(1);
            }
            for (const auto& item : data) {
                write_file << item << std::endl;
            }
            write_file.close();
        }
    };
}  // namespace interface

#endif
