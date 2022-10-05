#ifndef TOOL_PATH_FINDER_HPP
#define TOOL_PATH_FINDER_HPP

#include <cassert>
#include <cstdlib>
#include <map>
#include <vector>

#include "common_pmf_parser.hpp"

// find the optimal pathway connecting two points on a pmf
// Usage:
//   // run calculation
//   auto path = pathFinder(pmfData, initialPoint, endPoint, pbc)
//   path.Dijkstra()
//   // one may want to add external manhatton potential
//   path.setTargetedPoints(
//                          {{19.5,2.2},{20.0,2.5}},
//                          {{1.0,1.0},{1.0,1.0}}
//                         )
//   path.Dijkstra(pathFinder::manhattonPotential)
//   // get results
//   std::vector<std::vector<double> > trajectory
//   std::vector<double> energyResults
//   path.getResults(trajectory, energyResults)
//   // get all the points explored during the process
//   std::vector<std::vector<double> > pointList
//   path.getExploredPoints(pointList)
//   // get the number of points explored
//   auto num = path.getExploredPointNum()
//

namespace path_finder {

    // find the optimal pathway connecting two points on a PMF
    class PathFinder {
       public:
        // constructor
        // either explicitly provide pbc or read it from pmfData
        PathFinder(const pmf_parser::Pmf<double>& pmf_data,
                   const std::vector<double>& initial_point,
                   const std::vector<double>& end_point,
                   const std::vector<bool>& pbc = {}) {
            assert(initial_point.size() == end_point.size());
            assert(initial_point.size() == pbc.size() ||
                   pbc == std::vector<bool>{});
            assert(initial_point.size() == pmf_data.GetDimension());

            this->pmf_data_ = &pmf_data;
            // internally, all the analyses are performed in the internal RC
            // space
            this->lowerboundary_ =
                std::vector<int>(this->pmf_data_->GetDimension(), 0);
            this->upperboundary_ = this->pmf_data_->GetShape();
            // upperboundary is in fact shape - 1
            for (auto& b : this->upperboundary_) {
                b -= 1;
            }
            this->width_ = std::vector<int>(this->pmf_data_->GetDimension(), 1);
            this->initial_point_ = this->pmf_data_->RCToInternal(initial_point);
            this->end_point_ = this->pmf_data_->RCToInternal(end_point);

            // if pbc is not provided
            // then read it from pmfData
            if (pbc != std::vector<bool>{}) {
                this->pbc_ = pbc;
            } else {
                this->pbc_ = this->pmf_data_->GetPbc();
            }

            this->dimension_ = this->pmf_data_->GetDimension();

            this->open_list_.push_back(this->initial_point_);
        }

        // set targeted points and force constants
        // used in the A-star alg
        void SetTargetedPoints(
            const std::vector<std::vector<double> >& points,
            const std::vector<std::vector<double> >& force_const) {
            assert(points.size() == force_const.size());
            for (int i = 0; i < points.size(); i++) {
                this->targeted_points_.push_back(
                    this->pmf_data_->RCToInternal(points[i]));
                this->force_constants_.push_back(force_const[i]);
            }
        }

        // run the Dijkstra alg
        void Dijkstra(double (PathFinder::*func)(const std::vector<int>& point)
                          const = &PathFinder::DefaultFunc) {
            // just a simple translation of the classical dijkstera alg
            while (this->open_list_.size() != 0) {
                auto p = this->PopMin(this->open_list_, func);
                this->close_list_.push_back(p);

                // the end point is found
                if (p == this->end_point_) {
                    break;
                }

                for (const auto& q : this->FindAdjacentPoints(p)) {
                    if (common_tools::VectorInVectorOfVector(q,
                                                             this->open_list_)) {
                        ;  // pass
                    } else if (!common_tools::VectorInVectorOfVector(
                                   q, this->close_list_)) {
                        this->open_list_.push_back(q);
                        this->father_points_[q] = p;
                    }
                }
            }
        }

        // return the explored points of dijkstera calculation
        void GetExploredPoints(
            std::vector<std::vector<double> >& point_list) const {
            if (this->close_list_.size() == 0) {
                std::cerr << "Error, no information about results!\n";
                exit(1);
            }

            point_list = {};
            for (const auto& p : this->close_list_) {
                point_list.push_back(this->pmf_data_->InternalToRC(p));
            }
        }

        // how many points have been explored during the calculation
        int GetExploredPointNum() const {
            if (this->close_list_.size() == 0) {
                std::cerr << "Error, no information about results!\n";
                exit(1);
            }
            return this->close_list_.size();
        }

        // return the minimum energy pathway of dijkstera calculation
        void GetResults(std::vector<std::vector<double> >& trajectory,
                        std::vector<double>& energy_results) {
            if (this->close_list_.size() == 0) {
                std::cerr << "Error, no information about results!\n";
                exit(1);
            }

            std::vector<std::vector<int> > internal_trajectory = {};
            trajectory = {};
            energy_results = {};

            this->ConstructResults(this->end_point_, internal_trajectory);
            internal_trajectory.push_back(this->end_point_);

            for (const auto& p : internal_trajectory) {
                trajectory.push_back(this->pmf_data_->InternalToRC(p));
            }

            for (const auto& p : internal_trajectory) {
                energy_results.push_back((*(this->pmf_data_))[p]);
            }
        }

        // below are h(x) functions used in the A-star alg
        inline double DefaultFunc(const std::vector<int>& point) const {
            return 0;
        }

        // Manhatton potential
        double ManhattonPotential(const std::vector<int>& point) const {
            if (this->targeted_points_.size() == 0) {
                return 0;
            }
            double energy = 0;
            int distance = 0;
            int d1, d2, d3;
            for (int i = 0; i < this->targeted_points_.size(); i++) {
                for (int j = 0; j < this->targeted_points_[i].size(); j++) {
                    distance = 0;
                    if (this->pbc_[j] == false) {
                        distance = abs(point[j] - targeted_points_[i][j]);
                    } else {
                        d1 = abs(point[j] - targeted_points_[i][j]);
                        d2 = (abs(point[j] - lowerboundary_[j]) +
                              abs(targeted_points_[i][j] - upperboundary_[j]));
                        d3 = (abs(point[j] - upperboundary_[j]) +
                              abs(targeted_points_[i][j] - lowerboundary_[j]));
                        distance = ((d1 < d2 ? d1 : d2) < d3)
                                       ? (d1 < d2 ? d1 : d2)
                                       : d3;
                    }
                    energy += distance * this->force_constants_[i][j];
                }
            }
            return energy;
        }

       private:
        // used in getResults
        std::vector<int> ConstructResults(
            std::vector<int> point,
            std::vector<std::vector<int> >& internal_trajectory) {
            if (this->father_points_.find(point) != this->father_points_.end()) {
                internal_trajectory.push_back(ConstructResults(
                    this->father_points_[point], internal_trajectory));
            }
            return point;
        }

        // find the point that has the lowest energy in popList,
        // remove it from the list and return it
        std::vector<int> PopMin(
            std::vector<std::vector<int> >& pop_list,
            double (PathFinder::*func)(const std::vector<int>& point)
                const = &PathFinder::DefaultFunc) {
            assert(pop_list.size() != 0);

            // func(popList[0]) is h(x) in A-star alg
            // by default it is zero
            double min_energy =
                (*pmf_data_)[pop_list[0]] + (this->*func)(pop_list[0]);
            int min_pos = 0;

            for (int i = 0; i < pop_list.size(); i++) {
                if ((*pmf_data_)[pop_list[i]] + (this->*func)(pop_list[i]) <
                    min_energy) {
                    min_energy =
                        (*pmf_data_)[pop_list[i]] + (this->*func)(pop_list[i]);
                    min_pos = i;
                }
            }

            auto min_point = pop_list[min_pos];
            pop_list.erase(pop_list.begin() + min_pos);
            return min_point;
        }

        // find the adjacent points of the input point,
        // return them as a list
        std::vector<std::vector<int> > FindAdjacentPoints(
            const std::vector<int>& point) const {
            std::vector<std::vector<int> > adjacent_points;

            for (int i = 0; i < this->dimension_; i++) {
                // left side
                auto p = point;
                if (point[i] - this->width_[i] >=
                    this->lowerboundary_[i] - common_tools::kAccuracy) {
                    p[i] = point[i] - this->width_[i];
                } else {
                    if (this->pbc_[i]) {
                        p[i] = this->upperboundary_[i];
                    } else {
                        p = {};
                    }
                }

                // right side
                auto q = point;
                if (point[i] + this->width_[i] <=
                    this->upperboundary_[i] + common_tools::kAccuracy) {
                    q[i] = point[i] + this->width_[i];
                } else {
                    if (this->pbc_[i]) {
                        q[i] = this->lowerboundary_[i];
                    } else {
                        q = {};
                    }
                }

                if (p.size() != 0) {
                    adjacent_points.push_back(p);
                }
                if (q.size() != 0) {
                    adjacent_points.push_back(q);
                }
            }

            if (adjacent_points.size() != 0) {
                return adjacent_points;
            } else {
                std::cerr << "Error! No adjacent point is found!" << std::endl;
                exit(1);
            }
        }

        // the pmf data
        const pmf_parser::Pmf<double>* pmf_data_;
        std::vector<int> lowerboundary_;
        std::vector<int> upperboundary_;
        std::vector<int> width_;
        // initial and end point
        std::vector<int> initial_point_;
        std::vector<int> end_point_;
        // whether periodic for each dimension
        std::vector<bool> pbc_;
        // dimension of pmf
        int dimension_;

        // below are vars that will be used in dijkstra/A* algs
        // record points' father point
        std::map<std::vector<int>, std::vector<int> > father_points_ = {};
        // open and closeList in dijkstra/A* algs
        std::vector<std::vector<int> > open_list_ = {};
        std::vector<std::vector<int> > close_list_ = {};

        // in the A* alg, one can define manhatton potential
        // based on targeted points and force constants
        std::vector<std::vector<int> > targeted_points_;
        std::vector<std::vector<double> > force_constants_;
    };
}  // namespace path_finder

#endif  // PATHFINDER_HPP
