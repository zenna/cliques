/* @authors Zenna Tavares-zennatavares@gmail.com, Michael Schuab 2010-2011 */
#pragma once

namespace clq {

void dump_gaussian(std::vector<int> &comm, std::map<std::vector<int>, int> &community_to_dump_count) {


}

// TODO - Have a random walker which tends towards local maxima
// Every t_dump steps, the walker will release a pile of sand
// Repeat until full
// Estimate shape by sum of Gaussians

// Difficult things:
//a) Computing the escape probability - // A full basin will have a zero probability of going there from any point.
// If we add sand simply to the optima,it will no longer induce a basin, one or several of its immediate neighbours will.
// How can we differentiate between this case and really flattening the entire basin.
//b) choosing parameters
//c) Making it computationally, naively computing the current quality scales linearly with the number of dumnps
//feasible d)
template<typename G, typename M, typename QF>
void explore_landscape(G &graph, M &weights, QF &compute_quality, double const time,
        std::vector<std::vector<double> > all_stabilities,
        std::vector<std::vector<int> > communities) {

    std::map<std::vector<int>, int> community_to_dump_count;
    double gauss_sigma = 1.0;
    double gauss_mean = 1.0;
    int dump_frequency = 5;
    int num_iterations = 1000, current_timestep = 0;

    std::vector<int> comm = {0};

    for (int timestep = 0; timestep < num_iterations; ++timestep) {
        if (timestep % dump_frequency == 0) {
            community_to_dump_count[comm]++;
       }

        auto neighs = find_community_neighbours(graph, comm);
        for (auto neigh = neighs.begin(); neigh != neighs.end(); ++neigh) {
            clq::VectorPartition p = clq::community_to_partition(graph,
                    *neigh, 0);
            double neigh_quality = compute_quality(p, 1, time);
        }

        // Compute transition probabilities

        // Transition to neighbour with probability

    }
    // Until convergence
    // take a step
    // If it is time to dump, dump
    // Otherw
}

/**
 @brief  Optimise for maxima

 */
template<typename G, typename M, typename QF>
std::set<std::vector<int> > find_optimal_communities_huxley(G &graph,
        M &weights, QF &compute_quality, double time,
        std::vector<std::vector<double> > all_stabilities,
        std::vector<std::vector<int> > communities) {

    std::set<std::vector<int> > all_maxima;
    std::vector<int> comm = { 0 };
    int num_iterations = 1000;
    std::vector<int> buffer(lemon::countNodes(graph), 0);
    std::map<std::vector<int>, int> maxima_to_seen_count;
    int iterations_per_snapshot = 10;

    for (int i = 0; i < num_iterations; ++i) {
        //      clq::print_collection(comm);
        clq::VectorPartition p = clq::community_to_partition(graph,
                comm, 0);
        auto neighs = find_community_neighbours(graph, comm);
        double current_quality = compute_quality(p, 1, time);
        double real_quality = current_quality;
        for (auto maximum = maxima_to_seen_count.begin(); maximum
                != maxima_to_seen_count.end(); ++maximum) {
            int dist_to_maximum = clq::find_community_dist(graph, weights,
                    comm, maximum->first, buffer);
            //          clq::output("dist between following two is", dist_to_maximum);
            //          clq::print_collection(comm);
            //          clq::print_collection(maximum->first);
            current_quality -= discrete_gauss_kernel(dist_to_maximum, 0.5)
                    * maximum->second;
        }

        std::vector<double> neigh_qualities;
        //      clq::output("current quality:", current_quality);

        // Compute neighbour qualities then bias by history dependent maxima filling
        bool is_real_maximum = true;
        for (auto neigh = neighs.begin(); neigh != neighs.end(); ++neigh) {
            //          clq::print_collection(*neigh);

            clq::VectorPartition p = clq::community_to_partition(graph,
                    *neigh, 0);
            double neigh_quality = compute_quality(p, 1, time);

            if (neigh_quality > real_quality) {
                is_real_maximum = false;
            }
            for (auto maximum = maxima_to_seen_count.begin(); maximum
                    != maxima_to_seen_count.end(); ++maximum) {
                int dist_to_maximum = clq::find_community_dist(graph,
                        weights, *neigh, maximum->first, buffer);

                neigh_quality -= discrete_gauss_kernel(dist_to_maximum, 0.5)
                        * maximum->second;
            }
            //          clq::output("neigh quality:", neigh_quality);

            neigh_qualities.push_back(neigh_quality);
        }

        // Save if I am maxima on the unmodified landscape
        if (is_real_maximum == true) {
            all_maxima.insert(comm);
        }

        // Move to neighbour with probability dependent on difference in stability
        double best_quality_diff = -std::numeric_limits<double>::max();
        int j = 0, best_neighbour = -1;
        for (auto neigh_quality = neigh_qualities.begin(); neigh_quality
                != neigh_qualities.end(); ++neigh_quality) {
            double quality_diff = *neigh_quality - current_quality;
            if (quality_diff > best_quality_diff) {
                best_quality_diff = quality_diff;
                best_neighbour = j;
            }
            ++j;
        }
        // If not maxima in meta landscape
        if (best_quality_diff > 0.0) {
            comm = neighs[best_neighbour];
            clq::output("moving", i);
            //          clq::print_collection(comm);
        } else {
            clq::output("meta maxima", i);
            maxima_to_seen_count[comm]++;
        }

        if (i % iterations_per_snapshot == 0) {
            std::vector<double> stabilities;
            for (auto comm = communities.begin(); comm != communities.end(); ++comm) {
                clq::VectorPartition p = clq::community_to_partition(
                        graph, *comm, 0);
                double quality = compute_quality(p, 1, time);
                for (auto maximum = maxima_to_seen_count.begin(); maximum
                        != maxima_to_seen_count.end(); ++maximum) {
                    int dist_to_maximum = clq::find_community_dist(graph,
                            weights, *comm, maximum->first, buffer);

                    quality -= discrete_gauss_kernel(dist_to_maximum, 0.5)
                            * maximum->second;
                }
                stabilities.push_back(quality);
            }
            all_stabilities.push_back(stabilities);
        }
    }

    return all_maxima;
}

template<typename G, typename M, typename QF>
std::set<std::vector<int> > find_optimal_communities_huxley(G &graph,
        M &weights, QF &compute_quality, double time) {
    std::vector<std::vector<int> > communities;
    return find_optimal_communities_huxley(graph, weights, compute_quality,
            time, communities);
}

// Kernals
double discrete_gauss_kernel(int N, double T) {
    return 1.0 / (N + 1);
}

}
