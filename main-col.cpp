#include <mpi.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <unordered_map>
#include <string>
#include <chrono>
#include <cstdio>
#include "graph.hpp"
#include "col-queries.hpp"


long long N_SUBGRAPHS = 18980;
std::string SUBGRAPH_DIR = "./col/";
long long K = 1;

std::vector<long long> query_nodes;

int main(int argc, char** argv) {
    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point stop;
    if(world_rank == 0) {
        start = std::chrono::high_resolution_clock::now();
    }

    std::vector<Graph*> subgraphs;
    for(long long i = world_rank; i < N_SUBGRAPHS; i += world_size) {
        std::ifstream infile;
        infile.open(SUBGRAPH_DIR+"subgraph"+std::to_string(i));
        long long n;

        infile >> n;
        // cout<<n<<" Edges\n";
        std::vector<std::vector<long long> > edges;
        for(long long i = 0; i < n; i++) {
            //cout<<i<<" ";
            edges.push_back(std::vector<long long>(3, 0));
            infile>>edges[i][0];
            infile>>edges[i][1];
            infile>>edges[i][2];

        }

        long long m;

        infile >> m;

        std::vector<long long> boundary_vertices;

        for(long long i = 0; i < m; i++) {
            long long b;
            infile >> b;

            boundary_vertices.push_back(b);
        }

        Graph* subgraph = new Graph(edges, boundary_vertices);
        subgraph->initialise_bounding_paths(K);
        subgraphs.push_back(subgraph);
        infile.close();

        //std::cout << "Generated subgraph on processor rank - " << world_rank << std::endl;
    }

    // wait till all processors have created subgraphs
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<std::vector<long long> > edges;
    std::unordered_map<long long, std::unordered_map<long long, long long > > adj;
    // subgraph index to send for the current node
    long long si = 0;
    for(long long i = 0; i < N_SUBGRAPHS; i++) {
        long long pid = i%world_size;

        long long n;

        if(pid == world_rank) {
            n = 3*subgraphs[si]->lower_bounding_paths.size();
        }

        MPI_Bcast(&n, 1, MPI_LONG_LONG, pid, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        // now we create a data buffer
        long long data[n];
        //std::vector<int> data;
        long long di = 0;
        if(pid == world_rank) {
            // put data into buffer
            for(auto path : subgraphs[si]->lower_bounding_paths) {
                //if(di >= n) continue;
                data[di++] = path.front()->id;
                //if(di >= n) continue;
                data[di++] = path.back()->id;
                //if(di >= n) continue;
                data[di++] = path_distance(path, subgraphs[si]->adjacency_list);
            }
            si++;
        }
        MPI_Bcast(data, n, MPI_LONG_LONG, pid, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for(long long j = 0; j < n; j += 3) {
            //edges.push_back({data[j+0], data[j+1], data[j+2]});
            long long u = data[j+0];
            long long v = data[j+1];
            long long w = data[j+2];
            if(adj.find(u) == adj.end()) {
                adj[u] = std::unordered_map<long long, long long>();
            }
            if(adj[u].find(v) == adj[u].end() || w < adj[u][v]) {
                adj[u][v] = w;
            }
        }
    }
    for(auto v1: adj) {
        for(auto v2 : v1.second) {
            edges.push_back({v1.first, v2.first, v2.second});
        }
    }
    Graph skeleton_graph = Graph(edges, std::vector<long long>());
    
    if(world_rank == 0) {
        //std::cout <<  << std::endl;
        stop = std::chrono::high_resolution_clock::now();
        auto time_diff = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "#vertices in subgraph: " << skeleton_graph.vertices.size() << std::endl; 
        std::cout << "#Edges in subgraph: " << edges.size() << std::endl; 
        std::cout << "Time taken to generate skeleton graph on all nodes: " << (time_diff.count() / 1000000) << " seconds" << std::endl; 

        start = std::chrono::high_resolution_clock::now();
    }

    /*for(auto query : QUERIES) {
        long long source = query[0];
        long long dest = query[1];
    }*/
    for(long long i = 0; i < QUERIES.size(); i++) {
        long long pid = i%world_size;
        long long s = QUERIES[i][0];
        long long t = QUERIES[i][1];
        skeleton_graph.get_shortest_path(skeleton_graph.vertices[s], skeleton_graph.vertices[t]);
    }

    MPI_Barrier(MPI_Comm comm);

    if(world_rank == 0) {
        stop = std::chrono::high_resolution_clock::now();
        auto time_diff = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time taken to process 1000 queries: " << (time_diff.count() / 1000000) << " seconds" << std::endl; 

    }

    // Finalize the MPI environment.
    MPI_Finalize();
    
}