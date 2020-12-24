#include <iostream>
#include <cstdio>
#include "graph.hpp"

using namespace std;
int main() {
    ios::sync_with_stdio(true);
    long long n;
    
    cin >> n;
    // cout<<n<<" Edges\n";
    vector<vector<long long> > edges;
    for(long long i = 0; i < n; i++) {
        //cout<<i<<" ";
        edges.push_back(vector<long long>(3, 0));
        cin>>edges[i][0];
        cin>>edges[i][1];
        cin>>edges[i][2];

    }

    long long m;

    cin >> m;

    vector<long long> boundary_vertices;

    for(long long i = 0; i < m; i++) {
        long long b;
        cin >> b;

        boundary_vertices.push_back(b);
    }

    Graph subgraph = Graph(edges, boundary_vertices);
    subgraph.initialise_bounding_paths(1);
    cout<<"Graph generated successfully\n" << endl;
    cout << subgraph.edges.size() << endl;
    /*
    subgraph.print_neighbours(1230);
    subgraph.print_neighbours(955);
    // cout<<subgraph.vertices[0]<<" ";
    // subgraph.print_vertices();
    // return 0;
    vector<Vertex*> path = subgraph.get_shortest_path(subgraph.vertices[0], subgraph.vertices[1230]);
    cout<<path.size()<<" path size\n";
    cout<<"Start";
    for(Vertex* node: path) {
        cout<<" -> "<<node->id;
    }
    cout<<"\n\n";
    //subgraph.vertices[1362]->enabled = false;
    vector<vector<Vertex*>> paths = subgraph.get_k_bounding_paths(subgraph.vertices[0], subgraph.vertices[1230], 20);
    int i = 1;
    for(vector<Vertex*> path: paths) {
        cout<<i++<<" shortest path with distance - "<<path_distance(path, subgraph.adjacency_list)<<"\n";
        cout<<"Start";
        for(Vertex* node: path) {
            cout<<" -> "<<node->id;
        }
        cout<<"\n\n";
        
    }
    */
    //cout << subgraph.bp.size() << endl;
    //subgraph.initialise_bounding_paths(1);

    return 0;
}