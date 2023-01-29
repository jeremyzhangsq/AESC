

#ifndef GRAPH_H
#define GRAPH_H


#include "utils.h"


class Graph{
    public:
        std::vector<std::vector<uint>> m_edges; // edges of each node
        std::vector<uint> m_deg;     // degree of each node
    private:
        std::string m_folder;
        std::string m_graph;
        uint m_n; // the number of nodes
        uint64 m_m; // the number of edges
        double m_lambda; // the 2nd largest eigenvalue

        void readNM();
        void readGraph();
        void addEdge(uint u, uint v);
    public:
        uint getDeg(uint u) const;
        uint64 getM() const;
        uint getN() const;
        double getLambda() const;
        std::string getGraphFolder() const;

        Graph(const std::string& t_folder, const std::string& t_graph);
};
#endif