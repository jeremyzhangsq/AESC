#include "graph.h"





#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}

using namespace std;

void handle_error(const char* msg) {
	perror(msg);
	exit(255);
}

Graph::Graph(const string& t_folder, const string& t_graph): 
            m_folder(t_folder), m_graph(t_graph) {

    readNM();
    this->m_edges = std::vector<std::vector<uint>>(this->m_n, std::vector<uint>());
    this->m_deg = std::vector<uint>(this->m_n);

    readGraph();
}

void Graph::readNM(){
    ifstream fin((m_folder + "/" + m_graph + "/stat.txt").c_str());
    string s;
    if (fin.is_open()){
        while (fin >> s){
            if (s.substr(0, 2) == "n="){
                this->m_n = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "m="){
                this->m_m = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "l="){
                this->m_lambda = atof(s.substr(2).c_str());
                continue;
            }
        }
        fin.close();
    }
    else handle_error("Fail to open attribute file!");
}

void Graph::readGraph(){
    FILE *fin = fopen((m_folder + "/" + m_graph + "/graph.txt").c_str(), "r");
    uint64 readCnt = 0;
    uint u, v;
    while (fscanf(fin, "%d%d", &u, &v) != EOF) {
        readCnt++;
        ASSERT( u < this->m_n );
        ASSERT( v < this->m_n );
        addEdge(u, v);
    }
    fclose(fin);
//    cout << "read Graph Done!" << endl;
}

void Graph::addEdge(uint u, uint v){
    this->m_edges[u].push_back(v);
    this->m_edges[v].push_back(u);
    this->m_deg[u]+=1;
    this->m_deg[v]+=1;
}

uint Graph::getDeg(uint u) const{
    return this->m_deg[u];
}

uint Graph::getN() const{
    return this->m_n;
}

uint64 Graph::getM() const{
    return this->m_m;
}

double Graph::getLambda() const{
    return this->m_lambda;
}

std::string Graph::getGraphFolder() const{
    std::string folder(m_folder + "/" + m_graph + "/");
    return folder;
}
