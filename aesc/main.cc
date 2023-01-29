
#include <iostream>
#include <boost/program_options.hpp>
#include <random>
#include <string>
#include <iomanip>
#include <sstream>
#include "graph.h"
#include "algo.h"

using namespace std;
namespace po = boost::program_options;


namespace { 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace

Config parseParams(int argc, char** argv){
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("data-folder,f", po::value<string>()->required(), "graph data folder")
        ("graph-name,g", po::value<string>()->required(), "graph file name")
        ("algo,a", po::value<string>()->required(), "algorithm name")
//        ("edge-seed,s", po::value<int>()->default_value(0), "use edge seeds")
        ("epsilon,e", po::value<double>()->default_value(0.05), "epsilon")
        ("gamma,y", po::value<int>()->default_value(10), "gamma")
        ("omega,w", po::value<int>()->default_value(128), "omega largest eigen values")
        ("verbose,v", po::value<int>()->default_value(1), "evaluate error or not")
    ;

    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw 
    po::notify(vm);

    Config config;

    if (vm.count("help")){
        cout << desc << '\n';
        exit(0);
    }
    if (vm.count("data-folder")){
        config.strFolder = vm["data-folder"].as<string>();
    }
    if (vm.count("graph-name")){
        config.strGraph = vm["graph-name"].as<string>();
    }
    if (vm.count("algo")){
        config.strAlgo = vm["algo"].as<string>();
    }
    if (vm.count("epsilon")){
        config.epsilon = vm["epsilon"].as<double>();
    }
    if (vm.count("gamma")){
        config.gamma = vm["gamma"].as<int>();
    }
    if (vm.count("omega")){
        config.omega = vm["omega"].as<int>();
    }
    if (vm.count("verbose")){
        config.evaflag = vm["verbose"].as<int>();
    }
    return config;
}

void loadSeed(const string& folder, const string& file_name, const Graph& graph, vector<ipair>& vecSeeds, vector<double>& vecER){
    string strSeed = "/gt_sec.txt";

    FILE *fin = fopen((folder + "/" + file_name + strSeed).c_str(), "r");
    uint s, t;
    double st;
    // vector<ipair> seeds;
    int i=0;
    while (fscanf(fin, "%d %d %lf", &s, &t, &st) != EOF) {
        ipair p = MP(s,t);
        vecSeeds.push_back(p);
        vecER.push_back(st);
        i++;
    }
    fclose(fin);
}


void loadEigens(const string &folder, const string &file_name,const Graph& graph, vector<pair<double,vector<double>>> &Eigens, int omega) {
    string strEigen = "/sorted_eigens_" + to_string(omega) + ".txt";
    std::ifstream infile((folder + "/" + file_name + strEigen).c_str());
//    FILE *fin = fopen((folder + "/" + file_name + strEigen).c_str(), "r");
    double val;
    uint n = graph.getN();
    Eigens.reserve(omega);
    for (int i = 0; i < omega; ++i) {
        infile >> val;
        vector<double> vec(n,0);
        for (int j = 0; j < n; ++j) {
            infile >> vec[j];
        }
        pair<double,vector<double>> p = MP(val,vec);
        Eigens.push_back(p);
    }
    infile.close();
}



int main(int argc, char **argv){
    Config config;
    try{
        config = parseParams(argc, argv);
        config.check();
    }
    catch (const exception &ex){
        cerr << ex.what() << '\n';
        return ERROR_IN_COMMAND_LINE;
    }


    Graph graph(config.strFolder, config.strGraph);

    config.delta = 1.0/graph.getN(); // set failure probability to 1/n

    config.display();
//    int query_count = config.numQuery;

    vector<ipair> seeds;
    vector<double> exact_secs;
    vector<double> errors;

    // set the edge set as seed set
    if (config.evaflag){
        loadSeed(config.strFolder, config.strGraph,graph, seeds, exact_secs);
        errors.resize(exact_secs.size(),0);
    }
    else{
        for (uint u = 0; u < graph.getN(); ++u) {
            for (uint &v:graph.m_edges[u])
                seeds.emplace_back(make_pair(u,v));
        }
    }

    config.lambda = graph.getLambda();

    if(config.strAlgo==TGTP){
        vector<pair<double,vector<double>>> Eigens;
        map<ipair,int> taus;
        loadEigens(config.strFolder, config.strGraph, graph, Eigens,config.omega);
        calTau(Eigens, taus, graph, config.omega, config.epsilon/2.0);
        Timer tm(1, "TGT+");
        srand(time(NULL));
        map<ipair,double> pred_secs;
        truncatedGraphTravPlus(graph,pred_secs,taus,config);
        if (config.evaflag){
            for (int q = 0; q < graph.getM(); ++q) {
                ipair &edge = seeds[q];
                errors[q] = abs(exact_secs[q]-pred_secs[edge]);
            }
        }
    }
    else if(config.strAlgo==TGT){
        vector<pair<double,vector<double>>> Eigens;
        map<ipair,int> taus;
        loadEigens(config.strFolder, config.strGraph, graph, Eigens,config.omega);
        calTau(Eigens, taus, graph, config.omega,config.epsilon);
        Timer tm(1, "TGT");
        srand(time(NULL));
        map<ipair,double> pred_secs;
        truncatedGraphTrav(graph,pred_secs,taus);
        if (config.evaflag){
            for (int q = 0; q < graph.getM(); ++q) {
                ipair &edge = seeds[q];
                errors[q] = abs(exact_secs[q]-pred_secs[edge]);
            }
        }
    }
    else if(config.strAlgo==MC){
        vector<pair<double,vector<double>>> Eigens;
        map<ipair,int> taus;
        loadEigens(config.strFolder, config.strGraph, graph, Eigens,config.omega);
        calTau(Eigens, taus, graph, config.omega, config.epsilon/2.0);
        vector<int> edgeid(graph.getM());
        iota (begin(edgeid), end(edgeid), 0);
        shuffle(edgeid.begin(), edgeid.end(), std::mt19937(std::random_device()()));
        Timer tm(1, "MonteCarlo");
        srand(time(NULL));

	cout << "runnng.." << endl;
	vector<map<uint,double>> vecER;
	for(uint src=0; src<graph.getN(); src++){
	    int lenwalk = 0;
	    for(const auto& v: graph.m_edges[src]){
		ipair edge = MP(src,v);
	        if(lenwalk<taus[edge]){
		    lenwalk=taus[edge];
		}
	    }
	    auto numwalks = uint64(40*lenwalk*log(8*graph.getM()*lenwalk/config.delta)/config.epsilon/config.epsilon);
	    map<uint,double> vals;
	    monteCarlo(src, lenwalk, numwalks, vals, graph);
	    vecER.push_back(vals);
	}
        for (int i = 0; i < graph.getM(); ++i) {
            int q = edgeid[i];
            ipair &edge = seeds[q];
            uint s = edge.first;
            uint t = edge.second;
	        double er = vecER[s][s]+vecER[t][t]-vecER[s][t]-vecER[t][s];
            errors[q] = abs(exact_secs[q]-er);
        }

    }
    else if(config.strAlgo==MCC){
        vector<pair<double,vector<double>>> Eigens;
        map<ipair,int> taus;
        loadEigens(config.strFolder, config.strGraph, graph, Eigens,config.omega);
        calTau(Eigens, taus, graph, config.omega, config.epsilon/2.0);
        vector<int> edgeid(graph.getM());
        iota (begin(edgeid), end(edgeid), 0);
        shuffle(edgeid.begin(), edgeid.end(), std::mt19937(std::random_device()()));
        Timer tm(1, "MonteCarlo-C");
        srand(time(NULL));
	    vector<map<uint,double>> vecER;
        cout << "runnng.." << endl;
	    for(uint src=0; src<graph.getN(); src++){
            int lenwalk = 0;
            for(const auto& v: graph.m_edges[src]){
                ipair edge = MP(src,v);
                if(lenwalk<taus[edge]){
                    lenwalk=taus[edge];
                }
            }
            map<uint,double> vals;
            monteCarlo_C(src, lenwalk, config.epsilon, vals, graph);
            vecER.push_back(vals);
        }
        for (int i = 0; i < graph.getM(); ++i) {
            int q = edgeid[i];
            ipair &edge = seeds[q];
            uint s = edge.first;
            uint t = edge.second;

	    double er = vecER[s][s]+vecER[t][t]-vecER[s][t]-vecER[t][s];
            errors[q] = abs(exact_secs[q]-er);
        }
    }

    if (config.evaflag) {
        double sum_err = 0;
        for (int q = 0; q < errors.size(); q++) {
            sum_err += errors[q];
        }
            cout << "avg-error: " << sum_err / errors.size() << endl;
            Timer::show();
//        }
    } else
        Timer::show();



    return SUCCESS;
}
