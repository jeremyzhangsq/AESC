#include "algo.h"



using namespace std;

void monteCarlo(uint src, int len_walk, uint64 n_walk, map<uint,double>& mapVals, const Graph& graph){
	fastSrand();
	double x=0;
    vector<double> vec(graph.getN(),0);
    for(uint len=0; len<len_walk; len++){
	    for(uint64 i=0; i<n_walk; i++){
	    	uint cur = src;
	        for(uint j=0; j<len; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint k = fastRand()%deg;
	        	cur = graph.m_edges[cur][k];
	        }
		vec[cur]+=1.0/n_walk/graph.getDeg(cur);
	    }
	}

        mapVals[src]=vec[src];
	for(const auto& v: graph.m_edges[src]){
	   mapVals[v]= vec[v];
	}
}

void monteCarlo_C(uint src, int len_walk, double eps, map<uint,double>& mapVals, const Graph& graph){
    fastSrand();
    double x=0;
    for(uint len=1; len<len_walk; len++){
        double beta = pow(0.5,len); // this value is unknown
        uint64 n_walk = 20000*(pow(len_walk,1.5)*sqrt(beta)/eps + pow(len_walk,3)*pow(beta,1.5)/eps/eps)/len_walk;
        vector<double> xs(graph.getN(),0);
        for(uint64 i=0; i<n_walk; i++){
            uint cur = src;
            for(uint j=0; j<ceil(len/2); j++){
                uint deg = graph.m_deg[cur];
                uint k = fastRand()%deg;
                cur = graph.m_edges[cur][k];
                xs[cur]+=1.0/sqrt(graph.m_deg[cur])/n_walk;
            }
	}

        for(uint64 i=0; i<n_walk; i++){
            uint cur = src;
            for(uint j=0; j<floor(len/2); j++){
                uint deg = graph.m_deg[cur];
                uint k = fastRand()%deg;
                cur = graph.m_edges[cur][k];
                mapVals[src]+=xs[cur]*1.0/sqrt(graph.m_deg[cur])/n_walk;
            }
        }
        
	for(const auto& tgt: graph.m_edges[src]){
        
	    for(uint64 i=0; i<n_walk; i++){
                uint cur = tgt;
                for(uint j=0; j<floor(len/2); j++){
                    uint deg = graph.m_deg[cur];
                    uint k = fastRand()%deg;
                    cur = graph.m_edges[cur][k];
                    mapVals[tgt]+=xs[cur]*1.0/sqrt(graph.m_deg[cur])/n_walk;
                }
            }
	}
    }
}



double get_delta(uint u, uint v, const Graph &graph, int omega, vector<pair<double,vector<double>>> &Eigens, int t){
    double delta = 0;
    for (int i = 1; i < omega-1; ++i) {
        auto val = Eigens[i].first;
        auto & f = Eigens[i].second;
        delta += pow((f[u] - f[v]),2) * pow(val, (t + 1)) / (1 - val);
    }
    return delta / 2.0 / graph.getM();
}


double get_upsilon(uint u, uint v, const Graph &graph, int omega, vector<pair<double,vector<double>>> &Eigens){
    double upsilon = 0;
    for (int i = 1; i < omega-1; ++i) {
        auto val = Eigens[i].first;
        auto & f = Eigens[i].second;
        upsilon += pow((f[u] - f[v]),2) * (1 + val);
    }
    return upsilon / 2.0 / graph.getM();
}


int get_tau(uint u, uint v, const Graph &graph, double epsilon,double lamba,double delta,double upsilon){
    int du = graph.getDeg(u);
    int dv = graph.getDeg(v);
    double eps_delta = epsilon-delta>0?epsilon-delta:epsilon;
    double a = log(max((1.0/du+1.0/dv-2.0/du/dv-upsilon)/eps_delta/(1-pow(lamba,2)),1.0));
    double b = log(1.0/abs(lamba));
    int tau = max(ceil(a / b - 1),1.0);
    return tau % 2 == 0 ? int(tau + 1) : int(tau);
}


int calTau(uint u, uint v, vector<pair<double,vector<double>>> &Eigens,const Graph& graph, int omega, double epsilon){
    double delta= 0;
    double upsilon = 0;
    double lamba = Eigens[1].first;
    int tau = get_tau(u,v,graph,epsilon,lamba,delta,upsilon);
    int t=1;
    while(true){
        delta = get_delta(u,v,graph,omega,Eigens,t);
        upsilon = get_upsilon(u,v,graph,omega,Eigens);
        lamba = Eigens[omega-1].first;
        int tauprime = get_tau(u,v,graph,epsilon,lamba,delta,upsilon);
        if (t<=tauprime and tauprime < tau){
            tau = tauprime;
            t += 2;
        }
        else
            break;
    }

    return tau;
}

void calTau(vector<pair<double,vector<double>>> &Eigens, map<ipair,int>& Taus, const Graph& graph, int omega, double epsilon){
    for (uint u = 0; u < graph.getN(); ++u) {
        for(const auto& v: graph.m_edges[u]){
            int tau = calTau(u,v,Eigens,graph, omega, epsilon);
            ipair p = MP(u,v);
            Taus[p] = tau;
        }
    }
}



void push_single(uint src, const Graph &graph, vector<double> &pvec, vector<uint> &S, double ds, vector<double> &tmp_hvec) {
    vector<uint> tmpS;
    vector<bool> tmpSflag(graph.getN(),false);
    for(unsigned int v : S){
        double residue = pvec[v];
        pvec[v]=0;
        for(const auto& u: graph.m_edges[v]){
            if(!tmpSflag[u]){
                tmpS.push_back(u);
                tmpSflag[u] = true;
            }
            double update = residue/(double)graph.m_deg[u];
            pvec[u] += update;
        }
    }
    S = tmpS;
    for(const auto& vj: graph.m_edges[src]){
        double update = (pvec[src]-pvec[vj])/ds;
        tmp_hvec[vj] += update;
    }
}


void truncatedGraphTrav(const Graph& graph, map<ipair,double>& preds, map<ipair,int> &taus){
    std::vector<map<uint,double>> all_hvecs;
    all_hvecs.resize(graph.getN());
    double start = clock();
    int cnt = 0;
    for (int src = 0; src < graph.getN(); ++src) {
        int maxTau = 0;
        for(uint tgt: graph.m_edges[src]){
            maxTau = max(maxTau,taus[{src,tgt}]);
        }
        map<uint,double> hvec;
        std::vector<double> pvec(graph.getN(),0);
        pvec[src] = 1.0;
        vector<uint> S;
        auto ds = (double)graph.getDeg(src);
        S.push_back(src);
        std::vector<double> tmp_hvec(graph.getN(),0);
        for (const auto& v: graph.m_edges[src])
            tmp_hvec[v] = 1.0/ds;
        for(int i=0;i<maxTau;i++){
            push_single(src, graph, pvec, S,  ds, tmp_hvec);
        }
        for(const auto& vj: graph.m_edges[src]){
            hvec[vj]=tmp_hvec[vj];
        }
        all_hvecs[src] = hvec;
    }
    for (int u = 0; u < graph.getN(); ++u) {
        for (uint v: graph.m_edges[u]) {
            preds[{u,v}] = all_hvecs[u][v]+all_hvecs[v][u];
        }
    }
}


double calEdgeMax(const Graph &graph, vector<uint> &S, vector<double> &pvec, const double &global_max, int gamma) {

    int nnz_size = S.size();
    int real_gamma = min(gamma,nnz_size);

    nth_element(S.begin(), S.begin()+real_gamma-1, S.end(), [&](const uint A, const uint B) -> bool {
        return pvec[A] > pvec[B];});
    double gamma_max = 1;
    vector<bool> candidates_exist(graph.getN(), false);
    for (int i = 0; i < real_gamma; ++i){
        auto &node = S[i];
        auto &val = pvec[node];
        gamma_max = min(gamma_max,val);
        candidates_exist[node]=true;
    }

    gamma_max = gamma<nnz_size ? gamma_max : 0;
    double edge_max = global_max + gamma_max;

    for (int i = 0; i < real_gamma; ++i){
        auto &u = S[i];
        for (const auto& v:graph.m_edges[u]){
            if (!candidates_exist[v])
                continue;
            edge_max = max(edge_max, pvec[u] + pvec[v]);
        }
    }
    return edge_max;
}

double calChi(uint vi, uint vj, const Graph& graph, vector<double> &pvec, vector<uint> &S, int len_walk, double &global_min, double &global_max, double &edge_max){



    double src_local_max = 0;
    double src_local_min = MAXFLOAT;
    for (auto val:graph.m_edges[vi]){
        src_local_max = max(src_local_max, pvec[val]);
        src_local_min = min(src_local_min, pvec[val]);
    }
    double tgt_local_max = 0;
    double tgt_local_min = MAXFLOAT;
    for (auto val:graph.m_edges[vj]){
        tgt_local_max = max(tgt_local_max, pvec[val]);
        tgt_local_min = min(tgt_local_min, pvec[val]);
    }

    double chi = global_max + (src_local_max + tgt_local_max) / 2.0 + (len_walk - 1) * edge_max - src_local_min - tgt_local_min - 2 * (len_walk - 1) * global_min;
    return chi;
}

double RW_single(uint src,const Graph& graph, double num_walk, int len_walk, vector<double> &pvec){
    fastSrand();
    double x=0;
    for (int i = 0; i < num_walk; ++i) {
        uint cur = src;
        for(uint len=0; len<len_walk; len++){
            uint deg = graph.m_deg[cur];
            uint k = fastRand()%deg;
            cur = graph.m_edges[cur][k];
            x += pvec[cur];
        }
    }
    return x;
}

double RW_double(uint src, uint v, const Graph& graph, double num_walk, int len_walk, vector<double> &pvec){
    fastSrand();
    double x=0;
    for (int i = 0; i < num_walk; ++i) {
        uint cur = src;
        uint cur2 = v;
        for(uint len=0; len<len_walk; len++){
            uint k = fastRand()%graph.m_deg[cur];
            cur = graph.m_edges[cur][k];
            k = fastRand()%graph.m_deg[cur2];
            cur2 = graph.m_edges[cur2][k];
            x += (pvec[cur]-pvec[cur2]);
        }
    }
    return x;
}

void truncatedGraphTravPlus(const Graph& graph, map<ipair,double>& preds, map<ipair,int> &taus, Config &config){
    vector<map<uint,double>> all_hvecs;
    all_hvecs.resize(graph.getN());
//    double pi_time = 0, rw_time = 0, chi_time=0;
    int cnt = 0;
    double numrw = 0;
    for (int src = 0; src < graph.getN(); ++src) {

        int ell = 0;
        map<uint,double> hvec;
        vector<double> pvec(graph.getN(),0);
        pvec[src] = 1.0;
        vector<uint> S;
        auto ds = (double)graph.getDeg(src);
        S.push_back(src);
        vector<double> tmp_hvec(graph.getN(),0);
        for (const auto& v: graph.m_edges[src])
            tmp_hvec[v] = 1.0/ds;

        double global_min,global_max;
        while(true){
            global_min = 1;
            global_max = 0;
            double picost = 0,rwcost = 0;
//            double start = clock();
            push_single(src, graph, pvec, S,  ds, tmp_hvec);
//            pi_time += time_by(start);
            ell++;
            for(auto v: S){
                picost+=(double)graph.m_deg[v];
                auto &val = pvec[v];
                global_max = max(global_max,val);
                global_min = min(global_min,val);
            }
            global_min = S.size()<graph.getN()?0:global_min;
//            double edge_max = calEdgeMax(graph, S,pvec, global_max, config.gamma);
            for (const auto& v: graph.m_edges[src]){
                int rest_ell = taus[{src,v}]-ell;
                if (rest_ell<=0)
                    continue;
                double chi = 2.0*rest_ell*(global_max-global_min);
//                double chi = calChi(src, v, graph, pvec, S, rest_ell, global_min,global_max,edge_max, config.gamma);
                double nr = max(ceil(8*pow(chi,2)*log(2.0*graph.getM()/config.delta)/pow(ds*config.epsilon,2)),1.0);
                rwcost += nr;
            }
            if (picost>=25*rwcost)
                break;
        }

        double edge_max = calEdgeMax(graph, S,pvec, global_max, config.gamma);
        if (config.gamma==1)
            edge_max = 2*global_max;

        for (const auto& v: graph.m_edges[src]){
            int rest_ell = taus[{src,v}]-ell;
            if (rest_ell<=0)
                continue;
//            double start = clock();
            double chi = calChi(src, v, graph, pvec, S, rest_ell, global_min,global_max,edge_max);
            if (config.gamma<=0)
                chi = 2*rest_ell*global_max;
//            chi_time += time_by(start);
            double nr = max(ceil(8*pow(chi,2)*log(2.0*graph.getM()/config.delta)/pow(ds*config.epsilon,2)),1.0);
            numrw += nr;
//            start = clock();
            double xij = RW_double(src, v, graph, nr, rest_ell, pvec);
//            rw_time += time_by(start);
            tmp_hvec[v] += xij/(ds*nr); //(xi-xj)/(ds*nr);
//            tmp_hvec[v] += (xi-xj)/(ds*nr);
        }

        for(const auto& vj: graph.m_edges[src]){
            hvec[vj]=tmp_hvec[vj];
        }
        all_hvecs[src] = hvec;
    }

    for (int u = 0; u < graph.getN(); ++u) {
        for (uint v: graph.m_edges[u]) {
            preds[{u,v}] = all_hvecs[u][v]+all_hvecs[v][u];
        }
    }
    cout<<"num_rw: "<<numrw<<endl;

}
