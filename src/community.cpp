// File: community.h
// -- community detection source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <libgen.h>

#include "community.h"
#include "functional"
#include "algorithm"
#include "time.h"
using namespace std;

class LogMod{ public: int computeTimes; int passNum;      int nodeNum;       int neighNum; 
                      int src_deg;      int edgeWeight;   int dst;           double dst_tot; 
                      double dst_in;    int commSize;     double communityc; double increase; 
  public: LogMod(){} 
  public: LogMod(int ct, int pn, int nm, int nl, int sd, int ew, int d, int dd, int dc, int cs, double cc, double i){ 
    computeTimes = ct;
    passNum = pn;     nodeNum    = nm; 
    neighNum= nl;     src_deg = sd;  
    edgeWeight = ew;  dst     = d;   
    dst_tot    = dd;  dst_in  = dc;
    commSize   = cs;  
    communityc = cc;  increase   = i ;}};
class LogDst{ public: int computeTimes;  int passNum; 
                      int nodeNum;       int neighNum; 
                      int src_deg;       int edgeWeight; 
                      int dst;           double dst_tot; 
                      double dst_in;     int commSize;
                      double communityc; double increase;
  public: 
    LogDst(){} 
   public: LogDst(int ct, int pn, int nm, int nl, int sd, int ew, int d, int dd, int dc, int cs, double cc, double i){ 
    computeTimes = ct;
    passNum = pn;  nodeNum    = nm;
    neighNum= nl; 
    src_deg = sd;    edgeWeight = ew;   
    dst     = d;     dst_tot    = dd;   
    dst_in  = dc;    commSize   = cs;
    communityc = cc; increase   = i ;}};
// class LogCoe{ public: int computeTimes; int passNum; int commNum, double communityc;
//     public:
//       LogCoe(){}
//     public:LogCoe( int t, int pn, int c, double cc ){
//       computeTimes = t;  passNum = pn;
//       commNum = c;          communityc = cc;
//     }
// };
//class LogMod{ public: double increase; public: LogMod(){} public: LogMod(double i){ increase = i; }};
class LogBest{ public: int computeTimes; int passNum; int nodeNum; double increase; double deltaMod; public: LogBest(){} public: LogBest(int t, int pn, int nm, double i, double dq){ computeTimes = t; passNum = pn; nodeNum = nm; increase = i; deltaMod = dq; }};
class LogDeltaMod{ public: int computeTimes; double increase; public: LogDeltaMod(){} public: LogDeltaMod(int t, double i){ computeTimes = t; increase = i; }};
class LogMod1{ public: int computeTimes; double increase; public: LogMod1(){} public: LogMod1(int t, double i){ computeTimes = t; increase = i; }};

Community::Community(char * filename, char * filename_w, char * filename_coe, int type, int nbp, double minm) {
  g = Graph(filename, filename_w, type);
  size = g.nb_nodes;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  coefficient.resize(size);
  commCoeff.resize(size);
  n2c.resize(size);
  in.resize(size);
  tot.resize(size);
  isChanged.resize(size);
  changedComm.resize(size);

  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    tot[i] = g.weighted_degree(i);
    in[i]  = g.nb_selfloops(i);
    isChanged[i]   = false;
    changedComm[i] = -1;
  }

  nb_pass = nbp;
  min_modularity = minm;
  //read coefficient file 
  ifstream finput_c;
  finput_c.open(filename_coe,fstream::in | fstream::binary);
  for (int i = 0; i < size; ++i)
  {
    finput_c.read((char *)(&coefficient[i]), 8);
    commCoeff[i] = coefficient[i];
  }
}

Community::Community(Graph gc, int nbp, double minm) {
  g = gc;
  size = g.nb_nodes;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  n2c.resize(size);
  in.resize(size);
  tot.resize(size);
  isChanged.resize(size);
  changedComm.resize(size);

  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    tot[i] = g.weighted_degree(i);
    isChanged[i] = false;
    changedComm[i] = -1;
  }

  nb_pass = nbp;
  min_modularity = minm;
}

void
Community::init_partition(char * filename) {
  ifstream finput;
  finput.open(filename,fstream::in);

  // read partition
  while (!finput.eof()) {
    unsigned int node, comm;
    finput >> node >> comm;
    
    if (finput) {
      int old_comm = n2c[node];
      neigh_comm(node);

      remove(node, old_comm, neigh_weight[old_comm]);

      unsigned int i=0;
      for ( i=0 ; i<neigh_last ; i++) {
	unsigned int best_comm     = neigh_pos[i];
	float best_nblinks  = neigh_weight[neigh_pos[i]];
	if (best_comm==comm) {
	  insert(node, best_comm, best_nblinks);
	  break;
	}
      }
      if (i==neigh_last)
	insert(node, comm, 0);
    }
  }
  finput.close();
}

void
Community::display() {
  for (int i=0 ; i<size ; i++)
    cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i] ;
  cerr << endl;
}


double
Community::modularity() {
  double q  = 0.;
  double m2 = (double)g.total_weight;

  for (int i=0 ; i<size ; i++) {
    if (tot[i]>0)
      q += (double)in[i]/m2 - ((double)tot[i]/m2)*((double)tot[i]/m2);
  }

  return q;
}

inline int 
Community::getCommSize( int comm ){
  int commSize = 0;
  double  coef = 0.;      
  for (int i = 0; i < size; ++i)
  {
    int commNum = n2c[i];
    if ( commNum == comm ){
      coef += coefficient[i];
      commSize++;
    }
  }
  commCoeff[comm] = (double)coef / commSize;
  return commSize;
}

inline int
Community::neigh_comm( unsigned int node ) {
  for (unsigned int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1;
  neigh_last=0;

  pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(node);

  unsigned int deg = g.nb_neighbors(node);

  neigh_pos[0]=n2c[node];
  neigh_weight[neigh_pos[0]]=0;
  neigh_last=1;

  for (unsigned int i=0 ; i<deg ; i++) {
    unsigned int neigh        = *(p.first+i);
    unsigned int neigh_comm   = n2c[neigh];
    double neigh_w = (g.weights.size()==0)?1.:*(p.second+i);
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
        neigh_weight[neigh_comm]=0.;
        neigh_pos[neigh_last++]=neigh_comm;
    }
      neigh_weight[neigh_comm]+=neigh_w;
    }
  }
  return deg;
}

inline bool
Community::neigh_comm( unsigned int node, unsigned int nb_pass_done ) {
  bool isCompute = false;
  for (unsigned int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1;
  neigh_last=0;

  pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(node);

  unsigned int deg = g.nb_neighbors(node);

  neigh_pos[0]=n2c[node];
  neigh_weight[neigh_pos[0]]=0;
  neigh_last=1;

  for (unsigned int i=0 ; i<deg ; i++) {
    unsigned int neigh        = *(p.first+i);
    unsigned int neigh_comm   = n2c[neigh];
    double neigh_w = (g.weights.size()==0)?1.:*(p.second+i);
    if( nb_pass_done !=1 && isChanged[neigh_comm] ){
        isCompute = true;
        //continue;
    }
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
        neigh_weight[neigh_comm]=0.;
        neigh_pos[neigh_last++]=neigh_comm;
    }
      neigh_weight[neigh_comm]+=neigh_w;
    }
  }
  if(nb_pass_done == 1) return true;
  return isCompute;
  //return deg;
}

void
Community::partition2graph() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;


  for (int i=0 ; i<size ; i++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(i);

    int deg = g.nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
    }
  }
}

void
Community::display_partition() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  for (int i=0 ; i<size ; i++)
    cout << i << " " << renumber[n2c[i]] << endl;
}


Graph
Community::partition2graph_binary() {
  // Renumber communities
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  // Compute communities
  vector<vector<int> > comm_nodes(final);
  for (int node=0 ; node<size ; node++) {
    comm_nodes[renumber[n2c[node]]].push_back(node);
  }

  // Compute weighted graph
  Graph g2;
  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(comm_nodes.size());

  int comm_deg = comm_nodes.size();
  for (int comm=0 ; comm<comm_deg ; comm++) {
    map<int,float> m;
    map<int,float>::iterator it;

    int comm_size = comm_nodes[comm].size();
    for (int node=0 ; node<comm_size ; node++) {
      pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(comm_nodes[comm][node]);
      int deg = g.nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
	int neigh        = *(p.first+i);
	int neigh_comm   = renumber[n2c[neigh]];
	double neigh_weight = (g.weights.size()==0)?1.:*(p.second+i);

	it = m.find(neigh_comm);
	if (it==m.end())
	  m.insert(make_pair(neigh_comm, neigh_weight));
	else
	  it->second+=neigh_weight;
      }
    }
    g2.degrees[comm]=(comm==0)?m.size():g2.degrees[comm-1]+m.size();
    g2.nb_links+=m.size();

    
    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight  += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
  }

  return g2;
}


bool
Community::one_level( int level, char * fileName ) {
  vector<LogMod> logmod;
  vector<LogBest> logbest;
  vector<LogDeltaMod> logdelta;
  vector<LogDst> logdst;
  vector<LogMod1> logmod1;

  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  double new_mod   = modularity();
  double cur_mod   = new_mod;

  vector<int> random_order(size);
  for (int i=0 ; i<size ; i++)
    random_order[i]=i;
  for (int i=0 ; i< (size-1) ; i++) {
    int rand_pos = rand()%(size-i)+i;
    int tmp      = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // repeat while 
  //   there is an improvement of modularity
  //   or there is an improvement of modularity greater than a given epsilon 
  //   or a predefined number of pass have been done
  int nodeTimes = 0;
  int t = 0;
  int computeTimes = 0;
  double before = cur_mod;
  do {
    cur_mod = new_mod;
    nb_moves = 0;
    nb_pass_done++;
    //cerr<< "Pass::" << nb_pass_done << endl;
    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
      //int node = node_tmp;
      int node = random_order[node_tmp];
      int node_comm     = n2c[node];
      double w_degree = g.weighted_degree(node);
      // computation of all neighboring communities of current node
      bool isCompute = neigh_comm(node, nb_pass_done);
      if( !isCompute ) continue;
      //int deg = neigh_comm(node);
      // remove node from its current community
      remove(node, node_comm, neigh_weight[node_comm]);
      int    best_comm     = node_comm;
      double best_nblinks  = neigh_weight[best_comm];
      double best_increase = modularity_gain(node, best_comm, best_nblinks, w_degree);
      computeTimes++;
      //if (level == 0)logmod.push_back(LogMod( t, nb_pass_done, node, neigh_last, g.nb_neighbors(node) , best_nblinks, best_comm, tot[best_comm], in[best_comm], getCommSize(best_comm),commCoeff[best_comm] ,best_increase));
      nodeTimes++;
      bool isNum = false;
      // compute the nearest community for node
      // default choice for future insertion is the former community
      for (unsigned int i = 1 ; i < neigh_last; i++) {
          int    theN2C    = neigh_pos[i];
          double theWeight = neigh_weight[theN2C];
          double increase  = modularity_gain(node, theN2C, theWeight, w_degree);
          computeTimes++;
         // if (level == 0)logmod.push_back(LogMod( t, nb_pass_done, node, neigh_last,  g.nb_neighbors(node), theWeight, theN2C, tot[theN2C], in[theN2C], getCommSize(theN2C),commCoeff[theN2C] , increase));
          nodeTimes++;
          if(nodeTimes % 10000 == 0) isNum = true;
          //logmod.push_back(LogMod( nb_pass_done, node, deg, theWeight, increase));
          if (increase > best_increase) {
              best_comm     = theN2C;
              best_nblinks  = theWeight;
              best_increase = increase;
          }
      }

      // insert node in the nearest community
      insert(node, best_comm, best_nblinks);
      if (isNum && level == 0){
         double newMod = modularity();
         logdelta.push_back(LogDeltaMod( t, (newMod-before)));
         logmod1.push_back(LogMod1( t, (newMod)));
         t++;
         before = newMod;
         nodeTimes = 0;
         isNum = false;
      }

      if (best_comm!=node_comm){
          //if (level == 0)logdst.push_back(LogDst( t, nb_pass_done, node, neigh_last, g.nb_neighbors(node) , best_nblinks, best_comm, tot[best_comm], in[best_comm], getCommSize(best_comm), commCoeff[best_comm] , best_increase));
        changedComm[best_comm] = 1;
        changedComm[node_comm] = 1;
        nb_moves++;
      }
    }
    isChanged.clear();
    for (int i = 0; i < size; ++i)
    {
      int changed = changedComm[i];
      isChanged[i] = ( changed > 0 ) ? true : false;
      changedComm[i] = -1;
    }
    new_mod = modularity();
    //cerr<< " Pass" << nb_pass_done <<", computeTimes::"<< computeTimes <<endl;
    if (nb_moves>0)
      improvement=true;
    
  } while (nb_moves>0 && new_mod-cur_mod>min_modularity);
  if (level == 0)
  {
    ofstream csvdelta;//csvbest, csvdelta, sv1, csv2, csvDetail, csvmove;
    time_t  nowtime = time(NULL);  
    struct  tm  *p;  
    p = gmtime(&nowtime);  
    char    filename[256] = {0};
    // sprintf(filename, "/Users/donaldchi/Dropbox/x10/Data/5.8/%s-%d-deltaMod-original-time.csv", basename(fileName), level);  
    // csvdelta.open(filename);
    // for (vector<LogDeltaMod>::iterator iter1 = logdelta.begin(); iter1 != logdelta.end(); ++iter1) {
    //   LogDeltaMod entry1 = *iter1;
    //   csvdelta << entry1.increase << endl;
    // }
    // csvdelta.close();
    // sprintf(filename, "/Users/donaldchi/Dropbox/x10/Data/5.8/%s-%d-mod-original-time.csv", basename(fileName), level);  
    // csvdelta.open(filename);
    // for (vector<LogMod1>::iterator iter1 = logmod1.begin(); iter1 != logmod1.end(); ++iter1) {
    //   LogMod1 entry1 = *iter1;
    //   csvdelta << entry1.increase << endl;
    // }
    // csvdelta.close();
    // sprintf(filename, "/Users/donaldchi/Dropbox/x10/Data/5.8/%s-%d-all-original21.csv", basename(fileName), level);  
    // csvdelta.open(filename);
    // for (vector<LogMod>::iterator iter1 = logmod.begin(); iter1 != logmod.end(); ++iter1) {
    //   LogMod entry1 = *iter1;
    //   csvdelta << entry1.computeTimes << "," << entry1.passNum << "," << entry1.nodeNum    << "," 
    //            << entry1.neighNum << ","
    //            << entry1.src_deg << "," << entry1.edgeWeight << ","
    //            << entry1.dst     << "," << entry1.dst_tot    << "," 
    //            << entry1.dst_in  << "," << entry1.commSize << ","
    //            << entry1.communityc  << "," << entry1.increase <<endl;
    // }
    // csvdelta.close();
    // sprintf(filename, "/Users/donaldchi/Dropbox/x10/Data/5.8/%s-%d-moved-original21.csv", basename(fileName), level);  
    // csvdelta.open(filename);
    // for (vector<LogDst>::iterator iter1 = logdst.begin(); iter1 != logdst.end(); ++iter1) {
    //   LogDst entry1 = *iter1;
    //  // csvdelta << entry1.increase << endl;
    //   csvdelta << entry1.computeTimes << "," << entry1.passNum << "," << entry1.nodeNum    << "," 
    //            << entry1.neighNum << ","
    //            << entry1.src_deg << "," << entry1.edgeWeight << ","
    //            << entry1.dst     << "," << entry1.dst_tot    << "," 
    //            << entry1.dst_in  << "," << entry1.commSize << ","
    //            << entry1.communityc  << "," << entry1.increase <<endl;
    // }
    // csvdelta.close();
    //     sprintf(filename, "/Users/donaldchi/Dropbox/x10/Data/5.8/%s-%d-deltaQ-original2.csv", basename(fileName), level);  
    // csvdelta.open(filename);
    // for (vector<LogDeltaMod>::iterator iter1 = logdelta.begin(); iter1 != logdelta.end(); ++iter1) {
    //   LogDeltaMod entry1 = *iter1;
    //   csvdelta << entry1.computeTimes <<","<<entry1.increase << endl;
    // }
    // csvdelta.close();  
    //     sprintf(filename, "/Users/donaldchi/Dropbox/x10/Data/5.8/%s-%d-best-original2.csv", basename(fileName), level);  
    // csvdelta.open(filename);
    // for (vector<LogBest>::iterator iter1 = logbest.begin(); iter1 != logbest.end(); ++iter1) {
    //   LogBest entry1 = *iter1;
    //   csvdelta << entry1.computeTimes <<"," << entry1.passNum << "," << entry1.nodeNum << "," <<entry1.increase <<"," << (double)entry1.deltaMod << endl;
    // }
    // csvdelta.close();  
  }
  //cerr<<"neigh_num:: "<<neigh_num<<endl;
  return improvement;
}

