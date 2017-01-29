#include <vector>

#include "synarmosma/variety.h"
#include "synarmosma/lattice.h"
#include "synarmosma/polynomial.h"
#include "synarmosma/functional_equation.h"
#include "synarmosma/directed_graph.h"

extern SYNARMOSMA::Random RND;

int main(int argc,char** argv)
{
  RND.set_seed(12);
  int i,j,m,n = 25;
  double p = 0.15;

  SYNARMOSMA::Polynomial<SYNARMOSMA::Rational> p1(2),p2(3);
  SYNARMOSMA::Variety<SYNARMOSMA::Rational> v(2);
  SYNARMOSMA::Functional_Equation<SYNARMOSMA::Rational> feqn(3);

  std::cout << feqn << std::endl;
  std::cout << p1 << std::endl;
  std::cout << p2 << std::endl;
  std::cout << p1*p2 << std::endl;
  std::cout << p1 + p2 << std::endl;

  std::cout << v << std::endl;
  std::ofstream s1("polynomial.dat",std::ios::out | std::ios::binary);
  v.serialize(s1);
  s1.close();
  
  std::ifstream s2("polynomial.dat",std::ios::in | std::ios::binary);
  v.deserialize(s2);
  s2.close();
  std::cout << v << std::endl;
  std::system("rm -f polynomial.dat");
  
  SYNARMOSMA::Rational q1(1,5);
  SYNARMOSMA::Rational q2(2,3);
  std::cout << q1 + q2 << std::endl;

  SYNARMOSMA::Lattice L(4);
  SYNARMOSMA::Poset P;
  P.power_set(7);
  std::cout << P.consistent() << "  " << P.totality() << std::endl;
 
  // The Petersen graph
  SYNARMOSMA::Graph G(10);
  // Outer pentagon
  G.add_edge(0,1); G.add_edge(0,4); G.add_edge(0,5);
  G.add_edge(1,2); G.add_edge(1,6);
  G.add_edge(2,3); G.add_edge(2,7);
  G.add_edge(3,4); G.add_edge(3,8);
  G.add_edge(4,9);
  // Inner pentagram
  G.add_edge(5,7); G.add_edge(5,8);
  G.add_edge(6,8); G.add_edge(6,9);
  G.add_edge(7,9);

  std::vector<SYNARMOSMA::Monomial<int> > output;
  G.tutte_polynomial(output);
  n = (signed) output.size();
  for(i=0; i<n; ++i) {
    m = (signed) output[i].exponents.size();
    if (output[i].coefficient != 1) {
      std::cout << output[i].coefficient << "*";
    }
    if (m == 1) {
      if (output[i].exponents[0].first == 0) {
        std::cout << "x^(" << output[i].exponents[0].second << ")";
      }
      else {
        std::cout << "y^(" << output[i].exponents[0].second << ")";
      }
    }
    else {
      std::cout << "x^(" << output[i].exponents[0].second << ")*y^(" << output[i].exponents[1].second << ")";
    }
    std:: cout << " + " << std::endl;
  }


  SYNARMOSMA::Pseudograph CG1(4);
  CG1.add_edge(0,1); CG1.add_edge(0,2); CG1.add_edge(0,3);
  CG1.add_edge(1,2);
  CG1.add_edge(2,3);
  SYNARMOSMA::Pseudograph* CG2 = new SYNARMOSMA::Pseudograph(4);
  CG2->add_edge(0,1); CG2->add_edge(0,2); CG2->add_edge(0,3);
  CG2->add_edge(1,2);
  CG2->add_edge(2,3);
  CG1.contract(0,1,CG2);
  std::vector<int> cvector;
  int nb = CG2->get_candidates(cvector);
  std::cout << "There are " << nb << " bridges." << std::endl;
  std::cout << "There are " << cvector.size()/2 << " candidate edges." << std::endl;
  for(i=0; i<(signed) cvector.size()/2; ++i) {
    std::cout << cvector[2*i] << ":" << cvector[2*i+1] << std::endl;
  }

  std::string type = "ring";
  SYNARMOSMA::Graph G2(n,type);

  double c0 = G2.clustering_coefficient();
  double l0 = G2.mean_path_length();
  for(i=0; i<n; ++i) {
    if (RND.drandom() < p) {
      G2.drop_edge(i,(i+1)%n);
      do {
        j = RND.irandom(n);
        if (i == j) continue;
        if (G2.add_edge(i,j)) break;
      } while(true);
    }
  }
  std::cout << c0 << "  " << G2.clustering_coefficient() << std::endl; 
  std::cout << l0 << "  " << G2.mean_path_length() << "  " << G2.mean_path_length()/l0 << std::endl;

  SYNARMOSMA::Directed_Graph G3;
  for(int i=0; i<7; ++i) {
    G3.add_vertex();
  }
  G3.add_edge(0,1,SYNARMOSMA::OUTGOING,10.0);
  G3.add_edge(0,2,SYNARMOSMA::OUTGOING,5.0);
  G3.add_edge(0,3,SYNARMOSMA::OUTGOING,4.0);
  G3.add_edge(1,2,SYNARMOSMA::OUTGOING,3.0);
  G3.add_edge(1,4,SYNARMOSMA::OUTGOING,5.0);
  G3.add_edge(2,3,SYNARMOSMA::OUTGOING,5.0);
  G3.add_edge(2,4,SYNARMOSMA::OUTGOING,2.0);
  G3.add_edge(3,5,SYNARMOSMA::OUTGOING,8.0);
  G3.add_edge(4,6,SYNARMOSMA::OUTGOING,7.0);
  G3.add_edge(5,2,SYNARMOSMA::OUTGOING,3.0);
  G3.add_edge(5,6,SYNARMOSMA::OUTGOING,11.0);
  
  std::cout << G3.bipartite() << std::endl;
  std::cout << G3.entwinement() << std::endl;
  std::cout << G3.eccentricity(3) << std::endl;
  std::vector<double> foo;
  G3.vertex_centrality(foo,0.000001);
  for(int i=0; i<7; ++i) {
    std::cout << foo[i] << std::endl;
  }

  std::pair<int,int> pr;
  SYNARMOSMA::edge_hash dmap;
  SYNARMOSMA::edge_hash::const_iterator qt;
  G3.compute_distances(dmap);
  for(int i=0; i<7; ++i) {
    for(int j=0; j<7; ++j) {
      if (i == j) continue;
      pr.first = i; pr.second = j;
      qt = dmap.find(pr);
      std::cout << "Distance between " << i << " and " << j << " is " << qt->second << std::endl;
    }
  }
  std::cout << G3.compute_flow(0,6) << std::endl;

  SYNARMOSMA::Directed_Graph G4; 
  int d1 = SYNARMOSMA::OUTGOING,d2 = SYNARMOSMA::INCOMING,d3 = SYNARMOSMA::UNDIRECTED;
  for(int i=0; i<6; ++i) {
    G4.add_vertex();
  }
  G4.add_edge(0,1,d1);
  G4.add_edge(1,2,d1);
  G4.add_edge(1,3,d3);
  G4.add_edge(3,4,d3);
  G4.add_edge(1,4,d3);
  G4.add_edge(2,3,d2);
  G4.add_edge(3,5,d3);
  G4.add_edge(2,5,d3);
  std::cout << G4.distance(0,4) << std::endl; 
  std::cout << G4.distance(0,2) << std::endl;
  std::cout << G4.distance(2,3) << std::endl;
  std::cout << G4.distance(3,2) << std::endl;
  std::cout << G4.directedness() << "  " << double(G4.directedness())/double(G4.size()) << std::endl;

  
  SYNARMOSMA::Directed_Graph G5(n,0.25);
  std::cout << "Graph initialized" << std::endl;
  std::cout << G5.connected() << std::endl;
  
  std::set<int> S;
  G5.compute_sinks(S);
  std::cout << S.size() << std::endl;
  
  for(int i=0; i<n; ++i) {
    for(int j=1+i; j<n; ++j) {
      std::cout << i << ":" << j << "  " << G5.path_connected(i,j) << std::endl;
    }
  }
  
  int q = G5.directedness();
  std::cout << q << "  " << double(q)/double(G5.size()) << std::endl; 
  
  SYNARMOSMA::Graph G6(3);
  std::cout << G6.connected() << "  " << G6.circuit_rank() << std::endl;
  std::cout << G6.order() << "  " << G6.size() << std::endl;

  std::cout << G6.stellar_addition(0) << std::endl;
  std::cout << G6.connected() << "  " << G6.circuit_rank() << std::endl;
  std::cout << G6.order() << "  " << G6.size() << std::endl;
  std::cout << G6.stellar_deletion(n) << std::endl;
  std::cout << G6.connected() << "  " << G6.circuit_rank() << std::endl;
  std::cout << G6.order() << "  " << G6.size() << std::endl;
  return 0;
}

