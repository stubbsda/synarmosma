#include <vector>

#include "synarmosma/lattice.h"
#include "synarmosma/polynomial.h"
#include "synarmosma/directed_graph.h"

SYNARMOSMA::Random RND;

int main(int argc,char** argv)
{
  RND.set_seed(12);
  int i,j,n = 25;
  auto d1 = SYNARMOSMA::Relation::before,d2 = SYNARMOSMA::Relation::after,d3 = SYNARMOSMA::Relation::disparate;
  double p = 0.15;
  std::string type,name;
  std::vector<double> foo;
  std::pair<int,int> pr;
  SYNARMOSMA::pair_index dmap;
  SYNARMOSMA::pair_index::const_iterator qt;

  std::cout << "Testing rationals..." << std::endl;
  SYNARMOSMA::Rational q1(1,5);
  SYNARMOSMA::Rational q2(2,3);
  assert((q1 + q2) == SYNARMOSMA::Rational(13,15));

  std::cout << "Testing rational polynomials..." << std::endl;
  SYNARMOSMA::Polynomial<SYNARMOSMA::Rational> p1(2),p2(3);
  SYNARMOSMA::Polynomial<SYNARMOSMA::Rational> r = p1 + p2;
  assert(r.get_degree() == 3);
  r = p1*p2;
  assert(r.get_degree() == 5);
  std::ofstream s1("polynomial.dat",std::ios::out | std::ios::binary);
  p2.serialize(s1);
  s1.close();
  std::ifstream s2("polynomial.dat",std::ios::in | std::ios::binary);
  p1.deserialize(s2);
  s2.close();
  std::system("rm -f polynomial.dat");
  assert(p1 == p2);  

  std::cout << "Testing lattices and posets..." << std::endl;
  SYNARMOSMA::Lattice L(4);
  SYNARMOSMA::Poset P;
  P.power_set(7);
  assert(P.consistent());
 
  // The Petersen graph
  name = "petersen";
  std::cout << "Testing Tutte polynomial computation..." << std::endl;
  SYNARMOSMA::Graph G(name);
  G.entwinement();
  G.vertex_centrality(foo,0.000001);
  std::vector<SYNARMOSMA::Monomial<int> > output;
  G.tutte_polynomial(output);
  for(i=0; i<(signed) output.size(); ++i) {
    if (output[i].exponents.size() != 2) continue;
    if (output[i].exponents[0].second == 1 && output[i].exponents[1].second == 3) assert(output[i].coefficient == 65);
    if (output[i].exponents[0].second == 2 && output[i].exponents[1].second == 2) assert(output[i].coefficient == 105);
  }

  std::cout << "Testing pseudographs..." << std::endl;
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
  assert(nb == 0);
  assert(cvector.size() == 8);
  assert(cvector[0] == 0);
  assert(cvector[1] == 1);
  assert(cvector[2] == 0);
  assert(cvector[3] == 2);
  assert(cvector[4] == 0);
  assert(cvector[5] == 1);
  assert(cvector[6] == 1);
  assert(cvector[7] == 2);
  delete CG2;

  std::cout << "Testing clustering coefficients..." << std::endl;
  type = "cyclic";
  SYNARMOSMA::Graph G2(n,type);
  G2.clustering_coefficient();
  G2.mean_path_length();
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
  G2.clustering_coefficient(); 
  G2.mean_path_length();

  std::cout << "Testing directed graphs..." << std::endl;
  SYNARMOSMA::Directed_Graph G3;
  for(i=0; i<7; ++i) {
    G3.add_vertex();
  }
  G3.add_edge(0,1,d1,10.0);
  G3.add_edge(0,2,d1,5.0);
  G3.add_edge(0,3,d1,4.0);
  G3.add_edge(1,2,d1,3.0);
  G3.add_edge(1,4,d1,5.0);
  G3.add_edge(2,3,d1,5.0);
  G3.add_edge(2,4,d1,2.0);
  G3.add_edge(3,5,d1,8.0);
  G3.add_edge(4,6,d1,7.0);
  G3.add_edge(5,2,d1,3.0);
  G3.add_edge(5,6,d1,11.0);
  assert(G3.bipartite() == 1); 
  assert(G3.maximum_parents() == 3);
  assert(!G3.DAG());
  assert(!G3.singly_connected());
  assert(G3.eccentricity(3) == 2);
  G3.compute_distances(dmap);
  pr.first = 0; pr.second = 6;
  qt = dmap.find(pr);
  assert(qt->second == 3);
  pr.first = 3; pr.second = 1;
  qt = dmap.find(pr);
  assert(qt->second == -1);
  pr.first = 5; pr.second = 6;
  qt = dmap.find(pr);
  assert(qt->second == 1);
  assert(G3.compute_flow(0,6) == 15);

  SYNARMOSMA::Directed_Graph G4; 
  for(i=0; i<6; ++i) {
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
  G4.compute_directedness();
  assert(G4.distance(0,4) == -1);
  assert(G4.distance(0,2) == 2);
  assert(G4.distance(2,3) == -1);
  assert(G4.distance(3,2) == 1);
  assert(G4.directedness() == 3); 
  assert(G4.size() == 8);
  
  std::cout << "Testing undirected graphs..." << std::endl;
  type = "complete";
  SYNARMOSMA::Graph G6(3,type);
  assert(G6.connected());
  assert(G6.circuit_rank() == 1);
  assert(G6.order() == 3);
  assert(G6.size() == 3);
  assert(G6.stellar_addition(0) == 1);
  assert(G6.connected());
  assert(G6.circuit_rank() == 0);
  assert(G6.order() == 4);
  assert(G6.size() == 3);
  assert(G6.stellar_deletion(3) == 1);
  assert(G6.connected());
  assert(G6.circuit_rank() == 1);
  assert(G6.order() == 3);
  assert(G6.size() == 3);

  std::cout << "All tests passed successfully!" << std::endl;
  return 0;
}

