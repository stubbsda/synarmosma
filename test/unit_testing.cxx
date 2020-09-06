#include <synarmosma/lattice.h>
#include <synarmosma/polynomial.h>
#include <synarmosma/directed_graph.h>
#include <synarmosma/group.h>
#include <synarmosma/integer_matrix.h>

SYNARMOSMA::Random RND;

int main(int argc,char** argv)
{
  RND.set_seed(12);
  std::string test = std::string(argv[1]);

  if (test == "Rational") {
    SYNARMOSMA::Rational q1(1,5),q2(2,3),r(13,15),s(27,41);
    SYNARMOSMA::Rational q3 = q1 + q2;
    if (q3 != r) return 1;
    std::vector<SYNARMOSMA::Rational> qx;
    qx.push_back(SYNARMOSMA::Rational(1,2));
    qx.push_back(SYNARMOSMA::Rational(9,8));
    qx.push_back(SYNARMOSMA::Rational(3,5));
    if (SYNARMOSMA::compute_mean(qx,"Harmonic") != s) return 1;
  }
  else if (test == "Polynomial") {
    SYNARMOSMA::Polynomial<SYNARMOSMA::Rational> p(2),q(3);
    SYNARMOSMA::Polynomial<SYNARMOSMA::Rational> r = p + q;
    if (r.get_degree() != 3) return 1;
    r = p*q;
    if (r.get_degree() != 5) return 1;
    std::ofstream s1("polynomial.dat",std::ios::out | std::ios::binary);
    q.serialize(s1);
    s1.close();
    std::ifstream s2("polynomial.dat",std::ios::in | std::ios::binary);
    p.deserialize(s2);
    s2.close();
    std::system("rm -f polynomial.dat");
    if (p != q) return 1;
  }
  else if (test == "Poset") {
    SYNARMOSMA::Poset P;
    P.power_set(7);
    if (!P.consistent()) return 1;
  }
  else if (test == "Lattice") { 
    SYNARMOSMA::Lattice L(4);
    if (!L.consistent()) return 1;
  }
  else if (test == "Graph") {
    int i,j,n = 25;
    std::vector<double> vx;
    SYNARMOSMA::Graph G("Petersen"),H(n,"CYCLIC"),K(3,"Complete");
    std::vector<SYNARMOSMA::Monomial<int> > output;

    G.entwinement();
    G.vertex_centrality(vx,0.000001);
    G.tutte_polynomial(output);
    for(i=0; i<(signed) output.size(); ++i) {
      if (output[i].exponents.size() != 2) continue;
      if (output[i].exponents[0].second == 1 && output[i].exponents[1].second == 3) if (output[i].coefficient != 65) return 1;
      if (output[i].exponents[0].second == 2 && output[i].exponents[1].second == 2) if (output[i].coefficient != 105) return 1;
    }

    H.clustering_coefficient();
    H.mean_path_length();
    for(i=0; i<n; ++i) {
      if (RND.drandom() < 0.15) {
        H.drop_edge(i,(i+1)%n);
        do {
          j = RND.irandom(n);
          if (i == j) continue;
          if (H.add_edge(i,j)) break;
        } while(true);
      }
    }
    H.clustering_coefficient();
    H.mean_path_length();

    if (!K.connected()) return 1;
    if (K.circuit_rank() != 1) return 1;
    if (K.order() != 3) return 1;
    if (K.size() != 3) return 1;
    if (K.stellar_addition(0) != 1) return 1;
    if (!K.connected()) return 1;
    if (K.circuit_rank() != 0) return 1;
    if (K.order() != 4) return 1;
    if (K.size() != 3) return 1;
    if (K.stellar_deletion(3) != 1) return 1;
    if (!K.connected()) return 1;
    if (K.circuit_rank() != 1) return 1;
    if (K.order() != 3) return 1;
    if (K.size() != 3) return 1;
  }
  else if (test == "Pseudograph") {
    std::vector<int> vx;
    SYNARMOSMA::Pseudograph G(4),H(4);

    G.add_edge(0,1); G.add_edge(0,2); G.add_edge(0,3);
    G.add_edge(1,2);
    G.add_edge(2,3);

    H.add_edge(0,1); H.add_edge(0,2); H.add_edge(0,3);
    H.add_edge(1,2);
    H.add_edge(2,3);
    G.contract(0,1,&H);
    
    if (H.get_candidates(vx) != 0) return 1;
    if (vx.size() != 8) return 1;
    if (vx[0] != 0) return 1;
    if (vx[1] != 1) return 1;
    if (vx[2] != 0) return 1;
    if (vx[3] != 2) return 1;
    if (vx[4] != 0) return 1;
    if (vx[5] != 1) return 1;
    if (vx[6] != 1) return 1;
    if (vx[7] != 2) return 1;
  }
  else if (test == "Directed_Graph") {
    SYNARMOSMA::Directed_Graph G,H;
    std::pair<int,int> pr;
    SYNARMOSMA::pair_index dmap;
    SYNARMOSMA::pair_index::const_iterator qt;
    auto d1 = SYNARMOSMA::Relation::before,d2 = SYNARMOSMA::Relation::after,d3 = SYNARMOSMA::Relation::disparate;

    for(int i=0; i<7; ++i) {
      G.add_vertex();
    }
    G.add_edge(0,1,d1,10.0);
    G.add_edge(0,2,d1,5.0);
    G.add_edge(0,3,d1,4.0);
    G.add_edge(1,2,d1,3.0);
    G.add_edge(1,4,d1,5.0);
    G.add_edge(2,3,d1,5.0);
    G.add_edge(2,4,d1,2.0);
    G.add_edge(3,5,d1,8.0);
    G.add_edge(4,6,d1,7.0);
    G.add_edge(5,2,d1,3.0);
    G.add_edge(5,6,d1,11.0);
    if (G.bipartite() != 1) return 1;
    if (G.maximum_parents() != 3) return 1;
    if (G.DAG()) return 1;
    if (G.singly_connected()) return 1;
    if (G.eccentricity(3) != 2) return 1;
    G.compute_distances(dmap);
    pr.first = 0; pr.second = 6;
    qt = dmap.find(pr);
    if (qt->second != 3) return 1;
    pr.first = 3; pr.second = 1;
    qt = dmap.find(pr);
    if (qt->second != -1) return 1;
    pr.first = 5; pr.second = 6;
    qt = dmap.find(pr);
    if (qt->second != 1) return 1;
    if (G.compute_flow(0,6) != 15) return 1;

    for(int i=0; i<6; ++i) {
      H.add_vertex();
    }
    H.add_edge(0,1,d1);
    H.add_edge(1,2,d1);
    H.add_edge(1,3,d3);
    H.add_edge(3,4,d3);
    H.add_edge(1,4,d3);
    H.add_edge(2,3,d2);
    H.add_edge(3,5,d3);
    H.add_edge(2,5,d3);
    H.compute_directedness();
    if (H.distance(0,4) != -1) return 1;
    if (H.distance(0,2) != 2) return 1;
    if (H.distance(2,3) != -1) return 1;
    if (H.distance(3,2) != 1) return 1;
    if (H.directedness() != 3) return 1;
    if (H.size() != 8) return 1;
  }
  else if (test == "Group") {
    SYNARMOSMA::Group G("DIHEDRAL",5);
    std::set<unsigned int> S;

    if (!G.consistent(S)) return 1;
    if (G.order() != 10) return 1;
    if (!G.get_solvability()) return 1;
  }
  else if (test == "Integer_Matrix") {
    SYNARMOSMA::Integer_Matrix<int> M(4,true);
    M.set(0,2,-1);
    M.set(1,0,4);
    M.set(2,1,-3);
    M.set(3,0,1);
    M.set(3,1,-1);

    if (M.determinant() != 13) return 1;
    if (M.symmetric()) return 1;
  }
  else {
    std::cerr << "Unrecognized test type!" << std::endl;
    return 1;
  }

  return 0;
}

