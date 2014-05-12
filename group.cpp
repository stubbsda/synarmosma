#include "group.h"

extern Random RND;

Group::Group()
{
  clear();
}

Group::Group(const std::string& name,int n)
{
  assert(n > 2);
  if (name == "Dihedral") {
    // Dihedral group
    finite = true;
    abelian = (n > 2) ? false : true;
    free = false;
    solvable = true;
    cardinality = 2*n;
    ngenerator = 2;
    Word w1(ngenerator,0,2),w2(ngenerator,1,2),w3(ngenerator,0,n),w4(ngenerator,1,n);
    Word w5 = w3 * w4;
    relations.push_back(w1);
    relations.push_back(w2);
    relations.push_back(w5);
  }
  else if (name == "Braid") {
    // Braid group
    ngenerator = n - 1;
    finite = true;

  }
  else if (name == "Cyclic") {
    // Cyclic group
    finite = true;
    abelian = true;
    free = false;
    solvable = true;
    cardinality = n;
    ngenerator = 1;
    Word w(ngenerator,0,n);
    relations.push_back(w);
  }
  else if (name == "Symmetric") {
    // Symmetric group
    cardinality = factorial(n);
    finite = true;
    abelian = false;
    solvable = (n < 5) ? true : false;
  }
  else if (name == "Alternating") {
    // Alternating group
    finite = true;
    cardinality = factorial(n)/2;
    solvable = (n < 5) ? true : false;
  }
  else {
    std::cout << "Unknown group, creating a random group..." << std::endl;
    Group();
  }
}

Group::Group(unsigned int r,const std::vector<unsigned int>& torsion)
{
  initialize(r,torsion);
} 

Group::Group(int n)
{
  // The free group on n generators
  ngenerator = n;
  if (n == 0) {
    // The trivial group...
    cardinality = 1;
    finite = true;
    solvable = true;
    abelian = true;
    free = true;
  }
  else if (n == 1) {
    // This is just the infinite cyclic group...
    finite = false;
    abelian = true;
    free = true;
    solvable = true;
  }
} 

Group::Group(int n,const std::vector<Word>& R)
{
  initialize(n,R);
}

Group::Group(int n,int m)
{
  ngenerator = n;
  if (m > 0) {
    free = false;
    allocate(m);
  }
}

void Group::create_random()
{
  clear();
  // A random group...
  int m = 0;
  ngenerator = RND.irandom(3,15);
  double alpha = RND.drandom();
  if (alpha < 0.7) {
    m = RND.irandom(ngenerator/2);
    free = false;
  }
  else {
    free = true;
  }
  allocate(m);
}

unsigned int Group::implied_generators() const
{
  // This method computes the number of distinct 
  // generators in the various relations
  if (relations.empty()) return 0;
  unsigned int i,n;
  std::set<int> glabels;
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  for(i=0; i<relations.size(); ++i) {
    for(it=relations[i].content.begin(); it!=relations[i].content.end(); it++) {
      n = it->first;
      glabels.insert(n);
    }
  }
  n = glabels.size();
  return n;
}

bool Group::consistent() const
{
  unsigned int i,j,n;
  Word w(ngenerator,0);
  bool output = true;

  // Sanity check...
  for(i=0; i<relations.size(); ++i) {
    w = relations[i];
    for(j=0; j<w.content.size(); ++j) {
      n = w.content[j].first;
      if (n >= ngenerator) {
        std::cout << "Error in relation " << i << ": " << n << "  " << ngenerator << std::endl;
        output = false;
      }
    }
  }
  return output;
}

void Group::reduce()
{
  if (relations.empty()) return;

  unsigned int i,j,p,q,offset[ngenerator],n = 0;
  bool inv,reduction = false,change = false;
  std::set<unsigned int> trivial_generators;
  std::vector<Word> new_relations;
  Word w(ngenerator);

  for(i=0; i<relations.size(); ++i) {
    w = relations[i].normalize();
    if (w != relations[i]) change = true;
    if (w != relations[i]) change = true;
    if (w.empty()) continue;
    if (w.trivial()) {
      trivial_generators.insert(w.content[0].first);
    }
    else {
      new_relations.push_back(w);
    }
  }
  if (!trivial_generators.empty()) change = true;
  relations.clear();
  
  for(i=0; i<ngenerator; ++i) {
    if (trivial_generators.count(i) == 1) continue;
    offset[i] = n;
    n++;  
  }
  ngenerator -= trivial_generators.size();

  for(i=0; i<new_relations.size(); ++i) {
    w = new_relations[i].reduce(ngenerator,trivial_generators,offset);
    relations.push_back(w.normalize());
  }
  new_relations.clear();

  // The next step is to seek out quasi-trivial relations, 
  // such as x[i]^1*x[j]^(-1) = 1, which imply that x[i] = 
  // x[j]
  for(i=0; i<relations.size(); ++i) {
    w = relations[i];
    if (w.alias()) {
      p = w.content[0].first;
      q = w.content[1].first;
      inv = ((w.content[0].second + w.content[1].second) == 0) ? false : true;
      if (p < q) {
        // Convert x[q] into x[p]
        for(j=0; j<relations.size(); ++j) {
          if (j == i) continue;
          w = relations[j].swap(p,q,inv);
          new_relations.push_back(w.normalize());
        }
      }
      else {
        // Convert x[p] into x[q]
        for(j=0; j<relations.size(); ++j) {
          if (j == i) continue;
          w = relations[j].swap(q,p,inv);
          new_relations.push_back(w.normalize());
        }
      }
      ngenerator -= 1;
      change = true;
      reduction = true;
      break;
    }
  }
  if (reduction) relations = new_relations;

  if (ngenerator == 1 && change == false) {
    n = relations.size();
    p = 0;
    for(i=0; i<n; ++i) {
      if (relations[i].content[0].second < 0) p++;
    }
    if (p > n/2) {
      for(i=0; i<n; ++i) {
        relations[i].content[0].second *= -1;
      }
    }
  }
  if (change) {
    if (ngenerator == 0) {
      // The trivial group...
      cardinality = 1;
      finite = true;
      solvable = true;
      abelian = true;
      free = true;
      relations.clear();
      return;
    }
    reduce();
  }
  else {
    // If the presentation is simple enough, such as < {x} | {x^n} > then 
    // we can identify some properties of the group it represents, such as 
    // (in this case) it's cyclic and thus abelian, finite and of order n. 
    if (ngenerator == 1) {
       if (relations.empty()) {
         abelian = true;
         solvable = true;
         finite = false;
         rank = 1;
       }
       else if (relations.size() == 1) {
         abelian = true;
         solvable = true;
         finite = true;
         rank = 0;
         cardinality = relations[0].content[0].second;
       }
     }
  }
}

void Group::initialize(unsigned int r,const std::vector<unsigned int>& torsion)
{
  // The constructor for a finitely generated abelian group, whose presentation 
  // is thus of the form
  // \begin{equation*}
  // G = <g_1,\dotsc,g_{r+s};g_i g_j g_i^{-1}g_j^{-1}, 1 \le i < j \le r, g_{r+k}^{t_k}, 
  // 1 \le k \le s>
  // \end{equation*}
  // corresponding to
  // \begin{equation*} 
  // G = \mathbb{Z}^r \oplus \mathbb{Z}_{t_1} \oplus \dotsb \oplus \mathbb{Z}_{t_s}
  // \end{equation*}
  clear();

  unsigned int i,j;

  rank = r;
  abelian = true;
  solvable = true;
  free = false;
  finite = (rank > 0) ? false : true;
  ngenerator = rank + torsion.size();
  Word w(ngenerator);
  w.clear();
  for(i=0; i<rank; ++i) {
    for(j=i+1; j<rank; ++j) {
      w.content.push_back(std::pair<int,int>(i,1));
      w.content.push_back(std::pair<int,int>(j,1));
      w.content.push_back(std::pair<int,int>(i,-1));
      w.content.push_back(std::pair<int,int>(j,-1));
      relations.push_back(w);
      w.clear();
    }
  }
  if (rank == 0) cardinality = 0;
  for(i=0; i<torsion.size(); ++i) {
    assert(torsion[i] > 1);
    w.content.push_back(std::pair<int,int>(rank+i,torsion[i]));
    relations.push_back(w);
    w.clear();
    if (rank == 0) cardinality += torsion[i];
  }
}

void Group::initialize(unsigned int n,const std::vector<Word>& R)
{
  clear();

  ngenerator = n;
  relations = R;
  reduce();
  if (ngenerator == 0) {
    abelian = true;
    cardinality = 1;
    finite = true;
    solvable = true;
    free = true;
  }
}

void Group::compute_rank()
{
  if (!abelian) {
    std::cerr << "Error: Rank is undefined for a non-abelian group!" << std::endl;
    return;
  }
  if (finite) {
    rank = 0;
    return;
  }
  rank = ngenerator;
  unsigned int i;
  for(i=0; i<relations.size(); ++i) {
   
  }
}

void Group::clear()
{
  ngenerator = 0;
  cardinality = 0;
  relations.clear();
  abelian = false;
  finite = false;
  solvable = false;
  free = false;
  rank = 0;
  torsion.clear();
}

Group Group::abelianize() const
{
  // We need to add the word $aba^{-1}b^{-1}$ to the set of relations for each pair of generators
  // $a$ and $b$ in the group presentation.
  unsigned int i,j,k;
  bool found;
  Word w1(ngenerator),w2(ngenerator);
  Group output(*this);

  for(i=0; i<ngenerator; ++i) {
    for(j=1+i; j<ngenerator; ++j) {
      // Does this pair of generators already exist in the set of relations?
      w1.clear();

      // w1 = $aba^{-1}b^{-1}$
      w1.content.push_back(std::pair<unsigned int,int>(i,1));
      w1.content.push_back(std::pair<unsigned int,int>(j,1));
      w1.content.push_back(std::pair<unsigned int,int>(i,-1));
      w1.content.push_back(std::pair<unsigned int,int>(j,-1));

      // w2 = $bab^{-1}a^{-1}$
      w2 = !w1;

      found = false;
      for(k=0; k<relations.size(); ++k) {
        if (relations[k] == w1 || relations[k] == w2) {
          found = true;
          break;
        }
      }
      if (found) continue;
      output.relations.push_back(w1);
    }
  }
  return output;
}

void Group::allocate(unsigned int m)
{
  abelian = false;
  finite = false;
  solvable = false;
  cardinality = 0;

  if (m == 0) return;

  int e;
  unsigned int i,j,k,b,sum = 0;
  bool good;
  Word w(ngenerator);
  std::vector<unsigned int> length,base;
  std::vector<int> exponent,used;

  for(i=0; i<m; ++i) {
   relations.push_back(w);
  }
  for(i=0; i<m-1; ++i) {
    length.push_back(2 + RND.irandom(ngenerator/2));
    sum += length[i];
  }
  if (ngenerator > sum) {
    length.push_back(ngenerator - sum);
  }
  else {
    length.push_back(2 + RND.irandom(ngenerator/2));
  }
  for(i=0; i<ngenerator; ++i) {
    used.push_back(0);
  }

  i = 0;
  do {
    base.clear(); exponent.clear();
    for(j=0; j<length[i]; ++j) {
      e = RND.irandom(1,10);
      if (RND.drandom() < 0.5) e = -e;
      exponent.push_back(e);
      // How to ensure that every generator appears at least 
      // once in a relation?
      if (j > 0) {
        good = false;
        do {
          b = RND.irandom(used);
          if (b != base[j-1]) good = true;
        } while(!good);
        used[b] = 1;
        base.push_back(b);
        continue;
      }
      k = RND.irandom(used);
      used[k] = 1;
      base.push_back(k);
    }
    relations[i].initialize(base,exponent);
    // Check to make sure none of these relators are simply cyclic permutations 
    // of one of the others!
    good = true;
    for(j=0; j<i; ++j) {
      for(k=1; k<relations[j].size(); ++k) {
        relations[j].permute(k,w);
        if (w == relations[i]) {
          good = false;
          break;
        }
      }
      if (!good) break;
    }
    if (good) i += 1;
  } while(i < m);
}

Group::Group(const Group& source)
{
  relations = source.relations;
  ngenerator = source.ngenerator;
  abelian = source.abelian;
  finite = source.finite;
  solvable = source.solvable;
  free = source.free;
  cardinality = source.cardinality;
  rank = source.rank;
  torsion = source.torsion;
}

Group& Group::operator =(const Group& source)
{
  if (this == &source) return *this;

  relations = source.relations;
  ngenerator = source.ngenerator;
  abelian = source.abelian;
  finite = source.finite;
  solvable = source.solvable;
  free = source.free;
  cardinality = source.cardinality;
  rank = source.rank;
  torsion = source.torsion;

  return *this;
}

Group::~Group()
{

}

void Group::serialize(std::ofstream& s) const
{
  unsigned int i,n = relations.size();
  s.write((char*)(&ngenerator),sizeof(int));
  s.write((char*)(&cardinality),sizeof(int));
  s.write((char*)(&abelian),sizeof(bool));
  s.write((char*)(&finite),sizeof(bool));
  s.write((char*)(&solvable),sizeof(bool));
  s.write((char*)(&free),sizeof(bool));
  s.write((char*)(&rank),sizeof(int));
  n = torsion.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.write((char*)(&torsion[i]),sizeof(int));
  }
  n = relations.size();
  s.write((char*)(&n),sizeof(int));  
  for(i=0; i<n; ++i) {
    relations[i].serialize(s);
  }
}

void Group::deserialize(std::ifstream& s)
{
  unsigned int i,j,n;
  Word w(ngenerator);

  clear();

  s.read((char*)(&ngenerator),sizeof(int));
  s.read((char*)(&cardinality),sizeof(int));
  s.read((char*)(&abelian),sizeof(bool));
  s.read((char*)(&finite),sizeof(bool));
  s.read((char*)(&solvable),sizeof(bool));
  s.read((char*)(&free),sizeof(bool));
  s.read((char*)(&rank),sizeof(int));
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    torsion.push_back(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    w.deserialize(s);
    relations.push_back(w);
  }
}

std::string Group::compact_form() const
{
  unsigned int i,m;
  std::string output = "";

  if (abelian) {
    output += boost::lexical_cast<std::string>(rank);
    if (!torsion.empty()) {
      m = torsion.size();
      output += " | ";
      for(i=0; i<m-1; ++i) {
        output += boost::lexical_cast<std::string>(torsion[i]) + ", ";
      }
      output += boost::lexical_cast<std::string>(torsion[m-1]); 
    }
  }
  else {
    if (ngenerator == 0) {
      output += "{e}";
    }
    else {
      output += "{";
      for(i=0; i<ngenerator-1; ++i) {
        output += "x[" + boost::lexical_cast<std::string>(i+1) + "],";
      }
      output += "x[" + boost::lexical_cast<std::string>(ngenerator) + "]}";
      if (!relations.empty()) {
        m = relations.size();
        std::stringstream sstream;
        sstream << " | {";
        for(i=0; i<m-1; ++i) {
          sstream << relations[i] << ",";
        }
        sstream << relations[m-1] << "}";
        output += sstream.str();
      }
    }
  }
  return output;
}

std::ostream& operator <<(std::ostream& os,const Group& g)
{
  unsigned int i;
  os << g.finite << "  " << g.abelian << "  " << g.cardinality << "  " << g.solvable << "  " << g.free << std::endl;
  if (g.ngenerator == 0) {
    os << "< {e} | >";
    return os;
  }
  os << "< {";
  for(i=0; i<g.ngenerator-1; ++i) {
    os << "x[" << i+1 << "],";
  }
  os << "x[" << g.ngenerator << "]} | ";
  if (g.relations.empty()) {
    os << ">";
    return os;
  }
  os << "{";
  for(i=0; i<g.relations.size()-1; ++i) {
    std::cout << g.relations[i] << ",";
  }
  os << g.relations[g.relations.size()-1] << "} >";
  return os;
}

