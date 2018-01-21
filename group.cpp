#include "group.h"

using namespace SYNARMOSMA;

extern Random RND;

Group::Group()
{
  clear();
}

Group::Group(const std::string& name,unsigned int n)
{
#ifdef DEBUG
  assert(n > 0);
#endif
  if (name == "Order") {
    // Here the user has specified the desired order of the 
    // group that is to be constructed
    Word w;
    int alpha;
    cardinality = n;
    finite = true;
    switch (n) {
      case 1:
        // Fairly simple, it's just the trivial group
        abelian = true;
        free = true;
        solvable = true;
        ngenerator = 0;
        break;
      case 2:
        // Also simple, the only group of order two is Z/2
        abelian = true;
        free = false;
        solvable = true;
        ngenerator = 1;
        w.initialize(ngenerator,0,2);
        relations.push_back(w);
        break;
      case 3:
        // The only group here is Z/3
        abelian = true;
        free = false;
        solvable = true;
        ngenerator = 1;
        w.initialize(ngenerator,0,3);
        relations.push_back(w);
        break;
      case 4:
        abelian = true;
        free = false;
        solvable = true;
        if (RND.irandom(2) == 0) {
          // Z/4
          ngenerator = 1;
          // {e,a,a^2,a^3}
          w.initialize(ngenerator,0,4);
          relations.push_back(w);
        }
        else {
          // Z/2 x Z/2, i.e. the Klein Vierergruppe
          ngenerator = 2;
          // {e,a,b,ab}
          Word w1(ngenerator,0,2);
          relations.push_back(w1);
          Word w2(ngenerator,1,2);
          relations.push_back(w2);
          relations.push_back(w1*w2);
        }
        break;
      case 5:
        // The only group here is Z/5
        abelian = true;
        free = false;
        solvable = true;
        ngenerator = 1;
        w.initialize(ngenerator,0,5);
        relations.push_back(w);
        break;
      case 6:
        free = false;
        solvable = true;
        if (RND.irandom(2) == 0) {
          // This is Z/6
          abelian = true;
          ngenerator = 1;
          w.initialize(ngenerator,0,6);
          relations.push_back(w);
        }
        else {
          // The dihedral group D_3
          abelian = false;
          ngenerator = 2;
          Word w1(ngenerator,0,2);
          Word w2(ngenerator,1,2);
          Word w3(ngenerator,0,3);
          Word w4(ngenerator,1,3);
          Word w5 = w3*w4;
          relations.push_back(w1);
          relations.push_back(w2);
          relations.push_back(w5);
        } 
        break;
      case 7:
        // The only group here is Z/7
        abelian = true;
        free = false;
        solvable = true;
        ngenerator = 1;
        w.initialize(ngenerator,0,7);
        relations.push_back(w);
        break;
      case 8:
        // There are five possibilities here: Z/2 x Z/2 x Z/2, Z/4 x Z/2, Z/8, D_4 and H
        free = false;
        solvable = true;
        alpha = RND.irandom(5);
        if (alpha == 0) {
          // Z/2 x Z/2 x Z/2
          abelian = true;
          ngenerator = 3;
          // {e,a,b,c,ab,ac,bc,abc}
          Word w1(ngenerator,0,1);
          Word w2(ngenerator,1,1);
          Word w3(ngenerator,2,1);
          relations.push_back(w1*w1);
          relations.push_back(w2*w2);
          relations.push_back(w3*w3);
          relations.push_back(w1*w2*w1.invert()*w2.invert());
          relations.push_back(w1*w3*w1.invert()*w3.invert());
          relations.push_back(w2*w3*w2.invert()*w3.invert());
        }
        else if (alpha == 1) {
          // Z/4 x Z/2
          abelian = true;
          ngenerator = 2;
          // {e,b,b^2,b^3,a,ab,ab^2,ab^3}
          Word w1(ngenerator,0,1);
          Word w2(ngenerator,1,1);
          relations.push_back(w1*w1*w1*w1);
          relations.push_back(w2*w2);
          relations.push_back(w1*w2*w1.invert()*w2.invert());
        }
        else if (alpha == 2) {
          // Z/8
          abelian = true;
          ngenerator = 1;
          // {e,a,a^2,a^3,a^4,a^5,a^6,a^7}
          w.initialize(ngenerator,0,8);
          relations.push_back(w);
        }
        else if (alpha == 3) {
          // D_4
          abelian = false;
          ngenerator = 2;
          Word w1(ngenerator,0,2);
          Word w2(ngenerator,1,2);
          Word w3(ngenerator,0,4);
          Word w4(ngenerator,1,4);
          Word w5 = w3*w4;
          relations.push_back(w1);
          relations.push_back(w2);
          relations.push_back(w5);
        }
        else {
          // H, the quaternion group
          abelian = false;
          ngenerator = 2;
          Word w1(ngenerator,0,1);
          Word w2(ngenerator,1,1);
          Word w3(ngenerator,1,-2);
          relations.push_back(w1*w1*w1*w1);
          relations.push_back(w2.invert()*w2.invert()*w1*w1);
          relations.push_back(w1*w2.invert()*w1*w2);
        }  
        break;
      default:
        // Create a random group with this many elements...
        create_random();
    }
  }
  else if (name == "Dihedral") {
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
    // B_1 = {e}, B_2 = (Z,+)
    unsigned int i,j;
    braid = true;
    ngenerator = n - 1;
    finite = (n > 1) ? false : true;
    abelian = (n > 2) ? false : true;
    free = (n > 2) ? false : true;
    if (n == 2) solvable = true;
    for(i=0; i<ngenerator; ++i) {
      Word w1(ngenerator,i,1),w2(ngenerator,i,-1);
      for(j=i+2; j<ngenerator; ++j) {
        Word w3(ngenerator,j,1),w4(ngenerator,j,-1);
        relations.push_back(w2*w4*w1*w3);
      }
    }
    // The cubic relations
    for(i=0; i<ngenerator-1; ++i) {
      Word w1(ngenerator,i,1),w2(ngenerator,i,-1);
      Word w3(ngenerator,i+1,1),w4(ngenerator,i+1,-1);
      relations.push_back(w4*w2*w4*w1*w3*w1);
    }
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
    create_random();
  }
}

Group::Group(unsigned int r,const std::vector<unsigned int>& torsion)
{
  initialize(r,torsion);
} 

Group::Group(unsigned int n)
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

Group::Group(unsigned int n,const std::vector<Word>& R)
{
  initialize(n,R);
}

Group::Group(unsigned int n,unsigned int m)
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
    for(it=relations[i].content.begin(); it!=relations[i].content.end(); ++it) {
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

  // Sanity check...
  for(i=0; i<relations.size(); ++i) {
    w = relations[i];
    for(j=0; j<w.content.size(); ++j) {
      n = w.content[j].first;
      if (n >= ngenerator) return false;
    }
  }
  return true;
}

bool Group::equivalent(const Word& w1,const Word& w2) const
{
  // Make sure these are in fact words in this group...
#ifdef DEBUG
  assert(w1.legal() && w2.legal());
#endif
  // We only know how to solve the word problem in the braid group...
  if (!braid) return false;
  Word nword(ngenerator,0),cword = w1*w2.invert();
  cword.free_reduce();
  if (cword.empty()) return true;
  if (cword.homogeneous()) return false;
  int e1,e2,stype,dtype;
  unsigned int i,j,n,m,b1,b2,spoint,base;
  bool hfound,candidate,permissible,ppow,npow;
  std::vector<unsigned int> handle;
  do {
    n = cword.length();
    hfound = false;
    for(i=0; i<n; ++i) {
      spoint = i;
      b1 = cword.content[i].first;
      e1 = cword.content[i].second;
      handle.clear();
      candidate = false;
      for(j=1+i; j<n; ++j) {
        b2 = cword.content[j].first;
        e2 = cword.content[j].second;
        if (b2 == b1 && e1*e2 < 0 && !handle.empty()) {
          candidate = true;
          break;
        }
        handle.push_back(j); 
      }
      // See if this is a genuine handle...
      if (!candidate) continue;
      permissible = true;
      m = handle.size();
      for(j=0; j<m; ++j) {
        b2 = cword.content[handle[j]].first;
        if (b2 == b1) permissible = false;
        if (b1 > 0) {
          if (b2 == (b1-1)) permissible = false;
        }
      }
      ppow = false; npow = false;
      for(j=0; j<m; ++j) {
        b2 = cword.content[handle[j]].first;
        e2 = cword.content[handle[j]].second;
        if (b2 == (b1+1)) {
          if (e2 > 0) {
            ppow = true;
          }
          else if (e2 < 0) {
            npow = true;
          }
        }
      }
      if (ppow && npow) permissible = false; 
      if (permissible) {
        hfound = true;
        break;
      }
    }
    if (!hfound) break;
    // Do the handle reduction...
    base = cword.content[spoint].first;
    stype = (cword.content[spoint].second < 0) ? -1 : 1;
    m = handle.size();
    nword.clear();
    for(i=0; i<spoint; ++i) {
      nword.content.push_back(cword.content[i]);
    }
    for(i=0; i<m; ++i) {
      b1 = cword.content[handle[i]].first;
      if (b1 == (1+base)) {
        dtype = (cword.content[handle[i]].second < 0) ? -1 : 1;
        nword.content.push_back(std::pair<unsigned int,int>(1+base,-stype));
        nword.content.push_back(std::pair<unsigned int,int>(base,dtype));
        nword.content.push_back(std::pair<unsigned int,int>(1+base,stype));
      }
      else if (b1 == base) {
        continue;
      }
      else {
        nword.content.push_back(cword.content[handle[i]]);
      }
    }
    for(i=spoint+m+2; i<n; ++i) {
      nword.content.push_back(cword.content[i]);
    }
    nword.free_reduce();
    cword = nword;
    if (cword.homogeneous()) break;
  } while(true);
  if (cword.empty()) return true;
  return false;
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
#ifdef DEBUG
    assert(torsion[i] > 1);
#endif
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
  braid = false;
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
      for(k=1; k<relations[j].length(); ++k) {
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
  braid = source.braid;
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
  braid = source.braid;
  cardinality = source.cardinality;
  rank = source.rank;
  torsion = source.torsion;

  return *this;
}

Group::~Group()
{

}

int Group::serialize(std::ofstream& s) const
{
  unsigned int i,n = relations.size();
  int count = 0;

  s.write((char*)(&ngenerator),sizeof(int)); count += sizeof(int);
  s.write((char*)(&cardinality),sizeof(int)); count += sizeof(int);
  s.write((char*)(&abelian),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&finite),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&solvable),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&free),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&rank),sizeof(int)); count += sizeof(int);
  n = torsion.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&torsion[i]),sizeof(int)); count += sizeof(int);
  }
  n = relations.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int); 
  for(i=0; i<n; ++i) {
    count += relations[i].serialize(s);
  }
  return count;
}

int Group::deserialize(std::ifstream& s)
{
  unsigned int i,j,n;
  int count = 0;
  Word w(ngenerator);

  clear();

  s.read((char*)(&ngenerator),sizeof(int)); count += sizeof(int);
  s.read((char*)(&cardinality),sizeof(int)); count += sizeof(int);
  s.read((char*)(&abelian),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&finite),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&solvable),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&free),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&rank),sizeof(int)); count += sizeof(int);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    torsion.push_back(j);
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += w.deserialize(s);
    relations.push_back(w);
  }
  return count;
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

namespace SYNARMOSMA {
  std::ostream& operator <<(std::ostream& s,const Group& g)
  {
    unsigned int i;
    s << g.finite << "  " << g.abelian << "  " << g.cardinality << "  " << g.solvable << "  " << g.free << std::endl;
    if (g.ngenerator == 0) {
      s << "< {e} | >";
      return s;
    }
    s << "< {";
    for(i=0; i<g.ngenerator-1; ++i) {
      s << "x[" << i+1 << "],";
    }
    s << "x[" << g.ngenerator << "]} | ";
    if (g.relations.empty()) {
      s << ">";
      return s;
    }
    s << "{";
    for(i=0; i<g.relations.size()-1; ++i) {
      s << g.relations[i] << ",";
    }
    s << g.relations[g.relations.size()-1] << "} >";
    return s;
  }
}

