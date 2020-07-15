#include "group.h"

using namespace SYNARMOSMA;

extern Random RND;

Group::Group()
{

}

Group::Group(const std::string& type,unsigned int n)
{
  std::string utype = boost::to_upper_copy(type);
  if (utype == "ORDER") {
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
        w.initialize(0,2);
        relations.push_back(w);
        break;
      case 3:
        // The only group here is Z/3
        abelian = true;
        free = false;
        solvable = true;
        ngenerator = 1;
        w.initialize(0,3);
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
          w.initialize(0,4);
          relations.push_back(w);
        }
        else {
          // Z/2 x Z/2, i.e. the Klein Vierergruppe
          ngenerator = 2;
          // {e,a,b,ab}
          Word w1(0,2);
          relations.push_back(w1);
          Word w2(1,2);
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
        w.initialize(0,5);
        relations.push_back(w);
        break;
      case 6:
        free = false;
        solvable = true;
        if (RND.irandom(2) == 0) {
          // This is Z/6
          abelian = true;
          ngenerator = 1;
          w.initialize(0,6);
          relations.push_back(w);
        }
        else {
          // The dihedral group D_3
          abelian = false;
          ngenerator = 2;
          Word w1(0,2);
          Word w2(1,2);
          Word w3(0,3);
          Word w4(1,3);
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
        w.initialize(0,7);
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
          Word w1(0,1);
          Word w2(1,1);
          Word w3(2,1);
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
          Word w1(0,1);
          Word w2(1,1);
          relations.push_back(w1*w1*w1*w1);
          relations.push_back(w2*w2);
          relations.push_back(w1*w2*w1.invert()*w2.invert());
        }
        else if (alpha == 2) {
          // Z/8
          abelian = true;
          ngenerator = 1;
          // {e,a,a^2,a^3,a^4,a^5,a^6,a^7}
          w.initialize(0,8);
          relations.push_back(w);
        }
        else if (alpha == 3) {
          // D_4
          abelian = false;
          ngenerator = 2;
          Word w1(0,2);
          Word w2(1,2);
          Word w3(0,4);
          Word w4(1,4);
          Word w5 = w3*w4;
          relations.push_back(w1);
          relations.push_back(w2);
          relations.push_back(w5);
        }
        else {
          // H, the quaternion group
          abelian = false;
          ngenerator = 2;
          Word w1(0,1);
          Word w2(1,1);
          Word w3(1,-2);
          relations.push_back(w1*w1*w1*w1);
          relations.push_back(w2.invert()*w2.invert()*w1*w1);
          relations.push_back(w1*w2.invert()*w1*w2);
        }  
        break;
      default:
        throw std::invalid_argument("Group order > 8 is too high for constructor!");
    }
  }
  else if (utype == "DIHEDRAL") {
    // Dihedral group
    finite = true;
    abelian = (n > 2) ? false : true;
    free = false;
    solvable = true;
    cardinality = 2*n;
    ngenerator = 2;
    Word w1(0,2),w2(1,2),w3(0,n),w4(1,n);
    Word w5 = w3*w4;
    relations.push_back(w1);
    relations.push_back(w2);
    relations.push_back(w5);
  }
  else if (utype == "BRAID") {
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
      Word w1(i,1),w2(i,-1);
      for(j=i+2; j<ngenerator; ++j) {
        Word w3(j,1),w4(j,-1);
        relations.push_back(w2*w4*w1*w3);
      }
    }
    // The cubic relations
    for(i=0; i<ngenerator-1; ++i) {
      Word w1(i,1),w2(i,-1);
      Word w3(i+1,1),w4(i+1,-1);
      relations.push_back(w4*w2*w4*w1*w3*w1);
    }
  }
  else if (utype == "CYCLIC") {
    // Cyclic group
    finite = true;
    abelian = true;
    free = false;
    solvable = true;
    cardinality = n;
    ngenerator = 1;
    Word w(0,n);
    relations.push_back(w);
  }
  else if (utype == "SYMMETRIC") {
    // Symmetric group
    cardinality = factorial(n);
    finite = true;
    abelian = false;
    solvable = (n < 5) ? true : false;
  }
  else if (utype == "ALTERNATING") {
    // Alternating group
    finite = true;
    cardinality = factorial(n)/2;
    solvable = (n < 5) ? true : false;
  }
  else {
    throw std::invalid_argument("Unrecognized group type!");
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
  if (n == 0) {
    if (!R.empty()) throw std::invalid_argument("If there are no generators, then there should be no relations either!");
  }
  else {
    unsigned int i,m;
    std::set<unsigned int> S;

    free = false;
    // For each relation, get the maximum index of its letters and verify 
    // that this is less than value of ngenerator.
    for(i=0; i<R.size(); ++i) {
      R[i].get_alphabet(S);
      m = *S.rbegin();
      if (m >= n) throw std::invalid_argument("One of the relations in the Group constructor has an illegal generator!");
    }
  }
  initialize(n,R);
}

Group::Group(unsigned int n,unsigned int m)
{
  ngenerator = n;
  if (m > 0) {
    free = false;
    initialize(m);
  }
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

bool Group::consistent(std::set<unsigned int>& alphabet) const
{
  // This method computes the number of distinct 
  // generators in the various relations and verifies 
  // that there are no letters in the relation words 
  // that aren't generators.
  alphabet.clear();
  if (relations.empty()) return true;
  unsigned int i,n;
  bool output = true;
  std::set<unsigned int> S;
  std::set<unsigned int>::const_iterator it;

  for(i=0; i<relations.size(); ++i) {
    relations[i].get_alphabet(S);
    for(it=S.begin(); it!=S.end(); ++it) {
      n = *it;
      alphabet.insert(n);
      // Illegal letter in this relation!
      if (n >= ngenerator) output = false;
    }
  }
  return output;
}

bool Group::equivalent(const Word& w1,const Word& w2) const
{
  // Make sure these are in fact words in this group...
  if (!w1.legal() || !w2.legal()) throw std::invalid_argument("One or more of the words to be tested for equivalency is illegal!");

  // We only know how to solve the word problem in the braid group...
  if (!braid) return false;
  Word cword,nword = w1*w2.invert();
  cword = nword.reduce();
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
      nword.append(cword.content[i]);
    }
    for(i=0; i<m; ++i) {
      b1 = cword.content[handle[i]].first;
      if (b1 == (1+base)) {
        dtype = (cword.content[handle[i]].second < 0) ? -1 : 1;
        nword.append(std::pair<unsigned int,int>(1+base,-stype));
        nword.append(std::pair<unsigned int,int>(base,dtype));
        nword.append(std::pair<unsigned int,int>(1+base,stype));
      }
      else if (b1 == base) {
        continue;
      }
      else {
        nword.append(cword.content[handle[i]]);
      }
    }
    for(i=spoint+m+2; i<n; ++i) {
      nword.append(cword.content[i]);
    }
    cword = nword.reduce();
    if (cword.homogeneous()) break;
  } while(true);
  if (cword.empty()) return true;
  return false;
}

void Group::reduce()
{
  if (relations.empty()) return;
  if (ngenerator == 0) throw std::runtime_error("The number of generators must be greater than zero to reduce the group!");

  unsigned int i,j,p,q,offset[ngenerator],n = 0;
  bool inv,reduction = false,change = false;
  std::set<unsigned int> trivial_generators;
  std::vector<Word> new_relations;
  Word w;

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

  if (trivial_generators.empty()) {
    // Nothing to do in this case...
    relations = new_relations;
  }
  else {
    for(i=0; i<new_relations.size(); ++i) {
      w = new_relations[i].reduce(trivial_generators,offset);
      relations.push_back(w.normalize());
    }
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

void Group::initialize(unsigned int m)
{
  abelian = false;
  finite = false;
  solvable = false;
  cardinality = 0;

  if (m == 0) return;

  int e;
  unsigned int i,j,k,b,sum = 0;
  bool duplicate;
  Word w;
  std::vector<unsigned int> length,base;
  std::vector<int> exponent;

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

  i = 0;
  do {
    base.clear(); exponent.clear();
    j = RND.irandom(ngenerator);
    base.push_back(j);
    e = RND.irandom(1,10);
    if (RND.drandom() < 0.5) e = -e;
    exponent.push_back(e);
    for(j=1; j<length[i]; ++j) {
      e = RND.irandom(1,10);
      if (RND.drandom() < 0.5) e = -e;
      exponent.push_back(e);
      do {
        b = RND.irandom(ngenerator);
        if (b != base[j-1]) break;
      } while(true);
      base.push_back(b);
    }
    relations[i].initialize(base,exponent);
    // Check to make sure none of these relators are simply cyclic permutations 
    // of one of the others!
    duplicate = false;
    for(j=0; j<i; ++j) {
      for(k=1; k<relations[j].length(); ++k) {
        w = relations[j].permute(k);
        if (w == relations[i]) {
          duplicate = true;
          break;
        }
      }
      if (duplicate) break;
    }
    if (!duplicate) i++;
  } while(i < m);
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
  unsigned int i,j;

  rank = r;
  abelian = true;
  solvable = true;
  free = false;
  finite = (rank > 0) ? false : true;
  ngenerator = rank + torsion.size();
  Word w;
  
  for(i=0; i<rank; ++i) {
    for(j=i+1; j<rank; ++j) {
      w.append(std::pair<unsigned int,int>(i,1));
      w.append(std::pair<unsigned int,int>(j,1));
      w.append(std::pair<unsigned int,int>(i,-1));
      w.append(std::pair<unsigned int,int>(j,-1));
      relations.push_back(w);
      w.content.clear();
    }
  }
  if (rank == 0) cardinality = 0;
  for(i=0; i<torsion.size(); ++i) {
    if (torsion[i] <= 1) throw std::invalid_argument("The torsion elements for the initialization of an Abelian group must be greater than one!");
    w.append(std::pair<unsigned int,int>(rank+i,torsion[i]));
    relations.push_back(w);
    w.clear();
    if (rank == 0) cardinality += torsion[i];
  }
}

void Group::initialize(unsigned int n,const std::vector<Word>& R)
{
  if (n == 0 && !R.empty()) throw std::invalid_argument("If there are no generators then the group presentation cannot have any relations!");

  ngenerator = n;
  relations = R;
 
  if (ngenerator == 0) {
    abelian = true;
    cardinality = 1;
    finite = true;
    solvable = true;
    free = true;
  }
  else {
    reduce();
  }
}

void Group::compute_rank()
{
  if (!abelian) throw std::runtime_error("Rank is undefined for a non-abelian group!");
  if (finite) {
    rank = 0;
    return;
  }
  // The maximum possible value for the rank; this method needs to be completed by an 
  // algorithm which reduces the relations to their canonical form and enables the calculation 
  // of the rank and torsion coefficients.
  rank = ngenerator;
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
  Word w1,w2;
  Group output(*this);

  for(i=0; i<ngenerator; ++i) {
    for(j=1+i; j<ngenerator; ++j) {
      // Does this pair of generators already exist in the set of relations?
      w1.clear();

      // w1 = $aba^{-1}b^{-1}$
      w1.append(std::pair<unsigned int,int>(i,1));
      w1.append(std::pair<unsigned int,int>(j,1));
      w1.append(std::pair<unsigned int,int>(i,-1));
      w1.append(std::pair<unsigned int,int>(j,-1));

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
  Word w;

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
  if (ngenerator == 0) return std::string("{e}");

  unsigned int i,m;
  std::stringstream sstream;

  if (abelian) {
    sstream << rank;
    if (!torsion.empty()) {
      m = torsion.size();
      sstream << " | ";
      for(i=0; i<m-1; ++i) {
        sstream << torsion[i] << ", ";
      }
      sstream << torsion[m-1];
    }
  }
  else {
    sstream << "{";
    for(i=0; i<ngenerator-1; ++i) {
      sstream << "x[" << i+1 << "],";
    }
    sstream << "x[" << ngenerator << "]}";
    if (!relations.empty()) {
      m = relations.size();
      sstream << " | {";
      for(i=0; i<m-1; ++i) {
        if (relations[i].empty()) continue;
        sstream << relations[i] << ",";
      }
      sstream << relations[m-1] << "}";
    }
  }
  return sstream.str();
}

namespace SYNARMOSMA {
  std::ostream& operator <<(std::ostream& s,const Group& g)
  {
    unsigned int i;
    s << g.finite << "  " << g.abelian << "  " << g.cardinality << "  " << g.solvable << "  " << g.free << "  " << g.braid << std::endl;
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

