#include "cell.h"

Cell::Cell()
{

}

Cell::Cell(const std::string& vx)
{
  // Given a key, create the corresponding cell...
  std::string str;
  boost::char_separator<char> sep(":");
  boost::tokenizer<boost::char_separator<char> > tok(vx,sep);
  for(boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); beg++) {
    str = *beg;
    vertices.insert(boost::lexical_cast<int>(str));
  }
  string_assembly();
}

Cell::Cell(int n)
{
  for(int i=0; i<=n; ++i) {
    vertices.insert(i);
  }
  string_assembly();
}

Cell::Cell(int v1,int v2)
{
  // A specialized constructor for 0-cells and 1-cells
  std::stringstream s;
  if (v1 == -1) {
    vertices.insert(v2);
    s << v2;
  }
  else {
    std::stringstream s1,s2;

    vertices.insert(v1);
    vertices.insert(v2);
    if (v1 < v2) {
      s << v1 << ":" << v2;
    }
    else {
      s << v2 << ":" << v1;
    }
    s1 << v1;
    s2 << v2;
    faces.push_back(s1.str());
    faces.push_back(s2.str());
  }
  key = s.str();
}

Cell::Cell(const std::set<int>& v)
{
  vertices = v;
  string_assembly();
}

Cell::Cell(const Cell& source)
{
  vertices = source.vertices;
  key = source.key;
  entourage = source.entourage;
  faces = source.faces;
}

Cell& Cell::operator =(const Cell& source)
{
  if (this == &source) return *this;

  vertices = source.vertices;
  key = source.key;
  entourage = source.entourage;
  faces = source.faces;

  return *this;
}

void Cell::initialize(int v1,int v2)
{
  std::stringstream s;

  clear();

  if (v1 == -1) {
    vertices.insert(v2);
    s << v2;
  }
  else {
    std::stringstream s1,s2;

    vertices.insert(v1);
    vertices.insert(v2);
    if (v1 < v2) {
      s << v1 << ":" << v2;
    }
    else {
      s << v2 << ":" << v1;
    }
    s1 << v1;
    s2 << v2;
    faces.push_back(s1.str());
    faces.push_back(s2.str());
  }
  key = s.str();
}

Cell::~Cell()
{

}

void Cell::clear()
{
  vertices.clear();
  entourage.clear();
  faces.clear();
  key = "";
}

void Cell::string_assembly()
{
  int j,n = (signed) vertices.size(),i = 0;
  std::set<int>::const_iterator it;
  std::stringstream s;
  std::stringstream* sf = new std::stringstream[n];

  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    if (i < (n-1)) {
      s << *it << ":";
    }
    else {
      s << *it;
    }
    ++i;
  }
  key = s.str();

  faces.clear();
  for(i=0; i<n; ++i) {
    j = -1;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      ++j;
      if (j == i) continue;
      if (j == (n-1) || (i == (n-1) && j == (n-2))) {
        sf[i] << *it;
      }
      else {
        sf[i] << *it << ":";
      }
    }
    faces.push_back(sf[i].str());
  }
  delete[] sf;
}

void Cell::serialize(std::ofstream& s) const
{
  int n;
  std::set<int>::const_iterator it;

  n = (signed) vertices.size();
  s.write((char*)(&n),sizeof(int));
  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int));
  }
}

void Cell::deserialize(std::ifstream& s)
{
  int i,n,m;

  clear();

  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int));
    vertices.insert(m);
  }
  string_assembly();
}

bool Cell::face(const std::string& f) const
{
  for(int i=0; i<1+dimension(); ++i) {
    if (f == faces[i]) return true;
  }
  return false;
}

bool Cell::exchange(int p,int q)
{
  std::set<int> vx;
  std::set<int>::const_iterator it;
  bool found = false;

  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    if (*it == q) {
      found = true;
      continue;
    }
    vx.insert(*it);
  }
  if (found) {
    // We have to change several things in this case, this
    // simplex's key and also its faces
    vx.insert(p);
    vertices = vx;
    string_assembly();
  }
  return found;
}




