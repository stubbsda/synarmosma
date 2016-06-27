#include "cell.h"

using namespace SYNARMOSMA;

Cell::Cell()
{

}

Cell::Cell(int n)
{
  for(int i=0; i<=n; ++i) {
    vertices.insert(i);
  }
  calculate_faces();
}

Cell::Cell(int v1,int v2)
{
  // A specialized constructor for 1-cells
  vertices.insert(v1);
  faces.push_back(vertices);
  vertices.clear();
  vertices.insert(v2);
  faces.push_back(vertices);
  vertices.insert(v1);
}

Cell::Cell(const std::set<int>& v)
{
  vertices = v;
  calculate_faces();
}

Cell::Cell(const Cell& source)
{
  vertices = source.vertices;
  entourage = source.entourage;
  faces = source.faces;
}

Cell& Cell::operator =(const Cell& source)
{
  if (this == &source) return *this;

  vertices = source.vertices;
  entourage = source.entourage;
  faces = source.faces;

  return *this;
}

void Cell::initialize(int v1,int v2)
{
  clear();

  vertices.insert(v1);
  faces.push_back(vertices);
  vertices.clear();
  vertices.insert(v2);
  faces.push_back(vertices);
  vertices.insert(v1);
}

void Cell::initialize(const std::set<int>& vx)
{
  clear();
  vertices = vx;
  calculate_faces();
}

Cell::~Cell()
{

}

void Cell::clear()
{
  vertices.clear();
  entourage.clear();
  faces.clear();
}

void Cell::calculate_faces()
{
  int i,j;
  std::set<int> S;
  std::set<int>::const_iterator it;
  const int n = (signed) vertices.size(); 

  faces.clear();

  for(i=0; i<n; ++i) {
    j = -1;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      ++j;
      if (j == i) continue;
      S.insert(*it);
    }
    faces.push_back(S);
    S.clear();
  }
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
  calculate_faces();
}

bool Cell::face(const std::set<int>& f) const
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
    calculate_faces();
  }
  return found;
}




