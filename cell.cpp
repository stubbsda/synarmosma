#include "cell.h"

using namespace SYNARMOSMA;

Cell::Cell()
{

}

Cell::Cell(int n)
{
  if (n < 1) throw std::invalid_argument("The number of vertices in the Cell constructor must be greater than zero!");
  for(int i=0; i<=n; ++i) {
    vertices.insert(i);
  }
  calculate_faces();
}

Cell::Cell(int v1,int v2)
{
  if (v1 == v2) throw std::invalid_argument("The edge vertices in the Cell constructor must be distinct!");
  // A specialized constructor for 1-cells
  vertices.insert(v1);
  faces.push_back(vertices);
  vertices.clear();
  vertices.insert(v2);
  faces.push_back(vertices);
  vertices.insert(v1);
}

Cell::Cell(const std::set<int>& vx)
{
  if (vx.empty()) throw std::invalid_argument("The vertex set in the Cell constructor must not be empty!");
  vertices = vx;
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
  if (v1 == v2) throw std::invalid_argument("The edge vertices in the Cell constructor must be distinct!");

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
  if (vx.empty()) throw std::invalid_argument("The vertex set in the Cell constructor must not be empty!");

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

int Cell::serialize(std::ofstream& s) const
{
  int n,count = 0;
  std::set<int>::const_iterator it;

  n = (signed) vertices.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=vertices.begin(); it!=vertices.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  }
  return count;
}

int Cell::deserialize(std::ifstream& s)
{
  int i,n,m,count = 0;

  clear();

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    vertices.insert(m);
  }
  calculate_faces();
  return count;
}

bool Cell::exchange(int p,int q)
{
  if (vertices.count(q) > 0) {
    vertices.erase(q);    
    vertices.insert(p);
    calculate_faces();
    return true;
  }
  return false;
}




