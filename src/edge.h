#ifndef _EDGE_
#define _EDGE_
#include "head.h"
#include "eigen_head.h"
#include "myvector.h"
#include "vertex.h"

// a data structure for mass spring modeling 
template<typename T,size_t dim>
class edge
{
 public:
  size_t index_vertex[2]; //from zero to one
  T rest_length;
  T energy_now;
  T energy_maybe;
  
  edge(){}
  ~edge(){}
  edge(T rest_length,const size_t (&index_vertex)[2]);
  int calJacobianAndHessian(std::vector< vertex<T,dim > > &myvertexs,std::vector<Eigen::Triplet<double > > &tripletsForHessian,Eigen::VectorXd &Jacobian,T stiffness,const std::string &newton_fastMS);
  int calJacobian(std::vector< vertex<T,dim > > &myvertexs,Eigen::VectorXd &Jacobian,T stiffness,const std::string &newton_fastMS);
  bool checkInversion(std::vector< vertex<T,dim > > &myvertexs,T stiffness,const std::string &newton_fastMS);
  bool checkInversion(std::vector< vertex<T,dim > > &myvertexs,Eigen::VectorXd &Jacobian,T stiffness,const std::string &newton_fastMS);
  
};

template class edge<double,3>;
template class edge<float,3>;
#endif
