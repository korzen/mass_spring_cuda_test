#include "vertex.h"
#include "autodiff.h"
using namespace std;
using namespace Eigen;
//DECLARE_DIFFSCALAR_BASE();

template<typename T,size_t dim>
vertex<T,dim>::vertex()
{
  ;
}

template<typename T,size_t dim>
vertex<T,dim>::vertex(const myvector<T,dim> &location)
{
  this->location=location;
  this->location_original=location;
  this->location_maybe=location;
  this->velocity=myvector<T,dim>();
  this->force_external=myvector<T,dim>();
  isFixed=0; // all free when at beginning
}

template<typename T,size_t dim>
vertex<T,dim>::~vertex()
{
  ;
}

template<typename T,size_t dim>
int vertex<T,dim>::calPenetrationJacobianAndHessian(size_t index,std::vector<Eigen::Triplet<double > > &tripletsForHessian,Eigen::VectorXd &Jacobian,T intensity,T radius)
{
  {
    if(location.len()>=radius)
      {
	return 0;
      }
    typedef DScalar2<T,VectorXd, MatrixXd> DScalar; // use DScalar2 for calculating gradient and hessian and use DScalar1 for calculating gradient

    size_t i,j,row,col;
    VectorXd x(dim);
    for(j=0;j<dim;++j)
      {
	x(j)=this->location(j);
      }

    DiffScalarBase::setVariableCount(dim);
    DScalar x_d[dim];
    DScalar energy_d=DScalar(0);
    for(i=0;i<dim;i++)
      {
	x_d[i]=DScalar(i,x(i));
      }

    DScalar help_d(0);
    for(i=0;i<dim;i++)
      {
	help_d+=(x_d[i]*x_d[i]);
      }
    energy_d=(radius-sqrt(help_d))*(radius-sqrt(help_d))*intensity*0.5;
    
    MatrixXd grad(dim,1);
    MatrixXd hes(dim,dim);
    grad=energy_d.getGradient();
    hes=energy_d.getHessian();

    for(row=0;row<dim;++row)
      {
	Jacobian(index*dim+row)+=grad(row);

	for(col=0;col<dim;col++)
	  {
	      tripletsForHessian.emplace_back(index*dim+row,index*dim+col,hes(row,col));
	  }

      }
  }  
  return 0;
}
