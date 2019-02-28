#include "myvector.h"
#include "mass_spring_obj.h"
#include "io.h"
#include "para.h"

using namespace std;
using namespace Eigen;
using namespace wtyatzoo;

//using namespace __gnc_cxx;
template<typename T,size_t dim>
mass_spring_obj<T,dim >::mass_spring_obj()
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!mysimplexs.empty())
    {
      mysimplexs.pop_back();
    }
  while(!myedges.empty())
    {
      myedges.pop_back();
    }

  while(!time_norm_pair.empty())
    {
      time_norm_pair.pop_back();
    }
  time_all=norm_Jacobian_cal=0;
  this->mapIndexToLocInMartix=NULL;
   this->mapLocInMatrixToIndex=NULL;
}


template<typename T,size_t dim>
mass_spring_obj<T,dim>::~mass_spring_obj()
{
  // because the harmonic deformation do not need these two array, they can be null
  if(mapIndexToLocInMartix!=NULL)
    {
      delete[] mapIndexToLocInMartix;
    }
  if(mapLocInMatrixToIndex!=NULL)
    {
      delete[] mapLocInMatrixToIndex;
    }
  delete my_matrix_form_cuda_jacobi_solver;
}


template<typename T,size_t dim>
int mass_spring_obj<T,dim >::draw()
{
  myvector<T,dim> v1,v2,v3;
  int i,j;
  for(i=0;i<num_vertex;i++)
    {
      myvertexs[i].normal=myvector<T,dim>();
    }
  for(i=0;i<num_simplexs;i++)
    {
      switch(dim_simplex)
	{
	case 1:
	  {
	    break;
	  }
	case 2:
	  {
	    v1=myvertexs[mysimplexs[i].index_vertex[1]].location-myvertexs[mysimplexs[i].index_vertex[0]].location;
	    v2=myvertexs[mysimplexs[i].index_vertex[2]].location-myvertexs[mysimplexs[i].index_vertex[1]].location;
	    v3=v1.cross(v2);
	    v3.normalize();
	    for(j=0;j<3;j++)
	      {
		myvertexs[mysimplexs[i].index_vertex[j]].normal+=v3;
	      }
	    break;
	  }
	case 3:
	  {
	    break;
	  }
	}
      
    }
  glColor3d(0,0,1);
  //  num_surfaces=0;
  for(i=0;i<num_simplexs;i++)
    {
      switch(dim_simplex)
	{
	case 1:
	  {
	    break;
	  }
	case 2:
	  {
	    glBegin(GL_POLYGON);
	    {
	      for(j=0;j<3;j++)
		{
		  myvertexs[mysimplexs[i].index_vertex[j]].normal.normalize();
		  glNormal3d((double)myvertexs[mysimplexs[i].index_vertex[j]].normal.x[0],(double)myvertexs[mysimplexs[i].index_vertex[j]].normal.x[1],(double)myvertexs[mysimplexs[i].index_vertex[j]].normal.x[2]);
		  glVertex3d((double)myvertexs[mysimplexs[i].index_vertex[j]].location.x[0],(double)myvertexs[mysimplexs[i].index_vertex[j]].location.x[1],(double)myvertexs[mysimplexs[i].index_vertex[j]].location.x[2]);
		}
	    }
	    glEnd();

	    glBegin(GL_POLYGON);
	    {
	      for(j=2;j>=0;j--)
		{
		  myvertexs[mysimplexs[i].index_vertex[j]].normal.normalize();
		  glNormal3d((double)myvertexs[mysimplexs[i].index_vertex[j]].normal.x[0]*-1,(double)myvertexs[mysimplexs[i].index_vertex[j]].normal.x[1]*-1,(double)myvertexs[mysimplexs[i].index_vertex[j]].normal.x[2]*-1);
		  glVertex3d((double)myvertexs[mysimplexs[i].index_vertex[j]].location.x[0],(double)myvertexs[mysimplexs[i].index_vertex[j]].location.x[1],(double)myvertexs[mysimplexs[i].index_vertex[j]].location.x[2]);
		}
	    }
	    glEnd();
	    break;
	  }
	case 3:
	  {
	    break;
	  }	  
	}
    }
  return 0;
}


template<typename T,size_t dim >
mass_spring_obj<T,dim>::mass_spring_obj(std::string input_dir,T dt,T density,int line_search,T weight_line_search,T stiffness,std::string newton_fastMS)
{
  while(!myvertexs.empty())
    {
      myvertexs.pop_back();
    }
  while(!mysimplexs.empty())
    {
      mysimplexs.pop_back();
    }
  while(!myedges.empty())
    {
      myedges.pop_back();
    }
  while(!time_norm_pair.empty())
    {
      time_norm_pair.pop_back();
    }
  time_all=norm_Jacobian_cal=0;

  converge=0;

  pre_succeed=0;
  this->dt=dt;
  this->density=density;
  this->line_search=line_search;
  this->weight_line_search=weight_line_search;
  this->stiffness=stiffness;
  this->newton_fastMS=newton_fastMS;
  
  this->mapIndexToLocInMartix=NULL;
  this->mapLocInMatrixToIndex=NULL;
  size_t i,j,k;
  size_t x,y,z;

  io<T> myio=io<T>();

  myio.getVertexAndSimplex(myvertexs,mysimplexs,dim_simplex,input_dir);
  num_vertex=myvertexs.size();
  printf("num_vertex ::%u \n",num_vertex);
  num_simplexs=mysimplexs.size();
  printf("num_simplexs::%u\n",num_simplexs);
  
  num_all_dof=dim*num_vertex;

  prepare();
  
}

template<typename T,size_t dim >
int mass_spring_obj<T,dim>::prepare()
{
  getEdges();
  init_Energy_now_ForEdge();
  calMassForVertex();
  if(newton_fastMS=="newton")
    {
      
    }
  else if(newton_fastMS=="fastMS_original"||newton_fastMS=="fastMS_ChebyshevSIM")
    {
      max_iteration_num=50; // fast mass spring's max iteration number 
      // when we fix the max iteration number to compare the results of orginal fast mass spring method and the accelerated one by Chebyshev semi-iterative method, we let the both method to iterate to the same max_iteration_num times and compare the results' error to the newton method's result which is regarded as the groudtruth.

      if(newton_fastMS=="fastMS_ChebyshevSIM")
	{
	  rho=0.8; //need to be learned, here we hardcode it as a  constant for cloth simulation.
	  gamma=0.9;
	  start_iteration_num=5;
	  omega=1.0;
	}
    }

  // trick for no penetration of a sphere
  radius=para::radius;
  intensity=para::intensity;
  return 0;
}


template<typename T,size_t dim >
int mass_spring_obj<T,dim >::getEdges()
{
  std::map< std::pair<int,int > ,int> mpFromIndexVertexToIndexEdge;
  mpFromIndexVertexToIndexEdge.clear();

  size_t i;
  size_t index_vertex_for_edge[2];
  size_t index_edge;
  for(i=0;i<num_simplexs;i++)
    {
      if(dim_simplex==1)
	{
	  index_vertex_for_edge[0]=mysimplexs[i].index_vertex[0];
	  index_vertex_for_edge[1]=mysimplexs[i].index_vertex[1];
	  sort(index_vertex_for_edge,index_vertex_for_edge+2);

	  if(mpFromIndexVertexToIndexEdge.find(make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1]))==mpFromIndexVertexToIndexEdge.end())
	    {
	      index_edge=myedges.size();
	      mpFromIndexVertexToIndexEdge[make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1])]=index_edge;

	      myvector<T,dim> loc[2];
	      loc[0]=myvertexs[index_vertex_for_edge[0]].location_original;
	      loc[1]=myvertexs[index_vertex_for_edge[1]].location_original;
	      T rest_length=(loc[0]-loc[1]).len();
	      myedges.push_back(edge<T,dim>(rest_length,index_vertex_for_edge));
	    }
	}
      else if(dim_simplex==2)
	{
	  // edge 0
	  index_vertex_for_edge[0]=mysimplexs[i].index_vertex[0];
	  index_vertex_for_edge[1]=mysimplexs[i].index_vertex[1];
	  sort(index_vertex_for_edge,index_vertex_for_edge+2);

	  if(mpFromIndexVertexToIndexEdge.find(make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1]))==mpFromIndexVertexToIndexEdge.end())
	    {
	      index_edge=myedges.size();
	      mpFromIndexVertexToIndexEdge[make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1])]=index_edge;

	      myvector<T,dim> loc[2];
	      loc[0]=myvertexs[index_vertex_for_edge[0]].location_original;
	      loc[1]=myvertexs[index_vertex_for_edge[1]].location_original;
	      T rest_length=(loc[0]-loc[1]).len();
	      myedges.push_back(edge<T,dim>(rest_length,index_vertex_for_edge));
	    }

	  //edge 1
	  index_vertex_for_edge[0]=mysimplexs[i].index_vertex[0];
	  index_vertex_for_edge[1]=mysimplexs[i].index_vertex[2];
	  sort(index_vertex_for_edge,index_vertex_for_edge+2);

	  if(mpFromIndexVertexToIndexEdge.find(make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1]))==mpFromIndexVertexToIndexEdge.end())
	    {
	      index_edge=myedges.size();
	      mpFromIndexVertexToIndexEdge[make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1])]=index_edge;

	      myvector<T,dim> loc[2];
	      loc[0]=myvertexs[index_vertex_for_edge[0]].location_original;
	      loc[1]=myvertexs[index_vertex_for_edge[1]].location_original;
	      T rest_length=(loc[0]-loc[1]).len();
	      myedges.push_back(edge<T,dim>(rest_length,index_vertex_for_edge));
	    }

	  
	  // here  we use a triangle shape spring system for testing the chebyshev semi-iterative method to accelerate the convergence of the orginal fast mass spring method.
	  //do not use the long edge's spring
	  //edge 2
	  
	  index_vertex_for_edge[0]=mysimplexs[i].index_vertex[1];
	  index_vertex_for_edge[1]=mysimplexs[i].index_vertex[2];
	  sort(index_vertex_for_edge,index_vertex_for_edge+2);

	  if(mpFromIndexVertexToIndexEdge.find(make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1]))==mpFromIndexVertexToIndexEdge.end())
	    {
	      index_edge=myedges.size();
	      mpFromIndexVertexToIndexEdge[make_pair(index_vertex_for_edge[0],index_vertex_for_edge[1])]=index_edge;
	      
	      myvector<T,dim> loc[2];
	      loc[0]=myvertexs[index_vertex_for_edge[0]].location_original;
	      loc[1]=myvertexs[index_vertex_for_edge[1]].location_original;
	      T rest_length=(loc[0]-loc[1]).len();
	      myedges.push_back(edge<T,dim>(rest_length,index_vertex_for_edge));
	    }
	}
      else if(dim_simplex==3)
	{
	  ;
	}
    }
  num_edges=myedges.size();
  printf("num_edges ::%u \n",num_edges);

  printf("dim_simplex ::%d\n",dim_simplex);
  return 0;
}

template<typename T,size_t dim >
int mass_spring_obj<T,dim>::checkFixedOrFree()
{
  int i,j,k,a,b,c;  
  num_cal_dof=0;
  mapIndexToLocInMartix=new size_t[num_all_dof];
  mapLocInMatrixToIndex=new size_t[num_all_dof];
  for(i=0;i<num_vertex;++i)
    {
      if(myvertexs[i].isFixed==0)
	{
	  for(j=0;j<dim;++j)
	    {
	      //   printf("dim %u\n",dim);
	      mapIndexToLocInMartix[i*dim+j]=num_cal_dof;
	      mapLocInMatrixToIndex[num_cal_dof]=i*dim+j;
	      ++num_cal_dof;
	    }
	}
    }
  return 0;
}


template<typename T,size_t dim>
int mass_spring_obj<T,dim>::init_Energy_now_ForEdge()
{
  size_t i;
  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<num_edges;++i)
    {
      myedges[i].energy_now=0;
    }
  return 0;
}


template<typename T,size_t dim>
int mass_spring_obj<T,dim>::calMassForVertex()
{
  size_t i,j;
  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<num_vertex;++i)
    {
      myvertexs[i].mass=0;
    }

  printf("dim_simplex : %u\n",dim_simplex);
  T help=1.0/(dim_simplex+1);

  for(i=0;i<num_simplexs;++i)
    {
      for(j=0;j<dim_simplex+1;j++)
	{
	  myvertexs[mysimplexs[i].index_vertex[j]].mass+=(help*mysimplexs[i].vol*density);
	}
    }
  return 0;
}


template<typename T,size_t dim>
T mass_spring_obj<T,dim>::calElasticEnergy()
{
  size_t i;
  T elasticE=0;
  for(i=0;i<num_edges;++i)
    {
      elasticE+=myedges[i].energy_now;
    }
  return elasticE;
}


template<typename T,size_t dim>
int mass_spring_obj<T,dim>::dynamicSimulator()
{
  size_t i;
  //printf("num_fixed: %u\n",num_fixed);
  

  //  MatrixXd Hessian=MatrixXd::Random(num_all_dof,num_all_dof);
  // VectorXd Jacobian(num_all_dof);
  iteration_num=0;
  converge=0;

  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<num_vertex;++i)
    {
      myvertexs[i].location_lastFrame=myvertexs[i].location;
      myvertexs[i].velocity_lastFrame=myvertexs[i].velocity;
    }
  
  clock_t start,finish;
  T totaltime;
  while(!converge)
    { 
      vector<Triplet<double> > tripletsForHessian;
      VectorXd Jacobian(num_all_dof);
      while(!tripletsForHessian.empty())
	{
	  tripletsForHessian.pop_back();
	}
      Jacobian.fill(0);
      
      printf("--[Inf]:calJacobianAndHessianForEdge\n");
      start=clock();
      calJacobianAndHessianForEdge(tripletsForHessian,Jacobian);
      finish=clock();
      totaltime=(T)(finish-start)/CLOCKS_PER_SEC;
      time_all+=totaltime; 
      printf("Assemble Time Cost: %lf\n",totaltime);

      
      printf("--[Inf]:solve Matrix\n");
      start=clock();
      //  printf("here converge before solve: %d\n",converge);
      solve(tripletsForHessian,Jacobian);
      finish=clock();
      totaltime=(T)(finish-start)/CLOCKS_PER_SEC;
      time_all+=totaltime;
      printf("Solve Time Cost: %lf\n",totaltime);
      //  printf("here converge after solve: %d\n",converge);
    }
  return 0;
}



template<typename T,size_t dim>
int mass_spring_obj<T,dim>::calJacobianAndHessianForEdge(std::vector<Eigen::Triplet<double> > &tripletsForHessian,Eigen::VectorXd &Jacobian)
{
  size_t i,j;

  if(pre_succeed==0)
    {
      for(i=0;i<num_edges;++i)
	{
	  myedges[i].calJacobianAndHessian(myvertexs,tripletsForHessian,Jacobian,stiffness,newton_fastMS);
	}
      if(newton_fastMS=="newton")
	{
	  for(i=0;i<num_vertex;++i)
	    {
	       myvertexs[i].calPenetrationJacobianAndHessian(i,tripletsForHessian,Jacobian,intensity,radius);
	    }
	}
      
    }
  else if(pre_succeed==1)
    {
      for(i=0;i<num_edges;++i)
	{
	  myedges[i].calJacobian(myvertexs,Jacobian,stiffness,newton_fastMS);
	}
    }
  
  //auto diff
  Jacobian*=-1;
  
  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<num_vertex;++i)
    {
      for(j=0;j<dim;++j)
	{
	  Jacobian(i*dim+j)+=myvertexs[i].force_external(j); 
	  
	}
      // printf("force_external %lf %lf %lf\n",myvertexs[i].force_external(0),myvertexs[i].force_external(1),myvertexs[i].force_external(2));


      if(newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original")
	{
	  // trick for sphere non-penetration
	  // Jacobian is the same even though the penetration energy is different!!! 
	  T dis_now=myvertexs[i].location.len();
	  if(dis_now<radius)
	    {
	      myvector<T,dim> penetration_force=myvertexs[i].location*(T)((radius-dis_now)*intensity/dis_now);
	      for(j=0;j<dim;++j)
		{
		  Jacobian(i*dim+j)+=penetration_force(j);
		}
	    }
	}
      
    }

  T EPS=1e-5; // local constant to judge zero for dt-100
  T d1dt;
  if(fabs(dt-100)<EPS) 
    {
      printf("--[INF]:: No mass matrix\n");
      d1dt=0;
    }
  else
    {
      d1dt=1.0/dt;
    }
  T mdt,mdtdt;
  myvector<T,dim> vmdt;
  for(i=0;i<num_vertex;++i)
    {
      mdt=myvertexs[i].mass*d1dt;
      mdtdt=mdt*d1dt;
      vmdt=(myvertexs[i].velocity_lastFrame-(myvertexs[i].location-myvertexs[i].location_lastFrame)*d1dt)*mdt;
      for(j=0;j<dim;++j)
	{
	  if(pre_succeed==0)
	    {
	      tripletsForHessian.emplace_back(i*dim+j,i*dim+j,mdtdt);
	    }
	  Jacobian(i*dim+j)+=vmdt(j);
	}
    }
   
  // cout<< Jacobian<<endl;
  return 0;
}
 

template<typename T,size_t dim>
bool mass_spring_obj<T,dim>::checkInversion(Eigen::VectorXd &dx,Eigen::VectorXd &Jacobian)
{
  size_t i;
  size_t x,y;
  Jacobian.fill(0);
  for(i=0;i<num_cal_dof;++i)
    {
      y=mapLocInMatrixToIndex[i]%dim;
      x=mapLocInMatrixToIndex[i]/dim;
      myvertexs[x].location_maybe(y)=dx(i)+myvertexs[x].location(y);
    }
  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<num_edges;++i)
    {         
      myedges[i].checkInversion(myvertexs,Jacobian,stiffness,newton_fastMS);
    } 
  return false;
}


template<typename T,size_t dim>
bool mass_spring_obj<T,dim>::checkInversion(Eigen::VectorXd &dx)
{
  size_t i;
  size_t x,y;
  for(i=0;i<num_cal_dof;++i)
    {
      y=mapLocInMatrixToIndex[i]%dim;
      x=mapLocInMatrixToIndex[i]/dim;
      myvertexs[x].location_maybe(y)=dx(i)+myvertexs[x].location(y);
    }
  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<num_edges;++i)
    {         
      myedges[i].checkInversion(myvertexs,stiffness,newton_fastMS);
    }  
  return false;
}


template<typename T,size_t dim>
T mass_spring_obj<T,dim>::calEnergyDif()
{
  size_t i;
  T energy_old,energy_new;
  energy_new=energy_old=0;
  for(i=0;i<num_edges;++i)
    {
      energy_new+=myedges[i].energy_maybe;
      energy_old+=myedges[i].energy_now;
    }

  // printf("the first part old energy:%lf new energy:%lf  \n",energy_old,energy_new);
  T EPS=1e-5; // local constant to judge zero for dt-100
  T d1dt;
  if(fabs(dt-100)<EPS) 
    {  
      d1dt=0;
    }
  else
    {
      d1dt=1.0/dt;
    }
  
  myvector<T,dim> help1,help2;
  for(i=0;i<num_vertex;++i)
    {

      if(newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original")
	{
	  // trick for sphere non peneration
	  T dis_now_real=myvertexs[i].location.len();
	  myvector<T,dim> penetration_force=myvertexs[i].location*(T)((radius-dis_now_real)*intensity/dis_now_real);
	  if(dis_now_real<radius)
	    {	  
	      energy_old+=-1*myvertexs[i].location.dot(penetration_force);
	      energy_new+=-1*myvertexs[i].location_maybe.dot(penetration_force);
	    }
	}
      else if(newton_fastMS=="newton")
	{
	  T dis_now_real=myvertexs[i].location.len();
	  if(dis_now_real<radius)
	    {
	      energy_old+=0.5*intensity*(radius-dis_now_real)*(radius-dis_now_real);
	    }
	  T dis_now_maybe=myvertexs[i].location_maybe.len();
	  if(dis_now_maybe<radius)
	    {
	      energy_new+=0.5*intensity*(radius-dis_now_maybe)*(radius-dis_now_maybe);
	    }
	}
      
      
      help1=(myvertexs[i].location_maybe-myvertexs[i].location_lastFrame)*d1dt-myvertexs[i].velocity_lastFrame;
      help2=(myvertexs[i].location-myvertexs[i].location_lastFrame)*d1dt-myvertexs[i].velocity_lastFrame;
      
      energy_new+=help1.len_sq()*myvertexs[i].mass*0.5;
      energy_old+=help2.len_sq()*myvertexs[i].mass*0.5;
      
      energy_new+=-1*myvertexs[i].location_maybe.dot(myvertexs[i].force_external);
      energy_old+=-1*myvertexs[i].location.dot(myvertexs[i].force_external);
    }
  // printf("old energy:%lf new energy:%lf \n",energy_old,energy_new);
  return energy_new-energy_old;
}


template<typename T,size_t dim>
int mass_spring_obj<T,dim>::solve(std::vector<Eigen::Triplet<double> > &tripletsForHessian,Eigen::VectorXd &Jacobian)
{
  size_t i,j;
  size_t ii,jj;
  size_t row=num_cal_dof;
  size_t col=row;
  VectorXd Jacobian_cal(row),dx(row),dx_now(row);
  SparseMatrix<double,Eigen::RowMajor > Hessianspa(row,col);
  vector< Triplet<double > > tripletsForHessianspa;
  T EPS=1e-10; // local constant to judge zero for K(i,j)
  if(newton_fastMS=="newton"||((newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original")&&pre_succeed==0))
    {
      //  SimplicialLLT<SparseMatrix<double>> linearSolver;

      size_t size_Hessian=tripletsForHessian.size();
      
      for(i=0;i<size_Hessian;++i)
	{
	  int row=tripletsForHessian[i].row();
	  int col=tripletsForHessian[i].col();
	  double val=tripletsForHessian[i].value();
	  int map_row=mapIndexToLocInMartix[row];
	  int map_col=mapIndexToLocInMartix[col];
	  if(map_row>=0&&map_row<num_cal_dof&&map_col>=0&&map_col<num_cal_dof&&fabs(val)>=EPS)
	    {
	      //   printf("row:%d col:%d val:%lf\n",row,col,val);
	      tripletsForHessianspa.emplace_back(map_row,map_col,val);
	    }
	}
    }
  
  omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
  for(i=0;i<row;++i)
    {
      Jacobian_cal(i)=Jacobian(mapLocInMatrixToIndex[i]);
    }

  double norm_Jacobian_cal=0;
  for(i=0;i<row;++i)
    {
      norm_Jacobian_cal+=(Jacobian_cal(i)*Jacobian_cal(i));
    }
  norm_Jacobian_cal=sqrt(norm_Jacobian_cal);
  printf("norm_Jacobian_cal: %lf\n",norm_Jacobian_cal);  
  time_norm_pair.push_back(make_pair(time_all,norm_Jacobian_cal));

  SparseMatrix<double,Eigen::RowMajor > diagMatrix_spa(row,col);
  if(newton_fastMS=="newton"||((newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original")&&pre_succeed==0))
    {      
      vector< Triplet<double > > tripletsFordiagMatrix_spa;

      for(i=0;i<row;i++)
	{
	  tripletsFordiagMatrix_spa.emplace_back(i,i,1e-8);
	}
      diagMatrix_spa.setFromTriplets(tripletsFordiagMatrix_spa.begin(),tripletsFordiagMatrix_spa.end());
      diagMatrix_spa.makeCompressed();
      Hessianspa.setFromTriplets(tripletsForHessianspa.begin(),tripletsForHessianspa.end());
      Hessianspa.makeCompressed();
    }
  
  while(1)
    {
      if(newton_fastMS=="newton")
	{	  
	  linearSolver.compute(Hessianspa);
	  int info=(int)linearSolver.info();
	  //  cout<<info<<"info"<<endl;
	  if(info==0)
	    {
	      break;
	    }
	  else if(info==1)
	    {
	      Hessianspa=Hessianspa+diagMatrix_spa;
	      diagMatrix_spa=diagMatrix_spa*2;
	    }
	}
      else if(newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original")
	{
	  //  printf("[info]:: here is fast mass spring method\n");
	  if(pre_succeed==0)
	    {
	      // printf("[info]:: the first time to decompose the matrix\n");
	      //	      linearSolver.compute(Hessianspa);
	      my_matrix_form_cuda_jacobi_solver=new matrix_form_cuda_jacobi_solver<double>(Hessianspa.valuePtr(),Hessianspa.innerIndexPtr(),Hessianspa.outerIndexPtr(),Hessianspa.nonZeros(),num_cal_dof);
	      pre_succeed=1;
	      break;
	    }	  
	  else if(pre_succeed==1)
	    {
	      //printf("[info]:: no decomposition here!!!\n");
	      break;
	    }
	}
      
    }

  if(newton_fastMS=="newton")
    {
      //      printf("linearSolver\n");
      // cout<<Eigen::Success<<endl;
      dx=linearSolver.solve(Jacobian_cal);  
    }
  else if(newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original")
    {
      //printf("my_matrix_form_cuda_jacobi_solver\n");
      my_matrix_form_cuda_jacobi_solver->apply(Jacobian_cal.data(),dx.data(),1e-10,400,0);
    }

  T h=2;
  size_t max_lineSearch=20;
  T energy_dif,threshold;
  bool find=0;

  VectorXd Jacobian_loc_maybe(num_all_dof);
  VectorXd Jacobian_cal_loc_maybe(row);
  
  if(line_search==1)
    {
      printf("[control info]:: line_search open\n");
      for(i=0;i<max_lineSearch;++i)
	{
	  h*=0.5;
	  dx_now=dx*h;
	  if(checkInversion(dx_now,Jacobian_loc_maybe)==1)
	    {
	      continue;
	    }
	  else
	    {

	      /*
	      {
		Jacobian_loc_maybe*=-1;
	      }  
	      for(ii=0;ii<num_vertex;++ii)
		{
		  for(jj=0;jj<dim;++jj)
		    {
		      Jacobian_loc_maybe(ii*dim+jj)+=myvertexs[ii].force_external(jj); 
		    }

		  
		  // trick for sphere non-penetration
		  T dis_now=myvertexs[ii].location_maybe.len();
		  if(dis_now<radius)
		    {
		      myvector<T,dim> peneration_force=myvertexs[ii].location_maybe*(T)((radius-dis_now)*intensity/dis_now);
		      for(jj=0;jj<dim;++jj)
			{
			  Jacobian_loc_maybe(ii*dim+jj)+=peneration_force(jj);
			}
		     }
		  
		}
	      for(ii=0;ii<row;++ii)
		{
		  Jacobian_cal_loc_maybe(ii)=Jacobian_loc_maybe(mapLocInMatrixToIndex[ii]);
		}

	      double Jdx=dx.dot(Jacobian_cal_loc_maybe);
	      double Jdx_now=dx.dot(Jacobian_cal);

	      */
	      energy_dif=calEnergyDif();
	      threshold=weight_line_search*dx_now.dot(Jacobian_cal)*-1;
	      if((energy_dif<=threshold)/*&&(fabs(Jdx)<=fabs(Jdx_now))*/)
		{
		  find=1;
		  break;
		}
	    }
	}
    }
  else if(line_search==0)
    {
      printf("[control info]:: line_search closed\n");
      h=1;
      dx_now=dx*h;
      if(checkInversion(dx_now)==1)
	{
	  for(i=0;i<max_lineSearch;++i)
	    {
	      h*=0.5;
	      dx_now=dx*h;
	      if(checkInversion(dx_now)==1)
		{
		  continue;
		}
	      else
		{
		  energy_dif=calEnergyDif();
		  printf("energy_dif %lf\n",energy_dif);
		  threshold=weight_line_search*dx_now.dot(Jacobian_cal)*-1;
		  if(energy_dif<=threshold)
		    {
		      find=1;
		      break;
		    }
		}
	    }
	}
      find=1;  
      energy_dif=calEnergyDif();
    }  
      
  EPS=1e-4; // local constant to judge zero for energy_dif
  if(find==1)
    {
      if(fabs(norm_Jacobian_cal)<=EPS&&newton_fastMS=="newton")
	{
          printf("bingo\n");
	  converge=1;
	}

      omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
      for(i=0;i<num_vertex;++i)
	{
	  myvertexs[i].location=myvertexs[i].location_maybe;
	}
      omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
      for(i=0;i<num_edges;++i)
	{
	  myedges[i].energy_now=myedges[i].energy_maybe;
	}
    }
  else if(find==0)
    {
      converge=1;
    }

  ++iteration_num;
  /*
    if(iteration_num==max_iteration_num&&(newton_fastMS=="fastMS_ChebyshevSIM"||newton_fastMS=="fastMS_original"))
    {
    converge=1;
    }
  */
 
  if(newton_fastMS=="fastMS_original")
    {
      if(iteration_num==max_iteration_num)
	{
	  converge=1;
	}
    }
  else if(newton_fastMS=="fastMS_ChebyshevSIM")
    {
      if(iteration_num==1)
	{
	  ;
	}
      else if(iteration_num>=2)
	{
	  if(iteration_num<start_iteration_num)
	    {
	      omega=1;
	    }
	  else if(iteration_num==start_iteration_num)
	    {
	      omega=2.0/(2.0-rho*rho);
	    }
	  else if(iteration_num>start_iteration_num)
	    {
	      omega=4.0/(4.0-rho*rho*omega);
	    }
	  if(iteration_num>=3)
	    {
	      for(i=0;i<num_vertex;i++)
		{
		  myvertexs[i].location=omega*(gamma*(myvertexs[i].location-myvertexs[i].location_lastIteration)+myvertexs[i].location_lastIteration-myvertexs[i].location_lastlastIteration)+myvertexs[i].location_lastlastIteration;
		}
	    }	  	  
	}
      omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
      for(i=0;i<num_vertex;i++)
	{
	  myvertexs[i].location_lastlastIteration=myvertexs[i].location_lastIteration;
	}
      omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
      for(i=0;i<num_vertex;i++)
	{
	  myvertexs[i].location_lastIteration=myvertexs[i].location;
	}

      if(iteration_num==max_iteration_num)
	{
	  converge=1;
	}
    }
  
  if(converge==1)
    {
      T EPS=1e-5; // local constant to judge zero for dt-100
      T d1dt;
      if(fabs(dt-100)<EPS) 
	{
	  d1dt=0;
	}
      else
	{
	  d1dt=1.0/dt;
	}
      omp_set_num_threads(para::num_threads);  
  #pragma omp parallel for 
      for(i=0;i<num_vertex;++i)
	{      
	  myvertexs[i].velocity=(myvertexs[i].location-myvertexs[i].location_lastFrame)*d1dt;
	}
    }
  return 0;
}

