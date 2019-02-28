#include "camera_test.h"
#include "head.h"
//using namespace std;

template<typename T>
camera_test<T>::camera_test()
{
  lookpoint = myvector<T,3>(0, 0, 0);
  this->length = 5;
  eyepoint = myvector<T,3>(0, -5, 0);

  this->baseAngle = 2;
  this->baseLength = 0.1;

  vector_lr = myvector<T,3>(0, 0, 1);
  vector_ud = myvector<T,3>(1, 0, 0);

  calVectorTo();
  calVectorUp();

}

template<typename T>
camera_test<T>::~camera_test()
{
}

template<typename T>
int camera_test<T>::calVectorTo()
{
  vector_to=eyepoint-lookpoint;
  return 0;
}

template<typename T>
int camera_test<T>::calVectorUp()
{
  vectorUp.x[0]=vector_to.x[1]*vector_ud.x[2]-vector_ud.x[1]*vector_to.x[2];
  vectorUp.x[1]=vector_to.x[2]*vector_ud.x[0]-vector_to.x[0]*vector_ud.x[2];
  vectorUp.x[2]=vector_to.x[0]*vector_ud.x[1]-vector_ud.x[0]*vector_to.x[1];
  return 0;
}
template<typename T>
int camera_test<T>::calLength()
{
  length=sqrt(vector_to.len_sq());
  return 0;
}

template<typename T>
int camera_test<T>::calVector_UD()
{
  calVectorTo();

  vector_ud.x[0]=vectorUp.x[1]*vector_to.x[2]-vector_to.x[1]*vectorUp.x[2];
  vector_ud.x[1]=vector_to.x[0]*vectorUp.x[2]-vectorUp.x[0]*vector_to.x[2];
  vector_ud.x[2]=vectorUp.x[0]*vector_to.x[1]-vector_to.x[0]*vectorUp.x[1];
  return 0;
}

template<typename T>
int camera_test<T>::calEyePoint()
{
  eyepoint=vector_to+lookpoint;
  return 0;
}

template<typename T>
int camera_test<T>::matrixMult(double matrix[16],myvector<T,3>& vector_now)
{
  T vec[3];
  int i;
  for(i=0;i<3;i++)
    {
      vec[i]=vector_now.x[0]*matrix[i]+vector_now.x[1]*matrix[i+4]+vector_now.x[2]*matrix[i+8]+matrix[i+12]; 
    }
  vector_now.set(vec[0],vec[1],vec[2]);
  return 0;
}

template<typename T>
int  camera_test<T>::see()
{
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();
  gluLookAt(eyepoint.x[0],eyepoint.x[1],eyepoint.x[2],lookpoint.x[0],lookpoint.x[1],lookpoint.x[2],vectorUp.x[0],vectorUp.x[1],vectorUp.x[2]);
  return 0;
}

template<typename T>
int camera_test<T>::rotate_LR(int lr)
{
  double matrix_now[16];
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glRotated(lr*baseAngle,vector_lr.x[0],vector_lr.x[1],vector_lr.x[2]);
  glGetDoublev(GL_MODELVIEW_MATRIX,matrix_now);  
  glPopMatrix();
  matrixMult(matrix_now,vector_to); 
  calEyePoint();
  calVector_UD();

  see();
  return 0;
}

template<typename T>
int camera_test<T>::rotate_UD(int ud)
{
  double matrix_now[16]; 
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glRotated(ud*baseAngle,vector_ud.x[0],vector_ud.x[1],vector_ud.x[2]);
  glGetDoublev(GL_MODELVIEW_MATRIX,matrix_now);
  glPopMatrix();
  matrixMult(matrix_now,vector_to);
  calEyePoint();
  calVectorUp(); 

  see();
  return 0;
}

template<typename T>
int camera_test<T>::move(int cf) 
{
  double matrix_now[16];
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslated(vector_to.x[0]*cf*baseLength/length,vector_to.x[1]*cf*baseLength/length,vector_to.x[2]*cf*baseLength/length);
  glGetDoublev(GL_MODELVIEW_MATRIX,matrix_now);
  glPopMatrix();
  matrixMult(matrix_now,vector_to);
  calEyePoint();
  calVectorTo();
  calLength();
	
  see();
  return 0;
}
