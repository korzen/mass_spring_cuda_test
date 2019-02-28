#ifndef _CAMERA_TEST_
#define _CAMERA_TEST_
#include "myvector.h"

// for generating testing demo use,it will and can be replaced by the Physika's camera class 
template<typename T>
class camera_test
{
 public:
  myvector<T,3> eyepoint;
  myvector<T,3> lookpoint;
  myvector<T,3> vector_to;
  T length;
  T baseLength;
  myvector<T,3> vector_lr;
  myvector<T,3> vector_ud;
  T baseAngle;
  myvector<T,3> vectorUp;
  camera_test();
  ~camera_test();

  int calVectorTo();
  int calLength();
  int calVectorUp();
  int calVector_UD();
  int calEyePoint();
  int see();
  int rotate_LR(int lr);
  int rotate_UD(int ud);
  int move(int cf);
  int matrixMult(double matrix[16], myvector<T,3>& vector_now);
};

template class camera_test<double>;
template class camera_test<float>;
#endif
