//#include <windows.h>
#include "head.h"
#include "camera_test.h"
#include "object_creator.h"
#include "simulator.h"
#include "para.h"
/*
class camera_test
{
public:
  myvector eyepoint;
  myvector lookpoint;
  myvector vector_to;
  double length;
  double baseLength;
  myvector vector_lr;
  myvector vector_ud;
  double baseAngle;
  myvector vectorUp;
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
  int matrixMult(double matrix[16], myvector& vector_now);
};
*/
camera_test<double>* mycamera_test;
simulator<float> *mysimulator;

int kind;
int count_vtk;
int frame_now;
int save;

using namespace std;

void saveBmp(const char* name,int width,int height,unsigned char* data)
{
  
  /*
  BITMAPFILEHEADER hdr;
  BITMAPINFOHEADER infoHdr;
  infoHdr.biSize=40;
  infoHdr.biWidth=width;
  infoHdr.biHeight =height;
  infoHdr.biPlanes =1;
  infoHdr.biBitCount =24;
  infoHdr.biCompression =0;
  infoHdr.biSizeImage=width*height*3;
  infoHdr.biXPelsPerMeter=0;
  infoHdr.biYPelsPerMeter=0;
  infoHdr.biClrUsed=0;
  infoHdr.biClrImportant=0;
  
  hdr.bfType = 0x4D42;
  hdr.bfReserved1 = 0;
  hdr.bfReserved2 = 0;
  hdr.bfOffBits = 54;
  hdr.bfSize =(DWORD)(sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+width* height * 3);
  
  unsigned char x;
  for(int i=0;i<infoHdr.biSizeImage;i+=3)
    {
      x=data[i];
      data[i]=data[i+2];
      data[i+2]=x;
    }
  
  FILE *fd;
  fd=fopen(name,"wb");
  fwrite(&hdr,sizeof(BITMAPFILEHEADER),1,fd);
  fwrite(&infoHdr,sizeof(BITMAPINFOHEADER),1,fd);
  fwrite(data,width* height*3,1,fd);
  fclose(fd);*/
}


void saveSceneImage()
{
  /*
  int view[4];
  glGetIntegerv(GL_VIEWPORT,view);
  int bufferSize=view[2]*view[3]*3;
  void* color_buffer=malloc(bufferSize);
  glReadPixels(view[0],view[1],view[2],view[3],GL_RGB,GL_UNSIGNED_BYTE,color_buffer);
  
  string st,loc;
  loc="./img/";
  stringstream  ss;
  ss<<count_vtk;
  ss>>st;
  st+=".bmp";
  loc+=st;
  saveBmp(loc.c_str(),view[2],view[3],(unsigned char*)color_buffer);
  free(color_buffer);	 */
}


void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
  mysimulator->simulate(frame_now);
  frame_now++;

  if(frame_now==para::frame)
    {
      {
	delete mycamera_test;
	delete mysimulator;
	mycamera_test=NULL;
	mysimulator=NULL;
      }
      exit(0);
    }
  if(save==1)
    {
      saveSceneImage();
      //saveAsVTK(); 
      save=0;
      count_vtk++;
    }
  glutSwapBuffers();   glFlush();
  
  glutPostRedisplay();
}

void makeSmooth()
{
  glEnable(GL_LINE_SMOOTH); 
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
  glLineWidth(1.0);
  
  glPointSize(1.0);
  
  glEnable(GLUT_MULTISAMPLE);
  glEnable(GL_POLYGON_SMOOTH);
}

void lightControl()
{
  float light_position[] ={0,-6,0,0};
  float white_light[] ={1,1,1.5,1};
  float lmod_ambient[]={0.5,0.5,0.65,1}; 
  float spot_direction[]={0,1,0};
  glShadeModel(GL_SMOOTH);
  glLightfv(GL_LIGHT0,GL_POSITION,light_position);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,white_light);
  glLightfv(GL_LIGHT0,GL_SPECULAR,white_light);
  glLightf(GL_LIGHT0,GL_SPOT_CUTOFF,90);
  glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,spot_direction);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,lmod_ambient);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  
  glEnable(GL_NORMALIZE);
  
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST); 
  glClearDepth(1.0);
  
  glEnable(GL_TEXTURE_2D);
}

void init()
{
  glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA) ;  
  glutInitWindowSize(400,400) ; 
  glutInitWindowPosition(100,100) ;  
  glutCreateWindow("fast_mass_spring_test") ; 	
  glClearColor(1,1,1,1);
  makeSmooth(); 	
  
  
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  
  lightControl();
  
  mysimulator=new simulator<float>();
  mycamera_test=new camera_test<double>();
  kind=0;
  count_vtk=0;
  frame_now=0;
}

void reshape(int w,int h) 
{
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  
  gluPerspective(60.0,(double)w/(double)h,0.01,50000);
  //from the function:  void APIENTRY gluPerspective ( GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar) 
  glMatrixMode(GL_MODELVIEW);
  mycamera_test->see();
}


void keyboard(unsigned char key,int x,int y)
{
  switch(key)
    {
    case 'a':
      mycamera_test->rotate_LR(1);  
      break;
    case 'd':
      mycamera_test->rotate_LR(-1); 
      break;
    case 'w':
      mycamera_test->rotate_UD(1);
      break;
    case 's':
      mycamera_test->rotate_UD(-1);
      break;
    case 'k':
      mycamera_test->move(-1);
      break;
    case 'l':
      mycamera_test->move(1);
      break;
    case 'p':
      // myobject->signal=1;
      // myobject->add+=5;
      break;
    case 'o':
      // myobject->signal=2;
      break;
    case 'q':
      save=1;
      break;
    case 27:
      {
	delete mycamera_test;
	delete mysimulator;
	mycamera_test=NULL;
	mysimulator=NULL;
      }
      exit(0);
      break;
    }
  glutPostRedisplay();
}


void readcmdline(int argc, char* argv[],boost::property_tree::ptree &para_tree)
{
  size_t i;
  for(i=1;i<argc;++i)
    {
      string para_here=argv[i];
      size_t pos=para_here.find("=");
      if(pos!= string::npos)
	{
	  string key=para_here.substr(0,pos);
	  string value=para_here.substr(pos+1);
	  para_tree.put(key+".value",value);
	  printf("--[cmdline para] %s %s \n",key.c_str(),value.c_str());
	}
    }
  return;
}

int main(int argc,char**argv)
{
  /*
  // testing from here 
  object_creator<double>* myobject_creator=new object_creator<double>();
  delete myobject_creator;
  myobject_creator=NULL;
  
  */
  boost::property_tree::ptree para_tree;
  readcmdline(argc,argv,para_tree);

  para::out_dir_simulator=para_tree.get<string>("out_dir_simulator.value"); 
  para::simulation_type=para_tree.get<string>("simulation.value","static"); 
  para::newton_fastMS=para_tree.get<string>("newton_fastMS.value");
  para::dt=para_tree.get<double>("dt.value",0.01);
  para::density=para_tree.get<double>("density.value",10);
  para::stiffness=para_tree.get<double>("stiffness.value",8000);
  para::frame=para_tree.get<int>("frame.value",3000);
  para::line_search=para_tree.get<int>("line_search.value",1);
  para::weight_line_search=para_tree.get<double>("weight_line_search.value",1e-5);
  para::gravity=para_tree.get<double>("gravity.value",9.8); 
  para::object_name=para_tree.get<string>("object_name.value"); 
  para::input_object=para_tree.get<string>("input_object.value");
  para::input_constraint=para_tree.get<string>("input_constraint.value");
  para::force_function=para_tree.get<string>("force_function.value");
  para::radius=para_tree.get<double>("radius.value");
  para::intensity=para_tree.get<double>("intensity.value");

  /*
  glutInit(&argc,argv) ; 
  init();	
  
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMainLoop() ;*/

  mysimulator=new simulator<float>();
  frame_now=0;
  for(size_t i=0;i<para::frame;++i)
    {
      mysimulator->simulate(frame_now);
      frame_now++;
    }
  delete mysimulator;
  mysimulator=NULL;
  
  return 0 ;
}

