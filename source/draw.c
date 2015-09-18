#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "glut.h"
#include "path.h"



typedef struct atom {
    
    float      radius;
            
    vector     r;
    quaternion q;

    int        no;

} Particle;


Slice   *currentslice;
int     Npart;                  /* Number of particles */
int     Precision=10;              /* Number of sides per sphere, etc */
int     Quit = FALSE;           /* Quit flag */
int     Run  = FALSE;           /* Run flag */
int     Drawall =3;          /* Draw flag*/


Particle *pos;                   /* Storage of positions, etc */
vector  Boxl;                   /* Size of the periodic box */
vector  *Circle=NULL;


float	Far = 200.0;		/* The far clipping plane */
float   Boxdistance =1.25;
float	Distance = -10;		/* Distance, Azimuth, Inclination and Twist */
GLfloat	xRotation = -90,		/* define the viewing position */
        yRotation = 0,
        zRotation = 0;
char string[] = "A";


struct glstattype {
  
  int Animate,
      BackwardAnimation,
      OneUp            ,
      OneDown;

} glstat;



static float MaterialRed[]     =   { 0.7, 0.0, 0.0, 0.7 };
static float MaterialGreen[]   =   { 0.1, 0.5, 0.2, 1.0 };
static float MaterialBlue[]    =   { 0.0, 0.0, 0.7, 1.0 };
static float MaterialSome[]    =   { 0.0, 0.7, 0.0, 1.0 };
static float MaterialYellow[]  =   { 0.7, 0.7, 0.0, 1.0 };
static float MaterialCyan[]    =   { 0., 0.7, 0.7, 1.0 };
static float MaterialMagenta[]    =   { 0.7, 0., 0.7, 1.0 };


static float front_shininess[] =   {60.0};
static float front_specular[]  =   { 0.7, 0.7, 0.7, 1.0 };
static float ambient[]         =   { 0.0, 0.0, 0.0, 1.0 };
static float diffuse[]         =   { 1.0, 1.0, 1.0, 1.0 };
static float position0[]       =   { 1.0, 1.0, 1.0, 0.0 };
static float position1[]       =   {-1.0,-1.0, 1.0, 0.0 };
static float lmodel_ambient[]  =   { 0.5, 0.5, 0.5, 1.0 };
static float lmodel_twoside[]  =   {GL_TRUE};

static float vert[17][3] = {
    {-1.0,-1.0,-1.0},
    {-1.0,-1.0, 1.0},
    { 1.0,-1.0, 1.0},
    { 1.0,-1.0,-1.0},
    {-1.0,-1.0,-1.0},
    {-1.0, 1.0,-1.0},
    { 1.0, 1.0,-1.0},
    { 1.0,-1.0,-1.0},
    { 1.0,-1.0, 1.0},
    { 1.0, 1.0, 1.0},
    { 1.0, 1.0, 1.0},
    { 1.0, 1.0,-1.0},
    { 1.0, 1.0, 1.0},
    {-1.0, 1.0, 1.0},
    {-1.0, 1.0,-1.0},
    {-1.0, 1.0, 1.0},
    {-1.0,-1.0, 1.0},
};



static void DrawBitmapString(void *font, const char *string)
{
    int i;

    for (i = 0; string[i]; i++)
        glutBitmapCharacter(font, string[i]);
}

static void DrawStrokeString(void *font, const char *string)
{
    int i;

    for (i = 0; string[i]; i++)
        glutStrokeCharacter(font, string[i]);
}



//int readpath() {
//
//    int i,j,k;
//    char dum[20],dum2[20],nodum[20];
//    static int count;
//    int type,nrep,nslices,state_index;  
//    vector boxl;
//    Slice current;
//
//    fscanf(sys.fprc," ");
//    if(feof(sys.fprc)) return 1;
//    fscanf(sys.fprc,"%d %d %d %d %s",&nslices,&nrep,&state_index,&type,dum);
//    path.nslices=nslices;
//    path.type=nrep;
//    fscanf(sys.fprc,"%lf %lf %lf %s",&boxl.x,&boxl.y,&boxl.z,dum);
//    if(sys.sim_type==0) {
//        printf("reading path %d for %d slices\n",count,nslices);
//        printf("Initial state is %d and the path has type %d\n",state_index,type);
//        }
//    for(i=0; i<nslices; i++) {
//        readslice(&current);
//        reset_center(&current);
//        slice[i]=current;
//        if(nslices==1) {
//            printf("Size of the maximum cluster: %d\n",clustersize(&(slice[i])));
//            }
//        }
//    count++;
//    if(sys.sim_type==0) printf("Read complete path %d, now on display\n",count);
//
//    return 0;
//}
//
//
//void readslice(Slice *current) {
//
//    char dum[20];
//    int i;
//
//    fscanf(sys.fprc,"%s %lf",dum,&current->energy);
//    for(i=0; i<sys.npart; i++) {
//        fscanf(sys.fprc,"%lf %lf %lf %lf %lf %lf %lf", 
//         &current->pts[i].r.x,&current->pts[i].r.y,&current->pts[i].r.z,
//          &current->pts[i].q.r,&current->pts[i].q.u.x,&current->pts[i].q.u.y,&current->pts[i].q.u.z);
//        }
//
//    return;
//}



void convert_coordinates(Slice *psl )
{
    int ipart;

    //reset_center(psl);

    for(ipart=0; ipart<sys.npart; ipart++)    {
        pos[ipart].r = psl->pts[ipart].r;
        pos[ipart].no = ipart;
        pos[ipart].radius =0.5;
        pos[ipart].q = psl->pts[ipart].q;
    }


}


void reset_center(Slice *psl) {

    int ipart;
    vector dr,rcom;

    for( ipart=1; ipart<sys.npart; ipart++) {
        vector_minus(psl->pts[ipart].r,psl->pts[0].r,dr);
        pbc(dr,sys.boxl);
        vector_add(dr,psl->pts[0].r,psl->pts[ipart].r);
    }

    rcom=nulvec;
    for(ipart=0; ipart<sys.npart; ipart++) {
        vector_add(rcom,psl->pts[ipart].r,rcom);
    }

    scalar_divide(rcom,sys.npart,rcom);
    for( ipart=0; ipart<sys.npart; ipart++) {
        vector_minus(psl->pts[ipart].r,rcom,psl->pts[ipart].r);
    }

    return;
}


vector findnorm(vector vec1, vector vec2, vector vec3) {

	vector norm,edge1,edge2;

	vector_minus(vec2,vec1,edge1);
	vector_minus(vec3,vec1,edge2);

	scalar_divide(edge1,sqrt(vector_inp(edge1,edge1)),edge1);
	scalar_divide(edge2,sqrt(vector_inp(edge2,edge2)),edge2);

	vector_cross(edge1,edge2,norm);

	return norm;
}


void getdrawrotmatrix(quaternion q, float *mat) {

    float  q0,q1,q2,q3,q02,q12,q22,q32;

    q0 = q.q0;
    q1 = q.q1;
    q2 = q.q2;
    q3 = q.q3;
    
    q02 = q0*q0;
    q12 = q1*q1;
    q22 = q2*q2;
    q32 = q3*q3;
    
    mat[0] = q02 + q12 - q22 - q32;  mat[4] = 2*(q1*q2 - q0*q3)    ;  mat[8] = 2*(q1*q3 + q0*q2) ; mat[12]=0;
    mat[1] = 2*(q1*q2 + q0*q3)    ;  mat[5] = q02 - q12 + q22 - q32;  mat[9] = 2*(q2*q3 - q0*q1) ; mat[13]=0;
    mat[2] = 2*(q1*q3 - q0*q2)    ;  mat[6] = 2*(q2*q3 + q0*q1)    ;  mat[10]= q02-q12- q22 + q32; mat[14]=0;
    mat[3] = 0                    ;  mat[7] =0                     ;  mat[11]=0                  ; mat[15]=1;

    return;
}


quaternion quatVecToVec(vector vec1, vector vec2) {

	quaternion qrot;
	vector rotvec,xvec,yvec;
	double angle,s,invs,l;

	angle=vector_inp(vec1,vec2);
	if(angle>=1) {
		qrot.q0=1;
		qrot.q1=0;
		qrot.q2=0;
		qrot.q3=0;
	}

	else if(angle < -0.999999) {
		//xvec.x=1;
		//xvec.y=0;
		//xvec.z=0;
		//vector_cross(xvec,vec2,rotvec);	
		//l=sqrt(vector_inp(rotvec,rotvec));
		//if(l==0) {
		//	yvec.x=0;
		//	yvec.y=1;
		//	yvec.z=0;
		//	vector_cross(yvec,vec2,rotvec);
		//	l=sqrt(vector_inp(rotvec,rotvec));
		//}
		//scalar_divide(rotvec,l,rotvec);
		//qrot.r=PI;
		//qrot.u=rotvec;
		//scdivide_quat(qrot,sqrt(quat_inp(qrot,qrot)),qrot);
		qrot.q0=0;
		qrot.q1=0;
		qrot.q2=1;
		qrot.q3=0;
	}

	else {
		s=sqrt(2.*(1.+angle));	
		invs=1./s;
		vector_cross(vec1,vec2,rotvec);
		qrot.q0=0.5*s;
		scalar_times(rotvec,invs,rotvec);
		qrot.q1=rotvec.x;
		qrot.q2=rotvec.y;
		qrot.q3=rotvec.z;
		scdivide_quat(qrot,sqrt(quat_inp(qrot,qrot)),qrot);
	}
				

	return qrot;
}


void drawPatch(int ind,int isite)

{

	int i;
	float numTheta,numPhi;
	float maxPhi;
    static float maxTheta;
    static int initial;
	//maxTheta=sys.delta;
	maxPhi=2.*PI;
	numTheta=Precision+10;
	numPhi=Precision+10;
	int countfan;

    if(initial==0) {
        initial=1;
        maxTheta = (180./PI)*acos(0.5+0.5*cos((PI/180.)*sys.delta));
    }

	float theta = 0, phi = 0, deltaTheta = maxTheta / numTheta, deltaPhi = maxPhi / numPhi;
	float z1, x1, y1, z2, x2, y2, z3, x3, y3, z4, x4, y4;
	float anglerot;
	quaternion qrot,qrot2;
	vector zvec;
	vector vec1,vec2,vec3,vec4;
	vector norm,vecnorm1,vecnorm2;
	float r = 0.5+0.02;
	float matrix[16];
	float sina;
	vector vecfan[100];

	glPushMatrix();

	glDisable(GL_CULL_FACE);
	glTranslatef(pos[ind].r.x, pos[ind].r.y, pos[ind].r.z);

	getdrawrotmatrix(pos[ind].q,matrix);
	glMultMatrixf(matrix);

	zvec.x=0;
	zvec.y=0;
	zvec.z=1;
	qrot=quatVecToVec(zvec,sys.site[isite]);

	getdrawrotmatrix(qrot,matrix);
	glMultMatrixf(matrix);

	glBegin(GL_QUADS);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialRed);

	countfan=0;
	for(theta = 0; theta <= maxTheta ; theta += deltaTheta){
		for(phi = 0; phi <= maxPhi; phi += deltaPhi){
			x1 = r * sin(theta + deltaTheta) * cos(phi + deltaPhi);
			y1 = r * sin(theta + deltaTheta) * sin(phi + deltaPhi);
			z1 = r * cos(theta + deltaTheta);

			x2 = r * sin(theta) * cos(phi + deltaPhi);
		 	y2 = r * sin(theta) * sin(phi + deltaPhi);
			z2 = r * cos(theta);

			x3 = r * sin(theta) * cos(phi);
			y3 = r * sin(theta) * sin(phi);
			z3 = r * cos(theta);

			x4 = r * sin(theta + deltaTheta) * cos(phi);
			y4 = r * sin(theta + deltaTheta) * sin(phi);
			z4 = r * cos(theta + deltaTheta);

			if(theta>(maxTheta-deltaTheta)) {
				vecfan[countfan].x=x1;
				vecfan[countfan].y=y1;
				vecfan[countfan].z=z1;
				countfan++;
			}

			vec4.x=x4;
			vec4.y=y4;
			vec4.z=z4;

			vec1.x=x1;
			vec1.y=y1;
			vec1.z=z1;

			vec3.x=x3;
			vec3.y=y3;
			vec3.z=z3;

			norm = findnorm(vec4,vec1,vec3);

			glNormal3f(norm.x, norm.y, norm.z);
			glVertex3f(x4, y4, z4);
			glVertex3f(x1, y1, z1);
			glVertex3f(x2, y2, z2);
			glVertex3f(x3, y3, z3);
		}
	}

	glEnd();

	glBegin(GL_TRIANGLE_FAN);
	glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialRed);
	glVertex3f(0,0,0);
	countfan=0;
	vec1.x=vec1.y=vec1.z=0;
	for(i=(int)numPhi-1; i>=0; i--) {
		if(countfan==2) {
				countfan=0;
		}
		if(countfan==0) {
				vec2=vecfan[i];
				vec3=vecfan[i-1];
				norm=findnorm(vec1,vec2,vec3);
				glNormal3f(norm.x, norm.y, norm.z);
		}
		countfan++;
		glVertex3f(vecfan[i].x,vecfan[i].y,vecfan[i].z);

	}
	//glVertex3f(vecfan[0].x,vecfan[0].y,vecfan[0].z);
	glVertex3f(vecfan[(int)numPhi-1].x,vecfan[(int)numPhi-1].y,vecfan[(int)numPhi-1].z);
	glVertex3f(0,0,0);
	glEnd();


	glPopMatrix();

//	glFlush();

	return;

}



void display_sphere(int ind)

{
 
   glPushMatrix();
   glTranslatef(pos[ind].r.x,pos[ind].r.y,pos[ind].r.z);
   glScalef(pos[ind].radius,pos[ind].radius,pos[ind].radius);
   //   switch(pos[ind].no) {
   switch (0) {
   case 0:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialCyan)  ; break;
   case 2:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialRed)   ; break;
   case 4:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialYellow); break;
   case 3:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialGreen); break;
   case 1:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialBlue); break;
   case 5:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialMagenta); break;
   case 6:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialSome); break;
   case 7:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialCyan)  ; break;
   case 8:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialRed)   ; break;
   case 9:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialYellow); break;
   case 10:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialGreen); break;
   case 11:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialBlue); break;
   case 12:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialMagenta); break;
   case 13:  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialSome); break;
   }

   glutSolidSphere(1, Precision, Precision);
   
   glPopMatrix(); 
}



void create_circle( void )
{
    vector	  *ptch;
    float alpha, delta ;
    int	  i, j;

    delta = 2.0 * PI / (16.*Precision);
    if (Circle) free(Circle);
    Circle = (vector *)malloc((16*Precision+1)*sizeof(vector));
    ptch = Circle;
    
    alpha = 0.0;
    for (i = 0; i < 16*Precision; i++) 
    {  
      ptch[i].x = 0;
      ptch[i].y = sin(alpha);
      ptch[i].z = cos(alpha);
      alpha += delta;
    } 
  
}



void convert(vector *p, float f[3])
{
    f[0]=p->x;
    f[1]=p->y;
    f[2]=p->z;
    return;
}



void lighted_plate( void )
{

  vector	*ptch;
  int	i, j, k;
  float   f[3];
  
  ptch = Circle;
  {
    glBegin(GL_POLYGON);
    for (k = 0; k < Precision*16; k++) 
      {
	glVertex3f(ptch->x,ptch->y,ptch->z);
	glNormal3f(ptch->x,ptch->y,ptch->z);
	ptch++;
      }
    glEnd();

    }
 
  ptch--;
  {
    glBegin(GL_POLYGON);

    for (k = 0; k < Precision*16; k++) 
      {
	glVertex3f(ptch->x,ptch->y,ptch->z);
	glNormal3f(ptch->x,ptch->y,ptch->z);
	ptch--;
      }
    glEnd();
    }
}



void display_plates(int ind)

{
 
   glPushMatrix();
   glTranslatef(pos[ind].r.x,pos[ind].r.y,pos[ind].r.z);
   glScalef(1,pos[ind].radius,pos[ind].radius);
   glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialGreen);
   lighted_plate();
   glPopMatrix(); 
}



void draw_box( float *vec ,float boxdis)
{
  int i;
  
  
  glPushMatrix();
  glScalef(0.5*vec[0],0.5*vec[1],0.5*vec[2]);
  glTranslatef(0,boxdis,0);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,MaterialRed);
  glColor3f(1.0,0,0);
  glBegin(GL_LINE_STRIP);
  for(i=0;i<17;i++)	
    {
      glVertex3fv(vert[i]);
    }
  
  glEnd();
  glPopMatrix();
  
}

static void Reshape( int width, int height )
{
   glViewport( 0, 0, width, height );
   glMatrixMode( GL_PROJECTION );
   glLoadIdentity();
   glFrustum( -1.0, 1.0, -1.0, 1.0, 5.0, 100.0 );
   glMatrixMode( GL_MODELVIEW );
   glLoadIdentity();
   glTranslatef( 0.0, 0.0, -15.0 );
}

void Display( void )
{
    float vec[3];
    int ipart,isite;

    glLoadIdentity();
    glTranslatef( 0.0, 0.0, Distance );  
    glRotatef(xRotation, 1, 0, 0);
    glRotatef(yRotation, 0, 1, 0);
    glRotatef(zRotation, 0, 0, 1);

    /*  glPolarView(Distance,Azimuth,Inclination,Twist);    */
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
    vec[0]=sys.boxl.x;
    vec[1]=sys.boxl.y;
    vec[2]=sys.boxl.z;
    draw_box( vec, 0 );

    for(ipart=0; ipart<sys.npart; ipart++) {
        display_sphere(ipart);
        for( isite=0; isite<sys.nsites; isite++) {
            drawPatch(ipart,isite);
        }
    }
   
    glFlush();

    glutSwapBuffers();
}


static void CursorKeys(int key, int x, int y)
{


    switch (key) {
      case GLUT_KEY_LEFT:
        zRotation += 5;
        break;
      case GLUT_KEY_RIGHT:
        zRotation -= 5;
        break;
      case GLUT_KEY_UP:
        xRotation += 5;
        break;
      case GLUT_KEY_DOWN:
        xRotation -= 5;
        break;
      default:
        return;
    }
    glutPostRedisplay();

}



static void Key(unsigned char key, int x, int y)
{

    switch (key) {
      case 27:
        exit(1);

      case ',':
        yRotation += 5;
        break;
      case '.':
        yRotation -= 5;
        break;
      case 'x':
        Distance += 2;
        break;
      case 'z':
        Distance -= 2;
        break;
      case '1':
	Precision =2;
	create_circle();
	break;
      case '2':
        Precision =4;
	create_circle();
	break;
      case '3':
	Precision =8;
	create_circle();
	break; 
      case '4':
	Precision =16;
	create_circle();
	break;  
      case '5':
	Precision =32;
	create_circle();
	break;
      case 'l':
	Drawall++;
	if (Drawall>3) Drawall =0;
	break;
      case 'm':
	glstat.Animate = !glstat.Animate;
	break;
      case 'b':
	glstat.BackwardAnimation = !glstat.BackwardAnimation;
	break;
      case '=':
        glstat.OneUp = 1;
	break;
      case '-':	
	glstat.OneDown =1;
	break;
     default:
       return;
    } 
    glutPostRedisplay();

}


static void Init( void )
{
   static GLfloat mat[4] = { 0.8, 0.8, 0.0, 1.0 };
   static GLfloat pos[4] = { -1.0, 1.0, 1.0, 0.0 };

   /* setup lighting, etc */
   glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
   glLightfv(GL_LIGHT0, GL_POSITION, position0);
   glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
   glLightfv(GL_LIGHT1, GL_POSITION, position1);
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
   glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHT1);

   glEnable(GL_BLEND);
   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_shininess);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_specular);

   /*   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
	glLightfv(GL_LIGHT0, GL_POSITION, pos);*/

   glEnable(GL_CULL_FACE);

   glDisable(GL_RESCALE_NORMAL_EXT);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_NORMALIZE);
   glHint(GL_FOG_HINT, GL_FASTEST);
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
   glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

   create_circle();

}


#define UNSCALED  1
#define NORMALIZE 2
#define RESCALE   3
#define QUIT      4



static void ModeMenu(int entry)
{
   if (entry==UNSCALED) {
      glDisable(GL_RESCALE_NORMAL_EXT);
      glDisable(GL_NORMALIZE);
   }
   else if (entry==NORMALIZE) {
      glEnable(GL_NORMALIZE);
      glDisable(GL_RESCALE_NORMAL_EXT);
   }
   else if (entry==RESCALE) {
      glDisable(GL_NORMALIZE);
      glEnable(GL_RESCALE_NORMAL_EXT);
   }
   else if (entry==QUIT) {
      exit(0);
   }
   glutPostRedisplay();
}



static void Idle() {

    if  (glstat.Animate) {
        glstat.OneUp = glstat.Animate;
    }

    if( glstat.OneUp) {
        mainloop_for_graphics();
        currentslice = &slice[0];
        convert_coordinates(currentslice);
        glstat.OneUp =0;
        glutPostRedisplay();
    }

    return;

}



void Init_Graphics(int argc, char *argv[])
{

  printf("In graphics\n");


  glutInit( &argc, argv );
  glutInitWindowPosition(0,0);
  glutInitWindowSize(640,800);

  glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB );

  glutCreateWindow(argv[0]);
  glClearDepth(1.0);
  glClearColor( 1.0, 1.0, 1.0, 1.0 );
  glColor3f( 1.0, 1.0, 1.0 );

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glFlush();
  glutSwapBuffers();

  Init();

  pos = (Particle *) calloc( sys.npart , sizeof(Particle) );

  currentslice = &slice[0];
  convert_coordinates(currentslice);

  glutIdleFunc( Idle );
  glutReshapeFunc( Reshape );
  glutSpecialFunc( CursorKeys );
  glutKeyboardFunc(Key);
  glutDisplayFunc( Display );
  
  glutCreateMenu(ModeMenu);
  glutAddMenuEntry("Unscaled", UNSCALED);
  glutAddMenuEntry("Normalize", NORMALIZE);
  glutAddMenuEntry("Rescale EXT", RESCALE);
  glutAddMenuEntry("Quit", QUIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
   
  glutMainLoop();
  return ;
}

