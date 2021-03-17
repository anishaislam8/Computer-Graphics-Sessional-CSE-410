#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <iostream>
#include <stdlib.h>
#include <math.h>
#define pi (2*acos(0.0))

using namespace std;


double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double rotationAngle;

struct point
{
	double x,y,z;
};

struct point pos;
struct point u;
struct point r;
struct point l;

double degreeToradian(){
    return (rotationAngle * pi)/180;
}

double valueOfAVector(point a){
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

point crossProduct(point a, point b){
    /*cx = aybz − azby
    cy = azbx − axbz
    cz = axby − aybx*/
    point c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}


void drawAxes()
{
	if(drawaxes==1)
	{

		glBegin(GL_LINES);{
		    glColor3f(1,0,0);
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

            glColor3f(0,1,0);
			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glColor3f(0,0,1);
			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawUpperHemisphere(double radius,int slices,int stacks)
{
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

		for(j=0;j<slices;j++)
		{
		    int color = (j%2 == 0 ? 0 : 1);
            glColor3f(color,color,color);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                /*glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);*/
			}glEnd();
		}
	}
}

void drawLowerHemisphere(double radius,int slices,int stacks)
{
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

		for(j=0;j<slices;j++)
		{
		    int color = (j%2 == 0 ? 0 : 1);
            glColor3f(color,color,color);

			glBegin(GL_QUADS);{
			    //upper hemisphere
				/*glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);*/
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
		    int color = i%2 == 0 ? 0 : 1;
            glColor3f(color,color,color);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawCylinder(double radius,int slices,int stacks)
{
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

		for(j=0;j<slices;j++)
		{
		    int color = (j%2 == 0 ? 0 : 1);
            glColor3f(color,color,color);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                /*glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);*/
			}glEnd();
		}
	}
}

void drawOuterCircle(double radius,int slices,int stacks)
{
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		r = 2*radius - r;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

		for(j=0;j<slices;j++)
		{
		    int color = (j%2 == 0 ? 0 : 1);
            glColor3f(color,color,color);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                /*glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);*/
			}glEnd();
		}
	}
}
void drawSS()
{
    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			{
			    //look left
                //rotate l and r with respect to u by 5 degrees

                //rotate l
                //lrot = lcos(rotationAngle) + (uxl)sin(rotationAngle)
                point uCrossL = crossProduct(u,l);

                point lRot;
                lRot.x = l.x * cos(degreeToradian()) + uCrossL.x * sin(degreeToradian());
                lRot.y = l.y * cos(degreeToradian()) + uCrossL.y * sin(degreeToradian());
                lRot.z = l.z * cos(degreeToradian()) + uCrossL.z * sin(degreeToradian());

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate r
                //rRot = rcos(rotationAngle) + (uxr)sin(rotationAngle)
                point uCrossR = crossProduct(u,r);

                point rRot;
                rRot.x = r.x * cos(degreeToradian()) + uCrossR.x * sin(degreeToradian());
                rRot.y = r.y * cos(degreeToradian()) + uCrossR.y * sin(degreeToradian());
                rRot.z = r.z * cos(degreeToradian()) + uCrossR.z * sin(degreeToradian());

                val = valueOfAVector(rRot);
                rRot.x = rRot.x/val;
                rRot.y = rRot.y/val;
                rRot.z = rRot.z/val;

                r = rRot;

			}
			break;
        case '2':
            {
                //look right
                //rotate l and r with respect to u by 5 degrees clockwise

                //rotate l
                //lrot = lcos(rotationAngle) - (uxl)sin(rotationAngle)
                point uCrossL = crossProduct(u,l);

                point lRot;
                lRot.x = l.x * cos(degreeToradian()) - uCrossL.x * sin(degreeToradian());
                lRot.y = l.y * cos(degreeToradian()) - uCrossL.y * sin(degreeToradian());
                lRot.z = l.z * cos(degreeToradian()) - uCrossL.z * sin(degreeToradian());

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate r
                //rRot = rcos(rotationAngle) + (uxr)sin(rotationAngle)
                point uCrossR = crossProduct(u,r);

                point rRot;
                rRot.x = r.x * cos(degreeToradian()) - uCrossR.x * sin(degreeToradian());
                rRot.y = r.y * cos(degreeToradian()) - uCrossR.y * sin(degreeToradian());
                rRot.z = r.z * cos(degreeToradian()) - uCrossR.z * sin(degreeToradian());

                val = valueOfAVector(rRot);
                rRot.x = rRot.x/val;
                rRot.y = rRot.y/val;
                rRot.z = rRot.z/val;

                r = rRot;

            }
            break;
        case '3': {
                //look up
                //rotate l and u with respect to r by 5 degrees

                //rotate l
                //lrot = lcos(rotationAngle) + (rxl)sin(rotationAngle)
                point rCrossL = crossProduct(r,l);

                point lRot;
                lRot.x = l.x * cos(degreeToradian()) + rCrossL.x * sin(degreeToradian());
                lRot.y = l.y * cos(degreeToradian()) + rCrossL.y * sin(degreeToradian());
                lRot.z = l.z * cos(degreeToradian()) + rCrossL.z * sin(degreeToradian());

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate u
                //uRot = ucos(rotationAngle) + (rxu)sin(rotationAngle)
                point rCrossU = crossProduct(r,u);

                point uRot;
                uRot.x = u.x * cos(degreeToradian()) + rCrossU.x * sin(degreeToradian());
                uRot.y = u.y * cos(degreeToradian()) + rCrossU.y * sin(degreeToradian());
                uRot.z = u.z * cos(degreeToradian()) + rCrossU.z * sin(degreeToradian());

                val = valueOfAVector(uRot);
                uRot.x = uRot.x/val;
                uRot.y = uRot.y/val;
                uRot.z = uRot.z/val;

                u = uRot;

            }

            break;
        case '4':
            {
                //look down
                //rotate l and u with respect to r by 5 degrees clockwise

                //rotate l
                //lrot = lcos(rotationAngle) - (rxl)sin(rotationAngle)
                point rCrossL = crossProduct(r,l);

                point lRot;
                lRot.x = l.x * cos(degreeToradian()) - rCrossL.x * sin(degreeToradian());
                lRot.y = l.y * cos(degreeToradian()) - rCrossL.y * sin(degreeToradian());
                lRot.z = l.z * cos(degreeToradian()) - rCrossL.z * sin(degreeToradian());

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate u
                //uRot = ucos(rotationAngle) - (rxu)sin(rotationAngle)
                point rCrossU = crossProduct(r,u);

                point uRot;
                uRot.x = u.x * cos(degreeToradian()) - rCrossU.x * sin(degreeToradian());
                uRot.y = u.y * cos(degreeToradian()) - rCrossU.y * sin(degreeToradian());
                uRot.z = u.z * cos(degreeToradian()) - rCrossU.z * sin(degreeToradian());

                val = valueOfAVector(uRot);
                uRot.x = uRot.x/val;
                uRot.y = uRot.y/val;
                uRot.z = uRot.z/val;

                u = uRot;
            }
            break;
        case '5':
            {
                //tilt clockwise
                //rotate u and r with respect to l by 5 degrees

                //rotate u
                //uRot = ucos(rotationAngle) - (lxu)sin(rotationAngle)
                point lCrossU = crossProduct(l,u);

                point uRot;
                uRot.x = u.x * cos(degreeToradian()) - lCrossU.x * sin(degreeToradian());
                uRot.y = u.y * cos(degreeToradian()) - lCrossU.y * sin(degreeToradian());
                uRot.z = u.z * cos(degreeToradian()) - lCrossU.z * sin(degreeToradian());

                double val = valueOfAVector(uRot);
                uRot.x = uRot.x/val;
                uRot.y = uRot.y/val;
                uRot.z = uRot.z/val;

                u = uRot;

                //rotate r
                //rRot = rcos(rotationAngle) - (lxr)sin(rotationAngle)
                point lCrossR = crossProduct(l,r);

                point rRot;
                rRot.x = r.x * cos(degreeToradian()) - lCrossR.x * sin(degreeToradian());
                rRot.y = r.y * cos(degreeToradian()) - lCrossR.y * sin(degreeToradian());
                rRot.z = r.z * cos(degreeToradian()) - lCrossR.z * sin(degreeToradian());

                val = valueOfAVector(rRot);
                rRot.x = rRot.x/val;
                rRot.y = rRot.y/val;
                rRot.z = rRot.z/val;

                r = rRot;

            }
            break;
        case '6':
            {
                //tilt anticlockwise
                //rotate u and r with respect to l by 5 degrees

                //rotate u
                //uRot = ucos(rotationAngle) + (lxu)sin(rotationAngle)
                point lCrossU = crossProduct(l,u);

                point uRot;
                uRot.x = u.x * cos(degreeToradian()) + lCrossU.x * sin(degreeToradian());
                uRot.y = u.y * cos(degreeToradian()) + lCrossU.y * sin(degreeToradian());
                uRot.z = u.z * cos(degreeToradian()) + lCrossU.z * sin(degreeToradian());

                double val = valueOfAVector(uRot);
                uRot.x = uRot.x/val;
                uRot.y = uRot.y/val;
                uRot.z = uRot.z/val;

                u = uRot;

                //rotate r
                //rRot = rcos(rotationAngle) + (lxr)sin(rotationAngle)
                point lCrossR = crossProduct(l,r);

                point rRot;
                rRot.x = r.x * cos(degreeToradian()) + lCrossR.x * sin(degreeToradian());
                rRot.y = r.y * cos(degreeToradian()) + lCrossR.y * sin(degreeToradian());
                rRot.z = r.z * cos(degreeToradian()) + lCrossR.z * sin(degreeToradian());

                val = valueOfAVector(rRot);
                rRot.x = rRot.x/val;
                rRot.y = rRot.y/val;
                rRot.z = rRot.z/val;

                r = rRot;
            }
            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x -= 2*l.x;
			pos.y -= 2*l.y;
			pos.z -= 2*l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			pos.x += 2*l.x;
			pos.y += 2*l.y;
			pos.z += 2*l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x += 2*r.x;
			pos.y += 2*r.y;
			pos.z += 2*r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x -= 2*r.x;
			pos.y -= 2*r.y;
			pos.z -= 2*r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += 2*u.x;
			pos.y += 2*u.y;
			pos.z += 2*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x -= 2*u.x;
			pos.y -= 2*u.y;
			pos.z -= 2*u.z;
			break;

		/*case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;*/

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			//...
			break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);

	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);
    //
    glPushMatrix();
    {
        glRotatef(90,0,1,0);
        drawUpperHemisphere(12.5,50,25);
        drawLowerHemisphere(12.5,50,25);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(90,0,1,0);
        glTranslatef(0,0,-18.75);
        drawUpperHemisphere(6.25,50,25);
        glTranslatef(0,0,-6.25);
        drawCylinder(6.25,50,50);
        glTranslatef(0,0,-6.25);
        drawCylinder(6.25,50,50);
        glTranslatef(0,0,-6.25);
        drawCylinder(6.25,50,50);
        glTranslatef(0,0,-6.25);
        drawCylinder(6.25,50,50);
        glTranslatef(0,0,-6.25);
        drawCylinder(6.25,50,50);
        glTranslatef(0,0,-6.25);
        drawCylinder(6.25,50,50);

        glRotatef(180,0,1,0);
        drawOuterCircle(6.25,50,50);


    }
    glPopMatrix();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
    rotationAngle = 5.0;

    pos = {100,100,0};
    u = {0,0,1};
    r = {-1/sqrt(2),1/sqrt(2),0};
    l = {-1/sqrt(2),-1/sqrt(2),0};
	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
