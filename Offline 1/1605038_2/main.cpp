#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double boundaryLength;
double centerCircleRadius;
double bubbleVelocity;
double bubbleHighestVelocity;
double bubbleLowestVelocity;
bool bubbleMoving;
double saveBubbleVelocity;
int displayFunctionCalled;



struct point
{
	double x,y,z;
};

struct point bubbleDirection[5];
struct point position[5];


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

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
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}



void drawCircle(double radius,int segments, double translateX, double translateY)
{
    int i;
    struct point points[100];
    //glColor3f(0.7,0.7,0.7);
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
			glVertex3f(points[i].x + translateX,points[i].y + translateY,0);
			glVertex3f(points[i+1].x + translateX,points[i+1].y + translateY ,0);
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

void drawBoundary(double a)
{
    glBegin(GL_LINES);
    {
        glVertex3f( a, a,2);
        glVertex3f( a,-a,2);
        glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
    }glEnd();
     glBegin(GL_LINES);
    {
        glVertex3f( a, a,2);
        glVertex3f( -a,a,2);
        glVertex3f(-a,-a,2);
		glVertex3f(a, -a,2);
    }glEnd();
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
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
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

		case 'p':
			if(bubbleMoving){
                bubbleMoving = false;
                saveBubbleVelocity = bubbleVelocity;
                bubbleVelocity = 0;
			} else {
			    bubbleMoving = true;
                bubbleVelocity = saveBubbleVelocity;
			}
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			if (bubbleVelocity - 0.02 >= bubbleLowestVelocity){
                bubbleVelocity -= 0.02;
			}
			break;
		case GLUT_KEY_UP:		// up arrow key
			if (bubbleVelocity + 0.02 <= bubbleHighestVelocity){
                bubbleVelocity += 0.02;
			}
			break;

		/*case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
			break;

		case GLUT_KEY_INSERT:
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
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}


void display(){
    if(displayFunctionCalled > -1){
        displayFunctionCalled+=1;
    }


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
	gluLookAt(0,0,20,	0,0,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	//drawAxes();
	//drawGrid();


    glColor3f(0,1,0);
    drawBoundary(boundaryLength);

    //drawSS();

    glColor3f(1,0,0);
    drawCircle(centerCircleRadius, 50,0,0);


    // shoot bubbles :)
    /*if(displayFunctionCalled < 1500){
        drawCircle(1,50, position[0].x, position[0].y);
    }
    else if(displayFunctionCalled >= 1500 && displayFunctionCalled < 3000){
        drawCircle(1,50, position[1].x, position[1].y);
    }
    else if(displayFunctionCalled >= 3000 && displayFunctionCalled < 4500){
        drawCircle(1,50, position[2].x, position[2].y);
    }
    else if(displayFunctionCalled >= 4500 && displayFunctionCalled < 6000){
        drawCircle(1,50, position[3].x, position[3].y);
    }
    else if(displayFunctionCalled > 6000){
        drawCircle(1,50, position[4].x, position[4].y);
        displayFunctionCalled = -1;
    }*/

    for (int i = 0; i < 5; i++){
        drawCircle(1,50, position[i].x, position[i].y);
    }




    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){

	for(int i = 0; i < 5; i++){
        // reflection along x axis


        /*if((position[i].y  + (bubbleVelocity * bubbleDirection[i].y) - 1  <= -boundaryLength) || (position[i].y + (bubbleVelocity * bubbleDirection[i].y) + 1 >= boundaryLength) )
        {
            cout << "Inside y " << i << endl;
            cout << "Boundary " << boundaryLength << endl;
            cout << i << "  -> y: " << position[i].y << endl;
           bubbleDirection[i].y = bubbleDirection[i].y * -1;

        }
        // reflection along y axis
        else if((position[i].x + (bubbleVelocity * bubbleDirection[i].x) - 1 <= -boundaryLength) || (position[i].x + ((bubbleVelocity * bubbleDirection[i].x)) + 1 >= boundaryLength))
        {
            cout << "Inside x " << i << endl;
            cout << "Boundary " << boundaryLength << endl;
            cout << i << "  -> x: " << position[i].x << endl;
            bubbleDirection[i].x  = bubbleDirection[i].x * -1;
        }*/

        position[i].x += (bubbleVelocity * bubbleDirection[i].x);
        position[i].y += (bubbleVelocity * bubbleDirection[i].y);

	}


	glutPostRedisplay();
}

void init(){
	//codes for initialization
    srand (time(NULL) * 1000);
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	displayFunctionCalled = 0;

	bubbleVelocity = 0.1;
	bubbleHighestVelocity = 1.0;
	bubbleLowestVelocity = 0.01;
	bubbleMoving = true;
	boundaryLength = 10;
	centerCircleRadius = 7;

    double xVal = (rand() % 10 + 1) * 0.001;
    double yVal = (rand() % 10 + 1) * 0.002;


    position[0] = {-(boundaryLength), -(boundaryLength), 0};
    bubbleDirection[0] = {xVal, yVal,0};


    xVal = (rand() % 10 + 1) * 0.001;
    yVal = (rand() % 10 + 1) * 0.002;

    position[1] = {-(boundaryLength), -(boundaryLength), 0};
    bubbleDirection[1] = {xVal, yVal,0};

    xVal = (rand() % 10 + 1) * 0.001;
    yVal = (rand() % 10 + 1) * 0.002;

    position[2] = {-(boundaryLength), -(boundaryLength), 0};
    bubbleDirection[2] = {xVal, yVal,0};

    xVal = (rand() % 10 + 1) * 0.001;
    yVal = (rand() % 10 + 1) * 0.002;

    position[3] = {-(boundaryLength), -(boundaryLength), 0};
    bubbleDirection[3] = {xVal, yVal,0};

    xVal = (rand() % 10 + 1) * 0.001;
    yVal = (rand() % 10 + 1) * 0.002;

    position[4] = {-(boundaryLength), -(boundaryLength), 0};
    bubbleDirection[4] = {xVal, yVal,0};

    for(int i = 0 ; i < 5; i++){
        position[i].x += bubbleVelocity * bubbleDirection[i].x;
        position[i].y += bubbleVelocity * bubbleDirection[i].y;
    }


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
