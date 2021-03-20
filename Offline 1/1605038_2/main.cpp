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
double smallCircleRadius;
int displayFunctionCalled;
bool done;
double speed;
double maxSpeed;
double minSpeed;
double savedBubbleSpeed;



struct point
{
	double x,y,z;
};

struct point direction[5];
struct point position[5];
bool inside[5];
bool created[5];



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
        {
            double temp = speed;
            speed = savedBubbleSpeed;
            savedBubbleSpeed = temp;
        }
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:
		    if(speed - 0.001 >= minSpeed)
            {
                //down arrow key
                speed -= 0.001;
            }
			break;
		case GLUT_KEY_UP:		// up arrow key
            if(speed + 0.001 <= maxSpeed)
            {
                //down arrow key
                speed += 0.001;
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
    if(displayFunctionCalled <= 10001){
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


    /*for (int i = 0; i < nDone; i++){
        drawCircle(smallCircleRadius,50, position[i].x, position[i].y);
    }*/


   if(displayFunctionCalled >= 10000){
        for(int i = 0; i < 5; i++){
            created[i] = true;
        }
        drawCircle(smallCircleRadius,50, position[0].x, position[0].y);
        drawCircle(smallCircleRadius,50, position[1].x, position[1].y);
        drawCircle(smallCircleRadius,50, position[2].x, position[2].y);
        drawCircle(smallCircleRadius,50, position[3].x, position[3].y);
        drawCircle(smallCircleRadius,50, position[4].x, position[4].y);

    }
    else if(displayFunctionCalled >= 7500){
        for(int i = 0; i < 4; i++){
            created[i] = true;
        }
        drawCircle(smallCircleRadius,50, position[0].x, position[0].y);
        drawCircle(smallCircleRadius,50, position[1].x, position[1].y);
        drawCircle(smallCircleRadius,50, position[2].x, position[2].y);
        drawCircle(smallCircleRadius,50, position[3].x, position[3].y);

    }
    else if(displayFunctionCalled >= 5000){
        for(int i = 0; i < 3; i++){
            created[i] = true;
        }
        drawCircle(smallCircleRadius,50, position[0].x, position[0].y);
        drawCircle(smallCircleRadius,50, position[1].x, position[1].y);
        drawCircle(smallCircleRadius,50, position[2].x, position[2].y);

    }
    else if(displayFunctionCalled >= 2500){
        for(int i = 0; i < 2; i++){
            created[i] = true;
        }
        drawCircle(smallCircleRadius,50, position[0].x, position[0].y);
        drawCircle(smallCircleRadius,50, position[1].x, position[1].y);

    }

    else if(displayFunctionCalled >= 100){
        for(int i = 0; i < 1; i++){
            created[i] = true;
        }
        drawCircle(smallCircleRadius,50, position[0].x, position[0].y);

    }




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

bool isInsideTheCircle(point a)
{
    double c1c2 = sqrt((a.x - 0) * (a.x - 0) + (a.y - 0) * (a.y - 0) + (a.z - 0) * (a.z - 0));
    if(c1c2 + smallCircleRadius <= centerCircleRadius){
        return true;
    } else {
        return false;
    }
}

point normalVector(point a, point b){
    point normal = {a.x - b.x, a.y - b.y, a.z - b.z};
    double valueOfNormal = sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
    normal.x = normal.x / valueOfNormal;
    normal.y = normal.y / valueOfNormal;
    normal.z = normal.z / valueOfNormal;
    return normal;
}

double dotProduct(point a, point b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

void animate(){



    for(int i = 0; i < 5; i++){
        if(created[i]){

            //is within the circle?
            bool ans = isInsideTheCircle(position[i]);
            if( ans && !inside[i]){
                inside[i] = true;

                //overlap
                for(int k = 0; k < 5; k++){

                    if( k!= i && inside[k]){
                        double dist = sqrt(pow(position[i].x - position[k].x, 2) + pow(position[i].y - position[k].y, 2) + pow(position[i].z - position[k].z, 2) );
                        if(dist < 2*smallCircleRadius){

                            inside[i] = false;

                        }
                    }
                }

            }
            if(inside[i]){

                double c1c2 = sqrt(position[i].x * position[i].x + position[i].y * position[i].y + position[i].z * position[i].z);

                //reflection with circle
                if(c1c2 >= centerCircleRadius-smallCircleRadius){

                    //reflection
                    //step 1 : normal vector : c2 - c1, then normalize
                    point c1 = {0,0,0};
                    point normal = normalVector(position[i], c1);


                    //step 2 : r = d - 2(d.n)n
                    double val = 2 * dotProduct(direction[i], normal);
                    direction[i] = {direction[i].x - normal.x * val, direction[i].y - normal.y * val, direction[i].z - normal.z * val };

                    while(sqrt(position[i].x * position[i].x + position[i].y * position[i].y + position[i].z * position[i].z) >= centerCircleRadius-smallCircleRadius ){
                        position[i].x += (speed * 1.0 * direction[i].x);
                        position[i].y += (speed * 1.0 * direction[i].y);
                    }

                }



                //reflection with other bubbles
                for(int j = 0; j < 5; j++){
                    if (j != i) // other bubbles except this
                    {
                        if(inside[j] == true)
                        {
                            double c1c2 = sqrt(pow(position[i].x - position[j].x, 2) + pow(position[i].y - position[j].y, 2) + pow(position[i].z - position[j].z, 2) );

                            //reflection with bubble
                            if(c1c2 <= 2*smallCircleRadius){

                                //reflection
                                //step 1 : normal vector : c2 - c1, then normalize
                                point normalOne = normalVector(position[i], position[j]);

                                point normalTwo = normalVector(position[j], position[i]);



                                //step 2 : r = d - 2(d.n)n
                                double val = 2 * dotProduct(direction[i], normalOne);
                                direction[i] = {direction[i].x - normalOne.x * val, direction[i].y - normalOne.y * val, direction[i].z - normalOne.z * val };

                                val = 2 * dotProduct(direction[j], normalTwo);
                                direction[j] = {direction[j].x - normalTwo.x * val, direction[j].y - normalTwo.y * val, direction[j].z - normalTwo.z * val };

                                while(sqrt(pow(position[i].x - position[j].x, 2) + pow(position[i].y - position[j].y, 2) + pow(position[i].z - position[j].z, 2) ) <= 2 * smallCircleRadius){
                                    position[j].x += (speed * 1.0 * direction[j].x);
                                    position[j].y += (speed * 1.0 * direction[j].y);

                                }

                              }
                        }
                        else
                        {
                            //no reflection with bubbles outside the circle
                        }
                    }
                }
            } else {

                if(position[i].x + speed * direction[i].x >= boundaryLength || position[i].x + speed * direction[i].x <= -boundaryLength )
                {
                    direction[i].x *= -1;
                }
                else if (position[i].y + speed * direction[i].y >= boundaryLength || position[i].y + speed * direction[i].y <= -boundaryLength)
                {
                    direction[i].y *= -1;
                }



            }
            position[i].x += (speed * 1.0 * direction[i].x);
            position[i].y += (speed * 1.0 * direction[i].y);

        }

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
	speed = 0.001;
	maxSpeed = 0.05;
	minSpeed = 0.0005;
    savedBubbleSpeed = 0;


	displayFunctionCalled = 0;
	boundaryLength = 10;
	centerCircleRadius = 7;
	smallCircleRadius = 1;



	for(int i = 0; i < 5; i++){
        inside[i] = false;
        created[i] = false;
        position[i] = {-(boundaryLength), -(boundaryLength), 0};


        direction[i].x = (float)rand()/RAND_MAX + 0.001;
        direction[i].y= (float)rand()/RAND_MAX + 0.005;
        direction[i].z = 0;

        cout << direction[i].x << " " << direction[i].y << " " << direction[i].z << endl;


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
