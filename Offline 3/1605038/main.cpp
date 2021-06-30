#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#define pi (2*acos(0.0))
#include <vector>
#include "1605038_classes.h"
#include "bitmap_image.hpp"


using namespace std;


double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double rotationAngle;
double viewAngle = 80;
double windowHeight= 500.0;
double windowWidth = 500.0;



int pixels;
int numberOfObjects;
int numberOfLightSources;

ifstream inputFile;






struct point pos;
struct point u;
struct point r;
struct point l;



double degreeToradian(double r){
    return (r * pi)/180;
}



void capture(){ // hopefully done :3

    bitmap_image image(imageWidth,imageHeight);

    for(int i=0;i<imageWidth;i++){
        for(int j=0;j<imageHeight;j++){
            image.set_pixel(i,j,0,0,0); // 2nd ta row increment, amar row to j
        }
    }

    double planeDistance = (windowHeight/2) / (tan(degreeToradian(viewAngle/2)));

    struct point topleft = {
        pos.x + l.x*planeDistance - r.x*(windowWidth/2) + u.x * (windowHeight/2),
        pos.y + l.y*planeDistance - r.y*(windowWidth/2) + u.y * (windowHeight/2),
        pos.z + l.z*planeDistance - r.z*(windowWidth/2) + u.z * (windowHeight/2)
    };

    double du = (1.0 * windowWidth)/ (1.0 *imageWidth);

    double dv = (1.0 * windowHeight)/ (1.0 * imageHeight);

    topleft = {
        topleft.x + r.x * (0.5 * du) - u.x * (0.5 * dv),
        topleft.y + r.y * (0.5 * du) - u.y * (0.5 * dv),
        topleft.z + r.z * (0.5 * du) - u.z * (0.5 * dv)
    };

    double t;
    double tmin;
    for (int i = 0; i< imageWidth; i++ ){
        for(int j = 0; j < imageHeight; j++){
            // calculate curPixel using topleft,r,u,i,j,du,dv
            tmin = INT_MAX;
            struct point currentPixel = {
                topleft.x + r.x * i * du - u.x * j * dv,
                topleft.y + r.y * i * du - u.y * j * dv,
                topleft.z + r.z * i * du - u.z * j * dv
            };

            // cast ray from eye to (curPixel-eye) direction
            struct point direction = {
                currentPixel.x - pos.x,
                currentPixel.y - pos.y,
                currentPixel.z - pos.z,
            };

            double valueOfDirection = valueOfAVector(direction);

            direction ={
                direction.x/valueOfDirection,
                direction.y/valueOfDirection,
                direction.z/valueOfDirection
            };

            Ray * r;
            r = new Ray(pos, direction);

            double *color = new double[3];
            // for each object, o in objects
            Object *nearestObject = NULL;

            for(int k = 0; k < objects.size(); k++){
                t = objects[k]->intersect(r, color, 0);
                //update t so that it stores min +ve value
                //save the nearest object, on

                if(t >= 0 && t < tmin){
                    nearestObject = objects[k];
                    tmin = t;
                }


            }
            if(nearestObject){

                tmin = nearestObject->intersect(r, color, 1);
                //update image pixel (i,j)
                //image.set_pixel(i,j,nearestObject->color[0] * 255,nearestObject->color[1] * 255,nearestObject->color[2] * 255);
                image.set_pixel(i,j,color[0] * 255,color[1] * 255,color[2] * 255);
            }



        }
    }
    image.save_image("D:\\Academics\\4-1\\Lab\\Computer-Graphics-Sessional-CSE-410\\Offline 3\\1605038\\out.bmp");
    cout << "Here !!!" << endl;
    image.clear();
}





void drawAxes()
{
	if(drawaxes==1)
	{

		glBegin(GL_LINES);{
		    glColor3f(1,0,0);
			glVertex3f( 300,0,0);
			glVertex3f(-300,0,0);

            glColor3f(0,1,0);
			glVertex3f(0,-300,0);
			glVertex3f(0, 300,0);

			glColor3f(0,0,1);
			glVertex3f(0,0, 300);
			glVertex3f(0,0,-300);
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


void keyboardListener(unsigned char key, int x,int y){
    double cosTheta = cos(degreeToradian(rotationAngle));
    double sinTheta = sin(degreeToradian(rotationAngle));

	switch(key){
        case '0':{
            capture();
            break;
        }
		case '1':
			{
			    //look left
                //rotate l and r with respect to u by 5 degrees

                //rotate l
                //lrot = lcos(rotationAngle) + (uxl)sin(rotationAngle)
                point uCrossL = crossProduct(u,l);


                point lRot;
                lRot.x = l.x * cosTheta + uCrossL.x * sinTheta;
                lRot.y = l.y * cosTheta + uCrossL.y * sinTheta;
                lRot.z = l.z * cosTheta + uCrossL.z * sinTheta;

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate r
                //rRot = rcos(rotationAngle) + (uxr)sin(rotationAngle)
                point uCrossR = crossProduct(u,r);

                point rRot;
                rRot.x = r.x * cosTheta + uCrossR.x * sinTheta;
                rRot.y = r.y * cosTheta + uCrossR.y * sinTheta;
                rRot.z = r.z * cosTheta + uCrossR.z * sinTheta;

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
                lRot.x = l.x * cosTheta - uCrossL.x * sinTheta;
                lRot.y = l.y * cosTheta - uCrossL.y * sinTheta;
                lRot.z = l.z * cosTheta - uCrossL.z * sinTheta;

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate r
                //rRot = rcos(rotationAngle) + (uxr)sin(rotationAngle)
                point uCrossR = crossProduct(u,r);

                point rRot;
                rRot.x = r.x * cosTheta - uCrossR.x * sinTheta;
                rRot.y = r.y * cosTheta - uCrossR.y * sinTheta;
                rRot.z = r.z * cosTheta - uCrossR.z * sinTheta;

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
                lRot.x = l.x * cosTheta + rCrossL.x * sinTheta;
                lRot.y = l.y * cosTheta + rCrossL.y * sinTheta;
                lRot.z = l.z * cosTheta + rCrossL.z * sinTheta;

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate u
                //uRot = ucos(rotationAngle) + (rxu)sin(rotationAngle)
                point rCrossU = crossProduct(r,u);

                point uRot;
                uRot.x = u.x * cosTheta + rCrossU.x * sinTheta;
                uRot.y = u.y * cosTheta + rCrossU.y * sinTheta;
                uRot.z = u.z * cosTheta + rCrossU.z * sinTheta;

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
                lRot.x = l.x * cosTheta - rCrossL.x * sinTheta;
                lRot.y = l.y * cosTheta - rCrossL.y * sinTheta;
                lRot.z = l.z * cosTheta - rCrossL.z * sinTheta;

                double val = valueOfAVector(lRot);
                lRot.x = lRot.x/val;
                lRot.y = lRot.y/val;
                lRot.z = lRot.z/val;

                l = lRot;

                //rotate u
                //uRot = ucos(rotationAngle) - (rxu)sin(rotationAngle)
                point rCrossU = crossProduct(r,u);

                point uRot;
                uRot.x = u.x * cosTheta - rCrossU.x * sinTheta;
                uRot.y = u.y * cosTheta - rCrossU.y * sinTheta;
                uRot.z = u.z * cosTheta - rCrossU.z * sinTheta;

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
                uRot.x = u.x * cosTheta - lCrossU.x * sinTheta;
                uRot.y = u.y * cosTheta - lCrossU.y * sinTheta;
                uRot.z = u.z * cosTheta - lCrossU.z * sinTheta;

                double val = valueOfAVector(uRot);
                uRot.x = uRot.x/val;
                uRot.y = uRot.y/val;
                uRot.z = uRot.z/val;

                u = uRot;

                //rotate r
                //rRot = rcos(rotationAngle) - (lxr)sin(rotationAngle)
                point lCrossR = crossProduct(l,r);

                point rRot;
                rRot.x = r.x * cosTheta - lCrossR.x * sinTheta;
                rRot.y = r.y * cosTheta - lCrossR.y * sinTheta;
                rRot.z = r.z * cosTheta - lCrossR.z * sinTheta;

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
                uRot.x = u.x * cosTheta + lCrossU.x * sinTheta;
                uRot.y = u.y * cosTheta + lCrossU.y * sinTheta;
                uRot.z = u.z * cosTheta + lCrossU.z * sinTheta;

                double val = valueOfAVector(uRot);
                uRot.x = uRot.x/val;
                uRot.y = uRot.y/val;
                uRot.z = uRot.z/val;

                u = uRot;

                //rotate r
                //rRot = rcos(rotationAngle) + (lxr)sin(rotationAngle)
                point lCrossR = crossProduct(l,r);

                point rRot;
                rRot.x = r.x * cosTheta + lCrossR.x * sinTheta;
                rRot.y = r.y * cosTheta + lCrossR.y * sinTheta;
                rRot.z = r.z * cosTheta + lCrossR.z * sinTheta;

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
			// ...
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

	for(int i = 0; i < objects.size(); i++){
        //cout << objects[i].color[0] << " " << objects[i].color[1] << " " << objects[i].color[2] << endl;
        glPushMatrix();
        glColor3f(objects[i]->color[0], objects[i]->color[1], objects[i]->color[2]);
        objects[i]->draw();
        glPopMatrix();

	}

	for(int i = 0; i < lights.size(); i++){
        //cout << lights[i].color[0] << " " << lights[i].color[1] << " " << lights[i].color[2] << endl;
        glPushMatrix();
        glTranslatef(lights[i].light_pos.x, lights[i].light_pos.y, lights[i].light_pos.z);
        glColor3f(lights[i].color[0], lights[i].color[1], lights[i].color[2]);
        glRotatef(90,1,0,0);
        lights[i].draw();
        glTranslatef(-lights[i].light_pos.x, -lights[i].light_pos.y, -lights[i].light_pos.z);
        glPopMatrix();
	}



	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void loadSphereData(){
    //cout << "called loadspheredata" << endl;
    struct point center;
    double radius;
    double color[3];
    double coefficients[4];
    int shine;

    double val;

    //center
    for(int i = 0; i < 3; i++){
        inputFile >> val;
        if(i == 0){
            center.x = val;
        } else if(i == 1){
            center.y = val;
        } else {
            center.z= val;
        }
    }



    //radius
    inputFile >> val;
    radius = val;

    //color
    for(int m = 0; m < 3;  m++){
        inputFile >> val;
        color[m] = val;
    }

    /*cout << "color" <<endl;
    for(int j = 0; j < 3; j++){
        cout << color[j] << " ";
    }
    cout << endl;*/

    //coeffiecients
    for(int m = 0; m < 4;  m++){
        inputFile >> val;
        coefficients[m] = val;
    }

    //shininess
    inputFile >> val;
    shine = (int) val;



    Object *temp;
    temp = new Sphere(center, radius); // received as input
    temp->setColor(color);
    temp->setCoEfficients(coefficients);
    temp->setShine(shine);
    objects.push_back(temp);


    //capture();
   //cout << "????" << endl;
    //cout << temp->reference_point.x << " " << temp->reference_point.y <<" " << temp->reference_point.z << endl;

}

void loadLightData(){
    // lightsources
    inputFile >> numberOfLightSources;

    double val;

    for(int i = 0; i < numberOfLightSources; i++){
        Light light;
        //cout << "Light : " << endl;

        inputFile >> val;
        light.light_pos.x = val;
        inputFile >> val;
        light.light_pos.y = val;
        inputFile >> val;
        light.light_pos.z = val;
        //cout << light.light_pos.x << " " << light.light_pos.y << " " << light.light_pos.z << endl;

        double color[3];
        inputFile >> val;
        color[0] = val;
        inputFile >> val;
        color[1]= val;
        inputFile >> val;
        color[2] = val;
        light.setColor(color);

        /*cout << "color" <<endl;
        for(int j = 0; j < 3; j++){
            cout << color[j] << " ";
        }
        cout << endl;*/

        lights.push_back(light);


    }

    //lights[0].draw();

}

void loadTriangleData(){

    struct point points[3];

    int shine;
    double val;
    int m = 0;
    double arr[3];
    double color[3];
    double coefficients[4];

    // point 1
    //cout << "Triangle : " << endl;

    for(int m = 0; m < 3;  m++){
        inputFile >> val;
        arr[m] = val;
    }
    points[0] = {arr[0], arr[1], arr[2]};
    /*cout << "point 0" << endl;

    for(int j = 0; j < 3; j++){
        cout << arr[j] << " ";
    }
    cout << endl;*/


     // point 2
    for(int m = 0; m < 3;  m++){
        inputFile >> val;
        arr[m] = val;
    }
    points[1] = {arr[0], arr[1], arr[2]};

    /*cout << "point 1" << endl;

    for(int j = 0; j < 3; j++){
        cout << arr[j] << " ";
    }
    cout << endl;*/
     // point 3

    for(int m = 0; m < 3;  m++){
        inputFile >> val;
        arr[m] = val;
    }
    points[2] = {arr[0], arr[1], arr[2]};
    /*cout << "point 3" << endl;

    for(int j = 0; j < 3; j++){
        cout << arr[j] << " ";
    }
    cout << endl;*/


    //color

    for(int m = 0; m < 3;  m++){
        inputFile >> val;
        color[m] = val;
    }
    /*cout << "color" << endl;

    for(int j = 0; j < 3; j++){
        cout << color[j] << " ";
    }
    cout << endl;*/
    //coeffiecients
    for(int m = 0; m < 4;  m++){
        inputFile >> val;
        coefficients[m] = val;
    }

    /*cout << "coeff" << endl;

    for(int j = 0; j < 4; j++){
        cout << coefficients[j] << " ";
    }
    cout << endl;*/

    //shininess
    inputFile >> val;
    shine = (int) val;

    //cout << "Triangle shine : " << shine << endl;


    Object *temp;
    temp = new Triangle(points); // received as input
    temp->setColor(color);
    temp->setCoEfficients(coefficients);
    temp->setShine(shine);
    objects.push_back(temp);
    //cout << temp->reference_point.x << " " << temp->reference_point.y << " " << temp->reference_point.z << endl;
}

void loadGeneralShapeData(){


    double A,B,C,D,E,F,G,H,I,J;
    struct point ref_point;
    double length, width, height;
    int shine;
    double arr[4];
    double aToJ[10];
    double val;
    double color[3];
    double coefficients[4];


    // A - J
    for(int m= 0; m < 10; m++) {
        inputFile >> val;
        aToJ[m] = val;
    }
    A = aToJ[0]; B = aToJ[1]; C = aToJ[2]; D = aToJ[3]; E = aToJ[4]; F = aToJ[5]; G = aToJ[6]; H = aToJ[7]; I = aToJ[8]; J = aToJ[9];


    /*cout << "a-j" <<endl;
    for(int j = 0; j < 10; j++){
        cout << aToJ[j] << " ";
    }
    cout << endl;*/


    //reference Point + length, width, height
    double arr1[6];
    for(int m= 0; m < 6; m++) {
        inputFile >> val;
        arr1[m] = val;
    }

    ref_point = {arr1[0],arr1[1],arr1[2]};
    length = arr1[3];
    width = arr1[4];
    height = arr1[5];

    /*cout << "ref point length width height" <<endl;
    for(int j = 0; j < 6; j++){
        cout << arr1[j] << " ";
    }
    cout << endl;*/

    for(int m = 0; m < 3;  m++){
        inputFile >> val;
        color[m] = val;
    }
    /*cout << "color" <<endl;
    for(int j = 0; j < 3; j++){
        cout << color[j] << " ";
    }
    cout << endl;*/

    //coeffiecients
    for(int m = 0; m < 4; m++){
        inputFile >> val;
        coefficients[m] = val;
    }
    /*cout << "coeff" <<endl;
    for(int j = 0; j < 4; j++){
        cout << coefficients[j] << " ";
    }
    cout << endl;*/

    //shininess
    inputFile >> val;
    shine = (int) val;

    //cout<< "shine : " << shine << endl;

    Object *temp;
    temp = new General(ref_point, length, width, height); // received as input
    temp->A = A; temp->B = B; temp->C = C; temp->D = D; temp->E = E; temp->F = F; temp->G = G; temp->H = H; temp->I = I; temp->J = J;
    temp->setColor(color);
    temp->setCoEfficients(coefficients);
    temp->setShine(shine);
    objects.push_back(temp);
}

void loadData(){
    inputFile >> recursionLevel;
    inputFile >> pixels;
    inputFile >> numberOfObjects;
    imageHeight = (int) pixels;
    imageWidth = (int) pixels;



    string objectType;
    for(int i = 0;i < numberOfObjects; i++){
        inputFile >> objectType;

        if(objectType == "sphere"){
            loadSphereData();
        }
        else if(objectType == "triangle"){
            loadTriangleData();
        }
        else if(objectType == "general"){
            loadGeneralShapeData();
        }
    }

    loadLightData();

    //floor
    Object *temp;
    temp = new Floor(1000,20);
    objects.push_back(temp);

    //cout << objects.size() << endl;

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
	gluPerspective(viewAngle,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance


}

int main(int argc, char **argv){

    inputFile.open("D:\\Academics\\4-1\\Lab\\Computer-Graphics-Sessional-CSE-410\\Offline 3\\1605038\\scene.txt");
	if(!inputFile){
        cout << "Error Input" << endl;
        exit(1);
    }


	glutInit(&argc,argv);
	glutInitWindowSize(windowHeight, windowWidth);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");


	init();
	loadData();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL


	inputFile.close();
	return 0;
}
