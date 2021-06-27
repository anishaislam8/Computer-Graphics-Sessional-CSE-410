#include <GL/glut.h>
#define pi (2*acos(0.0))
#include <vector>
#include<math.h>


struct point
{
	double x,y,z;
};

class Object {
public:
    struct point reference_point; // should have x, y, z
    double height, width, length;
    double color[3];
    double coEfficients[4]; // reflection coefficients
    int shine; // exponent term of specular component
    Object(){};
    virtual void draw(){};
    void setColor(double newColor[]){
        for(int i = 0; i < 3; i++){
            color[i] = newColor[i];
        }
    };
    void setShine(int newShine){
        shine = newShine;
    };
    void setCoEfficients(double newCoEfficients[]){
        for(int i = 0; i < 4; i++){
            coEfficients[i] = newCoEfficients[i];
        }
    };
    double A,B,C,D,E,F,G,H,I,J;
    struct point trianglePoints[3];
};

class Triangle : public Object {
public:
    Triangle(struct point newTrianglePoints[3]){
        for(int i = 0; i < 3; i++){
            trianglePoints[i] = newTrianglePoints[i];
        }
    }
    void draw();
};

void Triangle::draw(){

    glBegin(GL_TRIANGLES);{
        glVertex3f(trianglePoints[0].x,trianglePoints[0].y,trianglePoints[0].z);
        glVertex3f(trianglePoints[1].x,trianglePoints[1].y,trianglePoints[1].z);
        glVertex3f(trianglePoints[2].x,trianglePoints[2].y,trianglePoints[2].z);
    }glEnd();

}

class General : public Object {
public:
    General(struct point ref_point, double length, double width, double height){
        this->length = length;
        this->height = height;
        this->width = width;
        reference_point= ref_point;
    }
    void draw();
};

void General::draw(){
}

class Sphere : public Object {
public:
    Sphere(struct point center, double radius){
        reference_point = center;
        length = radius;
    }
    void draw();
};

void Sphere::draw(){
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    struct point points[100][100];
	int i,j;
	double h,r;
	int slices = 50;
	int stacks = 25;
	double radius = length;


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

    glTranslatef(-reference_point.x, -reference_point.y, -reference_point.z);
}



class Light{
public:
    struct point light_pos;
    double color[3];
    void draw();
    void setColor(double newColor[]){
        for(int i = 0; i < 3; i++){
            color[i] = newColor[i];
        }
    };
};

void Light::draw(){

   int a = 5;
   glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();

}


class Floor: public Object{
public:

    Floor(int floorWidth, int tileWidth){

        reference_point={-floorWidth/2,-floorWidth/2,0};
        length=tileWidth;
    }
    void draw();
};

void Floor::draw(){
    // write codes for drawing a checkerboard-like
    // floor with alternate colors on adjacent tiles
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    for(int j = 0; j < abs(reference_point.y) * 2; j+=length){
        for(int i = 0; i < abs(reference_point.x) * 2; i+=length){
            if((i/20 + j/20) % 2 == 0){
                glColor3f(0,0,0);
            } else {
                glColor3f(1,1,1);
            }

            glBegin(GL_QUADS);{

                glVertex3f(i,j,0);
                glVertex3f(i + length,j,0);
                glVertex3f(i + length,j + length,0);
                glVertex3f(i,j + length,0);


            }glEnd();
        }
    }





	glTranslatef(-reference_point.x, -reference_point.y, -reference_point.z);
}

