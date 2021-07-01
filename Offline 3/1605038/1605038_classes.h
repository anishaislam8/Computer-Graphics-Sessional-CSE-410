#include <GL/glut.h>
#define pi (2*acos(0.0))
#include <vector>
#include<math.h>
#include <iostream>
using namespace std;

struct point
{
	double x,y,z;
};


int imageHeight, imageWidth;
int recursionLevel;

/***********************Ray********************/

double dotProduct(struct point a, struct point b){
    //cout << "dot : " << a.x*b.x + a.y*b.y + a.z*b.z << endl;
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

double valueOfAVector(point a){
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

void clipcolor(double *finalColor){
    if(finalColor[0] > 1) finalColor[0] = 1.0;
    if(finalColor[0] < 0) finalColor[0] = 0.0;
    if(finalColor[1] > 1) finalColor[1] = 1.0;
    if(finalColor[1] < 0) finalColor[1] = 0.0;
    if(finalColor[2] > 1) finalColor[2] = 1.0;
    if(finalColor[2] < 0) finalColor[2] = 0.0;
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

class Ray {
public:
    struct point start;
    struct point dir;

    Ray(struct point start, struct point dir){
        this->dir = dir; // already normalized
        this->start = start;
    }
};

/***********************Object********************/

class Object {
public:
    struct point reference_point; // should have x, y, z
    double height, width, length;
    double color[3];
    double coEfficients[4]; // reflection coefficients
    int shine; // exponent term of specular component
    Object(){};
    virtual void draw(){};
    virtual double intersect(Ray *r, double *color, int level){
        return -1.0;
    }

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


/***********************Light********************/

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



vector <Object*> objects;
vector <Light> lights;

/***********************Triangle********************/

class Triangle : public Object {
public:
    Triangle(struct point newTrianglePoints[3]){
        for(int i = 0; i < 3; i++){
            trianglePoints[i] = newTrianglePoints[i];
        }
    }
    void draw();
    double intersect(Ray *r, double *finalColor, int level);
};

void Triangle::draw(){

    glBegin(GL_TRIANGLES);{
        glVertex3f(trianglePoints[0].x,trianglePoints[0].y,trianglePoints[0].z);
        glVertex3f(trianglePoints[1].x,trianglePoints[1].y,trianglePoints[1].z);
        glVertex3f(trianglePoints[2].x,trianglePoints[2].y,trianglePoints[2].z);
    }glEnd();

}

double Triangle::intersect(Ray *r, double *finalColor, int level){



    struct point intersectionPoint;
    const double EPSILON = 0.0000001;
    struct point vertex0 = trianglePoints[0];
    struct point vertex1 = trianglePoints[1];
    struct point vertex2 = trianglePoints[2];

    struct point edge1, edge2, h, s, q;
    double t,a,f,u,v;


    edge1 = {
        vertex1.x - vertex0.x,
        vertex1.y - vertex0.y,
        vertex1.z - vertex0.z
    };
    edge2= {
        vertex2.x - vertex0.x,
        vertex2.y - vertex0.y,
        vertex2.z - vertex0.z
    };
    h = crossProduct(r->dir, edge2);
    a = dotProduct(edge1,h);
    if (a > -EPSILON && a < EPSILON){
        t =  -1.0;    // This ray is parallel to this triangle.
        return t;
    }
    else {
        f = 1.0/a;
        s = {
            r->start.x - vertex0.x,
            r->start.y - vertex0.y,
            r->start.z - vertex0.z
        };
        u = f * dotProduct(s,h);
        if (u < 0.0 || u > 1.0){
            t =  -1.0;
            return t;
        }
        else {
            q = crossProduct(s,edge1);
            v = f * dotProduct(r->dir,q);
            if (v < 0.0 || u + v > 1.0){
                t = -1.0;
                return t;
            }
            else {
                // At this stage we can compute t to find out where the intersection point is on the line.
                t = f * dotProduct(edge2,q);
                if (t > EPSILON) // ray intersection
                {
                    intersectionPoint = {
                        r->start.x + r->dir.x * t,
                        r->start.y + r->dir.y * t,
                        r->start.z + r->dir.z * t,

                    };
                }
                else {// This means that there is a line intersection but not a ray intersection.
                    t = -1.0;
                    return t;
                }
            }


        }


    }


    if(level == 0){
        //cout << "triangle" << t << endl;
        return t;
    } else {

        double *intersectionPointColor;
        intersectionPointColor = new double[3];
        for(int i = 0; i < 3; i++){
            intersectionPointColor[i] = color[i];
        }


        for(int i = 0; i < 3; i++){
             finalColor[i] = intersectionPointColor[i] * coEfficients[0];
        }
        clipcolor(finalColor);
        //cout << "final color : " << finalColor[0] << " " << finalColor[1] << " " << finalColor[2] << endl;
        //calculate normal at intersection point
        struct point normal = crossProduct(edge1, edge2);

        double normalVectorValue = valueOfAVector(normal);

        normal = {
            normal.x/normalVectorValue,
            normal.y/normalVectorValue,
            normal.z/normalVectorValue
        };

         for(int m = 0; m < lights.size(); m++){
            double tNew;
            double tMinNew;

            struct point direction = {
                lights[m].light_pos.x - intersectionPoint.x,
                lights[m].light_pos.y - intersectionPoint.y,
                lights[m].light_pos.z - intersectionPoint.z,
            };

            double valueOfDirection = valueOfAVector(direction);

            direction ={
                direction.x/valueOfDirection,
                direction.y/valueOfDirection,
                direction.z/valueOfDirection
            };
            //cout << "direction : " << direction.x << " " << direction.y << " " << direction.z << endl;
            Ray * rNew;
            struct point startingPosition = {
                intersectionPoint.x + direction.x * 0.0000000001,
                intersectionPoint.y + direction.y * 0.0000000001,
                intersectionPoint.z + direction.z * 0.0000000001
            };
            rNew = new Ray(startingPosition, direction);


            tMinNew = INT_MAX;
            double *colorNew = new double[3];
            // for each object, o in objects
            for(int k = 0; k < objects.size(); k++){
                tNew = objects[k]->intersect(rNew, colorNew, 0);
                if(tNew >= 0 && tNew < tMinNew){
                    tMinNew = tNew;
                }


            }

            if(tMinNew >= t){
                //calculate R  // r=d−2(d⋅n)n
                struct point d = {
                    -direction.x,
                    -direction.y,
                    -direction.z
                };

                double dDotN = dotProduct(d,normal);
                struct point R = {
                    d.x - 2 * dDotN * normal.x,
                    d.y - 2 * dDotN * normal.y,
                    d.z - 2 * dDotN * normal.z,

                };

                double valueOfR = valueOfAVector(R);
                R = {
                    R.x/valueOfR,
                    R.y/valueOfR,
                    R.z/valueOfR
                };

                 //cout << "R : " << R.x << " " << R.y << " " << R.z << endl;
                // calculate V

                struct point V = {
                    -(r->dir.x) ,
                    -(r->dir.y) ,
                    -(r->dir.z)
                };

                double valueOfV = valueOfAVector(V);
                V = {
                    V.x/valueOfV,
                    V.y/valueOfV,
                    V.z/valueOfV
                };


                 //cout << "V : " << V.x << " " << V.y << " " << V.z << endl;

                // calculate diffuse component
                double lambertComponent = max(dotProduct(direction, normal),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[1] * lambertComponent * intersectionPointColor[0];
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[1] * lambertComponent * intersectionPointColor[1];
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[1] * lambertComponent * intersectionPointColor[2];

                clipcolor(finalColor);
                // calculate specular component
                double phongComponent = max(pow(dotProduct(R,V),shine),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[2] * phongComponent;
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[2] * phongComponent;
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[2] * phongComponent;

                clipcolor(finalColor);

            }


        }

        if(level >= recursionLevel){
           return t;
        }

        struct point dNew = r->dir;

        double dNewDotNormal = dotProduct(dNew,normal);
        struct point reflectedDirection = {
            dNew.x - 2 * dNewDotNormal * normal.x,
            dNew.y - 2 * dNewDotNormal * normal.y,
            dNew.z - 2 * dNewDotNormal * normal.z,

        };

        double valueOfReflected = valueOfAVector(reflectedDirection);
        reflectedDirection = {
            reflectedDirection.x/valueOfReflected ,
            reflectedDirection.y/valueOfReflected ,
            reflectedDirection.z/valueOfReflected
        };


        Ray * rReflected;
        struct point startingPositionReflected = {
            intersectionPoint.x + reflectedDirection.x * 0.0000000001,
            intersectionPoint.y + reflectedDirection.y * 0.0000000001,
            intersectionPoint.z + reflectedDirection.z * 0.0000000001
        };
        rReflected = new Ray(startingPositionReflected, reflectedDirection);
        double tMinReflected = INT_MAX;
        double tReflected;
        double *colorReflected = new double[3];
        // for each object, o in objects
        Object *nearestIntersecting= NULL;
        for(int obj = 0;  obj < objects.size(); obj++){
            tReflected = objects[obj]->intersect(rReflected, colorReflected, 0);
            if(tReflected >= 0 && tReflected < tMinReflected){
                nearestIntersecting = objects[obj];
                tMinReflected = tReflected;
            }


        }
        if(nearestIntersecting){
            tMinReflected = nearestIntersecting->intersect(rReflected, colorReflected,level+1);
            finalColor[0] = finalColor[0] + colorReflected[0] * coEfficients[3] ;
            finalColor[1] = finalColor[1] + colorReflected[1] * coEfficients[3] ;
            finalColor[2] = finalColor[2] + colorReflected[2] * coEfficients[3] ;

            clipcolor(finalColor);
        }

         //cout << "final color : " << finalColor[0] << " " << finalColor[1] << " " << finalColor[2] << endl;
    }

}


/***********************General Shapes********************/

class General : public Object {
public:
    General(struct point ref_point, double length, double width, double height){
        this->length = length;
        this->height = height;
        this->width = width;
        reference_point= ref_point;
    }
    void draw();
    bool valid(double t, Ray *r);
    double intersect(Ray *r, double *color, int level);
};

void General::draw(){
}

bool General::valid(double t, Ray *r){
    struct point intersectionPoint= {
        r->start.x + r->dir.x * t,
        r->start.y + r->dir.y * t,
        r->start.z + r->dir.z * t
    };
    if(abs(length) > 0){
        if(intersectionPoint.x < reference_point.x ||  intersectionPoint.x > reference_point.x + length){
            return false;
        }
    }
    if(abs(width) > 0){
        if(intersectionPoint.y < reference_point.y ||  intersectionPoint.y > reference_point.y + width){
            return false;
        }
    }
    if(abs(height) > 0){
        if(intersectionPoint.z < reference_point.z ||  intersectionPoint.z > reference_point.z + height){
            return false;
        }
    }

    return true;
}

double General::intersect(Ray *r, double *finalColor, int level){
    double t;

    //struct point R0 = r->start;

    struct point R0 = {r->start.x, r->start.y, r->start.z};

    double a = A*pow(r->dir.x,2) + B*pow(r->dir.y,2) + C*pow(r->dir.z,2);
    a += D * r->dir.x * r->dir.y + F * r->dir.y * r->dir.z + E * r->dir.x * r->dir.z;
    double b = 2*A*R0.x*r->dir.x + 2*B*R0.y*r->dir.y + 2*C*R0.z*r->dir.z;
    b += D*(R0.x*r->dir.y + R0.y*r->dir.x) + F*(R0.y*r->dir.z + R0.z*r->dir.y) + E*(R0.z*r->dir.x + R0.x*r->dir.z);
    b += G*r->dir.x + H*r->dir.y + I*r->dir.z;
    double c = A*pow(R0.x,2) + B*pow(R0.y,2) + C*pow(R0.z,2);
    c+= D*R0.x*R0.y + F*R0.y*R0.z + E*R0.x*R0.z;
    c+= G*R0.x + H*R0.y + I*R0.z + J;


    double discriminant = b*b - 4*a*c;
    if(discriminant < 0){
        t = -1.0;
        return t;
    } else {
        double t1 = (-b + sqrt(discriminant))/(2.0 * a);
        double t2 = (-b - sqrt(discriminant))/(2.0 * a);
        double temp;
        if(t1 < t2){
            temp = t1;
            t1 = t2;
            t2 = temp;
        }

        if(t1 < 0 && t2 < 0){
            t= -1.0;
            return t;
        } else if(t2 >= 0 && valid(t2,r)){
            t = t2;

        } else if(t1 >= 0 && valid(t1,r)){
            t = t1;

        } else {
            return -1.0;
        }
    }

    if (level == 0){
        return t;
    }
    else {
        struct point intersectionPoint= {
            r->start.x + r->dir.x * t,
            r->start.y + r->dir.y * t,
            r->start.z + r->dir.z * t
        };


        //cout << "Intersection Point : " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << endl;
        double *intersectionPointColor;
        intersectionPointColor = new double[3];
        for(int i = 0; i < 3; i++){
            intersectionPointColor[i] = color[i];
        }


        for(int i = 0; i < 3; i++){
             finalColor[i] = intersectionPointColor[i] * coEfficients[0];
        }
        clipcolor(finalColor);
        //cout << "final color : " << finalColor[0] << " " << finalColor[1] << " " << finalColor[2] << endl;
        //calculate normal at intersection point
        struct point normal = {
            2*A*intersectionPoint.x + D*intersectionPoint.y + E*intersectionPoint.z + G,
            2*B*intersectionPoint.y + D*intersectionPoint.x + F*intersectionPoint.z + H,
            2*C*intersectionPoint.z + F*intersectionPoint.y + E*intersectionPoint.x + I
        };

        double normalVectorValue = valueOfAVector(normal);

        normal = {
            normal.x/normalVectorValue,
            normal.y/normalVectorValue,
            normal.z/normalVectorValue
        };

         for(int m = 0; m < lights.size(); m++){
            double tNew;
            double tMinNew;

            struct point direction = {
                lights[m].light_pos.x - intersectionPoint.x,
                lights[m].light_pos.y - intersectionPoint.y,
                lights[m].light_pos.z - intersectionPoint.z,
            };

            double valueOfDirection = valueOfAVector(direction);

            direction ={
                direction.x/valueOfDirection,
                direction.y/valueOfDirection,
                direction.z/valueOfDirection
            };
            //cout << "direction : " << direction.x << " " << direction.y << " " << direction.z << endl;
            Ray * rNew;
            struct point startingPosition = {
                intersectionPoint.x + direction.x * 0.0000000001,
                intersectionPoint.y + direction.y * 0.0000000001,
                intersectionPoint.z + direction.z * 0.0000000001
            };
            rNew = new Ray(startingPosition, direction);


            tMinNew = INT_MAX;
            double *colorNew = new double[3];
            // for each object, o in objects
            for(int k = 0; k < objects.size(); k++){
                tNew = objects[k]->intersect(rNew, colorNew, 0);
                if(tNew >= 0 && tNew < tMinNew){
                    tMinNew = tNew;
                }


            }

            if(tMinNew >= t){
                //calculate R  // r=d−2(d⋅n)n
                struct point d = {
                    -direction.x,
                    -direction.y,
                    -direction.z
                };

                double dDotN = dotProduct(d,normal);
                struct point R = {
                    d.x - 2 * dDotN * normal.x,
                    d.y - 2 * dDotN * normal.y,
                    d.z - 2 * dDotN * normal.z,

                };

                double valueOfR = valueOfAVector(R);
                R = {
                    R.x/valueOfR,
                    R.y/valueOfR,
                    R.z/valueOfR
                };

                 //cout << "R : " << R.x << " " << R.y << " " << R.z << endl;
                // calculate V

                struct point V = {
                    -(r->dir.x) ,
                    -(r->dir.y) ,
                    -(r->dir.z)
                };

                double valueOfV = valueOfAVector(V);
                V = {
                    V.x/valueOfV,
                    V.y/valueOfV,
                    V.z/valueOfV
                };


                 //cout << "V : " << V.x << " " << V.y << " " << V.z << endl;

                // calculate diffuse component
                double lambertComponent = max(dotProduct(direction, normal),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[1] * lambertComponent * intersectionPointColor[0];
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[1] * lambertComponent * intersectionPointColor[1];
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[1] * lambertComponent * intersectionPointColor[2];

                clipcolor(finalColor);
                // calculate specular component
                double phongComponent = max(pow(dotProduct(R,V),shine),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[2] * phongComponent;
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[2] * phongComponent;
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[2] * phongComponent;

                clipcolor(finalColor);
            }


        }

        if(level >= recursionLevel){
           return t;
        }

        struct point dNew = r->dir;

        double dNewDotNormal = dotProduct(dNew,normal);
        struct point reflectedDirection = {
            dNew.x - 2 * dNewDotNormal * normal.x,
            dNew.y - 2 * dNewDotNormal * normal.y,
            dNew.z - 2 * dNewDotNormal * normal.z,

        };

        double valueOfReflected = valueOfAVector(reflectedDirection);
        reflectedDirection = {
            reflectedDirection.x/valueOfReflected ,
            reflectedDirection.y/valueOfReflected ,
            reflectedDirection.z/valueOfReflected
        };


        Ray * rReflected;
        struct point startingPositionReflected = {
            intersectionPoint.x + reflectedDirection.x * 0.0000000001,
            intersectionPoint.y + reflectedDirection.y * 0.0000000001,
            intersectionPoint.z + reflectedDirection.z * 0.0000000001
        };
        rReflected = new Ray(startingPositionReflected, reflectedDirection);
        double tMinReflected = INT_MAX;
        double tReflected;
        double *colorReflected = new double[3];
        // for each object, o in objects
        Object *nearestIntersecting= NULL;
        for(int obj = 0;  obj < objects.size(); obj++){
            tReflected = objects[obj]->intersect(rReflected, colorReflected, 0);
            if(tReflected >= 0 && tReflected < tMinReflected){
                nearestIntersecting = objects[obj];
                tMinReflected = tReflected;
            }


        }
        if(nearestIntersecting){
            tMinReflected = nearestIntersecting->intersect(rReflected, colorReflected,level+1);
            finalColor[0] = finalColor[0] + colorReflected[0] * coEfficients[3] ;
            finalColor[1] = finalColor[1] + colorReflected[1] * coEfficients[3] ;
            finalColor[2] = finalColor[2] + colorReflected[2] * coEfficients[3] ;

            clipcolor(finalColor);
        }



    }


}




/***********************Sphere********************/

class Sphere : public Object {
public:
    Sphere(struct point center, double radius){
        reference_point = center;
        length = radius;
    }
    void draw();
    double intersect(Ray *r, double *finalcolor, int level);
};


double Sphere::intersect(Ray *r, double *finalColor, int level){
    double t;
    struct point R0 = {
        r->start.x - reference_point.x,
        r->start.y - reference_point.y,
        r->start.z - reference_point.z

    };

    //struct point R0 = {r->start.x, r->start.y, r->start.z};

    double a = 1.0;
    double b = 2.0 * dotProduct(R0, r->dir);
    double c = dotProduct(R0, R0) - length * length;


    double discriminant = b*b - 4*a*c;
    if(discriminant < 0){
        t = -1.0;
        return t;
    } else {
        double t1 = (-b + sqrt(discriminant))/(2.0 * a);
        double t2 = (-b - sqrt(discriminant))/(2.0 * a);
        if(t1 < 0 && t2 < 0){
            t= -1.0;
            return t;
        } else if(t1 >= 0 && t2 >= 0){
           t = t2;
        } else if(t2 < 0){
            t = t1;
        }
    }

    if (level == 0){
        return t;
    }
    else {
        struct point intersectionPoint= {
            r->start.x + r->dir.x * t,
            r->start.y + r->dir.y * t,
            r->start.z + r->dir.z * t
        };
        //cout << "Intersection Point : " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << endl;
        double *intersectionPointColor;
        intersectionPointColor = new double[3];
        for(int i = 0; i < 3; i++){
            intersectionPointColor[i] = color[i];
        }


        for(int i = 0; i < 3; i++){
             finalColor[i] = intersectionPointColor[i] * coEfficients[0];
        }
        clipcolor(finalColor);
        //cout << "final color : " << finalColor[0] << " " << finalColor[1] << " " << finalColor[2] << endl;
        //calculate normal at intersection point
        struct point normal = {
            intersectionPoint.x - reference_point.x,
            intersectionPoint.y - reference_point.y,
            intersectionPoint.z - reference_point.z
        };

        double normalVectorValue = valueOfAVector(normal);

        normal = {
            normal.x/normalVectorValue,
            normal.y/normalVectorValue,
            normal.z/normalVectorValue
        };

         for(int m = 0; m < lights.size(); m++){
            double tNew;
            double tMinNew;

            struct point direction = {
                lights[m].light_pos.x - intersectionPoint.x,
                lights[m].light_pos.y - intersectionPoint.y,
                lights[m].light_pos.z - intersectionPoint.z,
            };

            double valueOfDirection = valueOfAVector(direction);

            direction ={
                direction.x/valueOfDirection,
                direction.y/valueOfDirection,
                direction.z/valueOfDirection
            };
            //cout << "direction : " << direction.x << " " << direction.y << " " << direction.z << endl;
            Ray * rNew;
            struct point startingPosition = {
                intersectionPoint.x + direction.x * 0.0000000001,
                intersectionPoint.y + direction.y * 0.0000000001,
                intersectionPoint.z + direction.z * 0.0000000001
            };
            rNew = new Ray(startingPosition, direction);


            tMinNew = INT_MAX;
            double *colorNew = new double[3];
            // for each object, o in objects
            for(int k = 0; k < objects.size(); k++){
                tNew = objects[k]->intersect(rNew, colorNew, 0);
                if(tNew >= 0 && tNew < tMinNew){
                    tMinNew = tNew;
                }


            }

            if(tMinNew >= t){
                //calculate R  // r=d−2(d⋅n)n
                struct point d = {
                    -direction.x,
                    -direction.y,
                    -direction.z
                };

                double dDotN = dotProduct(d,normal);
                struct point R = {
                    d.x - 2 * dDotN * normal.x,
                    d.y - 2 * dDotN * normal.y,
                    d.z - 2 * dDotN * normal.z,

                };

                double valueOfR = valueOfAVector(R);
                R = {
                    R.x/valueOfR,
                    R.y/valueOfR,
                    R.z/valueOfR
                };

                 //cout << "R : " << R.x << " " << R.y << " " << R.z << endl;
                // calculate V

                struct point V = {
                    -(r->dir.x) ,
                    -(r->dir.y) ,
                    -(r->dir.z)
                };

                double valueOfV = valueOfAVector(V);
                V = {
                    V.x/valueOfV,
                    V.y/valueOfV,
                    V.z/valueOfV
                };


                 //cout << "V : " << V.x << " " << V.y << " " << V.z << endl;

                // calculate diffuse component
                double lambertComponent = max(dotProduct(direction, normal),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[1] * lambertComponent * intersectionPointColor[0];
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[1] * lambertComponent * intersectionPointColor[1];
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[1] * lambertComponent * intersectionPointColor[2];

                clipcolor(finalColor);
                // calculate specular component
                double phongComponent = max(pow(dotProduct(R,V),shine),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[2] * phongComponent;
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[2] * phongComponent;
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[2] * phongComponent;

                clipcolor(finalColor);
            }


        }

        if(level >= recursionLevel){
           return t;
        }

        struct point dNew = r->dir;

        double dNewDotNormal = dotProduct(dNew,normal);
        struct point reflectedDirection = {
            dNew.x - 2 * dNewDotNormal * normal.x,
            dNew.y - 2 * dNewDotNormal * normal.y,
            dNew.z - 2 * dNewDotNormal * normal.z,

        };

        double valueOfReflected = valueOfAVector(reflectedDirection);
        reflectedDirection = {
            reflectedDirection.x/valueOfReflected ,
            reflectedDirection.y/valueOfReflected ,
            reflectedDirection.z/valueOfReflected
        };


        Ray * rReflected;
        struct point startingPositionReflected = {
            intersectionPoint.x + reflectedDirection.x * 0.0000000001,
            intersectionPoint.y + reflectedDirection.y * 0.0000000001,
            intersectionPoint.z + reflectedDirection.z * 0.0000000001
        };
        rReflected = new Ray(startingPositionReflected, reflectedDirection);
        double tMinReflected = INT_MAX;
        double tReflected;
        double *colorReflected = new double[3];
        // for each object, o in objects
        Object *nearestIntersecting= NULL;
        for(int obj = 0;  obj < objects.size(); obj++){
            tReflected = objects[obj]->intersect(rReflected, colorReflected, 0);
            if(tReflected >= 0 && tReflected < tMinReflected){
                nearestIntersecting = objects[obj];
                tMinReflected = tReflected;
            }


        }
        if(nearestIntersecting){
            tMinReflected = nearestIntersecting->intersect(rReflected, colorReflected,level+1);
            finalColor[0] = finalColor[0] + colorReflected[0] * coEfficients[3] ;
            finalColor[1] = finalColor[1] + colorReflected[1] * coEfficients[3] ;
            finalColor[2] = finalColor[2] + colorReflected[2] * coEfficients[3] ;

            clipcolor(finalColor);
        }


    }


}



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




/***********************Floor********************/

class Floor: public Object{
public:

    Floor(int floorWidth, int tileWidth){

        reference_point={-floorWidth/2,-floorWidth/2,0};
        length=tileWidth;
    }
    void draw();
    double intersect(Ray *r, double *color, int level);
};

void Floor::draw(){
    // write codes for drawing a checkerboard-like
    // floor with alternate colors on adjacent tiles
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    for(int j = 0; j < abs(reference_point.y) * 2; j+=length){
        for(int i = 0; i < abs(reference_point.x) * 2; i+=length){
            if((i/20 + j/20) % 2 == 0){
                glColor3f(1,1,1);
            } else {
                glColor3f(0,0,0);
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

double Floor::intersect(Ray *r, double *finalColor, int level){
    //cout << "floor intersect called" << endl;
    double t;
    struct point R0 = {
        r->start.x - reference_point.x,
        r->start.y - reference_point.y,
        r->start.z - reference_point.z

    };
    struct point normal = {
        0,
        0,
        1
    };
    double denominator = dotProduct(normal, r->dir);
    if (denominator != 0) {

        t = - dotProduct(R0, normal) / denominator;

    } else {
        t = -1.0;
        return t;
    }


    if (level == 0){
        return t;
    }
    else {
            //cout << "called" << endl;
        struct point intersectionPoint= {
            r->start.x + r->dir.x * t,
            r->start.y + r->dir.y * t,
            r->start.z + r->dir.z * t
        };

        struct point newPoint = {
            intersectionPoint.x - reference_point.x,
            intersectionPoint.y - reference_point.y,
            intersectionPoint.z - reference_point.z
        };

       int newI = newPoint.x;
       int newJ = newPoint.y;
        //cout << newI << " " << newJ << endl;
        if(newI >= 1000 || newJ >= 1000 || newI < 0 || newJ < 0){
            return -1.0;
        }
        //cout << "Intersection Point : " << intersectionPoint.x << " " << intersectionPoint.y << " " << intersectionPoint.z << endl;
        for(int j = 0; j < abs(reference_point.y) * 2; j+=length){
            for(int i = 0; i < abs(reference_point.x) * 2; i+=length){

                if((newI>=i && newI<i+length) && (newJ >= j && newJ < j+length)){
                    if((i/20 + j/20) % 2 == 0){
                        color[0] = 1;
                        color[1] = 1;
                        color[2] = 1;
                    } else {
                        color[0] = 0;
                        color[1] = 0;
                        color[2] = 0;
                    }

                    break;
                }

            }
        }
        double *intersectionPointColor;
        intersectionPointColor = new double[3];
        for(int i = 0; i < 3; i++){
            intersectionPointColor[i] = color[i];
        }


        for(int i = 0; i < 3; i++){
             finalColor[i] = intersectionPointColor[i] * coEfficients[0];
        }
        clipcolor(finalColor);
        //cout << "final color : " << finalColor[0] << " " << finalColor[1] << " " << finalColor[2] << endl;
        //calculate normal at intersection point

         for(int m = 0; m < lights.size(); m++){
            double tNew;
            double tMinNew;

            struct point direction = {
                lights[m].light_pos.x - intersectionPoint.x,
                lights[m].light_pos.y - intersectionPoint.y,
                lights[m].light_pos.z - intersectionPoint.z,
            };

            double valueOfDirection = valueOfAVector(direction);

            direction ={
                direction.x/valueOfDirection,
                direction.y/valueOfDirection,
                direction.z/valueOfDirection
            };
            //cout << "direction : " << direction.x << " " << direction.y << " " << direction.z << endl;
            Ray * rNew;
            struct point startingPosition = {
                intersectionPoint.x + direction.x * 0.0000000001,
                intersectionPoint.y + direction.y * 0.0000000001,
                intersectionPoint.z + direction.z * 0.0000000001
            };
            rNew = new Ray(startingPosition, direction);


            tMinNew = INT_MAX;
            double *colorNew = new double[3];
            // for each object, o in objects
            for(int k = 0; k < objects.size(); k++){
                tNew = objects[k]->intersect(rNew, colorNew, 0);
                if(tNew >= 0 && tNew < tMinNew){
                    tMinNew = tNew;
                }


            }

            if(tMinNew >= t){
                //calculate R  // r=d−2(d⋅n)n
                struct point d = {
                    -direction.x,
                    -direction.y,
                    -direction.z
                };

                double dDotN = dotProduct(d,normal);
                struct point R = {
                    d.x - 2 * dDotN * normal.x,
                    d.y - 2 * dDotN * normal.y,
                    d.z - 2 * dDotN * normal.z,

                };

                double valueOfR = valueOfAVector(R);
                R = {
                    R.x/valueOfR,
                    R.y/valueOfR,
                    R.z/valueOfR
                };

                 //cout << "R : " << R.x << " " << R.y << " " << R.z << endl;
                // calculate V

                struct point V = {
                    -(r->dir.x) ,
                    -(r->dir.y) ,
                    -(r->dir.z)
                };

                double valueOfV = valueOfAVector(V);
                V = {
                    V.x/valueOfV,
                    V.y/valueOfV,
                    V.z/valueOfV
                };


                 //cout << "V : " << V.x << " " << V.y << " " << V.z << endl;

                // calculate diffuse component
                double lambertComponent = max(dotProduct(direction, normal),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[1] * lambertComponent * intersectionPointColor[0];
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[1] * lambertComponent * intersectionPointColor[1];
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[1] * lambertComponent * intersectionPointColor[2];

                clipcolor(finalColor);
                // calculate specular component
                double phongComponent = max(pow(dotProduct(R,V),shine),0.0);
                finalColor[0] = finalColor[0] + lights[m].color[0] * coEfficients[2] * phongComponent;
                finalColor[1] = finalColor[1] + lights[m].color[1] * coEfficients[2] * phongComponent;
                finalColor[2] = finalColor[2] + lights[m].color[2] * coEfficients[2] * phongComponent;

                clipcolor(finalColor);
            }


        }

        if(level >= recursionLevel){
           return t;
        }

        struct point dNew = r->dir;

        double dNewDotNormal = dotProduct(dNew,normal);
        struct point reflectedDirection = {
            dNew.x - 2 * dNewDotNormal * normal.x,
            dNew.y - 2 * dNewDotNormal * normal.y,
            dNew.z - 2 * dNewDotNormal * normal.z,

        };

        double valueOfReflected = valueOfAVector(reflectedDirection);
        reflectedDirection = {
            reflectedDirection.x/valueOfReflected ,
            reflectedDirection.y/valueOfReflected ,
            reflectedDirection.z/valueOfReflected
        };


        Ray * rReflected;
        struct point startingPositionReflected = {
            intersectionPoint.x + reflectedDirection.x * 0.0000000001,
            intersectionPoint.y + reflectedDirection.y * 0.0000000001,
            intersectionPoint.z + reflectedDirection.z * 0.0000000001
        };
        rReflected = new Ray(startingPositionReflected, reflectedDirection);
        double tMinReflected = INT_MAX;
        double tReflected;
        double *colorReflected = new double[3];
        // for each object, o in objects
        Object *nearestIntersecting= NULL;
        for(int obj = 0;  obj < objects.size(); obj++){
            tReflected = objects[obj]->intersect(rReflected, colorReflected, 0);
            if(tReflected >= 0 && tReflected < tMinReflected){
                nearestIntersecting = objects[obj];
                tMinReflected = tReflected;
            }


        }
        if(nearestIntersecting){
            tMinReflected = nearestIntersecting->intersect(rReflected, colorReflected,level+1);
            finalColor[0] = finalColor[0] + colorReflected[0] * coEfficients[3] ;
            finalColor[1] = finalColor[1] + colorReflected[1] * coEfficients[3] ;
            finalColor[2] = finalColor[2] + colorReflected[2] * coEfficients[3] ;

            clipcolor(finalColor);
        }


    }


}




