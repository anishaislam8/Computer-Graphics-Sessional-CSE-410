#include<bits/stdc++.h>
using namespace std;

#define pi (2*acos(0.0))

struct Matrix{
    double matrix[4][4];
};

struct point{
    double x,y,z;
};

struct point eye,look,up;
struct point translate,scale,rotation;
struct point triangle1,triangle2,triangle3;
double fovY, aspectRatio, near, far, rotationAngle;

double identity[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

ifstream inputFile;
ofstream outputFile;

stack<Matrix> myStack;
stack<int> stackLength;



void readInitialInfo(){
    //gluLookAt
    inputFile >> eye.x;
    inputFile >> eye.y;
    inputFile >> eye.z;

    inputFile >> look.x;
    inputFile >> look.y;
    inputFile >> look.z;

    inputFile >> up.x;
    inputFile >> up.y;
    inputFile >> up.z;

    // gluPerspective
    inputFile >> fovY;
    inputFile >> aspectRatio;
    inputFile >> near;
    inputFile >> far;
}

void readTriangle(){

    inputFile >> triangle1.x;
    inputFile >> triangle1.y;
    inputFile >> triangle1.z;

    inputFile >> triangle2.x;
    inputFile >> triangle2.y;
    inputFile >> triangle2.z;

    inputFile >> triangle3.x;
    inputFile >> triangle3.y;
    inputFile >> triangle3.z;
}

void readTranslate(){
    inputFile >> translate.x;
    inputFile >> translate.y;
    inputFile >> translate.z;
}

void readRotate(){
    inputFile >> rotationAngle;
    inputFile >> rotation.x;
    inputFile >> rotation.y;
    inputFile >> rotation.z;
}

void readScale(){
    inputFile >> scale.x;
    inputFile >> scale.y;
    inputFile >> scale.z;
}

void handlePush(){
    stackLength.push(myStack.size());
}

void handlePop(){
    int upto = stackLength.top();
    stackLength.pop();
    while(myStack.size() != upto){
        myStack.pop();
    }
}

double degreeToradian(double r){
    return (r * pi)/180;
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

double dotProduct(point a, point b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

point normalize(point a){
    double val = valueOfAVector(a);
    a.x = a.x / val;
    a.y = a.y / val;
    a.z = a.z / val;
    return a;
}

void multiplyMatrices(double** matrix, int arrSize){

    double mult[4][arrSize];
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<arrSize;j++){
            mult[i][j] = 0.0;
        }
    }


    for(int i = 0; i < 4; i++){
        for(int j = 0; j < arrSize; j++){
            for(int k = 0; k < 4; k++)
            {
                mult[i][j] += myStack.top().matrix[i][k] * matrix[k][j];
            }
        }
    }

    if(arrSize == 1){
        outputFile << fixed << mult[0][0]/mult[3][0] << " " << mult[1][0]/mult[3][0]  << " " << mult[2][0]/mult[3][0] << endl;
    } else {
        Matrix mat;
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                mat.matrix[i][j] = mult[i][j];
            }
        }
        myStack.push(mat);
    }
}

void openFilesStage1(){
    inputFile.open("4/scene.txt");
    outputFile.open("output/stage1.txt");
    if(!inputFile){
        cout << "Error Input" << endl;
        exit(1);
    }
    if(!outputFile){
        cout << "Error Output" << endl;
        exit(1);
    }
}

void openFilesStage2(){
    inputFile.open("output/stage1.txt");
    outputFile.open("output/stage2.txt");
    if(!inputFile){
        cout << "Error Input" << endl;
        exit(1);
    }
    if(!outputFile){
        cout << "Error Output" << endl;
        exit(1);
    }
}

void openFilesStage3(){
    inputFile.open("output/stage2.txt");
    outputFile.open("output/stage3.txt");
    if(!inputFile){
        cout << "Error Input" << endl;
        exit(1);
    }
    if(!outputFile){
        cout << "Error Output" << endl;
        exit(1);
    }
}

void processPoint1(){

    double** matrix;
    matrix = new double*[4];
    matrix[0] = new double[1];
    matrix[0][0] = triangle1.x;
    matrix[1] = new double[1];
    matrix[1][0] = triangle1.y;
    matrix[2] = new double[1];
    matrix[2][0] = triangle1.z;
    matrix[3] = new double[1];
    matrix[3][0] = 1;

    multiplyMatrices(matrix, 1);
}

void processPoint2(){

    double** matrix;
    matrix = new double*[4];
    matrix[0] = new double[1];
    matrix[0][0] = triangle2.x;
    matrix[1] = new double[1];
    matrix[1][0] = triangle2.y;
    matrix[2] = new double[1];
    matrix[2][0] = triangle2.z;
    matrix[3] = new double[1];
    matrix[3][0] = 1;

    multiplyMatrices(matrix, 1);
}

void processPoint3(){

    double** matrix;
    matrix = new double*[4];
    matrix[0] = new double[1];
    matrix[0][0] = triangle3.x;
    matrix[1] = new double[1];
    matrix[1][0] = triangle3.y;
    matrix[2] = new double[1];
    matrix[2][0] = triangle3.z;
    matrix[3] = new double[1];
    matrix[3][0] = 1;

    multiplyMatrices(matrix , 1);
    outputFile << endl;

}

void processTranslate(){

    double **matrix;
    matrix = new double*[4];
    for(int i=0;i<4;i++)
    {
        matrix[i] = new double[4];
        for(int j=0;j<4;j++)
            matrix[i][j] = 0; //initialize the matrix cells to 0
    }
    matrix[0][0] = 1;
    matrix[1][1] = 1;
    matrix[2][2] = 1;
    matrix[3][3] = 1;
    matrix[0][3] = translate.x;
    matrix[1][3] = translate.y;
    matrix[2][3] = translate.z;

    multiplyMatrices(matrix,4);
}

void processScale(){

    double **matrix;
    matrix = new double*[4];
    for(int i=0;i<4;i++)
    {
        matrix[i] = new double[4];
        for(int j=0;j<4;j++)
            matrix[i][j] = 0; //initialize the matrix cells to 0
    }
    matrix[0][0] = scale.x;
    matrix[1][1] = scale.y;
    matrix[2][2] = scale.z;
    matrix[3][3] = 1;

    multiplyMatrices(matrix,4);
}

point rodrigues(point axis){
    // R(x) = cos(theta) * x + (1-cos(theta))*(a dot x)*a + sin(theta) * (a cross x)
    point rotationVector = normalize(rotation);

    double angle = degreeToradian(rotationAngle);
    double cosTheta = cos(angle);
    double sinTheta = sin(angle);

    point one = {cosTheta * axis.x , cosTheta * axis.y , cosTheta * axis.z };

    double multiplyBy = (1-cosTheta) * dotProduct(axis,rotationVector);
    point two = {multiplyBy * rotationVector.x, multiplyBy * rotationVector.y, multiplyBy * rotationVector.z};

    point aCrossx = crossProduct(rotationVector,axis);
    point three = {sinTheta * aCrossx.x, sinTheta * aCrossx.y, sinTheta * aCrossx.z};

    point resultVal = {one.x + two.x + three.x, one.y + two.y + three.y, one.z + two.z + three.z};
    return resultVal;

}

void processRotate(){

    point i = {1.0,0.0,0.0};
    point j = {0.0,1.0,0.0};
    point k = {0.0,0.0,1.0};

    point c1 = rodrigues(i);
    point c2 = rodrigues(j);
    point c3 = rodrigues(k);

    double **matrix;
    matrix = new double*[4];
    for(int i=0;i<4;i++)
    {
        matrix[i] = new double[4];
        for(int j=0;j<4;j++)
            matrix[i][j] = 0; //initialize the matrix cells to 0
    }
    matrix[0][0] = c1.x;
    matrix[1][0] = c1.y;
    matrix[2][0] = c1.z;
    matrix[3][0] = 0;

    matrix[0][1] = c2.x;
    matrix[1][1] = c2.y;
    matrix[2][1] = c2.z;
    matrix[3][1] = 0;

    matrix[0][2] = c3.x;
    matrix[1][2] = c3.y;
    matrix[2][2] = c3.z;
    matrix[3][2] = 0;

    matrix[0][3] = 0;
    matrix[1][3] = 0;
    matrix[2][3] = 0;
    matrix[3][3] = 1;

    multiplyMatrices(matrix,4);

}

void insertIdentity(){
    Matrix mat;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            mat.matrix[i][j] = identity[i][j];
        }
    }
    myStack.push(mat);

    stackLength.push(1);

}

void handleStage2(){
    point l = {look.x-eye.x, look.y-eye.y, look.z-eye.z};
    l = normalize(l);

    point r = crossProduct(l, up);
    r = normalize(r);

    point u = crossProduct(r,l);

    double T[4][4] = {{1.0,0.0,0.0,-eye.x},{0.0,1.0,0.0,-eye.y},{0.0,0.0,1.0,-eye.z},{0.0,0.0,0.0,1.0}};
    double R[4][4] = {{r.x,r.y,r.z,0.0},{u.x,u.y,u.z,0.0},{-l.x,-l.y,-l.z,0.0},{0.0,0.0,0.0,1.0}};
    double V[4][4];

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            V[i][j] = 0;
        }
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            for(int k = 0; k < 4; k++)
            {
                V[i][j] += R[i][k] * T[k][j];
            }
        }
    }

    int m = 0;
    int totalLine = 0;

    double triangleStage1[3];
    double val;

    while(inputFile >> val){

        triangleStage1[m] = val;
        m+=1;


        if(m == 3){
            m = 0;
            totalLine+=1;
            double trianglePoint[4][1]= {{triangleStage1[0]},{triangleStage1[1]},{triangleStage1[2]},{1}};
            double result[4][1];
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 1; j++){
                    result[i][j] = 0;
                }
            }
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 1; j++){
                    for(int k = 0; k < 4; k++)
                    {
                        result[i][j] += V[i][k] * trianglePoint[k][j];
                    }
                }
            }
            outputFile << fixed << result[0][0]/result[3][0] << " " << result[1][0]/result[3][0]  << " " << result[2][0]/result[3][0] << endl;
        }
        if(totalLine == 3){
            outputFile << endl;
            totalLine = 0;
        }

    }
}

void handleStage3(){
    double fovX = fovY * aspectRatio;
    double t = near * tan(degreeToradian(fovY/2));
    double r = near * tan(degreeToradian(fovX/2));

    double P[4][4] = {{near/r, 0.0,0.0,0.0}, {0.0, near/t, 0.0,0.0}, {0.0,0.0,-((far + near)/(far - near)), -((2*far*near)/(far-near))}, {0.0, 0.0, -1.0, 0.0}};

    int m = 0;
    int totalLine = 0;

    double triangleStage2[3];
    double val;

    while(inputFile >> val){

        triangleStage2[m] = val;
        m+=1;


        if(m == 3){
            m = 0;
            totalLine+=1;

            double trianglePoint2[4][1]= {{triangleStage2[0]},{triangleStage2[1]},{triangleStage2[2]},{1}};

            double result2[4][1];
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 1; j++){
                    result2[i][j] = 0.0;
                }
            }
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 1; j++){
                    for(int k = 0; k < 4; k++)
                    {
                        result2[i][j] += (P[i][k] * trianglePoint2[k][j]);
                    }
                }
            }
            outputFile << fixed << result2[0][0]/result2[3][0] << " " << result2[1][0]/result2[3][0]  << " " << result2[2][0]/result2[3][0] << endl;
        }
        if(totalLine == 3){
            outputFile << endl;
            totalLine = 0;
        }

    }

}

int main(){

    /*********************************Stage 1*************************/

    openFilesStage1();
    insertIdentity();
    readInitialInfo();

    string command;

    while (inputFile >> command) {

      if(command.compare("triangle") == 0){

        readTriangle();
        processPoint1();
        processPoint2();
        processPoint3();

      }
      else if(command == "translate"){
        readTranslate();
        processTranslate();
      }
      else if(command == "rotate"){
        readRotate();
        processRotate();
      }
      else if(command == "scale"){
        readScale();
        processScale();
      }
      else if(command == "push"){
        handlePush();
      }
      else if(command == "pop"){
        handlePop();
      }
      else if(command == "end"){
        break;
      }
    }

    inputFile.close();
    outputFile.close();

    /******************************************************************/



    /************************Stage2***********************************/
    openFilesStage2();
    handleStage2();



    inputFile.close();
    outputFile.close();

    /*****************************************************************/


    /************************Stage3***********************************/
    openFilesStage3();
    handleStage3();

    inputFile.close();
    outputFile.close();

    /*****************************************************************/


    return 0;
}
