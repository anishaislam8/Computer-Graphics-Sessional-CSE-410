#include<bits/stdc++.h>
using namespace std;

struct point{
    double x,y,z;
};
struct point eye;
struct point look;
struct point up;
double fovY, aspectRatio, near, far;

struct point triangle1;
struct point triangle2;
struct point triangle3;

double identity[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

ifstream inputFile;

stack<double (*)[4]> myStack;


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
}

void readRotate(){
}

void readScale(){
}

void handlePush(){
}

void handlePop(){
}

int main(){

    myStack.push(identity);
    string command;
    string blank;

    inputFile.open("1/scene.txt");
    if(!inputFile){
        cout << "Error" << endl;
        exit(1);
    }

    readInitialInfo();

    while (inputFile >> command) {

      if(command.compare("triangle") == 0){
        readTriangle();
      }
      else if(command == "translate"){
        readTranslate();
      }
      else if(command == "rotate"){
        readRotate();
      }
      else if(command == "scale"){
        readScale();
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
    return 0;
}
