#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <windows.h>
#include <string>

using namespace std;

//Defining the functions u(x,y,t) and v(x,y,t)

double pi = 4*atan(1); // Value of the constant pi

//Defining parameters

//The values of the constants a1, a2 and phi would be accessed from an array for the cases explored in the problem 
double a1;
double a2;
double phi;

//Universal parameters for all the cases;
double lambda  =  2*pi;
double x_0     =  0;
double y_0     =  0;
double T       =  10;

double alpha = 2*pi/lambda;
double c = lambda;


double ubar(double y){
    return 1+ tanh(y);
}

double uprime(double x, double y, double t, double phi){

    double theta1 = alpha*(x - c*t);
    double theta2 = 0.5*alpha*(x - c*t + phi);

    return 2*(tanh(y)/cosh(y))*(a1*sin(theta1) + a2*sin(theta2));
}


double u(double x, double y, double t, double phi){
    return ubar(y) + uprime(x,y,t,phi);
}


double v(double x, double y, double t, double phi){

    double theta1 = alpha*(x - c*t);
    double theta2 = 0.5*alpha*(x - c*t + phi);

    return (alpha/cosh(y))*(2*a1*cos(theta1) + a2*cos(theta2));
}


/*  Code to solve the set of equations:
*   dx/dt = u(x,y,t) 
*   dy/dt = v(x,y,t)
*/

//Timestep
double dt = 0.01;

//No. of particles
double part = 1000;

int main(){
    // Matrix containing the values of coordinates X and Y for
    vector<vector<double>> X_all;
    vector<vector<double>> Y_all;
    

    //Names of Cases
    vector<string>CASE = {"a","b","c", "d", "e", "f","g","h"};

    // Vectors containing the values of a1, a2 and phi for all cases
    vector<double> A1  = {0.01, 0, 0.02, 0.05, 0.01, 0.02, 0.02, 0.04};
    vector<double> A2  = {0, 0.01, 0.02, 0.05, 0.02, 0.02, 0.04, 0.02};
    vector<double> PHI = {0, 0, 0, 0, 0, 0.25*pi, 0.5*pi, -0.25*pi};

    vector<double> X;
    vector<double> Y;

    for(int k = 0; k < 8; k++){
        
        X.clear(); 
        Y.clear();

        //Assigning values for each particular case
        a1 = A1[k];
        a2 = A2[k];
        phi = PHI[k];

        for(int j = 0; j < part; j++){
            
            double t = j*dt;

            //Normalized values
            double x = x_0;
            double y = y_0;

            while(t < T){

                //Functions u and v applied after normalizing
                double x2 = u(x*lambda,y*lambda,t,phi)*dt;
                double y2 = v(x*lambda,y*lambda,t,phi)*dt;
    
                x += x2;
                y += y2;
                t += dt;
            }
            
            X.push_back(x);
            Y.push_back(y);
        }
        X_all.push_back(X);
        Y_all.push_back(Y);
        cout << "Progress :  " << (k/8.0)*100 << " %" <<endl;
    }

    //Outputting data obtained to a .csv file
    ofstream outfile;
    outfile.open("Fig4-DataPoints.csv");

    //The Case id (a-h)
    for(int i = 0; i < 8; i++){
        outfile << CASE[i] <<", , , ,";
    }
    outfile << endl;

    //Column titles X and Y
    for(int i = 0; i < 8; i++){
        outfile << "X" <<", Y , , ,";
    }
    outfile << endl;
    for(int i = 0; i <part; i++){
        for(int j = 0; j < 8; j++){
        outfile << X_all[j][i] << "," << Y_all[j][i] << ", , ,"; 
        }
        outfile << endl;
    }
    cout << "Progress :  " << 100 << " %" <<endl;
    Beep(523.2511,500);
}
