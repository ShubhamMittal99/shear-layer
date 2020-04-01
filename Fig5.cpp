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
double ypt_0     =  0;
double xpt_0     =  0;
double T       =  15;

double alpha = 2*pi/lambda;
double c = lambda;

// Initial x-values for the streamlines
vector<double> x_0 = {};
vector<double> y_0 = {-0.1,-0.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,0.08,0.1};
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


double f(double x, double y, double t, double phi){
    return v(x,y,t,phi)/u(x,y,t,phi);
}
/*  Code to solve the equation:
*   dy/dx = v/u = f
*/

//Calculation step
double dx = 0.1;

//Names of Cases
vector<string>CASE = {"a","b","c"};

// Vectors containing the values of a1, a2 and phi for all cases
vector<double> A1  = {0.02, 0.02, 0.01};
vector<double> A2  = {0.02, 0.02, 0.02};
vector<double> PHI = {0, 0.25*pi, 0};

int main(){
    vector<vector<double>> X_all, Y_all;
    vector<double> X,Y;

    //.csv file to output data
    ofstream outfile;
    outfile.open("Fig5-DataPoints.csv");

    
    for(int k = 0; k < CASE.size(); k++){

        X_all.clear(); Y_all.clear();

        //Assigning values for each particular case
        a1 = A1[k];
        a2 = A2[k];
        phi = PHI[k];

        //Initial value of both coordinates taken as xpt__0 and ypt_0
        X.clear(); Y.clear();
        double x = xpt_0;
        double y = ypt_0;
        while(x < T){
            //if(y>=0){
                X.push_back(x);
                Y.push_back(y);
            //}
            y += f(x,y,T,phi)*(dx);
            x += dx;
        }
        X_all.push_back(X);
        Y_all.push_back(Y);

        //Different initial values of y coordinate
        for(int i = 0; i < y_0.size(); i++){
            X.clear(); Y.clear();
            double x = xpt_0;
            double y = y_0[i];
            while(x < T){
                //if(y>=0){
                    X.push_back(x);
                    Y.push_back(y); 
                //}
                y += f(x,y,T,phi)*(dx);
                x += dx;
            }
            X_all.push_back(X);
            Y_all.push_back(Y);
        }
        
        //Case ID (a-c)
        outfile << CASE[k] << endl;
        outfile << endl;

        //Column headings
        outfile << "(x0 | y0) : " << ",("<< xpt_0 <<" | "<< ypt_0 << ") , ,";

        for(int i = 0; i < y_0.size(); i++){
            outfile << "(x0 | y0) : " << ",("<< xpt_0 <<" | "<< y_0[i] << ") , ,";
        }
        outfile << endl;

        for(int i = 0; i < X_all.size(); i++){
            outfile << "X,Y, ,";
        }
        outfile << endl;

        //Output data to CSV file
        for(int i = 0; i < 155; i++){
            for(int j = 0; j < X_all.size(); j++){
                if(i < X_all[j].size()){
                    outfile << X_all[j][i] << "," << Y_all[j][i] << " , ,"; 
                }
                else{
                    outfile << ",,,";
                }
            }
            outfile << endl;
        }
    }
    Beep(523.2511,500);
}