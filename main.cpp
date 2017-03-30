/* 
 * File:   main.cpp
 * Author: almaredan
 *
 * Created on 30 марта 2017 г., 13:14
 * 
 * Реализация метода Галёркина для неявной схемы:
 * Отрезок: [0, 1]
 * f(0,t) = 3
 * f(x,0) = 1
 * Функция в правой части: x^3
 * 
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "vector3d.h"
#include <vector>
#include <math.h>
//#include <bits/c++config.h>
#include <algorithm>
#include <functional>
//#include <complex>

typedef double real;
#define tol12 1.e-12

using namespace std;
/*
 * 
 */

void scheme(vector<double>&, vector<double>&, double, double, double);
double integrate(std::function< double (double y) > f, double x0, double h);
double f_x(double x);

int main(int argc, char** argv) {
    double x1 = 0;
    double x2 = 1;
    double vel_l = 1;                                                           //Значение f(0,t)
    int N_u = 50;                                                               //Количество точек для грубой сетки
    int N_v = 100;
//    int N_w = 200;
    vector<double> mesh_u, mesh_v; //, mesh_w;
    vector<double> u, v;//, w;
    u.assign(N_u, 1);
    v.assign(N_v, 1);
//    w.assign(N_w, u_l);
        
    scheme(u, mesh_u, x1, x2, vel_l);
    cout << "Первая схема отработала" << endl;
    double err_u = 0;
//    double x_u = 0;
    for (int i = 0; i <  N_u; i++) {
        double err_tmp = abs(u[i] - exp(mesh_u[i]));               //Не забудь заменить здесь и ниже
        err_u = (err_u > err_tmp) ? err_u : err_tmp;                            //константу рядом с f_x(), если она изменится
//        x_u = (err_u > err_tmp) ? x_u : mesh_u[i];
        cout << mesh_u[i] << " " << u[i] << endl;
    } cout << endl;
    
    
    scheme(v, mesh_v, x1, x2, vel_l);
    cout << "Вторая схема отработала" << endl;
    double err_v = 0;
    for (int i = 0; i <  N_v; i++) {
        double err_tmp = abs(v[i]  - exp(mesh_v[i]));            //Не забудь заменить здесь и ниже
        err_v = (err_v > err_tmp) ? err_v : err_tmp;                    //константу рядом с f_x(), если она изменится
        cout << mesh_v[i] << " " << v[i] << endl;
    } cout << endl;
    
//    scheme(w, mesh_w, x1, x2);
//    double err_w = 0;
//    for (int i = 0; i <  N_w; i++) {
//        double err_tmp = abs(w[i] - pow(mesh_w[i],4)/4 - u_l);            //Не забудь заменить здесь и ниже
//        err_w = (err_w > err_tmp) ? err_w : err_tmp;                    //константу рядом с f_x(), если она изменится
//        cout << w[i] << " ";
//    } cout << endl;
    
    double p = abs( logf( err_u / err_v ) / logf( 0.5 ) );
    cout << "Errors: " << err_u << " " << err_v << endl;
    cout << "Порядок точности: " << p << endl;
    
    return 0;
}

double integrate(double (*f)(double), double x0, double h) {
    return f(x0)*h  + ( f(x0+h/2) + f(x0-h/2) - 2*f(x0) )*h/6;
}

double integrate(std::function< double (double y) > f, double x0, double h) {
    return f(x0)*h  + ( f(x0+h/2) + f(x0-h/2) - 2*f(x0) )*h/6;
}

double f_x(double x) {
    return exp(x);
}

Vector3D discrepancy(const Vector3D & vel, double & uLeft,
        double h, double x0) {
    /*
     * Функция невязки в ячейки:
     * Inputs: Vector3D vel - вектор скорости на предыдущем временном слое в данной ячейке
     *         double uLeft - нулевое приближение скорости в левой ячейки
     *         double h - шаг ячейки
     *         double x0 - координата центра ячейки
     * Output: Vector3D - значение невязки в ячейке    
     */
    
    Vector3D phi_i (1., .5, .0);
    real uRight =  phi_i * vel;
    Vector3D firstVec ( -uRight + uLeft, - 0.5*(uRight+uLeft), .0 );
    Vector3D secondVec (.0, vel[0], .0 );
    std::function< double (double y) > f = [&](double y) -> double
    {
        return f_x(y)*(y - x0)/h;
    };
    Vector3D thirdVec ( integrate(f_x, x0, h),
            integrate(f, x0, h), .0 );
    
    Vector3D tempVec = firstVec + secondVec + thirdVec;
    
    return tempVec;
}

Vector3D step(const Vector3D & vel, double tau, double h, double x0,
                    double& uLeft) {
    Vector3D phi_i (1., .5, .0);
    real uRight =  phi_i * vel;
    Vector3D tempVec = discrepancy(vel, uLeft, h, x0);
    
    std::function< Vector3D (Vector3D u) > R = [&](Vector3D u) -> Vector3D {    //Приведение невязки к функции от одной переменной
        return discrepancy(u, uLeft, h, x0);                                    //для более удобного дифференцирования в дальнейшем
    };
        
    double eps = h/2;
    Vector3D diffVel1 (eps, .0, .0);                                            //Отклонения скоростей
    Vector3D diffVel2 (.0, eps, .0);                                            //для поиска производных
    
    real R11 = (tempVec.x - R(vel + diffVel1).x)/eps;                            //первая строка по первой компоненте
    real R12 = (tempVec.x - R(vel + diffVel2).x)/eps;                            //первая строка по второй компоненте
    real R21 = (tempVec.y - R(vel + diffVel1).y)/eps;                            //вторая --//-- по первой
    real R22 = (tempVec.y - R(vel + diffVel2).y)/eps;                            //вторая --//-- по второй
    
    real a = h/tau + R11;
    real b = R12;
    real c = R21;
    real d = h/(12*tau) + R22;
    
    
    Vector3D rVec (( d*tempVec[0] - b*tempVec[1] ),
            ( -c*tempVec[0] + a*tempVec[1] ), .0);
    rVec *= 1/abs(a*d - b*c);
    
    Vector3D velNew = vel + rVec;         //дописал
    
    uLeft = uRight;
    return velNew;
}



void scheme(vector<double>& u, vector<double>& mesh, double x1,
        double x2, double vel_l) {
    vector< Vector3D > vel;
    double Cu = 0.3;
    double h = (x2 - x1)/u.size();
    double tau = h*Cu;
    for (int i = 0; i < u.size(); i++) mesh.push_back(x1 + i*h + h/2);
    for ( double value : u) {
        vel.push_back( Vector3D ( value, .0, .0 ) );
    }
    double epsilon = 100;
    while (epsilon > 1.e-20){
        double uLeft = vel_l;//u.at(0);
        vector< Vector3D > velTmp;
        for (int i = 0; i < vel.size(); i++) {
            velTmp.push_back(step(vel.at(i), tau, h, mesh[i], uLeft));
        }
        
        epsilon = 0;
        
        for (int i = 0; i < vel.size(); i++) {
            double diff = abs(vel.at(i)[0] - velTmp.at(i)[0]);
            epsilon = ( epsilon < diff ) ? diff : epsilon;
            vel.at(i) = velTmp.at(i);
        }
    }
    for (int i = 0; i < u.size(); i++) {
        u.at(i) = vel.at(i)[0];
    }
}