#ifndef _FUNCTION_H
#include <string>
#include <fstream>
#include <iostream>
#include "Dense"
#include "Sparse"

class Function {
  private:
    Eigen::RowVectorXd _rho, _p, _u, _c, _u_m, _lambda;
    double _n, _x0, _gamma, _lambdaMAX, _M, _CFL, _rho_0, _u_0, _p_0, _rho_1, _u_1, _p_1, _dx, _dt;
    double _lambda_1_moy, _lambda_2_moy, _lambda_3_moy;
    Eigen::MatrixXd _Ue;
    int _test;
    std::string _Methode, _Correctif;
    std::ofstream _solution_lambda;
    std::ofstream _solution_f;


	public: // Méthodes et opérateurs de la classe
    Function(int test, double x0, double M, double rho_0, double u_0, double p_0, double rho_1, double u_1, double p_1, double CFL, std::string Methode, std::string Correctif);
    void InitialCondition();
    Eigen::RowVector3d LeftBoundary();
    Eigen::RowVector3d RightBoundary();
    void Celerity();
    void Lambda();
    double CFL();
    Eigen::RowVector3d MethodeLaxFriedrichs(int i, int j);
    Eigen::RowVector3d MethodeRoe(int i);
    Eigen::RowVector3d Flow(int i);
    void Update();
    void SaveSol();
};

#define _FUNCTION_H
#endif
