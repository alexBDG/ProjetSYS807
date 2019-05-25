#ifndef _FUNCTION_CPP

#include "Function.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Eigen;

Function::Function(int test, double x0, double M, double rho_0, double u_0, double p_0, double rho_1, double u_1, double p_1, double CFL, string Methode, string Correctif) :
_x0(x0), _M(M), _gamma(1.4), _CFL(CFL), _dx(1./M),
_rho_0(rho_0), _u_0(u_0), _p_0(p_0),
_rho_1(rho_1), _u_1(u_1), _p_1(p_1),
_lambdaMAX(0.), _dt(0.),
_test(test),
_Methode(Methode), _Correctif(Correctif)
{
  _rho.resize(_M+1), _u.resize(_M+1), _p.resize(_M+1);
  _c.resize(_M+1), _u_m.resize(_M+1),_lambda.resize(_M+1);
  _Ue.resize(_M+1,3);
  _n = 0;

  if (_test == 1)
  {
    _rho_0 = 1.0;
    _u_0 = 0.75;
    _p_0 = 1.0;

    _rho_1 = 0.125;
    _u_1 = 0.0;
    _p_1 = 0.1;
  }
  else if (_test == 2)
  {
    _rho_0 = 1.0;
    _u_0 = -2.0;
    _p_0 = 0.4;

    _rho_1 = 1.0;
    _u_1 = 2.0;
    _p_1 = 0.4;
  }
  else if (_test == 3)
  {
    _rho_0 = 1.0;
    _u_0 = 0.0;
    _p_0 = 1000.0;

    _rho_1 = 1.0;
    _u_1 = 0.0;
    _p_1 = 0.01;
  }
  else if (_test == 4)
  {
    _rho_0 = 5.99924;
    _u_0 = 19.5975;
    _p_0 = 460.894;

    _rho_1 = 5.99242;
    _u_1 = -6.19633;
    _p_1 = 46.0950;
  }
  else
  {
    _rho_0 = 1.0;
    _u_0 = -19.59745;
    _p_0 = 1000.0;

    _rho_1 = 1.0;
    _u_1 = -19.59745;
    _p_1 = 0.01;
  }
}


void Function::InitialCondition()
{
  for (int i=0; i<=_M; i++)
  {
    if (i*_dx < _x0)
    {
      _rho[i] = _rho_0;
      _u[i] = _u_0;
      _p[i] = _p_0;
    }
    else
    {
      _rho[i] = _rho_1;
      _u[i] = _u_1;
      _p[i] = _p_1;
    }
    _Ue(i,0) = _rho[i];
    _Ue(i,1) = _rho[i] * _u[i];
    _Ue(i,2) = _p[i]/(_gamma-1) + _rho[i] * _u[i]*_u[i] / 2.;
  }
}


RowVector3d Function::LeftBoundary()
{
  RowVector3d U_cal = Vector3d::Zero(3);
  RowVector3d F_plus_demi = Vector3d::Zero(3);
  RowVector3d F_moins_demi = Vector3d::Zero(3);

  if (_Methode == "Roe")
  {
    F_plus_demi = MethodeRoe(0);
    F_moins_demi = Flow(0);
  }
  else if (_Methode == "LaxFriedrichs")
  {
    F_plus_demi = MethodeLaxFriedrichs(0,0);
    F_moins_demi = Flow(0);
  }

  U_cal = _Ue.row(0) - (_dt/_dx) * (F_plus_demi - F_moins_demi);

  _p[0] = (_gamma-1) * (U_cal(2) - U_cal(1)*U_cal(1)/(2*U_cal(0)));
  _u[0] = U_cal(1)/U_cal(0);
  _rho[0] =  U_cal(0);

  return U_cal;
}


RowVector3d Function::RightBoundary()
{
  RowVector3d U_cal = Vector3d::Zero(3);
  RowVector3d F_plus_demi = Vector3d::Zero(3);
  RowVector3d F_moins_demi = Vector3d::Zero(3);

  if (_Methode == "Roe")
  {
    F_plus_demi = Flow(_M);
    F_moins_demi = MethodeRoe(_M-1);
  }
  else if (_Methode == "LaxFriedrichs")
  {
    F_plus_demi = Flow(_M);
    F_moins_demi = MethodeLaxFriedrichs(_M-1,_M);
  }

  U_cal = _Ue.row(_M) - (_dt/_dx) * (F_plus_demi - F_moins_demi);

  _p[_M] = (_gamma-1) * (U_cal(2) - U_cal(1)*U_cal(1)/(2*U_cal(0)));
  _u[_M] = U_cal(1)/U_cal(0);
  _rho[_M] =  U_cal(0);

  return U_cal;
}


void Function::Celerity()
{
  if (_Methode == "Roe")
  {
    double H_R, H_L, H_m;
    for (int i=0; i<_M; i++)
    {
      H_R = (_Ue(i+1,2) + _p[i+1])/_rho[i+1];
      H_L = (_Ue(i,2) + _p[i])/_rho[i];

      // Calculer les valeurs moyennes _u_moy, _H_moy et _a_moy
      _u_m[i] = (sqrt(_rho[i]) * _u[i] + sqrt(_rho[i+1]) * _u[i+1]) / (sqrt(_rho[i]) + sqrt(_rho[i+1]));
      H_m = (sqrt(_rho[i]) * H_L + sqrt(_rho[i+1]) * H_R) / (sqrt(_rho[i]) + sqrt(_rho[i+1]));
      _c[i] = sqrt((_gamma-1)* (H_m - 0.5 * _u_m[i] * _u_m[i]) );
    }
    H_L = (_Ue(_M,2) + _p[_M])/_rho[_M];
    _u_m[_M] = sqrt(_rho[_M]) * _u[_M] / sqrt(_rho[_M]);
    H_m = sqrt(_rho[_M]) * H_L / sqrt(_rho[_M]);
    _c[_M] =sqrt( (_gamma-1)*(H_m - 0.5 * _u_m[_M] * _u_m[_M]) );
  }
  else if (_Methode == "LaxFriedrichs")
  {
    for (int i=0; i<=_M; i++)
    {
      _c[i] = sqrt( _gamma * _p[i] / _rho[i] );
    }
  }
}

void Function::Lambda()
{
  if (_Methode == "Roe")
  {
    _lambda[0] = max(abs(_u_m[0] + _c[0]),abs(_u_m[0] - _c[0]));
    _lambdaMAX = abs( _lambda[0] );
    for (int i=1; i<=_M; i++)
    {
      _lambda[i] = max(abs(_u_m[i] + _c[i]),abs(_u_m[i] - _c[i]));
      _lambdaMAX = max( _lambdaMAX , abs(_lambda[i]) );
    }
  }
  else if (_Methode == "LaxFriedrichs")
  {
    _lambda[0] = abs( _Ue(0,1) / _Ue(0,0) ) + _c[0];
    _lambdaMAX = abs( _lambda[0] );
    for (int i=1; i<=_M; i++)
    {
      _lambda[i] = abs( _u[i] ) + _c[i];
      _lambdaMAX = max( _lambdaMAX , abs(_lambda[i]) );
    }
  }

}

double Function::CFL()
{
  if (_n == 5)
  {
    _dt = _CFL * _dx / _lambdaMAX;
  }
  else
  {
    _dt = 0.2 * _CFL * _dx / _lambdaMAX;
    _n += 1;
  }
  return _dt;
}

RowVector3d Function::Flow(int i)
{
  RowVector3d F = RowVector3d::Zero(3);
  F[0] = _Ue(i,1);
  F[1] = _Ue(i,1)*_Ue(i,1)/_Ue(i,0) + (_gamma-1) * (_Ue(i,2) - _Ue(i,1)*_Ue(i,1) / (2*_Ue(i,0)) );
  F[2] = _Ue(i,2)*_Ue(i,1)/_Ue(i,0) + (_gamma-1) * (_Ue(i,1) / _Ue(i,0)) * (_Ue(i,2) - _Ue(i,1)*_Ue(i,1) / (2*_Ue(i,0)) );
  return F;
}

RowVector3d Function::MethodeRoe(int i)
{
  double rho_R = _rho[i+1];
  double rho_L = _rho[i];

  double u_R = _u[i+1];
  double u_L = _u[i];

  double H_R = (_Ue(i+1,2) + _p[i+1])/_rho[i+1];
  double H_L = (_Ue(i,2) + _p[i])/_rho[i];

  // Calculer les valeurs moyennes _u_moy, _H_moy et _a_moy
  double _u_moy = (sqrt(rho_L) * u_L + sqrt(rho_R) * u_R) / (sqrt(rho_L) + sqrt(rho_R));
  double _H_moy = (sqrt(rho_L) * H_L + sqrt(rho_R) * H_R) / (sqrt(rho_L) + sqrt(rho_R));
  double _a_moy = sqrt( (_gamma-1)*(_H_moy - 0.5 * _u_moy * _u_moy) );

  // Calculer les valeurs propres moyennes _lambda_i_moy
  _lambda_1_moy = _u_moy - _a_moy;
  _lambda_2_moy = _u_moy;
  _lambda_3_moy = _u_moy + _a_moy;

  //Eventuellement, on applique un correctif d'entropie
  if (_Correctif == "1")
  {
    double lambda_L, lambda_R, epsilon;
    double a_L = sqrt( (_gamma-1)*(H_L - 0.5 * u_L * u_L) );
    double a_R = sqrt( (_gamma-1)*(H_R - 0.5 * u_R * u_R) );
    // 1ère valeur propre
    lambda_L = u_L - a_L;
    lambda_R = u_R - a_R;
    epsilon = max(0., max(_lambda_1_moy - lambda_L, lambda_R - _lambda_1_moy) );
    if (abs(_lambda_1_moy)<epsilon)
    {
      _lambda_1_moy = (_lambda_1_moy*_lambda_1_moy + epsilon*epsilon)/(2*epsilon);
    }
    // 2ème valeur propre
    lambda_L = u_L;
    lambda_R = u_R;
    epsilon = max(0., max(_lambda_2_moy - lambda_L, lambda_R - _lambda_2_moy) );
    if (abs(_lambda_2_moy)<epsilon)
    {
      _lambda_2_moy = (_lambda_2_moy*_lambda_2_moy + epsilon*epsilon)/(2*epsilon);
    }
    // 3ème valeur propre
    lambda_L = u_L + a_L;
    lambda_R = u_R + a_R;
    epsilon = max(0., max(_lambda_2_moy - lambda_L, lambda_R - _lambda_2_moy) );
    if (abs(_lambda_3_moy)<epsilon)
    {
      _lambda_3_moy = (_lambda_3_moy*_lambda_3_moy + epsilon*epsilon)/(2*epsilon);
    }
  }

  if (_Correctif == "2")
  {
    double lambda_L, lambda_R, epsilon;
    double a_L = sqrt( (_gamma-1)*( H_L - 0.5 * u_L * u_L) );
    double a_R = sqrt( (_gamma-1)*( H_R - 0.5 * u_R * u_R) );
    // 1ère valeur propre
    lambda_L = u_L - a_L;
    lambda_R = u_R - a_R;
    epsilon = max(0., max(_lambda_1_moy - lambda_L, lambda_R - _lambda_1_moy) );
    if (abs(_lambda_1_moy)<epsilon)
    {
      _lambda_1_moy = epsilon;
    }
    // 2ème valeur propre
    lambda_L = u_L;
    lambda_R = u_R;
    epsilon = max(0., max(_lambda_2_moy - lambda_L, lambda_R - _lambda_2_moy) );
    if (abs(_lambda_2_moy)<epsilon)
    {
      _lambda_2_moy = epsilon;
    }
    // 3ème valeur propre
    lambda_L = u_L + a_L;
    lambda_R = u_R + a_R;
    epsilon = max(0., max(_lambda_3_moy - lambda_L, lambda_R - _lambda_3_moy) );
    if (abs(_lambda_3_moy)<epsilon)
    {
      _lambda_3_moy = epsilon;
    }
  }

  if (_Correctif == "3")
  {
    double epsilon = 0.3;
    if (abs(_lambda_1_moy)<epsilon)
    {
      _lambda_1_moy = (_lambda_1_moy*_lambda_1_moy + epsilon*epsilon)/(2*epsilon);
    }
    if (abs(_lambda_2_moy)<epsilon)
    {
      _lambda_2_moy = (_lambda_2_moy*_lambda_2_moy + epsilon*epsilon)/(2*epsilon);
    }
    if (abs(_lambda_3_moy)<epsilon)
    {
      _lambda_3_moy = (_lambda_3_moy*_lambda_3_moy + epsilon*epsilon)/(2*epsilon);
    }
  }

  // Calculer les vecteurs propres moyens _K_i_moy
  RowVector3d _K_1_moy = RowVector3d::Zero(3);
  RowVector3d _K_2_moy = RowVector3d::Zero(3);
  RowVector3d _K_3_moy = RowVector3d::Zero(3);

  _K_1_moy[0] = 1.;
  _K_1_moy[1] = _u_moy - _a_moy;
  _K_1_moy[2] = _H_moy - _u_moy * _a_moy;

  _K_2_moy[0] = 1.;
  _K_2_moy[1] = _u_moy;
  _K_2_moy[2] = 0.5 * _u_moy * _u_moy;

  _K_3_moy[0] = 1.;
  _K_3_moy[1] = _u_moy + _a_moy;
  _K_3_moy[2] = _H_moy + _u_moy * _a_moy;

  // Calculer les amplitude d'ondes _alpha_i_moy
  double _alpha_2_moy = (_gamma-1)/(_a_moy*_a_moy)*( (_Ue(i+1,0)-_Ue(i,0)) * (_H_moy - _u_moy * _u_moy) + _u_moy *(_Ue(i+1,1)-_Ue(i,1)) - (_Ue(i+1,2) - _Ue(i,2)) );
  double _alpha_1_moy = 1./(2*_a_moy) * ( (_Ue(i+1,0)-_Ue(i,0)) * (_u_moy + _a_moy) - (_Ue(i+1,1)-_Ue(i,1)) - _a_moy * _alpha_2_moy );
  double _alpha_3_moy = (_Ue(i+1,0)-_Ue(i,0)) - (_alpha_1_moy + _alpha_2_moy);

  // Calculer _F_plus_demi
  RowVector3d F = (Flow(i) + Flow(i+1))/2. - 0.5 * (_alpha_1_moy * abs(_lambda_1_moy) * _K_1_moy + _alpha_2_moy * abs(_lambda_2_moy) * _K_2_moy + _alpha_3_moy * abs(_lambda_3_moy) * _K_3_moy);
  return F;
}


RowVector3d Function::MethodeLaxFriedrichs(int i, int j)
{
  RowVector3d F = (Flow(i) + Flow(i+1))/2. - (_lambda(j)/2.) * (_Ue.row(i+1) - _Ue.row(i));
  return F;
}


void Function::Update()
{
  RowVector3d F_plus_demi = RowVector3d::Zero(3);
  RowVector3d F_moins_demi = RowVector3d::Zero(3);
  MatrixXd Ue_temp = MatrixXd::Zero(_M+1,3);

  Ue_temp.row(0) = LeftBoundary();
  for (int i=1; i<_M; i++)
  {

    if (_Methode == "Roe")
    {
      F_plus_demi = MethodeRoe(i);
      F_moins_demi = MethodeRoe(i-1);
    }
    else if (_Methode == "LaxFriedrichs")
    {
      F_plus_demi = MethodeLaxFriedrichs(i,i);
      F_moins_demi = MethodeLaxFriedrichs(i-1,i);
    }

    Ue_temp.row(i) = _Ue.row(i) - (_dt/_dx) * (F_plus_demi - F_moins_demi);
  }
  Ue_temp.row(_M) = RightBoundary();

  //Mise à jour de _p, _u, _rho
  for (int i=0; i<=_M; i++)
  {
    //Mise à jour de _Ue
    _Ue(i,0) = Ue_temp(i,0);
    _Ue(i,1) = Ue_temp(i,1);
    _Ue(i,2) = Ue_temp(i,2);
    _p[i] = (_gamma-1.) * (_Ue(i,2) - _Ue(i,1)*_Ue(i,1)/(2*_Ue(i,0)));
    _u[i] = _Ue(i,1)/_Ue(i,0);
    _rho[i] =  _Ue(i,0);
  }
}


void Function::SaveSol()
{
	string name_file = "Methode_" + _Methode + "_Correctif_" + _Correctif + "_M_" + std::to_string((int)_M) + "_Test_" + std::to_string((int)_test) + ".dat";

  ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  for (int i = 0 ; i <= _M ; ++i)
  {
    solution << i*_dx  << " " << float(_rho[i])  << " "  << float(_u[i])  << " "  << float(_p[i]) << " " << float(_Ue(i,2)/_Ue(i,0)-0.5*(_Ue(i,1)/_Ue(i,0))*(_Ue(i,1)/_Ue(i,0))) << endl;
  }
  solution << endl;

  _solution_lambda << endl;
  _solution_f << endl;

	solution.close();
  _solution_lambda.close();
  _solution_f.close();
}


#define _FUNCTION_CPP
#endif
