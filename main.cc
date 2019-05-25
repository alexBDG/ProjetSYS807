#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "Function.h"

using namespace std;
using namespace Eigen;

int main()
{
  int M, Nb, test, Nb_Methode, Nb_Correctif;
  double dt, t, x0, rho_0, u_0, p_0, rho_1, u_1, p_1, CFL;
  string Methode, Correctif;



////////////////////////////////////////////////////////////////////////////////

  //Séparation
  x0 = 0.5;

  //paramètres physiques
  test = 2;

  //discrétisation du maillage
  M = 100;

  //Durée désirée
  t = 0.15;

  //condition de CFL
  CFL = 0.01;

////////////////////////////////////////////////////////////////////////////////

  cout << "Numéro du test : " << endl;
  cin >> test;

  if (test == 1)
  {
    x0 = 0.3;
    t = 0.2;
    cout << "Nombre de maille : " << endl;
    cin >> M;
  }
  else if (test == 2)
  {
    x0 = 0.5;
    t = 0.15;
    cout << "Nombre de maille : " << endl;
    cin >> M;
  }
  else if (test == 3)
  {
    x0 = 0.5;
    t = 0.012;
    cout << "Nombre de maille : " << endl;
    cin >> M;
  }
  else if (test == 4)
  {
    x0 = 0.4;
    t = 0.035;
    cout << "Nombre de maille : " << endl;
    cin >> M;
  }

  cout << "Méthode utilisée : " << endl;
  cout << "1 - Lax-Friedrichs" << endl;
  cout << "2 - Roe" << endl;
  cin >> Nb_Methode;

  if (Nb_Methode == 1)
  {
    Methode = "LaxFriedrichs";
    Correctif = "0";
  }
  else if (Nb_Methode == 2)
  {
    Methode = "Roe";
    cout << "Correction d'entropie : " << endl;
    cout << "1 - Non" << endl;
    cout << "2 - HH2" << endl;
    cout << "3 - HH1" << endl;
    cout << "4 - HC" << endl;
    cin >> Nb_Correctif;

    if (Nb_Correctif == 1)
    {
        Correctif = "0";
    }
    else if (Nb_Correctif == 2)
    {
        Correctif = "1";
    }
    else if (Nb_Correctif == 3)
    {
        Correctif = "2";
    }
    else if (Nb_Correctif == 4)
    {
        Correctif = "3";
    }
  }
  else
  {
    cout << "La méthode est inconnue ou n'est pas implémentée...";
  }

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();


  // Conditions initiales
  Function* fct = new Function(test,x0,M,rho_0,u_0,p_0,rho_1,u_1,p_1,CFL,Methode,Correctif);
  fct->InitialCondition();
  cout << "-------------------------------------------------" << endl;
  cout << "Calcul des solutions dans le temps : " << endl;


  // Boucle sur le temps
  double D_t = 0;
  while (t > (D_t))
  {
    //Avancement
    cout.flush();
    cout << "Progression : " << D_t/t*100  << " %        \r";

    // 1) Calcul des célérités
    fct->Celerity();

    // 2) Calcul des lambda i
    // 3) Calcul du maximum des lambda i
    fct->Lambda();

    // 4) Calcul du pas de temps
    D_t += fct->CFL();
//    cout << "D_t = " << fct->CFL() << endl;

    // 5) Mise à jour des variables
    fct->Update();

  }

  //Enregistrement de la solution finale
  fct->SaveSol();

  cout << "-------------------------------------------------" << endl;

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double chrono_temps = chrono::duration_cast<chrono::milliseconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< chrono_temps/1000. << " secondes" << endl;

  return 0;
}
