#include <cstdlib>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include "qsim.h"

using namespace std;

int main(int argc, char *argv[])
{
  qubit psi(1, 2);
  psi.display();
  
  vector<qubit> q;
  q.push_back(qubit::ZERO());
  q.push_back(qubit::ZERO());
  qreg reg(q);
  
  // Quantum Teleportation
  // http://en.wikipedia.org/wiki/Quantum_teleportation
  reg.feed(qopr::H(), qopr::ID());   // Generate EPR Pair
  reg.feed(qopr::CNOT());            // Generate EPR Pair
  reg.display();
  
  reg.push_front(psi);
  reg.display();
  reg.feed(qopr::CNOT(), qopr::ID());
  reg.display();
  reg.feed(qopr::H(), qopr::ID(), qopr::ID());
  reg.display();
  int M0 = reg.MES(0);
  cout<<M0<<"  ";
  reg.display();
  int M1 = reg.MES(0);
  cout<<M1<<"  ";
  reg.display();
  if(M1)
    reg.feed(qopr::X());
  if(M0)
    reg.feed(qopr::Z());
  reg.display(); 
  
  system("PAUSE");
  return EXIT_SUCCESS;
}
