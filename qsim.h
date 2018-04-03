
#include <cstdlib>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

class qubit
{
  public:
    std::complex<float> a[2];
    qubit(float, float);
    qubit(std::complex<float>, std::complex<float>);
    void display();
    static qubit ZERO();
    static qubit ONE();
  private:
    void normalize();
};

qubit::qubit(float a0, float a1)
{
  a[0] = a0;
  a[1] = a1;
  normalize();
}

qubit::qubit(std::complex<float> a0, std::complex<float> a1)
{
  a[0] = a0;
  a[1] = a1;
  normalize();
}

void qubit::display()
{
  for(int i = 0; i < 2; i++)
  {
    if(std::abs(std::imag(a[i])) > 0)
      std::cout<<a[i];
    else
      std::cout<<std::real(a[i]);
    std::cout<<" |"<<i<<">";
    if(i==0)
      std::cout<<" + ";
  }
  std::cout<<std::endl;
}

qubit qubit::ZERO()
{
  return qubit(1, 0);
}

qubit qubit::ONE()
{
  return qubit(0, 1);
}

void qubit::normalize()
{
  std::complex<float> norm = sqrt(a[0]*a[0] + a[1]*a[1]);
  a[0] = a[0] / norm;
  a[1] = a[1] / norm;
}

class qopr
{
  public:
    std::vector<std::vector< std::complex<float> > > u;
    bool IsControl;
    int ControlType;  // 1 or 0
    qopr(int dim);
    qopr(float, float, float, float); // 2x2 matrix
    qopr(std::complex<float>, std::complex<float>, std::complex<float>, std::complex<float>); // 2x2 matrix
    qopr(std::vector<qopr>);
    qopr(qopr, qopr);
    qopr(qopr, qopr, qopr);
    static qopr ID();
    static qopr NOT();
    static qopr CTRL();
    static qopr CTRL0();
    static qopr H();  // Hadamard
    static qopr X();  // Pauli
    static qopr Y();  // Pauli
    static qopr Z();  // Pauli
    static qopr S();  // Phase
    static qopr T();  // Pi/8
    static qopr CNOT();
    static qopr SWAP();
    static qopr CZ();
    static qopr CS();
    static qopr TOFFOLI();
    static qopr FREDKIN();
    void display();
  private:
    void Init(int);
};

void qopr::Init(int dim)
{
//  this -> dim = dim;
  this -> ControlType = 1;
  this -> IsControl = false;
  u.resize(dim);
  for(int i=0; i<u.size(); i++)
    u[i].resize(dim);
  for(int i=0; i<u.size(); i++)
  {
    for(int j=0; j<u[i].size(); j++)
    {
      u[i][j] = 0;
    }
  }
}

qopr::qopr(int dim = 2)
{
  Init(dim);
}

qopr::qopr(float u00, float u01, float u10, float u11)
{
  Init(2);
  u[0][0] = u00;
  u[0][1] = u01;
  u[1][0] = u10;
  u[1][1] = u11;
}

qopr::qopr(std::complex<float> u00, std::complex<float> u01, std::complex<float> u10, std::complex<float> u11)
{
  Init(2);
  u[0][0] = u00;
  u[0][1] = u01;
  u[1][0] = u10;
  u[1][1] = u11;
}

qopr::qopr(std::vector<qopr> U)
{
  int dim = 1;
  for(int k=0; k<U.size(); k++)
  {
    dim = dim * U[k].u.size();
  }
  Init(dim);
  for(int i=0; i<u.size(); i++)
  {
    for(int j=0; j<u[i].size(); j++)
    {
      u[i][j] = 1;
    }
  }
  int len = dim;
  for(int k=0; k<U.size(); k++)
  {
    len = len/U[k].u.size();
    for(int i=0; i<u.size(); i++)
    {
      for(int j=0; j<u[i].size(); j++)
      {
        u[i][j] = u[i][j] * U[k].u[ (i/len) % U[k].u.size() ][ (j/len) % U[k].u.size() ];
      }
    }
  }
  
  len = dim;
  for(int k=0; k<U.size(); k++)
  {
    len = len/U[k].u.size();
    if(U[k].IsControl)
    {
      for(int i=0; i<u.size(); i++)
      {
        if( (i/len) % 2 == 1 - U[k].ControlType)
        {
          for(int j=0; j<u[i].size(); j++)
          {
            u[i][j] = (int)(i==j);
          }
        }
      }
    }
  }    
}


qopr::qopr(qopr U, qopr V)
{
  std::vector<qopr> UV;
  UV.push_back(U); 
  UV.push_back(V);
  *this = qopr(UV);
}

qopr::qopr(qopr U, qopr V, qopr W)
{
  std::vector<qopr> UVW;
  UVW.push_back(U);
  UVW.push_back(V);
  UVW.push_back(W);
  *this = qopr(UVW);
}

void qopr::display()
{
  for(int i=0; i<u.size(); i++)
  {
    for(int j=0; j<u[i].size(); j++)
    {
      if(std::abs(std::imag(u[i][j]))<0.01)
        std::cout<<std::real(u[i][j])<<" ";
      else
        std::cout<<u[i][j]<<" ";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

qopr qopr::ID()
{
  return qopr(1, 0, 0, 1);
}

qopr qopr::NOT()
{
  return qopr(0, 1, 1, 0);
}

qopr qopr::CTRL()
{
  qopr CTRL = qopr::ID();
  CTRL.ControlType = 1;
  CTRL.IsControl = true;
  return CTRL;
}

qopr qopr::CTRL0()
{
  qopr CTRL0 = qopr::ID();
  CTRL0.ControlType = 0;
  CTRL0.IsControl = true;
  return CTRL0;
}

qopr qopr::H()
{
  float s2 = 1/sqrt(2.0);
  return qopr(s2, s2, s2, -s2);
}

qopr qopr::X()
{
  return qopr(0, 1, 1, 0);
}

qopr qopr::Y()
{
  return qopr(std::complex<float>(0, 0), std::complex<float>(0, -1),
                   std::complex<float>(0, 1), std::complex<float>(0, 0));
}

qopr qopr::Z()
{
  return qopr(1, 0, 0, -1);
}

qopr qopr::S()
{
  return qopr(std::complex<float>(1, 0), std::complex<float>(0, 0),
                   std::complex<float>(0, 0), std::complex<float>(0, 1));
}

qopr qopr::T()
{
  return qopr(std::complex<float>(1, 0), std::complex<float>(0, 0),
                   std::complex<float>(0, 0), exp(std::complex<float>(0, 3.14/4)));
}

qopr qopr::CNOT()
{
  return qopr(qopr::CTRL(), qopr::NOT());
}

qopr qopr::SWAP()
{
  qopr SWAP = qopr(4);
  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      SWAP.u[i][j] = 0;
    }
  }
  SWAP.u[0][0] = 1;
  SWAP.u[1][2] = 1;
  SWAP.u[2][1] = 1;
  SWAP.u[3][3] = 1;
  return SWAP;
}

qopr qopr::CZ()
{
  return qopr(qopr::CTRL(), qopr::Z());
}

qopr qopr::CS()
{
  return qopr(qopr::CTRL(), qopr::S());
}

qopr qopr::TOFFOLI()
{
  return qopr(qopr::CTRL(), qopr::CTRL(), qopr::NOT());
}

qopr qopr::FREDKIN()
{
  return qopr(qopr::CTRL(), qopr::SWAP());
}

class qreg
{
  public:
    qreg(std::vector<qubit>);
    void push_back(qubit);
    void push_front(qubit);
    void feed(qopr U);
    void feed(qopr U, qopr V);
    void feed(qopr U, qopr V, qopr W);
    int  MES(int, int);
    void display();
    
  private:
    int HowManyqubits;
    std::vector< std::complex<float> > a;  // amplitude
    std::vector< qopr > CurrentOperatorArray;

    void normalize();
};

qreg::qreg(std::vector<qubit> q)
{
  srand(time(NULL));
  HowManyqubits = q.size();
  a.resize((int)pow(2.0, HowManyqubits));
  for(int j = 0; j<a.size(); j++)
  {
    a[j] = 1;
  }
  int len = a.size();
  for(int i = 0; i<HowManyqubits; i++)
  {
    len = len/2;
    for(int j = 0; j<a.size(); j++)
    {
      a[j] = a[j] * q[i].a[ (j/len) % 2 ];
    }
  }
  normalize();
  CurrentOperatorArray.clear();
}

void qreg::push_back(qubit q)
{
  a.resize(a.size() * 2);
  for(long int i = a.size()-1; i>=0; i--)
  {
    a[i] = a[i/2];
  }  
  HowManyqubits = HowManyqubits + 1;
  for(long int i = 0; i < a.size(); i++)
  {
    a[i] = a[i] * q.a[i%2];
  }
  normalize();
}

void qreg::push_front(qubit q)
{
  a.resize(a.size() * 2);
  for(long int i = a.size()/2; i < a.size(); i++)
  {
    a[i] = a[i - a.size()/2];
  }  
  HowManyqubits = HowManyqubits + 1;
  for(long int i = 0; i < a.size()/2; i++)
  {
    a[i] = a[i] * q.a[0];
  }
  for(long int i = a.size()/2; i < a.size(); i++)
  {
    a[i] = a[i] * q.a[1];
  }
  normalize();
}

void qreg::feed(qopr U)
{
  long int CurrentOperatorSize = 1;
  for(long int i = 0; i<CurrentOperatorArray.size(); i++)
    CurrentOperatorSize *= CurrentOperatorArray[i].u.size();
  CurrentOperatorSize *= U.u.size();
  if( CurrentOperatorSize <= a.size() )
  {
    CurrentOperatorArray.push_back(U);
  }
  if( CurrentOperatorSize == a.size() )
  {
    // matrix multiplication
    qopr CurrentOperator = qopr(CurrentOperatorArray);
    std::vector< std::complex<float> > result(a.size());
    for(long int i=0; i<a.size(); i++)
    {
      result[i] = 0;
      for(long int j = 0; j<a.size(); j++)
      {
        result[i] += CurrentOperator.u[i][j] * a[j];
      }
    }
    a = result;
    CurrentOperatorArray.clear();
  }
  if( CurrentOperatorSize > a.size() )
  {
    // Do nothing, Throw an exceptoion
  }
}

void qreg::feed(qopr U, qopr V)
{
  feed(U);
  feed(V);
}

void qreg::feed(qopr U, qopr V, qopr W)
{
  feed(U);
  feed(V);
  feed(W);
}

int qreg::MES(int qubitIndex, int DebugValue = -1)
{
  // qubitIndex'th qubit collapses into pure |0> or }1>
  // The MEaSurement makes the quantum register one qubit smaller
  if ( qubitIndex >= 0 && qubitIndex < HowManyqubits )
  {
    HowManyqubits = HowManyqubits - 1;
    long int len = (int)pow(2.0, HowManyqubits - qubitIndex);
    float prob0 = 0;
    for(long int i = 0; i<a.size()/2 /len; i++)
    {
      for(long int j = 0; j<len; j++)
      {
        prob0 = prob0 + real( a[ i*2 * len + j ] *
                              a[ i*2 * len + j ] );
      }
    }
    int X = rand()%100;
    int result = (int) (X >= 50); // prob0 * RAND_MAX
    if(DebugValue == 0 || DebugValue == 1)
    {
      result = DebugValue;
    }
    
    for(long int i = 0; i < a.size()/2 /len; i++)
    {
      for(long int j = 0; j < len; j++)
      {
        a[i*len + j] = a[ (i*2 + result) * len + j];
      }
    }
    a.resize(a.size()/2);
    normalize();
    return result;
  }
  else
  {
    // Throw exception
    std::cout<<"qubitIndex out of range"<<std::endl;
    return -1;
  }
}

void qreg::display()
{
  for(long int i = 0; i < a.size(); i++)
  {
    if(std::abs(std::imag(a[i]))<0.01)
      std::cout<<std::real(a[i])<<" ";
    else
      std::cout<<a[i]<<" ";
    std::cout << "|";
    int Currentqubit = i;
    int ToDisplay[HowManyqubits];
    for(int j = 0; j < HowManyqubits; j++)
    {
      ToDisplay[j] = Currentqubit%2;
      Currentqubit = Currentqubit / 2;
    }
    for(int j = HowManyqubits-1; j >=0; j--)
    {
      std::cout<<ToDisplay[j];
    }
    std::cout << ">";
    if(i < a.size() - 1)
      std::cout << " + ";
  }
  std::cout<<std::endl<<std::endl;
}

void qreg::normalize()
{
  std::complex<float> norm = 0;
  for(long int i = 0; i < a.size(); i++)
  {
    norm += a[i] * a[i];
  }
  if(abs(norm)>0)
  {
    for(long int i = 0; i < a.size(); i++)
    {
      a[i] = a[i] / sqrt(norm);
    }  
  }
  else
  {
    // Throw exception
  }
}
