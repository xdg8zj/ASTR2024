
// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

#include <iostream>
using namespace std;

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

Real HistFunc(MeshBlock *pmb, int iout);

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real gm = pin->GetOrAddReal("problem","GM",0.0);
  Real pressb   = pin->GetOrAddReal("problem", "pressb", 1.e6);
  Real rhob   = pin->GetOrAddReal("problem", "rhob", 1.e-4);
  Real gamma = peos->GetGamma();

  Real rb = pcoord->x1f(0);

  // Real K = pressb / pow(rhob,gamma);
  // Real Kcr = gm*(gamma-1.0)/(gamma*rb*pow(rhob,gamma-1.0)) ;
  // Real Kratio = K/Kcr;
  // Real rratio;
  // cout << "K= " << K << endl;
  // cout << "K= " << K << endl;
  // cout << "Kcr= " << Kcr << endl;
  // cout << "K/Kcr= " << Kratio << endl;

  Real asq = pressb/rhob;
  cout << "sqrt(GM/a^2 rb)= " << sqrt(gm/asq/rb) << endl;

  Real rho, press;
  Real r;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
        // rratio = r/rb;
        // rho = rhob * pow( (1.0/rratio + Kratio - 1.0)/Kratio, 1.0/(gamma-1.0) );
        // rho = pow( (gm*(gamma-1.0))/(K*gamma*r), 1.0/(gamma-1.0) );
        // press = K * pow(rho,gamma);
        rho = rhob * exp( (gm/asq/rb)*(rb/r-1.0) );
        press = asq * rho;
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = press/(gamma-1.0);
      }
    }
  }

  Real mass,Qxx,Qxy,Qyy,Qxz,Qyz,Qzz,dvol,theta,phi,x,y,z;
  Qxx=0.0;
  Qyy=0.0;
  Qxy=0.0;
  Qxz=0.0;
  Qyz=0.0;
  Qzz=0.0;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
	theta = pcoord->x2v(j);
        phi = pcoord->x3v(k);
        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);
        rho = phydro->u(IDN,k,j,i);
        dvol = (pow(pcoord->x1f(i+1),3)-pow(pcoord->x1f(i),3))/3.0 * (cos(pcoord->x2f(j))-cos(pcoord->x2f(j+1))) * (pcoord->x3f(k+1)-pcoord->x3f(k));
        mass = mass + rho * dvol;
        Qxx = Qxx + rho * dvol * pow(x,2);
        Qxy = Qxy + rho * dvol * x*y;
        Qyy = Qyy + rho * dvol * pow(y,2);
        Qxz = Qxz + rho * dvol * x*z;
        Qyz = Qyz + rho * dvol * y*z;
        Qzz = Qzz + rho * dvol * pow(z,2);
      }
    }
  }

  cout << "mass= " << mass << endl;
  cout << "Qxx= " << Qxx << endl;
  cout << "Qxy= " << Qxy << endl;
  cout << "Qyy= " << Qyy << endl;
  cout << "Qxz= " << Qxz << endl;
  cout << "Qyz= " << Qyz << endl;
  cout << "Qzz= " << Qzz << endl;

}

void Mesh::InitUserMeshData(ParameterInput *pin)
    {
      AllocateUserHistoryOutput(10);
      EnrollUserHistoryOutput(0, HistFunc, "mass");
      EnrollUserHistoryOutput(1, HistFunc, "Qxx");
      EnrollUserHistoryOutput(2, HistFunc, "Qxy");
      EnrollUserHistoryOutput(3, HistFunc, "Qxz");
      EnrollUserHistoryOutput(4, HistFunc, "Qyx");
      EnrollUserHistoryOutput(5, HistFunc, "Qyy");
      EnrollUserHistoryOutput(6, HistFunc, "Qyz");
      EnrollUserHistoryOutput(7, HistFunc, "Qzx");
      EnrollUserHistoryOutput(8, HistFunc, "Qzy");
      EnrollUserHistoryOutput(9, HistFunc, "Qzz");
      return;
    }

Real HistFunc(MeshBlock *pmb, int iout)
{
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real mass,Qxx,Qxy,Qxz,Qyx,Qyy,Qyz,Qzx,Qzy,Qzz;
  mass=0.0;
  Qxx=0.0;
  Qxy=0.0;
  Qxz=0.0;
  Qyx=0.0;
  Qyy=0.0;
  Qyz=0.0;
  Qzx=0.0;
  Qzy=0.0;
  Qzz=0.0;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real r = pmb->pcoord->x1v(i);
        Real theta = pmb->pcoord->x2v(j);
        Real phi = pmb->pcoord->x3v(k);
        Real x = r*sin(theta)*cos(phi);
        Real y = r*sin(theta)*sin(phi);
        Real z = r*cos(theta);
        Real rho = pmb->phydro->u(IDN,k,j,i);
        Real dvol = (pow(pmb->pcoord->x1f(i+1),3)-pow(pmb->pcoord->x1f(i),3))/3.0 \
        * (cos(pmb->pcoord->x2f(j))-cos(pmb->pcoord->x2f(j+1))) \
        * (pmb->pcoord->x3f(k+1)-pmb->pcoord->x3f(k));
        mass = mass + rho * dvol;
        Qxx = Qxx + rho * dvol * pow(x,2);
        Qxy = Qxy + rho * dvol * x*y;
        Qxz = Qxz + rho * dvol * x*z;
        Qyx = Qyx + rho * dvol * y*x;
        Qyy = Qyy + rho * dvol * pow(y,2);
        Qyz = Qyz + rho * dvol * y*z;
        Qzx = Qzx + rho * dvol * z*x;
        Qzy = Qzy + rho * dvol * z*y;
        Qzz = Qzz + rho * dvol * pow(z,2);
      }
    }
  }
  if (iout==0){
    return mass;
  } else if (iout==1){
    return Qxx;
  } else if (iout==2){
    return Qxy;
  } else if (iout==3){
    return Qxz;
  } else if (iout==4){
    return Qyx;
  } else if (iout==5){
    return Qyy;
  } else if (iout==6){
    return Qyz;
  } else if (iout==7){
    return Qzx;
  } else if (iout==8){
    return Qzy;
  } else if (iout==9){
    return Qzz;
  } else {
    return 0.0;
  }
}

