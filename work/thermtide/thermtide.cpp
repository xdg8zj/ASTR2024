
// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
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

//void Cooling(MeshBlock *pmb, const Real itme, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
//  Real gamma = pmb->peos->GetGamma();
//  Real temp_goal = 10.0;
//  Real tau = pin->GetOrAddReal("problem", "rotation_period", 8.64e4);
//    for (int k = pmb->ks; k<= pmb->ke; ++k) {
//      for (int j = pmb->js; j <= pmb->je; ++j{
//        for (int i = pmb->is; i <=pmb->ie; ++i){
//          Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i);
//            if (temp > temp_goal) {
//                cons(IEN,k,j,i) -= dt/tau*prim(IDN,k,j,i)*(temp-temp_goal)/(gamma - 1.0)
//            }
//        }
//    }
//    
//}

//Real HistFunc(MeshBlock *pmb, int iout);

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real gm = pin->GetOrAddReal("problem","GM",0.0);
  Real press0   = pin->GetOrAddReal("problem", "press0", 1.e6);
  Real rho0   = pin->GetOrAddReal("problem", "rho0", 1.e-4);
  Real gamma = peos->GetGamma();
  Real Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",7.27220521664303958333e-5);
  Real tau = pin->GetOrAddReal("problem", "rotation_period", 8.64e4);
//
  Real rb = pcoord-> x1f(0);
//  Real K = press0/pow(rho0,gamma)
//  Real Kcrit = gm*(gamma-1)/(gamma*pow(rho0,gamma-1))
    
  Real asq = press0/rho0;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real r = pcoord->x1v(i);
          // rratio = r/rb;
          // rho = rhob * pow( (1.0/rratio + Kratio - 1.0)/Kratio, 1.0/(gamma-1.0) );
          // rho = pow( (gm*(gamma-1.0))/(K*gamma*r), 1.0/(gamma-1.0) );
          // press = K * pow(rho,gamma);
        Real rho = rho0 * exp( (gm/asq/r)*(rb/r-1.0) );
        Real press = asq * rho;
        phydro->u(IDN,k,j,i) = rho0; //density
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = press0/(gamma-1.0);
        Real rad;


      }
    }
  }
   
  
  //Real exact_mass = 4.0*M_PI/3.0*((pow(pcoord-> x1f(ie+1),3))-(pow(pcoord-> x1f(is),3)))*rho0;
    
  Real Qxx_calculation = (4.0*M_PI/15.0)*rho0*((pow(pcoord-> x1f(ie+1),5))-(pow(pcoord-> x1f(is),5)));
    
  Real total_mass = 0.0;
    
    
    //    Quadripole moments:
  Real Qxx = 0.0;
  Real Qxy = 0.0;
  Real Qxz = 0.0;
  Real Qyy = 0.0;
  Real Qyx = 0.0;
  Real Qyz = 0.0;
  Real Qzz = 0.0;
  Real Qzx = 0.0;
  Real Qzy = 0.0;
  Real Qtotal = 0.0;
  
    
    
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x1volr = 1.0/3.0 * (pow(pcoord-> x1f(i+1),3.0) - pow(pcoord-> x1f(i),3.0));
        Real x2volt = cos(pcoord-> x2f(j)) - cos(pcoord-> x2f(j+1.0));
        Real x3volp = pcoord-> x3f(k+1.0) - pcoord-> x3f(k);
        Real volume = x1volr*x2volt*x3volp;
        Real calc_mass = phydro->u(IDN,k,j,i)*volume;
//        cout << calc_mass << \n;
        total_mass = total_mass + calc_mass;
          
          
          
          
        // Calculating quadripole moments
        
          
        //x cartesian to polar: x = rsin(phi)cos(theta)
        //y cartesian to polar: y = rsin(phi)sin(theta)
        //z cartesian to polar: z = rcos(phi)
        Real x_polar = pcoord-> x1f(i) * sin(pcoord-> x2f(j)) * cos(pcoord-> x3f(k));
        Real y_polar = pcoord-> x1f(i) * sin(pcoord-> x2f(j)) * sin(pcoord-> x3f(k));
        Real z_polar = pcoord-> x1f(i) * cos(pcoord-> x2f(j));
          
          
        // Qxx:
        Real calc_Qxx = calc_mass*pow(x_polar,2.0);
        Qxx = Qxx + calc_Qxx;
          
        // Qxy:
        Real calc_Qxy = calc_mass*x_polar*y_polar;
        Qxy = Qxy + calc_Qxy;
        
       // Qxz
        Real calc_Qxz = calc_mass*x_polar*z_polar;
        Qxz = Qxz + calc_Qxz;
    
       //Qyy:
        Real calc_Qyy = calc_mass*pow(y_polar,2.0);
        Qyy = Qyy + calc_Qyy;
          
      // Qyx
        Real calc_Qyx = calc_mass*y_polar*x_polar;
        Qyx = Qyx + calc_Qyx;
          
      // Qyz
        Real calc_Qyz = calc_mass*y_polar*z_polar;
        Qyz = Qyz + calc_Qyz;
          
      // Qzz
        Real calc_Qzz = calc_mass*pow(z_polar,2.0);
        Qzz = Qzz + calc_Qzz;
          
      // Qzx
        Real calc_Qzx = calc_mass*z_polar*x_polar;
        Qzx = Qzx + calc_Qzx;
          
      //Qzy
        Real calc_Qzy = calc_mass*z_polar*y_polar;
        Qzy = Qzy + calc_Qzy;
            
      }
    }
  }
    
    
  Real Spin_z = 0.0;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x1volr = 1.0/3.0 * (pow(pcoord-> x1f(i+1),3.0) - pow(pcoord-> x1f(i),3.0));
        Real x2volt = cos(pcoord-> x2f(j)) - cos(pcoord-> x2f(j+1.0));
        Real x3volp = pcoord-> x3f(k+1.0) - pcoord-> x3f(k);
        Real volume = x1volr*x2volt*x3volp;
        Real calc_mass = phydro->u(IDN,k,j,i)*volume;
        Real calc_spin = calc_mass*pcoord-> x1f(i)*sin(pcoord-> x2f(j))*x3volp;
        Spin_z = Spin_z+ calc_spin;
           
      }
    }
  }
          
          
    
    

    
  // Calc quadrapole moments 256 nodes
    cout << "mass =" << total_mass << "\n";
    //mass =6.24591e+24
    
    Real exact_mass = 4.0*M_PI/3.0*((pow(pcoord-> x1f(ie+1),3))-(pow(pcoord-> x1f(is),3)))*rho0;
    cout << "Exact mass = " << exact_mass << "\n";
    //Exact mass = 6.24591e+24
    
    cout << "Qxx = " << Qxx << "\n";
    //Qxx = 1.03479e+44
//
    cout << "Qxy = " << Qxy << "\n";
//    
    cout << "Qxz = " << Qxz << "\n";
    
    cout << "Qyy = " << Qyy << "\n" ;
    //Qyy = 1.03479e+44
//
    cout << "Qyx = " << Qyx << "\n";
//    
    cout << "Qyz = " << Qyz << "\n";
    
    cout << "Qzz = " << Qzz << "\n";
    //Qzz = 1.03487e+44
    
    cout << "Qzx = " << Qzx << "\n";
//    
    cout << "Qzy = " << Qzy << "\n";
  
    
    cout << "Spin(z) = " << Spin_z << "\n";
    //Spin(z) = 8.48801e+32

    
    
    
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

