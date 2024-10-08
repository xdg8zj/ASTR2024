
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

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real gm = pin->GetOrAddReal("problem","GM",0.0);
  Real press0   = pin->GetOrAddReal("problem", "press0", 1.e6);
  Real rho0   = pin->GetOrAddReal("problem", "rho0", 1.e-4);
  Real gamma = peos->GetGamma();
  Real rotation_period = pin->GetOrAddReal("problem", "rotation_period", 8.64e4);
//  Real rb = 7.0e9
//  Real K = press0/pow(rho0,gamma)
//  Real Kcrit = gm*(gamma-1)/(gamma*pow(rho0,gamma-1))
  

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real r = pcoord->x1v(i);
          // rratio = r/rb;
          // rho = rhob * pow( (1.0/rratio + Kratio - 1.0)/Kratio, 1.0/(gamma-1.0) );
          // rho = pow( (gm*(gamma-1.0))/(K*gamma*r), 1.0/(gamma-1.0) );
          // press = K * pow(rho,gamma);
//        Real rho = rhob * exp( (gm/asq/rb)*(rb/r-1.0) );
//        Real press = asq * rho;
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
          
          
    
    

    
  // Calc quadrapole moments
    cout << "mass =" << total_mass << "\n";
    
    Real exact_mass = 4.0*M_PI/3.0*((pow(pcoord-> x1f(ie+1),3))-(pow(pcoord-> x1f(is),3)))*rho0;
    cout << "Exact mass = " << exact_mass << "\n";
    
    cout << "Qxx = " << Qxx << "\n";
//    
//    cout << "Qxy = " << Qxy << "\n";
//    
//    cout << "Qxz = " << Qxz << "\n";
    
    cout << "Qyy = " << Qyy << "\n" ;
//    
//    cout << "Qyx = " << Qyx << "\n";
//    
//    cout << "Qyz = " << Qyz << "\n";
    
    cout << "Qzz = " << Qzz << "\n";
    
//    cout << "Qzx = " << Qzx << "\n";
//    
//    cout << "Qzy = " << Qzy << "\n";
  
    
    cout << "Spin(z) = " << Spin_z << "\n";

    
    
    
}

