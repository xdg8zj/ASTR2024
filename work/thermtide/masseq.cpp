
// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

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
    
    void calc_mass(double total_mass){
        Real total_mass = 0;
        for (int k=ks; k<=ke; k++) {
            for (int j=js; j<=je; j++) {
                for (int i=is; i<=ie; i++) {
                    
                    
                    phydro->u(IDN,k,j,i) = rho0;
                    phydro->u(IM1,k,j,i) = 0.0;
                    phydro->u(IM2,k,j,i) = 0.0;
                    phydro->u(IM3,k,j,i) = 0.0;
                    phydro->u(IEN,k,j,i) = press0/(gamma-1.0);
                    Real rad;
                    
                    Real x1rho = 1/3 * (pcoord-> x1f(i+1) - pcoord-> x1f(i));
                    Real x2rho = cos(pcoord-> x2f(j)) - cos(pcoord-> x2f(j));
                    Real x3rho = pcoord-> x3f(k+1) - pcoord-> x3f(k)
                    
                    Real calc_mass = x1c+x2c+x3c
                    
                }
            }
        }
        total_mass = total_mass + calc_mass;
        for (double dt = 0; double dt = 1e3; dt++){
            cout << dt << " " << "mass =" << " " << calc_mass
        }
    }
}
    
    
//    for (int k=ks; k<=ke; k++) {
//
//        for (int j=js; j<=je; j++) {
//
//            for (int i=is; i<=ie; i++) {
//
//
//                phydro->u(IDN,k,j,i) = rho0;
//                phydro->u(IM1,k,j,i) = 0.0;
//                phydro->u(IM2,k,j,i) = 0.0;
//                phydro->u(IM3,k,j,i) = 0.0;
//                phydro->u(IEN,k,j,i) = press0/(gamma-1.0);
//                Real rad;
//
//                Real x1rho = 1/3 * (pcoord-> x1f(i+1) - pcoord-> x1f(i));
//
//                cout << pcoord-> x1f(i) << " " << pcoord-> x1f(i+1) << " " << x1rho
//            }
//        }
//    }
//}
    

