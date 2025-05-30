/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * 2023: Cheng Li, University of Michigan
 * 2025: Francisco Spaulding-Astudillo, UCLA
 * Contact: chengcli@umich.edu, fspauldinga@ucla.edu
 * Reference: Bryan and Fritsch, 2002
 * -------------------------------------------------------------------------------------
 */

// athena
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>

// climath
#include <climath/core.h>  // sqr

#include <climath/root.hpp>

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// basic calculations
#include <cmath>

// Global variables
int iH2O, iH2Oc;
Real p0, grav;
Real btau, btem;
Real bu1, bu2, bu3;
Real dTdt_body;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(12);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "tempv");
  SetUserOutputVariableName(2, "enthalpy");
  SetUserOutputVariableName(3, "entropy");
  SetUserOutputVariableName(4, "intEng");

  SetUserOutputVariableName(5, "theta");
  SetUserOutputVariableName(6, "theta_v");
  SetUserOutputVariableName(7, "mse");

  SetUserOutputVariableName(8, "rh_H2O");
  SetUserOutputVariableName(9, "theta_e");
  SetUserOutputVariableName(10, "qtol");
 
  SetUserOutputVariableName(11, "sqrt(u^2+v^2)");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
        user_out_var(1, k, j, i) =
            user_out_var(0, k, j, i) * pthermo->RovRd(w.at(k, j, i));
        user_out_var(2, k, j, i) = pthermo->GetEnthalpy(w.at(k, j, i));
        user_out_var(3, k, j, i) = pthermo->GetEntropy(w.at(k, j, i));
        user_out_var(4, k, j, i) = pthermo->GetInternalEnergy(w.at(k, j, i));

        // theta
        user_out_var(5, k, j, i) = potential_temp(pthermo, w.at(k, j, i), p0);
        // theta_v
        user_out_var(6, k, j, i) =
            user_out_var(5, k, j, i) * pthermo->RovRd(w.at(k, j, i));

        // mse
        user_out_var(7, k, j, i) =
            moist_static_energy(pthermo, w.at(k, j, i), grav * pcoord->x1v(i));

        // rh
        user_out_var(8, k, j, i) = relative_humidity(pthermo, w.at(k, j, i))[iH2O];
        // theta_e
        user_out_var(9, k, j, i) = equivalent_potential_temp(
            pthermo, w.at(k, j, i), user_out_var(8, k, j, i), p0);

        // total mixing ratio
        user_out_var(10, k, j, i) = w(iH2O, k, j, i) + w(iH2Oc, k, j, i);

        // wind speed
        Real U = phydro->w(IVX,k,j,i);
        Real V = phydro->w(IVY,k,j,i);
        user_out_var(11, k, j, i) = std::sqrt(U*U + V*V);

      }
}

// Surface forcing function (Newtonian relaxation)
// The first step is to define a function that takes a pointer to
// the MeshBlock as an argument such that we can access all physics via this
// pointer. The name of this function is not important. It can be anything. But
// the order and types of the arguments must be <code>(MeshBlock *, Real const,
// Real const, AthenaArray<Real> const&, AthenaArray<Real> const&,
// AthenaArray<Real> &)</code>. They are called the signature of the function.
// src/forcing has examples of relaxation schemes
void RelaxBotTemp(MeshBlock *pmb, Real const time, Real const dt,
               AthenaArray<Real> const &w, AthenaArray<Real> const &r,
               AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
               AthenaArray<Real> &s) {
  // <code>pcoord</code> is a pointer to the Coordinates class and it is a
  // member of the MeshBlock class. We use the pointer to the MeshBlock class,
  // <code>pmb</code> to access <code>pcoord</code> and use its member function
  // to get the spacing of the grid.
  Real dx = pmb->pcoord->dx1f(pmb->is);
  Real dy = pmb->pcoord->dx2f(pmb->js);

  // Similarly, we use <code>pmb</code> to find the pointer to the
  // Thermodynamics class, <code>pthermo</code>.
  auto pthermo = Thermodynamics::GetInstance();

  // Loop over all grids at the "surface"
  // solve for u(prognostic) using w(diagnostic)
  int is = pmb->is, js = pmb->js, ks = pmb->ks; // bot vertical index (is)
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      Real temp = pthermo->GetTemp(w.at(k,j,is)); // (pmb,k,j,i) or (w.at(k,j,i))
      Real dTdt = -(temp-btem)/btau; // (units of K/s)
      Real cv   = pthermo->GetCv(w.at(k,j,is)); // src/snap/thermodynamics/ (J/kg/K)
      Real rho  = w(IDN,k,j,is);
      // update the internal energy
      u(IEN,k,j,is)  += dt*rho*cv*dTdt; // J/m3
    }
  }
}

void RelaxBotVelo(MeshBlock *pmb, Real const time, Real const dt,
               AthenaArray<Real> const &w, AthenaArray<Real> const &r,
               AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
               AthenaArray<Real> &s) {
  // Loop over all grids at the bottom boundary
  // solve for u(prognostic) using w(diagnostic)
  int is = pmb->is, js = pmb->js, ks = pmb->ks; // bot vertical index (is)
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      Real rho  = w(IDN,k,j,is);
      Real u1   = w(IVX,k,j,is);
      Real u2   = w(IVY,k,j,is);
      Real u3   = w(IVZ,k,j,is);
      // apply dissipation as a sink to the momentum budget: d(rho*u)
      u(IM1,k,j,is) += dt*rho*(bu1-u1)/btau;
      u(IM2,k,j,is) += dt*rho*(bu2-u2)/btau;
      u(IM3,k,j,is) += dt*rho*(bu3-u3)/btau;
    }
  }
}

void RelaxBotComp(MeshBlock *pmb, Real const time, Real const dt,
               AthenaArray<Real> const &w, AthenaArray<Real> const &r,
               AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
               AthenaArray<Real> &s) {
  // Similarly, we use <code>pmb</code> to find the pointer to the
  // Thermodynamics class, <code>pthermo</code>.
  auto pthermo = Thermodynamics::GetInstance();

  // Loop over all grids at the "surface"
  // solve for u(prognostic) using w(diagnostic)
  int is = pmb->is, js = pmb->js, ks = pmb->ks; // bot vertical index (is)
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  // specific heat of dry air
  Real gammad = pthermo->GetGammad();
  Real cvd = pthermo->GetRd() / (gammad - 1.);

  // Loop over all vapors in the atmosphere
  for (int n=0; n<NVAPOR; ++n) { 
    if (n!=iH2O) continue; // only adjust water vapor composition
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        // perturb the atmosphere towards saturation at the current T
        Real rh   = relative_humidity(pthermo,w.at(k,j,is))[n];
        // Real rh   = pthermo->RelativeHumidity(pmb,n,k,j,is); //pmb is a pointer to MeshBlock
        Real rho  = w(IDN,k,j,is); 
        Real temp = pthermo->GetTemp(w.at(k,j,is));
        // retreive the vapor mass fraction (or use AirParcelHelper)
        Real qv = w(n,k,j,is); 
        // solve for the saturation mass fraction from known quantities
        Real qsat = qv/rh;
        // relax to saturation over characteristic timescale (kg/kg/s)
        Real dqdt = (qsat-qv)/btau;
        // total change in mixing ratio over model timestep
        Real dq   = dt*dqdt;
        // Step 1: change in density (check 1-qv)
        u(IDN,k,j,is) += dq*rho;
        // Step 2: change in momentum
        Real u1 = w(IVX,k,j,is);
        Real u2 = w(IVY,k,j,is);
        Real u3 = w(IVZ,k,j,is);        

        u(IM1,k,j,is) += dq*u1*rho;
        u(IM2,k,j,is) += dq*u2*rho;
        u(IM3,k,j,is) += dq*u3*rho;
        // Step 3: change in total energy
        // effective specific heat of air
        Real cv   = pthermo->GetCv(w.at(k,j,is)); // src/snap/thermodynamics/ (J/kg/K)
        // get c_{v,v}/c_{v,d}
        Real cv_ratio = pthermo->GetCvRatio(n); // check with Cheng -> how do I get cvv??
        Real dEN = u(IEN,k,j,is)*dq + rho*(cvd*cv_ratio-cv)*temp*(dq/(1+dq));
        u(IEN,k,j,is) += dEN;
      }
    }
  }
}

void BodyHeating(MeshBlock *pmb, Real const time, Real const dt,
               AthenaArray<Real> const &w, AthenaArray<Real> const &r,
               AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
               AthenaArray<Real> &s) {
  // Similarly, we use <code>pmb</code> to find the pointer to the
  // Thermodynamics class, <code>pthermo</code>.
  auto pthermo = Thermodynamics::GetInstance();

  // Loop over all grids
  int is = pmb->is, js = pmb->js, ks = pmb->ks; 
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i){
        //if (w(IPR,k,j,i)<1.E4 || i==is) continue; // skip the bottommost layer 
        if (w(IPR,k,j,i)<1.E4 || i==is) continue; // skip the bottommost layer
        Real cv   = pthermo->GetCv(w.at(k,j,is)); // src/snap/thermodynamics/ (J/kg/K)
        Real rho  = w(IDN,k,j,is);
        // update the total energy (J/m3)
        u(IEN,k,j,i) += dt*rho*cv*dTdt_body; 
      }
    }
  }
}

// Initialize surface pressure from input file.
// This is the place where program specific variables are initialized.
// Note that the function is a member function of the Mesh class rather than the
// MeshBlock class we have been working with. The difference between class Mesh
// and class MeshBlock is that class Mesh is an all-encompassing class that
// manages multiple MeshBlocks while class MeshBlock manages all physics
// modules. During the instantiation of the classes. class Mesh is instantiated
// first and then it instantiates all MeshBlocks inside it. Therefore, this
// subroutine runs before any MeshBlock.
void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  grav = -pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");

  // index
  iH2O = pthermo->SpeciesIndex("H2O");
  iH2Oc = pthermo->SpeciesIndex("H2O(l)");
    
  // The program specific forcing parameters are set here.
  // They come from the input file, which is parased by the ParameterInput class.
  Real Ts   = pin->GetReal("problem", "Ts");
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd; // NOTE: does not equal cp at runtime!!
  Real cv = Rd/(gamma-1.); // NOTE: does not equal cv at runtime!!

  // parameters associated with thermal relaxation scheme at bottom boundary
  btau = pin->GetReal("forcing", "btau"); // relaxation timescale
  btem = pin->GetReal("forcing", "btem"); // relaxation temperature

  // parameters associated with velocity relaxation at bottom boundary
  bu1   = pin->GetReal("forcing", "bu1"); // relaxation velocity in x-dir
  bu2   = pin->GetReal("forcing", "bu2"); // relaxation velocity in y-dir
  bu3   = pin->GetReal("forcing", "bu3"); // relaxation velocity in z-dir

  // parameters associated with bulk tropospheric heating/cooling
  dTdt_body = pin->GetReal("forcing", "dTdt_body"); // (K/s)

  // This line code enrolls any forcing functions
  //EnrollUserExplicitSourceFunction(RelaxBotTemp);
  //EnrollUserExplicitSourceFunction(RelaxBotVelo);
  //EnrollUserExplicitSourceFunction(RelaxBotComp);
  //EnrollUserExplicitSourceFunction(BodyHeating);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd;  
  Real grav = -phydro->hsrc.GetG1();

  Real Ps = p0;
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");

  // construct a reversible adiabat
  std::vector<Real> yfrac({1. - qt, qt, 0.});
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pthermo->SetMassFractions<Real>(yfrac.data());
      pthermo->EquilibrateTP(Ts, Ps);

      // half a grid to cell center
      pthermo->Extrapolate_inplace(pcoord->dx1f(is) / 2., "reversible", grav);

      for (int i = is; i <= ie; ++i) {
        pthermo->GetPrimitive(w.at(k, j, i));
        pthermo->Extrapolate_inplace(pcoord->dx1f(i), "reversible", grav);
      }
    }

  // add temperature anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real L = sqrt(sqr((x2 - xc) / xr) + sqr((x1 - zc) / zr));

        if (L < 1.) {
          pthermo->SetPrimitive(w.at(k, j, i));
          Real temp = pthermo->GetTemp();
          Real pres = pthermo->GetPres();

          Real temp_v = temp * pthermo->RovRd();
          temp_v *= 1. + dT * sqr(cos(M_PI * L / 2.)) / 300.;

          int err = root(temp, temp + dT * Ts / 300., 1.E-8, &temp,
                         [&pthermo, temp_v, pres](Real temp) {
                           pthermo->EquilibrateTP(temp, pres);
                           return temp * pthermo->RovRd() - temp_v;
                         });

          //   if (err) throw RuntimeError("pgen", "TVSolver doesn't converge");

          pthermo->SetPrimitive(w.at(k, j, i));
          pthermo->EquilibrateTP(temp, pres);
          pthermo->GetPrimitive(w.at(k, j, i));
        }
      }

  peos->PrimitiveToConserved(w, pfield->bcc, phydro->u, pcoord, is, ie, js, je,
                             ks, ke);
}
