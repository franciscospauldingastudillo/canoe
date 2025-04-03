/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 * Francisco Spaulding-Astudillo, UCLA
 *
 * Year: 2025
 * Contact: chengcli@umich.edu, fspauldinga@gmail.com
 * Reference: Robert et al., 1992
 * -------------------------------------------------------------------------------------
 */

// @sect3{Include files}

// These input files are just like those in the @ref straka, so additional
// comments are not required.
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <index_map.hpp>

// canoe
#include <impl.hpp>

// snap
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// calculations
#include <cmath>

// @sect3{Preamble}

// We need 3 global variables here
// for communication between InitUserMeshData and forcing functions
Real p0;
Real btau, btem;
Real bu1, bu2, bu3;
int iH2O;
Real dTdt_body;

// Same as that in @ref straka, make outputs of temperature and potential
// temperature.
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "rh_H2O");
  SetUserOutputVariableName(3, "sqrt(u^2+v^2)");
}

// Set temperature and potential temperature.
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;
  Real gamma = peos->GetGamma();
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k,j,i));
        user_out_var(1, k, j, i) = potential_temp(pthermo,w.at(k,j,i), p0);
        user_out_var(2, k, j, i) = relative_humidity(pthermo,w.at(k,j,i))[iH2O];
        Real U = phydro->w(IVX,k,j,i);
        Real V = phydro->w(IVY,k,j,i);
        user_out_var(3, k, j, i) = std::sqrt(U*U + V*V); 
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


// In subsequent versions, implement a viscous dissipation and heating at the surface
// horizontal momentum dissipation is dt*(K*rho*lap(IVX,IVY))
// mom ~ (kg/m2 s), dmom/dt ~ (kg/m2s2), so dt*(dmom/dt)~flux
// horizontal dissipative heating  is du = dt*flux*Area_horiz/Vol

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
  // The program specific forcing parameters are set here.
  // They come from the input file, which is parased by the ParameterInput class.
  Real p0   = pin->GetReal("problem", "p0"); 
  Real Ts   = pin->GetReal("problem", "Ts");
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd; // NOTE: does not equal cp at runtime!!
  Real cv = Rd/(gamma-1.); // NOTE: does not equal cv at runtime!!

  // determine the index of water vapor
  auto pindex = IndexMap::GetInstance();
  iH2O = pindex->GetVaporId("H2O"); // from earth.yaml

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

// @sect3{Initial condition}

// We do not need forcings other than gravity in this problem,
// so we go directly to the initial condition.
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Similar to @ref straka, read variables in the input file
  Real gamma = pin->GetReal("hydro", "gamma");
  Real grav = -phydro->hsrc.GetG1();
  Real Ts = pin->GetReal("problem", "Ts");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd;

  Real xc = pin->GetReal("problem", "xc");
  Real yc = pin->GetReal("problem", "yc");
  Real zc = pin->GetReal("problem", "zc");
  Real s = pin->GetReal("problem", "s");
  Real a = pin->GetReal("problem", "a");
  Real dT = pin->GetReal("problem", "dT");

  // Whether to do a uniform bubble or not.
  bool uniform_bubble =
      pin->GetOrAddBoolean("problem", "uniform_bubble", false);

  // Loop over the grids and set initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real r = sqrt((x3 - yc) * (x3 - yc) + (x2 - xc) * (x2 - xc) +
                      (x1 - zc) * (x1 - zc));
        Real temp = Ts - grav * x1 / cp;
        phydro->w(IPR, k, j, i) = p0 * pow(temp / Ts, cp / Rd);
        if (r <= a)
          temp += dT * pow(phydro->w(IPR, k, j, i) / p0, Rd / cp);
        else if (!uniform_bubble)
          temp += dT * exp(-(r - a) * (r - a) / (s * s)) *
                  pow(phydro->w(IPR, k, j, i) / p0, Rd / cp);
        phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * temp);
        phydro->w(IVX, k, j, i) = phydro->w(IVY, k, j, i) = 0.;
      }

  // Change primitive variables to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
