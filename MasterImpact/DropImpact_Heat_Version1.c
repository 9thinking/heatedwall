// Author: Radu Cimpeanu (base)
// Modified by: ChatGPT
// Date: 2026-03-09
//
// Goal:
// - Add heat transfer (advection-diffusion of temperature) to droplet impact case
// - Set hot wall at 300 C (573.15 K)
// - No phase change in this version (pure heat transfer)
//
// Notes:
// 1) This is a practical "first working version" for Basilisk.
// 2) Temperature is solved as a tracer with explicit diffusion.
// 3) Thermal properties are kept constant per phase and blended by VOF fraction.
// 4) If your local Basilisk setup uses slightly different diffusion API, adjust the
//    diffusion call inside event thermal_diffusion accordingly.

#include "axi.h"                     // axisymmetric geometry
#include "navier-stokes/centered.h"  // solve NS equations
#define FILTERED 1
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#include "vof.h"
#include "fractions.h"
#include "view.h"
#include "draw.h"
#include "tag.h"
#include "tracer.h"                  // tracer advection support
#include "diffusion.h"               // diffusion solver for scalar

// -------------------- dimensional quantities (input) --------------------
double rhoLiquid; // kg/m^3
double rhoGas;    // kg/m^3
double muLiquid;  // kg/(m s)
double muGas;     // kg/(m s)
double sig;       // N/m
double g_accel;   // m/s^2
double dRadius;   // m
double v_init;    // m/s

// -------------------- dimensionless groups ------------------------------
#define rho_ratio   (rhoGas/rhoLiquid)
#define mu_ratio    (muGas/muLiquid)

#define poolHeight 0.0
#define domainSize 8.0

face vector av[];

FILE * fp_stats;
FILE * fp_droplets;
FILE * fp_thermal;

double ND_Weber;
double ND_Reynolds;
double ND_Froude;
double ND_Bond;
double ND_Ohnesorge;

double filmHeight;

int minLevel = 4;
int maxLevel;
double tEnd;

// -------------------- thermal model -------------------------------------
// Temperature scalar (dimensionless):
// Theta = (T - T_ref)/(T_wall - T_ref), so wall Theta=1.
scalar T[];

// Physical temperatures (for output only)
double T_ref   = 293.15; // K (20 C)
double T_wallK = 573.15; // K (300 C)

// Thermal diffusivities [m^2/s] (typical values, edit as needed)
double alphaLiquid = 1.4e-7; // water-like order
double alphaGas    = 2.0e-5; // air-like order

// Non-dimensional thermal diffusivities based on:
// x* = x/R, u* = u/U, t* = t U/R  -> alpha* = alpha/(U R)
double alphaL_nd, alphaG_nd;

// Cell-centered blended diffusivity and face counterpart
scalar alphaC[];
face vector alphaF[];

// -------------------- BCs -----------------------------------------------
// Bottom of domain = LEFT boundary in this setup
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
f[left]   = dirichlet(0.);

// Lateral
u.n[top]  = dirichlet(0.);
u.t[top]  = neumann(0.);

// Top/outlet
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// Thermal BCs (dimensionless Theta)
// hot wall: Theta=1
T[left]  = dirichlet(1.0);
// outlet / lateral: zero normal gradient
T[right] = neumann(0.);
T[top]   = neumann(0.);

// tracer list so T gets advected with flow
scalar * tracers = {T};

int main(int argc, char * argv[]) {

  if (argc < 11) {
    fprintf(stderr,
      "Usage:\n"
      "./a.out rhoL rhoG muL muG sigma g R U tEnd maxLevel\n");
    return 1;
  }

  rhoLiquid = atof(argv[1]);
  rhoGas    = atof(argv[2]);
  muLiquid  = atof(argv[3]);
  muGas     = atof(argv[4]);
  sig       = atof(argv[5]);
  g_accel   = atof(argv[6]);
  dRadius   = atof(argv[7]);
  v_init    = atof(argv[8]);
  tEnd      = atof(argv[9]);
  maxLevel  = atoi(argv[10]);

  ND_Weber     = (rhoLiquid*sq(v_init)*dRadius)/sig;
  ND_Reynolds  = (rhoLiquid*v_init*dRadius)/muLiquid;
  ND_Froude    = v_init/sqrt(dRadius*g_accel);
  ND_Bond      = rhoLiquid*g_accel*sq(dRadius)/sig;
  ND_Ohnesorge = muLiquid/sqrt(rhoLiquid*sig*dRadius);

  // thermal non-dimensionalization
  alphaL_nd = alphaLiquid/(v_init*dRadius);
  alphaG_nd = alphaGas/(v_init*dRadius);

  init_grid(1 << 6);
  size(domainSize);
  origin(-0.5*domainSize, 0.0);

  mkdir("Slices", 0700);
  mkdir("Animations", 0700);
  mkdir("Interfaces", 0700);

  fprintf(stdout, "Reynolds number   = %0.6f\n", ND_Reynolds);
  fprintf(stdout, "Weber number      = %0.6f\n", ND_Weber);
  fprintf(stdout, "Froude number     = %0.6f\n", ND_Froude);
  fprintf(stdout, "Bond number       = %0.6f\n", ND_Bond);
  fprintf(stdout, "Ohnesorge number  = %0.6f\n", ND_Ohnesorge);
  fprintf(stdout, "alphaL*           = %0.6e\n", alphaL_nd);
  fprintf(stdout, "alphaG*           = %0.6e\n", alphaG_nd);
  fflush(stdout);

  rho1 = 1.;
  rho2 = rho_ratio;

  mu1 = 1./ND_Reynolds;
  mu2 = mu_ratio*mu1;

  f.sigma = 1./ND_Weber;
  a = av;

  fp_stats    = fopen("logstats.dat", "w");
  fp_droplets = fopen("logdroplets.dat", "w");
  fp_thermal  = fopen("logthermal.dat", "w");

  DT = 5e-4;         // smaller than original due to thermal diffusion stability
  NITERMIN = 1;
  NITERMAX = 300;
  TOLERANCE = 1e-4;

  run();

  fclose(fp_stats);
  fclose(fp_droplets);
  fclose(fp_thermal);
}

// gravity
event acceleration (i++) {
  foreach_face(x)
    av.x[] -= 1./sq(ND_Froude);
  foreach_face(y)
    av.y[] += 0.0;
}

scalar omega[], viewingfield[], mylevel[];

event init (t = 0.0) {

  filmHeight = -domainSize/2. + poolHeight;

  refine (((sq(x - (filmHeight + 2.0)) + sq(y) < sq(1.0*1.05) &&
            sq(x - (filmHeight + 2.0)) + sq(y) > sq(1.0*0.95)) ||
           fabs(x - filmHeight) <= 0.005) && level < maxLevel);

  // liquid = union of drop + pool film
  fraction (f, sq(1.0) - sq(x - (filmHeight + 2.0)) - sq(y));

  foreach() {
    u.x[] = -1.0*f[];
    u.y[] = 0.0;
    p[]   = 0.0;
    omega[] = 0.0;

    // initial dimensionless temperature: ambient = 0
    T[] = 0.0;
  }
}

// build blended diffusivity each step
event set_thermal_diffusivity (i++) {
  foreach() {
    double ff = clamp(f[], 0., 1.);
    alphaC[] = ff*alphaL_nd + (1. - ff)*alphaG_nd;      
  }
  boundary ({alphaC});

  foreach_face() {
    alphaF.x[] = (alphaC[] + alphaC[-1,0])/2.;
  }
  boundary ((scalar *){alphaF});
}

// explicit diffusion for temperature
event thermal_diffusion (i++) {
  // diffusion(dt, scalar, face diffusivity)
  diffusion (T, dt, alphaF);
}

event adapt (i++) {
  // include T in adaptation to keep thermal boundary layer resolved
  adapt_wavelet ((scalar *){f, u.x, u.y, T},
                 (double[]){1e-6, 4e-3, 4e-3, 2e-3},
                 maxLevel, minLevel);
}

event gfsview (t = 0.0; t += 1.0; t <= tEnd) {
  char name_gfs[200];
  sprintf(name_gfs,"Slices/DropImpact-%0.1f.gfs",t);
  FILE* fp_gfs = fopen (name_gfs, "w");
  output_gfs(fp_gfs);
  fclose(fp_gfs);
}

event saveInterfaces (t += 0.01) {
  char nameInterfaces1[200];
  sprintf(nameInterfaces1,"Interfaces/interfaceDrop-%0.2f.dat",t);
  FILE * fp1 = fopen(nameInterfaces1, "w");
  output_facets (f, fp1);
  fclose(fp1);
}

event small_droplet_removal (i++) {
  remove_droplets(f, 8);
  remove_droplets(f, 8, true);
}

event droplets (t += 0.01) {
  scalar m[];
  foreach()
    m[] = f[] > 1e-2;
  int n = tag (m);

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;

  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord pnt = {x,y,z};
      foreach_dimension()
        b[j].x += dv()*f[]*pnt.x;
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int j = 0; j < n; j++)
    fprintf (fp_droplets, "%d %g %d %g %g %g\n", i, t,
             j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fp_droplets);
}

// thermal diagnostics:
// 1) average liquid temperature
// 2) wall heat flux proxy q* = -alpha * dT/dn at wall (n = +x direction here)
event thermal_stats (t += 0.01; t <= tEnd) {
  double Vliq = 0., Tint = 0.;
  foreach(reduction(+:Vliq) reduction(+:Tint)) {
    Vliq += dv()*f[];
    Tint += dv()*f[]*T[];
  }

  double Tavg_nd = (Vliq > 1e-14) ? Tint/Vliq : 0.;
  double Tavg_K  = T_ref + Tavg_nd*(T_wallK - T_ref);

  // integrate wall heat flux at left boundary from first interior cells
  double qwall_int = 0.;
  foreach_boundary(left, reduction(+:qwall_int)) {
    // one-sided gradient approx: dT/dx ~ (T[] - Twall)/Delta
    double dTdx = (T[] - 1.0)/Delta;
    // local blended alpha at adjacent cell
    double ff = clamp(f[], 0., 1.);
    double a_loc = ff*alphaL_nd + (1. - ff)*alphaG_nd;
    double qn = -a_loc*dTdx; // dimensionless flux
    qwall_int += qn*Delta;   // line integral in axi slice
  }

  fprintf(fp_thermal, "%d %.8g %.8g %.8g\n", i, t, Tavg_nd, Tavg_K);
  fflush(fp_thermal);

  // optional screen output
  fprintf(stdout, "t=%g Tavg_liq(K)=%.4f qwall_int*=%.6e\n", t, Tavg_K, qwall_int);
  fflush(stdout);
}

event movies (t += 0.01) {
  char timestring[100];

  foreach() {
    viewingfield[] = T[];  // show temperature instead of 1-f
    mylevel[] = level;
  }

  view(width=1200, height=800, fov=30.0, ty = 0.0,
       quat = { 0, 0, -0.707, 0.707 });
  clear();

  draw_vof("f", lw=2);
  squares("viewingfield", map = cool_warm, min = 0.0, max = 1.0);

  mirror({0,1}) {
    draw_vof("f", lw=2);
    cells(lw=0.5);
    squares("mylevel", map = cool_warm, min = minLevel, max = maxLevel);
  }

  sprintf(timestring, "t=%2.02f", t);
  draw_string(timestring, pos=1, lc={0,0,0}, lw=2);

  save ("Animations/ImpactSummary.mp4");
}

event logstats (t += 0.01; t <= tEnd) {
  timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
  fprintf(fp_stats,
          "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g\n",
          i, t, dt, grid->n, perf.t, s.cpu);
  fflush(fp_stats);
}