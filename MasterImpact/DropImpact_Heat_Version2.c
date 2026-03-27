#include "axi.h"                     // 轴对称几何
#include "navier-stokes/centered.h"  // NS 方程
#define FILTERED 1                   // 必须给值，避免 #if FILTERED 预处理报错
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#include "vof.h"
#include "fractions.h"
#include "view.h"
#include "draw.h"
#include "tag.h"
#include "tracer.h"                  // 标量对流
#include "diffusion.h"               // 扩散求解

// ======================== 输入的有量纲参数 ========================
double rhoLiquid; // 液体密度 kg/m^3
double rhoGas;    // 气体密度 kg/m^3
double muLiquid;  // 液体动力黏度 kg/(m s)
double muGas;     // 气体动力黏度 kg/(m s)
double sig;       // 表面张力 N/m
double g_accel;   // 重力加速度 m/s^2
double dRadius;   // 液滴半径 m
double v_init;    // 初始速度 m/s

// ======================== 无量纲组合 ========================
#define rho_ratio   (rhoGas/rhoLiquid)
#define mu_ratio    (muGas/muLiquid)

#define poolHeight 0.0
#define domainSize 8.0

face vector av[];

// 文件指针
FILE * fp_stats;
FILE * fp_droplets;
FILE * fp_thermal;

// 常用无量纲数
double ND_Weber;
double ND_Reynolds;
double ND_Froude;
double ND_Bond;
double ND_Ohnesorge;

double filmHeight;
int minLevel = 4;
int maxLevel;
double tEnd;

// ======================== 温度模型 ========================
// 无量纲温度 T：Theta = (T_phys - T_ref)/(T_wall - T_ref)
scalar T[];

// 物理温度（K）
double T_ref   = 293.15; // 20℃
double T_wallK = 573.15; // 300℃

// 热扩散率（m^2/s）
double alphaLiquid = 1.4e-7;
double alphaGas    = 2.0e-5;

// 无量纲热扩散率 alpha* = alpha/(U R)
double alphaL_nd, alphaG_nd;

// 扩散系数场
scalar alphaC[];
face vector alphaF[];

// ======================== 变物性模型参数 ========================
// 参考温度
double T0_refK = 293.15;

// 液体密度线性膨胀 rho_l(T)=rho_l0*(1-beta*(T-T0))
double rhoL0_dim  = 998.2;
double betaL_th   = 2.57e-4;

// 液体黏度 Andrade/Arrhenius: mu_l = A*exp(B/T)
double A_muL      = 2.414e-5;
double B_muL      = 247.8;

// 表面张力线性：sigma = sigma0 - k*(T-T0)
double sigma0_dim = 0.0728;
double k_sigma    = 1.5e-4;

// 气体密度：理想气体近似 rho_g(T)=rho_g0*T0/T
double rhoG0_dim  = 1.204;

// 气体黏度 Sutherland
double muG0_dim   = 1.81e-5;
double S_suth     = 110.4;

// 稳定性截断
double TminK = 273.15, TmaxK = 900.0;
double muMin_dim = 1e-6, muMax_dim = 5e-2;
double sigMin_dim = 1e-4;

// 局部物性场
scalar TLK[];
scalar rhoL_loc[], rhoG_loc[];
scalar muL_loc[],  muG_loc[];
scalar sigma_loc[];

// NS 使用的“面黏度”场（无量纲动力黏度）
face vector muv[];

// ======================== 边界条件 ========================
// 左边界：壁面（无滑移，不可穿透）
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
f[left]   = neumann(0.);//这里设置为0，表示壁面不含液相（即完全干燥），如果需要考虑液膜润湿，可以改为 Neumann 或 Dirichlet 非零值。

// 上边界：侧边（不可穿透，切向滑移）
u.n[top]  = dirichlet(0.);
u.t[top]  = neumann(0.);

// 右边界：开口/出流
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// 温度边界（无量纲）
T[left]  = dirichlet(1.0); // 热壁 300℃
T[right] = neumann(0.);    // 零梯度
T[top]   = neumann(0.);

// 把温度加入 tracer，对流项由框架自动处理
scalar * tracers = {T};

// ======================== 工具函数 ========================
static inline double clamp01 (double v) {
  return clamp(v, 0., 1.);
}

static inline double T_to_K (double Tnd) {
  return T_ref + Tnd*(T_wallK - T_ref);
}

static inline double alpha_blend_nd (double ff) {
  return ff*alphaL_nd + (1. - ff)*alphaG_nd;
}

static inline double mu_to_nd (double mu_dim) {
  return mu_dim/(rhoL0_dim*v_init*dRadius);
}

static inline double sigma_to_nd (double sigma_dim) {
  return sigma_dim/(rhoL0_dim*sq(v_init)*dRadius);
}

int main(int argc, char * argv[]) {

  // 参数：rhoL rhoG muL muG sigma g R U tEnd maxLevel [TwallK可选]
  if (argc < 11) {
    fprintf(stderr,
      "用法:\n"
      "./a.out rhoL rhoG muL muG sigma g R U tEnd maxLevel [TwallK]\n");
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

  if (argc >= 12)
    T_wallK = atof(argv[11]); // 可选输入壁温

  ND_Weber     = (rhoLiquid*sq(v_init)*dRadius)/sig;
  ND_Reynolds  = (rhoLiquid*v_init*dRadius)/muLiquid;
  ND_Froude    = v_init/sqrt(dRadius*g_accel);
  ND_Bond      = rhoLiquid*g_accel*sq(dRadius)/sig;
  ND_Ohnesorge = muLiquid/sqrt(rhoLiquid*sig*dRadius);

  alphaL_nd = alphaLiquid/(v_init*dRadius);
  alphaG_nd = alphaGas/(v_init*dRadius);

  init_grid(1 << 6);
  size(domainSize);
  origin(-0.5*domainSize, 0.0);

  mkdir("Slices", 0700);
  mkdir("Animations", 0700);
  mkdir("Interfaces", 0700);

  fprintf(stdout, "Re = %0.6f\n", ND_Reynolds);
  fprintf(stdout, "We = %0.6f\n", ND_Weber);
  fprintf(stdout, "Fr = %0.6f\n", ND_Froude);
  fprintf(stdout, "Bo = %0.6f\n", ND_Bond);
  fprintf(stdout, "Oh = %0.6f\n", ND_Ohnesorge);
  fprintf(stdout, "alphaL* = %0.6e, alphaG* = %0.6e\n", alphaL_nd, alphaG_nd);
  fprintf(stdout, "T_wall = %.2f K\n", T_wallK);
  fflush(stdout);

  // 两相基准参数（仍保留；后续由 mu=muv 做局部替代）
  rho1 = 1.;
  rho2 = rho_ratio;
  mu1 = 1./ND_Reynolds;
  mu2 = mu_ratio*mu1;

  // 表面张力初值（后面会在 event 中按液滴平均温度更新）
  f.sigma = sigma_to_nd(sigma0_dim);

  a = av;

  fp_stats    = fopen("logstats.dat", "w");
  fp_droplets = fopen("logdroplets.dat", "w");
  fp_thermal  = fopen("logthermal.dat", "w");

  DT = 5e-4;
  NITERMIN = 1;
  NITERMAX = 300;
  TOLERANCE = 1e-4;

  run();

  fclose(fp_stats);
  fclose(fp_droplets);
  fclose(fp_thermal);
  return 0;
}

// ======================== 事件：重力（含简化热浮力） ========================
event acceleration (i++) {
  foreach_face(x) {
    double gnd = 1./sq(ND_Froude);

    // 使用局部温度做简单 Boussinesq 修正（液相主导）
    double Tnd_f = (T[] + T[-1,0])/2.;
    double Tk_f  = clamp(T_to_K(Tnd_f), TminK, TmaxK);
    double buoy  = betaL_th*(Tk_f - T0_refK);

    av.x[] = -gnd*(1.0 - buoy);
  }

  foreach_face(y)
    av.y[] = 0.0;
}

scalar omega[], viewingfield[], mylevel[];

// ======================== 事件：初始条件 ========================
event init (t = 0.0) {

  filmHeight = -domainSize/2. + poolHeight;

  refine (((sq(x - (filmHeight + 2.0)) + sq(y) < sq(1.0*1.05) &&
            sq(x - (filmHeight + 2.0)) + sq(y) > sq(1.0*0.95)) ||
           fabs(x - filmHeight) <= 0.005) && level < maxLevel);

  // 初始化液滴（+液膜）
  fraction (f, sq(1.0) - sq(x - (filmHeight + 2.0)) - sq(y));

  foreach() {
    u.x[] = -1.0*f[];
    u.y[] = 0.0;
    p[]   = 0.0;
    omega[] = 0.0;

    // 初始环境温度：无量纲 0
    T[] = 0.0;
  }
}

// ======================== 事件：更新热扩散率场 ========================
event set_thermal_diffusivity (i++) {
  foreach() {
    double ff = clamp01(f[]);
    alphaC[] = alpha_blend_nd(ff);
  }
  boundary ({alphaC});

  foreach_face()
    alphaF.x[] = (alphaC[] + alphaC[-1,0])/2.;

  boundary ((scalar *){alphaF});
}

// ======================== 事件：温度扩散 ========================
event thermal_diffusion (i++) {
  diffusion (T, dt, alphaF);
}

// ======================== 事件：更新温度变物性 ========================
event update_variable_properties (i++) {

  foreach() {
    double Tk = clamp(T_to_K(T[]), TminK, TmaxK);
    TLK[] = Tk;

    // 液体密度
    double rhoL = rhoL0_dim*(1.0 - betaL_th*(Tk - T0_refK));
    if (rhoL < 0.2*rhoL0_dim) rhoL = 0.2*rhoL0_dim;
    rhoL_loc[] = rhoL;

    // 气体密度
    double rhoG = rhoG0_dim*(T0_refK/Tk);
    rhoG_loc[] = rhoG;

    // 液体黏度
    double muL = A_muL*exp(B_muL/Tk);
    muL = clamp(muL, muMin_dim, muMax_dim);
    muL_loc[] = muL;

    // 气体黏度
    double muG = muG0_dim*pow(Tk/T0_refK, 1.5)*(T0_refK + S_suth)/(Tk + S_suth);
    muG = clamp(muG, 5e-6, 2e-4);
    muG_loc[] = muG;

    // 局部表面张力（这里用于统计/可视化；实际求解用下面的等效 sigma）
    double sigL = sigma0_dim - k_sigma*(Tk - T0_refK);
    sigma_loc[] = max(sigL, sigMin_dim);
  }
  boundary ({TLK, rhoL_loc, rhoG_loc, muL_loc, muG_loc, sigma_loc});

  // 组装面黏度给 NS 求解器
  foreach_face() {
    double ff = clamp01((f[] + f[-1,0])/2.);

    double muL_f = (muL_loc[] + muL_loc[-1,0])/2.;
    double muG_f = (muG_loc[] + muG_loc[-1,0])/2.;

    // 谐和平均更稳定
    double mu_dim = 1.0/(ff/max(muL_f,1e-30) + (1.0 - ff)/max(muG_f,1e-30));

    // 无量纲化：mu* = mu/(rhoL0 U R)
    muv.x[] = mu_to_nd(mu_dim);
  }
  boundary ((scalar *){muv});

  // 让 NS 使用变黏度场
  mu = muv;
}

// ======================== 事件：更新等效表面张力 ========================
event update_effective_sigma (i++) {
  double Vliq = 0., TliqK_int = 0.;

  foreach(reduction(+:Vliq) reduction(+:TliqK_int)) {
    Vliq     += dv()*f[];
    TliqK_int += dv()*f[]*TLK[];
  }

  double TliqK_avg = (Vliq > 1e-14) ? TliqK_int/Vliq : T0_refK;

  double sig_eff = sigma0_dim - k_sigma*(TliqK_avg - T0_refK);
  sig_eff = max(sig_eff, sigMin_dim);

  // 无量纲 sigma
  f.sigma = sigma_to_nd(sig_eff);
}

// ======================== 事件：自适应网格 ========================
event adapt (i++) {
  adapt_wavelet ((scalar *){f, u.x, u.y, T},
                 (double[]){1e-6, 4e-3, 4e-3, 2e-3},
                 maxLevel, minLevel);
}

// ======================== 输出事件 ========================
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
    if (v[j] > 1e-30)
      fprintf (fp_droplets, "%d %g %d %g %g %g\n", i, t,
               j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fp_droplets);
}

// 热统计：液滴平均温度 + 壁面热流积分代理
event thermal_stats (t += 0.01; t <= tEnd) {
  double Vliq = 0., Tint = 0.;
  foreach(reduction(+:Vliq) reduction(+:Tint)) {
    Vliq += dv()*f[];
    Tint += dv()*f[]*T[];
  }

  double Tavg_nd = (Vliq > 1e-14) ? Tint/Vliq : 0.;
  double Tavg_K  = T_ref + Tavg_nd*(T_wallK - T_ref);

  double qwall_int = 0.;
  foreach_boundary(left, reduction(+:qwall_int)) {
    double dTdx = (T[] - 1.0)/Delta;
    double ff = clamp01(f[]);
    double a_loc = alpha_blend_nd(ff);
    double qn = -a_loc*dTdx;
    qwall_int += qn*Delta;
  }

  fprintf(fp_thermal, "%d %.8g %.8g %.8g %.8g\n", i, t, Tavg_nd, Tavg_K, qwall_int);
  fflush(fp_thermal);

  fprintf(stdout, "t=%g Tavg_liq(K)=%.4f qwall_int*=%.6e sigma*=%.6e\n",
          t, Tavg_K, qwall_int, f.sigma);
  fflush(stdout);
}

event movies (t += 0.01) {
  char timestring[100];

  foreach() {
    viewingfield[] = T[];  // 绘制全域温度（空气和液体都会显示）
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
