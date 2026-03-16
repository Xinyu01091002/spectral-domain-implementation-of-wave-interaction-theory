#include "mf12_cpp/spectral.hpp"

#include <filesystem>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace mf12_cpp {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct SpectralCoefficients {
  int order = 0;
  int n_comp = 0;
  double g = 0.0;
  double h = 0.0;
  double Ux = 0.0;
  double Uy = 0.0;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> kx;
  std::vector<double> ky;
  std::vector<double> kappa;
  std::vector<double> omega1;
  std::vector<double> omega;
  std::vector<double> mu;
  std::vector<double> mu_star;
  std::vector<double> f;
  std::vector<double> g2;
  std::vector<double> f2;
  std::vector<double> gamma2v;
  std::vector<double> a2;
  std::vector<double> b2;
  std::vector<double> mu2;
  std::vector<double> kx2;
  std::vector<double> ky2;
  std::vector<double> omega2;
  std::vector<double> f_npm;
  std::vector<double> g_npm;
  std::vector<double> a_npm;
  std::vector<double> b_npm;
  std::vector<double> mu_npm;
  std::vector<double> gamma_npm;
  std::vector<double> omega_npm;
  std::vector<double> kx_npm;
  std::vector<double> ky_npm;
  std::vector<double> c;
  std::vector<double> f13;
  std::vector<double> a3;
  std::vector<double> b3;
  std::vector<double> f3;
  std::vector<double> g3;
  std::vector<double> mu3;
  std::vector<double> kx3;
  std::vector<double> ky3;
  std::vector<double> omega3v;
  std::vector<double> a_np2m;
  std::vector<double> b_np2m;
  std::vector<double> f_np2m;
  std::vector<double> g_np2m;
  std::vector<double> mu_np2m;
  std::vector<double> omega_np2m;
  std::vector<double> kx_np2m;
  std::vector<double> ky_np2m;
  std::vector<double> a_2npm;
  std::vector<double> b_2npm;
  std::vector<double> f_2npm;
  std::vector<double> g_2npm;
  std::vector<double> mu_2npm;
  std::vector<double> omega_2npm;
  std::vector<double> kx_2npm;
  std::vector<double> ky_2npm;
  std::vector<double> a_npmpp;
  std::vector<double> b_npmpp;
  std::vector<double> f_npmpp;
  std::vector<double> g_npmpp;
  std::vector<double> mu_npmpp;
  std::vector<double> omega_npmpp;
  std::vector<double> kx_npmpp;
  std::vector<double> ky_npmpp;
  std::vector<double> idx_npmpp_nm;
  std::vector<double> idx_npmpp_np;
  std::vector<double> idx_npmpp_mp;
};

struct SurfaceFields {
  Matrix eta;
  Matrix phi;
  Matrix x;
  Matrix y;
};

double coth(double x) {
  return 1.0 / std::tanh(x);
}

double gamma2(double omega1n, double knx, double kny, double kappan, double omega1m, double kmx, double kmy, double kappam, double omega_npm, double beta_npm, double g, double h) {
  const double knkm = knx * kmx + kny * kmy;
  return h / (2.0 * omega1n * omega1m * beta_npm) *
         (omega1n * omega1m * omega_npm * (omega_npm * omega_npm - omega1n * omega1m) -
          g * g * omega1n * (kappam * kappam + 2.0 * knkm) -
          g * g * omega1m * (kappan * kappan + 2.0 * knkm));
}

double lambda2(double omega1n, double knx, double kny, double kappan, double omega1m, double kmx, double kmy, double kappam, double omega_npm, double alpha_npm, double gamma_npm, double beta_npm, double g, double h) {
  const double knkm = knx * kmx + kny * kmy;
  return h / (2.0 * omega1n * omega1m * beta_npm) *
         (g * alpha_npm * (omega1n * (kappam * kappam + knkm) + omega1m * (kappan * kappan + knkm)) +
          gamma_npm * (g * g * knkm + omega1n * omega1n * omega1m * omega1m - omega1n * omega1m * omega_npm * omega_npm));
}

double theta_a(double an, double bn, double am, double bm, double ap, double bp, double h) {
  return (an * am * ap - bn * bm * ap - bn * am * bp - an * bm * bp) / (h * h);
}

double theta_b(double an, double bn, double am, double bm, double ap, double bp, double h) {
  return (bn * am * ap + an * bm * ap + an * am * bp - bn * bm * bp) / (h * h);
}

double upsilon_nm(double omega1n, double knx, double kny, double kappan, double omega1m, double kmx, double kmy, double kappam, double fnpm, double fnmm, double gnpm, double gnmm, double kappanpm, double kappanmm, double g, double h) {
  const double knkm = knx * kmx + kny * kmy;
  return g / (4.0 * omega1n * omega1m * std::cosh(h * kappan)) * (omega1m * (kappan * kappan - kappam * kappam) - omega1n * knkm) +
         (gnpm + gnmm) / (4.0 * h * omega1n * omega1n * omega1m * std::cosh(h * kappan)) * (g * g * knkm + omega1m * omega1m * omega1m * omega1n) -
         1.0 / (4.0 * h * std::cosh(h * kappan)) * (fnpm * kappanpm * std::sinh(h * kappanpm) + fnmm * kappanmm * std::sinh(h * kappanmm)) +
         g * fnpm * std::cosh(h * kappanpm) / (4.0 * h * omega1n * omega1n * omega1m * std::cosh(h * kappan)) *
             ((omega1n + omega1m) * (knkm + kappam * kappam) - omega1m * kappanpm * kappanpm) +
         g * fnmm * std::cosh(h * kappanmm) / (4.0 * h * omega1n * omega1n * omega1m * std::cosh(h * kappan)) *
             ((omega1n - omega1m) * (knkm - kappam * kappam) - omega1m * kappanmm * kappanmm);
}

double xi_nm(double omega1n, double kappan, double omega1m, double gnpm, double fnpm, double gamma_npm, double gnmm, double fnmm, double gamma_nmm, double h, double g) {
  return 1.0 / (2.0 * h) * (omega1m * (gnpm - gnmm) + fnpm * gamma_npm + fnmm * gamma_nmm - g * h * kappan * kappan / (2.0 * omega1n));
}

double omega_nm_fn(double omega1n, double knx, double kny, double omega1m, double kmx, double kmy, double kappam, double fnpm, double fnmm, double gnpm, double gnmm, double kappanpm, double kappanmm, double g, double h) {
  const double knkm = knx * kmx + kny * kmy;
  return 1.0 / (kappam * kappam) * (((2.0 * omega1m * omega1m + omega1n * omega1n) / (4.0 * omega1n * omega1m)) * knkm + 0.25 * kappam * kappam) +
         (gnpm + gnmm) / (kappam * kappam) * (g * knkm / (4.0 * h * omega1n * omega1m) - omega1m * omega1m / (4.0 * g * h)) +
         omega1n / (4.0 * g * h * kappam * kappam) * (fnpm * kappanpm * std::sinh(h * kappanpm) + fnmm * kappanmm * std::sinh(h * kappanmm)) -
         fnpm * std::cosh(h * kappanpm) / (4.0 * h * omega1n * omega1m * kappam * kappam) * ((omega1n - omega1m) * (kappam * kappam + knkm) + omega1m * kappanpm * kappanpm) +
         fnmm * std::cosh(h * kappanmm) / (4.0 * h * omega1n * omega1m * kappam * kappam) * ((omega1n + omega1m) * (kappam * kappam - knkm) - omega1m * kappanmm * kappanmm);
}

double lambda3(
    double omega1n, double knx, double kny, double kappan,
    double omega1m, double kmx, double kmy, double kappam,
    double omega1p, double kpx, double kpy, double kappap,
    double kappanpm, double gammanpm, double gnpm, double fnpm,
    double kappanpp, double gammanpp, double gnpp, double fnpp,
    double kappampp, double gammampp, double gmpp, double fmpp,
    double omega_npmpp, double alpha_npmpp, double gamma_npmpp, double beta_npmpp,
    double g, double h) {
  const double knkm = knx * kmx + kny * kmy;
  const double knkp = knx * kpx + kny * kpy;
  const double kmkp = kmx * kpx + kmy * kpy;
  return (
      h * h / (4.0 * beta_npmpp) *
          (alpha_npmpp *
               (omega1n * (knkm + knkp + kappan * kappan) +
                omega1m * (knkm + kmkp + kappam * kappam) +
                omega1p * (knkp + kmkp + kappap * kappap)) +
           gamma_npmpp *
               (g / omega1n * (omega1m * knkm + omega1p * knkp - omega_npmpp * kappan * kappan) +
                g / omega1m * (omega1n * knkm + omega1p * kmkp - omega_npmpp * kappam * kappam) +
                g / omega1p * (omega1n * knkp + omega1m * kmkp - omega_npmpp * kappap * kappap))) -
      h * fnpm / (2.0 * beta_npmpp) *
          (alpha_npmpp * std::cosh(h * kappanpm) * (knkp + kmkp + kappanpm * kappanpm) +
           gamma_npmpp * (g / omega1p * (knkp + kmkp) * std::cosh(h * kappanpm) - gammanpm * omega_npmpp)) -
      h * fnpp / (2.0 * beta_npmpp) *
          (alpha_npmpp * std::cosh(h * kappanpp) * (knkm + kmkp + kappanpp * kappanpp) +
           gamma_npmpp * (g / omega1m * (knkm + kmkp) * std::cosh(h * kappanpp) - gammanpp * omega_npmpp)) -
      h * fmpp / (2.0 * beta_npmpp) *
          (alpha_npmpp * std::cosh(h * kappampp) * (knkm + knkp + kappampp * kappampp) +
           gamma_npmpp * (g / omega1n * (knkm + knkp) * std::cosh(h * kappampp) - gammampp * omega_npmpp)) +
      h * gnpm / (2.0 * beta_npmpp) * (alpha_npmpp * g / omega1p * (knkp + kmkp + kappap * kappap) - gamma_npmpp * omega1p * omega1p) +
      h * gnpp / (2.0 * beta_npmpp) * (alpha_npmpp * g / omega1m * (knkm + kmkp + kappam * kappam) - gamma_npmpp * omega1m * omega1m) +
      h * gmpp / (2.0 * beta_npmpp) * (alpha_npmpp * g / omega1n * (knkm + knkp + kappan * kappan) - gamma_npmpp * omega1n * omega1n));
}

double gamma3(
    double omega1n, double knx, double kny, double kappan,
    double omega1m, double kmx, double kmy, double kappam,
    double omega1p, double kpx, double kpy, double kappap,
    double kappanpm, double gammanpm, double gnpm, double fnpm,
    double kappanpp, double gammanpp, double gnpp, double fnpp,
    double kappampp, double gammampp, double gmpp, double fmpp,
    double omega_npmpp, double beta_npmpp,
    double g, double h) {
  const double knkm = knx * kmx + kny * kmy;
  const double knkp = knx * kpx + kny * kpy;
  const double kmkp = kmx * kpx + kmy * kpy;
  return (
      -g * h * h / (4.0 * beta_npmpp) *
          (omega1n * (knkm + knkp + kappan * kappan) +
           omega1m * (knkm + kmkp + kappam * kappam) +
           omega1p * (knkp + kmkp + kappap * kappap) +
           omega_npmpp / omega1n * (omega1m * knkm + omega1p * knkp - omega_npmpp * kappan * kappan) +
           omega_npmpp / omega1m * (omega1n * knkm + omega1p * kmkp - omega_npmpp * kappam * kappam) +
           omega_npmpp / omega1p * (omega1n * knkp + omega1m * kmkp - omega_npmpp * kappap * kappap)) +
      h * fnpm / (2.0 * beta_npmpp) *
          (g * std::cosh(h * kappanpm) * ((knkp + kmkp + kappanpm * kappanpm) + omega_npmpp / omega1p * (knkp + kmkp)) -
           gammanpm * omega_npmpp * omega_npmpp) +
      h * fnpp / (2.0 * beta_npmpp) *
          (g * std::cosh(h * kappanpp) * ((knkm + kmkp + kappanpp * kappanpp) + omega_npmpp / omega1m * (knkm + kmkp)) -
           gammanpp * omega_npmpp * omega_npmpp) +
      h * fmpp / (2.0 * beta_npmpp) *
          (g * std::cosh(h * kappampp) * ((knkm + knkp + kappampp * kappampp) + omega_npmpp / omega1n * (knkm + knkp)) -
           gammampp * omega_npmpp * omega_npmpp) +
      h * gnpm / (2.0 * beta_npmpp) * (omega1p * omega1p * omega_npmpp - g * g / omega1p * (knkp + kmkp + kappap * kappap)) +
      h * gnpp / (2.0 * beta_npmpp) * (omega1m * omega1m * omega_npmpp - g * g / omega1m * (knkm + kmkp + kappam * kappam)) +
      h * gmpp / (2.0 * beta_npmpp) * (omega1n * omega1n * omega_npmpp - g * g / omega1n * (knkm + knkp + kappan * kappan)));
}

double pi_fn(
    double omega1n, double kappan,
    double omega1m, double kappam,
    double omega1p, double kappap,
    double gamma_npm, double gnpm, double fnpm,
    double gamma_npp, double gnpp, double fnpp,
    double gamma_mpp, double gmpp, double fmpp,
    double fnpmpp, double kappa_npmpp,
    double g, double h) {
  return fnpmpp * std::cosh(h * kappa_npmpp) -
         g * h * h / 4.0 * (kappan * kappan / omega1n + kappam * kappam / omega1m + kappap * kappap / omega1p) -
         h / 2.0 * (omega1n * gmpp + omega1m * gnpp + omega1p * gnpm) +
         h / 2.0 * (fnpm * gamma_npm + fnpp * gamma_npp + fmpp * gamma_mpp);
}

void add_third_order_terms(SpectralCoefficients& coeffs, const CaseInputs& inp) {
  const int n_comp = coeffs.n_comp;
  coeffs.c.resize(n_comp);
  coeffs.f13.resize(n_comp);
  coeffs.a3.resize(n_comp);
  coeffs.b3.resize(n_comp);
  coeffs.f3.resize(n_comp);
  coeffs.g3.resize(n_comp);
  coeffs.mu3.resize(n_comp);
  coeffs.kx3.resize(n_comp);
  coeffs.ky3.resize(n_comp);
  coeffs.omega3v.resize(n_comp);

  std::vector<double> gamma2v_local(n_comp);
  for (int i = 0; i < n_comp; ++i) {
    coeffs.c[i] = std::hypot(coeffs.a[i], coeffs.b[i]);
    gamma2v_local[i] = 2.0 * coeffs.kappa[i] * std::sinh(inp.h * 2.0 * coeffs.kappa[i]);
    const double hk = inp.h * coeffs.kappa[i];
    const double upsilon = coeffs.omega1[i] * coeffs.kappa[i] * (-13.0 + 24.0 * std::cosh(2.0 * hk) + std::cosh(4.0 * hk)) / (64.0 * std::pow(std::sinh(hk), 5.0));
    const double xi = (coeffs.omega1[i] * coeffs.g2[i] + coeffs.f2[i] * 2.0 * coeffs.kappa[i] * std::sinh(inp.h * 2.0 * coeffs.kappa[i]) - inp.g * inp.h * coeffs.kappa[i] * coeffs.kappa[i] / (2.0 * coeffs.omega1[i])) / (4.0 * inp.h);
    coeffs.f13[i] = coeffs.c[i] * coeffs.c[i] * upsilon;
    coeffs.mu_star[i] = coeffs.c[i] * coeffs.c[i] * xi;
    const double omega_cap = (8.0 + std::cosh(4.0 * hk)) / (16.0 * std::pow(std::sinh(hk), 4.0));
    coeffs.omega3v[i] = coeffs.c[i] * coeffs.c[i] * coeffs.kappa[i] * coeffs.kappa[i] * omega_cap;
  }

  auto eval_pair = [&](int n, int m, double pm, double& omega_out, double& kappa_out, double& gamma_out, double& f_out, double& g_out) {
    omega_out = coeffs.omega1[n] + pm * coeffs.omega1[m];
    const double kx_out = coeffs.kx[n] + pm * coeffs.kx[m];
    const double ky_out = coeffs.ky[n] + pm * coeffs.ky[m];
    kappa_out = std::hypot(kx_out, ky_out);
    gamma_out = kappa_out * std::sinh(inp.h * kappa_out);
    const double beta_out = omega_out * omega_out * std::cosh(inp.h * kappa_out) - inp.g * kappa_out * std::sinh(inp.h * kappa_out);
    f_out = gamma2(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], pm * coeffs.omega1[m], pm * coeffs.kx[m], pm * coeffs.ky[m], coeffs.kappa[m], omega_out, beta_out, inp.g, inp.h);
    g_out = lambda2(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], pm * coeffs.omega1[m], pm * coeffs.kx[m], pm * coeffs.ky[m], coeffs.kappa[m], omega_out, omega_out * std::cosh(inp.h * kappa_out), gamma_out, beta_out, inp.g, inp.h);
  };

  for (int n = 0; n < n_comp; ++n) {
    for (int m = 0; m < n_comp; ++m) {
      if (m == n) {
        continue;
      }
      double omega_np = 0.0, kappa_np = 0.0, gamma_np = 0.0, f_np = 0.0, g_np = 0.0;
      double omega_nm = 0.0, kappa_nm = 0.0, gamma_nm = 0.0, f_nm = 0.0, g_nm = 0.0;
      eval_pair(n, m, 1.0, omega_np, kappa_np, gamma_np, f_np, g_np);
      eval_pair(n, m, -1.0, omega_nm, kappa_nm, gamma_nm, f_nm, g_nm);
      coeffs.f13[n] += coeffs.c[m] * coeffs.c[m] * upsilon_nm(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], f_np, f_nm, g_np, g_nm, kappa_np, kappa_nm, inp.g, inp.h);
      coeffs.mu_star[n] += coeffs.c[m] * coeffs.c[m] * xi_nm(coeffs.omega1[n], coeffs.kappa[n], coeffs.omega1[m], g_np, f_np, gamma_np, g_nm, f_nm, gamma_nm, inp.h, inp.g);
      coeffs.omega3v[n] += coeffs.c[m] * coeffs.c[m] * coeffs.kappa[m] * coeffs.kappa[m] * omega_nm_fn(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], f_np, f_nm, g_np, g_nm, kappa_np, kappa_nm, inp.g, inp.h);
    }
  }

  for (int n = 0; n < n_comp; ++n) {
    coeffs.omega[n] += coeffs.omega3v[n] * coeffs.omega1[n];
    coeffs.mu_star[n] += coeffs.f13[n] * std::cosh(inp.h * coeffs.kappa[n]);
    coeffs.a3[n] = 0.5 * theta_a(coeffs.a[n], coeffs.b[n], coeffs.a[n], coeffs.b[n], coeffs.a[n], coeffs.b[n], inp.h);
    coeffs.b3[n] = 0.5 * theta_b(coeffs.a[n], coeffs.b[n], coeffs.a[n], coeffs.b[n], coeffs.a[n], coeffs.b[n], inp.h);
    coeffs.f3[n] = (inp.h * inp.h * coeffs.kappa[n] * coeffs.omega1[n] / (32.0 * std::pow(std::sinh(inp.h * coeffs.kappa[n]), 7.0))) * (-11.0 + 2.0 * std::cosh(2.0 * inp.h * coeffs.kappa[n]));
    coeffs.g3[n] = (3.0 * inp.h * inp.h * coeffs.kappa[n] * coeffs.kappa[n] / (128.0 * std::pow(std::sinh(inp.h * coeffs.kappa[n]), 6.0))) *
                   (14.0 + 15.0 * std::cosh(2.0 * inp.h * coeffs.kappa[n]) + 6.0 * std::cosh(4.0 * inp.h * coeffs.kappa[n]) + std::cosh(6.0 * inp.h * coeffs.kappa[n]));
    coeffs.mu3[n] = coeffs.f3[n] * std::cosh(inp.h * 3.0 * coeffs.kappa[n]) - inp.g * inp.h * inp.h * coeffs.kappa[n] * coeffs.kappa[n] / (4.0 * coeffs.omega1[n]) + 0.5 * inp.h * (coeffs.f2[n] * gamma2v_local[n] - coeffs.omega1[n] * coeffs.g2[n]);
    coeffs.kx3[n] = 3.0 * coeffs.kx[n];
    coeffs.ky3[n] = 3.0 * coeffs.ky[n];
    coeffs.omega3v[n] = 3.0 * coeffs.omega[n];
  }

  std::vector<int> m_nm_row_odd(n_comp * n_comp, 0);
  std::vector<int> m_nm_col_odd(n_comp * n_comp, 0);
  int pair_count = 0;
  for (int n = 0; n < n_comp; ++n) {
    for (int m = n + 1; m < n_comp; ++m) {
      ++pair_count;
      m_nm_row_odd[n * n_comp + m] = 2 * pair_count - 1;
    }
  }
  pair_count = 0;
  for (int col = 1; col < n_comp; ++col) {
    for (int row = 0; row < col; ++row) {
      ++pair_count;
      m_nm_col_odd[row * n_comp + col] = 2 * pair_count - 1;
    }
  }

  const int num_pairs = n_comp * (n_comp - 1) / 2;
  coeffs.a_np2m.assign(num_pairs, 0.0);
  coeffs.b_np2m.assign(num_pairs, 0.0);
  coeffs.f_np2m.assign(num_pairs, 0.0);
  coeffs.g_np2m.assign(num_pairs, 0.0);
  coeffs.mu_np2m.assign(num_pairs, 0.0);
  coeffs.omega_np2m.assign(num_pairs, 0.0);
  coeffs.kx_np2m.assign(num_pairs, 0.0);
  coeffs.ky_np2m.assign(num_pairs, 0.0);
  coeffs.a_2npm.assign(num_pairs, 0.0);
  coeffs.b_2npm.assign(num_pairs, 0.0);
  coeffs.f_2npm.assign(num_pairs, 0.0);
  coeffs.g_2npm.assign(num_pairs, 0.0);
  coeffs.mu_2npm.assign(num_pairs, 0.0);
  coeffs.omega_2npm.assign(num_pairs, 0.0);
  coeffs.kx_2npm.assign(num_pairs, 0.0);
  coeffs.ky_2npm.assign(num_pairs, 0.0);

  int pair_idx = 0;
  for (int n = 0; n < n_comp; ++n) {
    for (int m = n + 1; m < n_comp; ++m) {
      const int idx_sum_nm = m_nm_row_odd[n * n_comp + m] - 1;
      coeffs.omega_np2m[pair_idx] = coeffs.omega1[n] + 2.0 * coeffs.omega1[m];
      coeffs.kx_np2m[pair_idx] = coeffs.kx[n] + 2.0 * coeffs.kx[m];
      coeffs.ky_np2m[pair_idx] = coeffs.ky[n] + 2.0 * coeffs.ky[m];
      const double kappa_np2m = std::hypot(coeffs.kx_np2m[pair_idx], coeffs.ky_np2m[pair_idx]);
      const double alpha_np2m = coeffs.omega_np2m[pair_idx] * std::cosh(inp.h * kappa_np2m);
      const double gamma_np2m = kappa_np2m * std::sinh(inp.h * kappa_np2m);
      const double beta_np2m = coeffs.omega_np2m[pair_idx] * coeffs.omega_np2m[pair_idx] * std::cosh(inp.h * kappa_np2m) - inp.g * kappa_np2m * std::sinh(inp.h * kappa_np2m);
      coeffs.a_np2m[pair_idx] = 0.5 * theta_a(coeffs.a[n], coeffs.b[n], coeffs.a[m], coeffs.b[m], coeffs.a[m], coeffs.b[m], inp.h);
      coeffs.b_np2m[pair_idx] = 0.5 * theta_b(coeffs.a[n], coeffs.b[n], coeffs.a[m], coeffs.b[m], coeffs.a[m], coeffs.b[m], inp.h);
      coeffs.g_np2m[pair_idx] = lambda3(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], 2.0 * coeffs.kappa[m], gamma2v_local[m], coeffs.g2[m], coeffs.f2[m], coeffs.omega_np2m[pair_idx], alpha_np2m, gamma_np2m, beta_np2m, inp.g, inp.h);
      coeffs.f_np2m[pair_idx] = gamma3(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], 2.0 * coeffs.kappa[m], gamma2v_local[m], coeffs.g2[m], coeffs.f2[m], coeffs.omega_np2m[pair_idx], beta_np2m, inp.g, inp.h);
      coeffs.mu_np2m[pair_idx] = pi_fn(coeffs.omega1[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kappa[m], coeffs.omega1[m], coeffs.kappa[m], coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], gamma2v_local[m], coeffs.g2[m], coeffs.f2[m], coeffs.f_np2m[pair_idx], kappa_np2m, inp.g, inp.h);

      coeffs.omega_2npm[pair_idx] = 2.0 * coeffs.omega1[n] + coeffs.omega1[m];
      coeffs.kx_2npm[pair_idx] = 2.0 * coeffs.kx[n] + coeffs.kx[m];
      coeffs.ky_2npm[pair_idx] = 2.0 * coeffs.ky[n] + coeffs.ky[m];
      const double kappa_2npm = std::hypot(coeffs.kx_2npm[pair_idx], coeffs.ky_2npm[pair_idx]);
      const double alpha_2npm = coeffs.omega_2npm[pair_idx] * std::cosh(inp.h * kappa_2npm);
      const double gamma_2npm = kappa_2npm * std::sinh(inp.h * kappa_2npm);
      const double beta_2npm = coeffs.omega_2npm[pair_idx] * coeffs.omega_2npm[pair_idx] * std::cosh(inp.h * kappa_2npm) - inp.g * kappa_2npm * std::sinh(inp.h * kappa_2npm);
      coeffs.a_2npm[pair_idx] = 0.5 * theta_a(coeffs.a[n], coeffs.b[n], coeffs.a[n], coeffs.b[n], coeffs.a[m], coeffs.b[m], inp.h);
      coeffs.b_2npm[pair_idx] = 0.5 * theta_b(coeffs.a[n], coeffs.b[n], coeffs.a[n], coeffs.b[n], coeffs.a[m], coeffs.b[m], inp.h);
      coeffs.g_2npm[pair_idx] = lambda3(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], 2.0 * coeffs.kappa[n], gamma2v_local[n], coeffs.g2[n], coeffs.f2[n], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], coeffs.omega_2npm[pair_idx], alpha_2npm, gamma_2npm, beta_2npm, inp.g, inp.h);
      coeffs.f_2npm[pair_idx] = gamma3(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kx[m], coeffs.ky[m], coeffs.kappa[m], 2.0 * coeffs.kappa[n], gamma2v_local[n], coeffs.g2[n], coeffs.f2[n], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], std::hypot(coeffs.kx[n] + coeffs.kx[m], coeffs.ky[n] + coeffs.ky[m]), coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], coeffs.omega_2npm[pair_idx], beta_2npm, inp.g, inp.h);
      coeffs.mu_2npm[pair_idx] = pi_fn(coeffs.omega1[n], coeffs.kappa[n], coeffs.omega1[n], coeffs.kappa[n], coeffs.omega1[m], coeffs.kappa[m], gamma2v_local[n], coeffs.g2[n], coeffs.f2[n], coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], coeffs.gamma_npm[idx_sum_nm], coeffs.g_npm[idx_sum_nm], coeffs.f_npm[idx_sum_nm], coeffs.f_2npm[pair_idx], kappa_2npm, inp.g, inp.h);
      ++pair_idx;
    }
  }

  const int num_triplets = n_comp * (n_comp - 1) * (n_comp - 2) / 6;
  coeffs.a_npmpp.assign(num_triplets, 0.0);
  coeffs.b_npmpp.assign(num_triplets, 0.0);
  coeffs.f_npmpp.assign(num_triplets, 0.0);
  coeffs.g_npmpp.assign(num_triplets, 0.0);
  coeffs.mu_npmpp.assign(num_triplets, 0.0);
  coeffs.omega_npmpp.assign(num_triplets, 0.0);
  coeffs.kx_npmpp.assign(num_triplets, 0.0);
  coeffs.ky_npmpp.assign(num_triplets, 0.0);
  coeffs.idx_npmpp_nm.assign(num_triplets, 0.0);
  coeffs.idx_npmpp_np.assign(num_triplets, 0.0);
  coeffs.idx_npmpp_mp.assign(num_triplets, 0.0);

  int c3 = 0;
  auto pair_index_col = [](int row, int col) {
    return 2 * ((col * (col - 1)) / 2 + row);
  };
  for (int n = 0; n < n_comp; ++n) {
    for (int m = n + 1; m < n_comp; ++m) {
      const int idx_sum_nm = m_nm_row_odd[n * n_comp + m] - 1;
      for (int p = m + 1; p < n_comp; ++p) {
        const int idx_sum_np = pair_index_col(n, p);
        const int idx_sum_mp = pair_index_col(m, p);
        coeffs.idx_npmpp_nm[c3] = static_cast<double>(idx_sum_nm);
        coeffs.idx_npmpp_np[c3] = static_cast<double>(idx_sum_np);
        coeffs.idx_npmpp_mp[c3] = static_cast<double>(idx_sum_mp);
        coeffs.omega_npmpp[c3] = coeffs.omega1[n] + coeffs.omega1[m] + coeffs.omega1[p];
        coeffs.kx_npmpp[c3] = coeffs.kx[n] + coeffs.kx[m] + coeffs.kx[p];
        coeffs.ky_npmpp[c3] = coeffs.ky[n] + coeffs.ky[m] + coeffs.ky[p];
        const double kappa_npmpp = std::hypot(coeffs.kx_npmpp[c3], coeffs.ky_npmpp[c3]);
        const double alpha_npmpp = coeffs.omega_npmpp[c3] * std::cosh(inp.h * kappa_npmpp);
        const double gamma_npmpp = kappa_npmpp * std::sinh(inp.h * kappa_npmpp);
        const double beta_npmpp = coeffs.omega_npmpp[c3] * coeffs.omega_npmpp[c3] * std::cosh(inp.h * kappa_npmpp) - inp.g * kappa_npmpp * std::sinh(inp.h * kappa_npmpp);
        const double kappanpm = std::hypot(coeffs.kx_npm[idx_sum_nm], coeffs.ky_npm[idx_sum_nm]);
        const double kappanpp = std::hypot(coeffs.kx_npm[idx_sum_np], coeffs.ky_npm[idx_sum_np]);
        const double kappampp = std::hypot(coeffs.kx_npm[idx_sum_mp], coeffs.ky_npm[idx_sum_mp]);
        const double gammanpm = coeffs.gamma_npm[idx_sum_nm];
        const double gammanpp = coeffs.gamma_npm[idx_sum_np];
        const double gammampp = coeffs.gamma_npm[idx_sum_mp];
        const double gnpm = coeffs.g_npm[idx_sum_nm];
        const double gnpp = coeffs.g_npm[idx_sum_np];
        const double gmpp = coeffs.g_npm[idx_sum_mp];
        const double fnpm = coeffs.f_npm[idx_sum_nm];
        const double fnpp = coeffs.f_npm[idx_sum_np];
        const double fmpp = coeffs.f_npm[idx_sum_mp];
        const double knkm = coeffs.kx[n] * coeffs.kx[m] + coeffs.ky[n] * coeffs.ky[m];
        const double knkp = coeffs.kx[n] * coeffs.kx[p] + coeffs.ky[n] * coeffs.ky[p];
        const double kmkp = coeffs.kx[m] * coeffs.kx[p] + coeffs.ky[m] * coeffs.ky[p];
        coeffs.a_npmpp[c3] = 0.5 * theta_a(coeffs.a[n], coeffs.b[n], coeffs.a[m], coeffs.b[m], coeffs.a[p], coeffs.b[p], inp.h);
        coeffs.b_npmpp[c3] = 0.5 * theta_b(coeffs.a[n], coeffs.b[n], coeffs.a[m], coeffs.b[m], coeffs.a[p], coeffs.b[p], inp.h);
        coeffs.g_npmpp[c3] =
            inp.h * inp.h / (4.0 * beta_npmpp) *
                (alpha_npmpp *
                     (coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]) +
                      coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]) +
                      coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p])) +
                 gamma_npmpp *
                     (inp.g / coeffs.omega1[n] * (coeffs.omega1[m] * knkm + coeffs.omega1[p] * knkp - coeffs.omega_npmpp[c3] * coeffs.kappa[n] * coeffs.kappa[n]) +
                      inp.g / coeffs.omega1[m] * (coeffs.omega1[n] * knkm + coeffs.omega1[p] * kmkp - coeffs.omega_npmpp[c3] * coeffs.kappa[m] * coeffs.kappa[m]) +
                      inp.g / coeffs.omega1[p] * (coeffs.omega1[n] * knkp + coeffs.omega1[m] * kmkp - coeffs.omega_npmpp[c3] * coeffs.kappa[p] * coeffs.kappa[p]))) -
            inp.h * fnpm / (2.0 * beta_npmpp) *
                (alpha_npmpp * std::cosh(inp.h * kappanpm) * (knkp + kmkp + kappanpm * kappanpm) +
                 gamma_npmpp * (inp.g / coeffs.omega1[p] * (knkp + kmkp) * std::cosh(inp.h * kappanpm) - gammanpm * coeffs.omega_npmpp[c3])) -
            inp.h * fnpp / (2.0 * beta_npmpp) *
                (alpha_npmpp * std::cosh(inp.h * kappanpp) * (knkm + kmkp + kappanpp * kappanpp) +
                 gamma_npmpp * (inp.g / coeffs.omega1[m] * (knkm + kmkp) * std::cosh(inp.h * kappanpp) - gammanpp * coeffs.omega_npmpp[c3])) -
            inp.h * fmpp / (2.0 * beta_npmpp) *
                (alpha_npmpp * std::cosh(inp.h * kappampp) * (knkm + knkp + kappampp * kappampp) +
                 gamma_npmpp * (inp.g / coeffs.omega1[n] * (knkm + knkp) * std::cosh(inp.h * kappampp) - gammampp * coeffs.omega_npmpp[c3])) +
            inp.h * gnpm / (2.0 * beta_npmpp) * (alpha_npmpp * inp.g / coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p]) - gamma_npmpp * coeffs.omega1[p] * coeffs.omega1[p]) +
            inp.h * gnpp / (2.0 * beta_npmpp) * (alpha_npmpp * inp.g / coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]) - gamma_npmpp * coeffs.omega1[m] * coeffs.omega1[m]) +
            inp.h * gmpp / (2.0 * beta_npmpp) * (alpha_npmpp * inp.g / coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]) - gamma_npmpp * coeffs.omega1[n] * coeffs.omega1[n]);
        coeffs.f_npmpp[c3] =
            -inp.g * inp.h * inp.h / (4.0 * beta_npmpp) *
                (coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]) +
                 coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]) +
                 coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p]) +
                 coeffs.omega_npmpp[c3] / coeffs.omega1[n] * (coeffs.omega1[m] * knkm + coeffs.omega1[p] * knkp - coeffs.omega_npmpp[c3] * coeffs.kappa[n] * coeffs.kappa[n]) +
                 coeffs.omega_npmpp[c3] / coeffs.omega1[m] * (coeffs.omega1[n] * knkm + coeffs.omega1[p] * kmkp - coeffs.omega_npmpp[c3] * coeffs.kappa[m] * coeffs.kappa[m]) +
                 coeffs.omega_npmpp[c3] / coeffs.omega1[p] * (coeffs.omega1[n] * knkp + coeffs.omega1[m] * kmkp - coeffs.omega_npmpp[c3] * coeffs.kappa[p] * coeffs.kappa[p])) +
            inp.h * fnpm / (2.0 * beta_npmpp) *
                (inp.g * std::cosh(inp.h * kappanpm) * ((knkp + kmkp + kappanpm * kappanpm) + coeffs.omega_npmpp[c3] / coeffs.omega1[p] * (knkp + kmkp)) -
                 gammanpm * coeffs.omega_npmpp[c3] * coeffs.omega_npmpp[c3]) +
            inp.h * fnpp / (2.0 * beta_npmpp) *
                (inp.g * std::cosh(inp.h * kappanpp) * ((knkm + kmkp + kappanpp * kappanpp) + coeffs.omega_npmpp[c3] / coeffs.omega1[m] * (knkm + kmkp)) -
                 gammanpp * coeffs.omega_npmpp[c3] * coeffs.omega_npmpp[c3]) +
            inp.h * fmpp / (2.0 * beta_npmpp) *
                (inp.g * std::cosh(inp.h * kappampp) * ((knkm + knkp + kappampp * kappampp) + coeffs.omega_npmpp[c3] / coeffs.omega1[n] * (knkm + knkp)) -
                 gammampp * coeffs.omega_npmpp[c3] * coeffs.omega_npmpp[c3]) +
            inp.h * gnpm / (2.0 * beta_npmpp) * (coeffs.omega1[p] * coeffs.omega1[p] * coeffs.omega_npmpp[c3] - inp.g * inp.g / coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p])) +
            inp.h * gnpp / (2.0 * beta_npmpp) * (coeffs.omega1[m] * coeffs.omega1[m] * coeffs.omega_npmpp[c3] - inp.g * inp.g / coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m])) +
            inp.h * gmpp / (2.0 * beta_npmpp) * (coeffs.omega1[n] * coeffs.omega1[n] * coeffs.omega_npmpp[c3] - inp.g * inp.g / coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]));
        coeffs.mu_npmpp[c3] =
            coeffs.f_npmpp[c3] * std::cosh(inp.h * kappa_npmpp) -
            inp.g * inp.h * inp.h / 4.0 * (coeffs.kappa[n] * coeffs.kappa[n] / coeffs.omega1[n] + coeffs.kappa[m] * coeffs.kappa[m] / coeffs.omega1[m] + coeffs.kappa[p] * coeffs.kappa[p] / coeffs.omega1[p]) -
            inp.h / 2.0 * (coeffs.omega1[n] * gmpp + coeffs.omega1[m] * gnpp + coeffs.omega1[p] * gnpm) +
            inp.h / 2.0 * (fnpm * gammanpm + fnpp * gammanpp + fmpp * gammampp);
        ++c3;
      }
    }
  }

  for (int n = 0; n < n_comp; ++n) {
    coeffs.omega2[n] = 2.0 * coeffs.omega[n];
  }

  pair_idx = 0;
  for (int n = 0; n < n_comp; ++n) {
    for (int m = n + 1; m < n_comp; ++m) {
      const int idx_sum = 2 * pair_idx;
      const int idx_diff = idx_sum + 1;
      coeffs.omega_npm[idx_sum] = coeffs.omega[n] + coeffs.omega[m];
      coeffs.omega_npm[idx_diff] = coeffs.omega[n] - coeffs.omega[m];
      coeffs.omega_np2m[pair_idx] = coeffs.omega[n] + 2.0 * coeffs.omega[m];
      coeffs.omega_2npm[pair_idx] = 2.0 * coeffs.omega[n] + coeffs.omega[m];
      ++pair_idx;
    }
  }
  c3 = 0;
  for (int n = 0; n < n_comp; ++n) {
    for (int m = n + 1; m < n_comp; ++m) {
      for (int p = m + 1; p < n_comp; ++p) {
        coeffs.omega_npmpp[c3] = coeffs.omega[n] + coeffs.omega[m] + coeffs.omega[p];
        ++c3;
      }
    }
  }
}

SpectralCoefficients compute_coefficients(const LoadedCase& loaded) {
  const auto& inp = loaded.manifest.inputs;
  SpectralCoefficients coeffs;
  coeffs.order = inp.order;
  coeffs.g = inp.g;
  coeffs.h = inp.h;
  coeffs.Ux = inp.Ux;
  coeffs.Uy = inp.Uy;

  auto copy_vector = [](const Matrix& mat) {
    return mat.values;
  };

  coeffs.a = copy_vector(loaded.arrays.at("a"));
  coeffs.b = copy_vector(loaded.arrays.at("b"));
  coeffs.kx = copy_vector(loaded.arrays.at("kx"));
  coeffs.ky = copy_vector(loaded.arrays.at("ky"));
  coeffs.n_comp = static_cast<int>(coeffs.a.size());

  coeffs.kappa.resize(coeffs.n_comp);
  coeffs.omega1.resize(coeffs.n_comp);
  coeffs.omega.resize(coeffs.n_comp);
  coeffs.mu.resize(coeffs.n_comp);
  coeffs.mu_star.assign(coeffs.n_comp, 0.0);
  coeffs.f.resize(coeffs.n_comp);
  coeffs.kx2.resize(coeffs.n_comp);
  coeffs.ky2.resize(coeffs.n_comp);
  coeffs.omega2.resize(coeffs.n_comp);

  for (int i = 0; i < coeffs.n_comp; ++i) {
    coeffs.kappa[i] = std::hypot(coeffs.kx[i], coeffs.ky[i]);
    coeffs.omega1[i] = std::sqrt(inp.g * coeffs.kappa[i] * std::tanh(inp.h * coeffs.kappa[i]));
    coeffs.omega[i] = coeffs.kx[i] * inp.Ux + coeffs.ky[i] * inp.Uy + coeffs.omega1[i];
    coeffs.f[i] = -coeffs.omega1[i] / (coeffs.kappa[i] * std::sinh(inp.h * coeffs.kappa[i]));
    coeffs.mu[i] = coeffs.f[i] * std::cosh(inp.h * coeffs.kappa[i]);
    coeffs.kx2[i] = 2.0 * coeffs.kx[i];
    coeffs.ky2[i] = 2.0 * coeffs.ky[i];
    coeffs.omega2[i] = 2.0 * coeffs.omega[i];
  }

  if (inp.order >= 2) {
    coeffs.g2.resize(coeffs.n_comp);
    coeffs.f2.resize(coeffs.n_comp);
    coeffs.gamma2v.resize(coeffs.n_comp);
    coeffs.a2.resize(coeffs.n_comp);
    coeffs.b2.resize(coeffs.n_comp);
    coeffs.mu2.resize(coeffs.n_comp);

    for (int i = 0; i < coeffs.n_comp; ++i) {
      const double hk = inp.h * coeffs.kappa[i];
      coeffs.g2[i] = 0.5 * inp.h * coeffs.kappa[i] * (2.0 + std::cosh(2.0 * hk)) * coth(hk) / std::pow(std::sinh(hk), 2.0);
      coeffs.f2[i] = -0.75 * inp.h * coeffs.omega1[i] / std::pow(std::sinh(hk), 4.0);
      coeffs.gamma2v[i] = 2.0 * coeffs.kappa[i] * std::sinh(inp.h * 2.0 * coeffs.kappa[i]);
      coeffs.a2[i] = (coeffs.a[i] * coeffs.a[i] - coeffs.b[i] * coeffs.b[i]) / (2.0 * inp.h);
      coeffs.b2[i] = (coeffs.a[i] * coeffs.b[i]) / inp.h;
      coeffs.mu2[i] = coeffs.f2[i] * std::cosh(inp.h * 2.0 * coeffs.kappa[i]) - inp.h * coeffs.omega1[i];
    }

    const int num_pairs = coeffs.n_comp * (coeffs.n_comp - 1) / 2;
    const int len2 = 2 * num_pairs;
    coeffs.f_npm.assign(len2, 0.0);
    coeffs.g_npm.assign(len2, 0.0);
    coeffs.gamma_npm.assign(len2, 0.0);
    coeffs.a_npm.assign(len2, 0.0);
    coeffs.b_npm.assign(len2, 0.0);
    coeffs.mu_npm.assign(len2, 0.0);
    coeffs.omega_npm.assign(len2, 0.0);
    coeffs.kx_npm.assign(len2, 0.0);
    coeffs.ky_npm.assign(len2, 0.0);

    int pair_count = 0;
    for (int n = 0; n < coeffs.n_comp; ++n) {
      for (int m = n + 1; m < coeffs.n_comp; ++m) {
        const int idx_plus = 2 * pair_count;
        const int idx_minus = idx_plus + 1;
        ++pair_count;
        for (int sign_index = 0; sign_index < 2; ++sign_index) {
          const double pm = sign_index == 0 ? 1.0 : -1.0;
          const int idx = sign_index == 0 ? idx_plus : idx_minus;
          const double omega_out = coeffs.omega1[n] + pm * coeffs.omega1[m];
          const double kx_out = coeffs.kx[n] + pm * coeffs.kx[m];
          const double ky_out = coeffs.ky[n] + pm * coeffs.ky[m];
          const double kappa_out = std::hypot(kx_out, ky_out);
          const double gamma_out = kappa_out * std::sinh(inp.h * kappa_out);
          const double beta_out = omega_out * omega_out * std::cosh(inp.h * kappa_out) - inp.g * kappa_out * std::sinh(inp.h * kappa_out);
          const double f_out = gamma2(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], pm * coeffs.omega1[m], pm * coeffs.kx[m], pm * coeffs.ky[m], coeffs.kappa[m], omega_out, beta_out, inp.g, inp.h);
          const double g_out = lambda2(coeffs.omega1[n], coeffs.kx[n], coeffs.ky[n], coeffs.kappa[n], pm * coeffs.omega1[m], pm * coeffs.kx[m], pm * coeffs.ky[m], coeffs.kappa[m], omega_out, omega_out * std::cosh(inp.h * kappa_out), gamma_out, beta_out, inp.g, inp.h);
          coeffs.omega_npm[idx] = omega_out;
          coeffs.kx_npm[idx] = kx_out;
          coeffs.ky_npm[idx] = ky_out;
          coeffs.f_npm[idx] = f_out;
          coeffs.g_npm[idx] = g_out;
          coeffs.gamma_npm[idx] = gamma_out;
          coeffs.mu_npm[idx] = f_out * std::cosh(inp.h * kappa_out) - 0.5 * inp.h * (coeffs.omega1[n] + pm * coeffs.omega1[m]);
        }
        coeffs.a_npm[idx_plus] = (coeffs.a[n] * coeffs.a[m] - coeffs.b[n] * coeffs.b[m]) / inp.h;
        coeffs.b_npm[idx_plus] = (coeffs.a[m] * coeffs.b[n] + coeffs.a[n] * coeffs.b[m]) / inp.h;
        coeffs.a_npm[idx_minus] = (coeffs.a[n] * coeffs.a[m] + coeffs.b[n] * coeffs.b[m]) / inp.h;
        coeffs.b_npm[idx_minus] = (coeffs.a[m] * coeffs.b[n] - coeffs.a[n] * coeffs.b[m]) / inp.h;
      }
    }
  }

  if (inp.order == 3) {
    add_third_order_terms(coeffs, inp);
  }

  return coeffs;
}

std::vector<std::complex<double>> make_complex(const std::vector<double>& real_part, const std::vector<double>& imag_part) {
  std::vector<std::complex<double>> out(real_part.size());
  for (std::size_t i = 0; i < real_part.size(); ++i) {
    out[i] = std::complex<double>(real_part[i], imag_part[i]);
  }
  return out;
}

std::vector<std::complex<double>> exp_phase(const std::vector<double>& omega, double t) {
  std::vector<std::complex<double>> out(omega.size());
  for (std::size_t i = 0; i < omega.size(); ++i) {
    out[i] = std::exp(std::complex<double>(0.0, -omega[i] * t));
  }
  return out;
}

void accumulate_spectrum(
    std::vector<std::complex<double>>& spec,
    int nx,
    int ny,
    double dkx,
    double dky,
    const std::vector<double>& k_in_x,
    const std::vector<double>& k_in_y,
    const std::vector<std::complex<double>>& values) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    const double ux = k_in_x[i] / dkx;
    const double uy = k_in_y[i] / dky;
    if (!std::isfinite(ux) || !std::isfinite(uy)) {
      continue;
    }
    const long long ix0 = static_cast<long long>(std::floor(ux));
    const long long iy0 = static_cast<long long>(std::floor(uy));
    const double fx = ux - static_cast<double>(ix0);
    const double fy = uy - static_cast<double>(iy0);
    const auto index = [nx, ny](long long iy, long long ix) {
      const long long wrap_y = ((iy % ny) + ny) % ny;
      const long long wrap_x = ((ix % nx) + nx) % nx;
      return static_cast<std::size_t>(wrap_y * nx + wrap_x);
    };
    spec[index(iy0, ix0)] += values[i] * (1.0 - fx) * (1.0 - fy);
    spec[index(iy0, ix0 + 1)] += values[i] * fx * (1.0 - fy);
    spec[index(iy0 + 1, ix0)] += values[i] * (1.0 - fx) * fy;
    spec[index(iy0 + 1, ix0 + 1)] += values[i] * fx * fy;
  }
}

std::vector<std::complex<double>> inverse_dft_1d(const std::vector<std::complex<double>>& in) {
  const std::size_t n = in.size();
  std::vector<std::complex<double>> out(n, std::complex<double>(0.0, 0.0));
  for (std::size_t x = 0; x < n; ++x) {
    std::complex<double> sum(0.0, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
      const double angle = 2.0 * kPi * static_cast<double>(x * k) / static_cast<double>(n);
      sum += in[k] * std::exp(std::complex<double>(0.0, angle));
    }
    out[x] = sum;
  }
  return out;
}

Matrix inverse_dft2_unscaled(const std::vector<std::complex<double>>& spec, int nx, int ny) {
  std::vector<std::complex<double>> temp(static_cast<std::size_t>(nx * ny));
  for (int row = 0; row < ny; ++row) {
    std::vector<std::complex<double>> row_in(nx);
    for (int col = 0; col < nx; ++col) {
      row_in[col] = spec[static_cast<std::size_t>(row * nx + col)];
    }
    const auto row_out = inverse_dft_1d(row_in);
    for (int col = 0; col < nx; ++col) {
      temp[static_cast<std::size_t>(row * nx + col)] = row_out[col];
    }
  }
  Matrix out;
  out.rows = static_cast<std::size_t>(ny);
  out.cols = static_cast<std::size_t>(nx);
  out.values.assign(static_cast<std::size_t>(nx * ny), 0.0);
  for (int col = 0; col < nx; ++col) {
    std::vector<std::complex<double>> col_in(ny);
    for (int row = 0; row < ny; ++row) {
      col_in[row] = temp[static_cast<std::size_t>(row * nx + col)];
    }
    const auto col_out = inverse_dft_1d(col_in);
    for (int row = 0; row < ny; ++row) {
      out.values[static_cast<std::size_t>(row * nx + col)] = col_out[row].real();
    }
  }
  return out;
}

SurfaceFields reconstruct_surface(const SpectralCoefficients& coeffs, const CaseInputs& inp) {
  const int nx = inp.Nx;
  const int ny = inp.Ny;
  const double dx = inp.Lx / static_cast<double>(nx);
  const double dy = inp.Ly / static_cast<double>(ny);
  const double dkx = 2.0 * kPi / inp.Lx;
  const double dky = 2.0 * kPi / inp.Ly;

  std::vector<std::complex<double>> spec_eta(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> spec_phi(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));

  const auto phase_lin = exp_phase(coeffs.omega, inp.t);
  std::vector<std::complex<double>> z_lin(coeffs.n_comp);
  std::vector<std::complex<double>> phi_lin(coeffs.n_comp);
  for (int i = 0; i < coeffs.n_comp; ++i) {
    const std::complex<double> amp(coeffs.a[i], coeffs.b[i]);
    z_lin[i] = amp * phase_lin[i];
    phi_lin[i] = z_lin[i] * std::complex<double>(0.0, coeffs.mu[i] + coeffs.mu_star[i]);
  }
  accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx, coeffs.ky, z_lin);
  accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx, coeffs.ky, phi_lin);

  if (coeffs.order >= 2) {
    std::vector<std::complex<double>> z2(coeffs.n_comp);
    std::vector<std::complex<double>> eta2(coeffs.n_comp);
    std::vector<std::complex<double>> phi2(coeffs.n_comp);
    for (int i = 0; i < coeffs.n_comp; ++i) {
      const std::complex<double> amp(coeffs.a2[i], coeffs.b2[i]);
      z2[i] = amp * std::exp(std::complex<double>(0.0, -coeffs.omega2[i] * inp.t));
      eta2[i] = z2[i] * coeffs.g2[i];
      phi2[i] = z2[i] * std::complex<double>(0.0, coeffs.mu2[i]);
    }
    accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx2, coeffs.ky2, eta2);
    accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx2, coeffs.ky2, phi2);

    const std::size_t pair_len = coeffs.a_npm.size();
    std::vector<std::complex<double>> znpm(pair_len);
    std::vector<std::complex<double>> etanpm(pair_len);
    std::vector<std::complex<double>> phinpm(pair_len);
    for (std::size_t i = 0; i < pair_len; ++i) {
      const std::complex<double> amp(coeffs.a_npm[i], coeffs.b_npm[i]);
      znpm[i] = amp * std::exp(std::complex<double>(0.0, -coeffs.omega_npm[i] * inp.t));
      etanpm[i] = znpm[i] * coeffs.g_npm[i];
      phinpm[i] = znpm[i] * std::complex<double>(0.0, coeffs.mu_npm[i]);
    }
    accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx_npm, coeffs.ky_npm, etanpm);
    accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx_npm, coeffs.ky_npm, phinpm);
  }

  if (coeffs.order >= 3) {
    std::vector<std::complex<double>> z3(coeffs.n_comp);
    std::vector<std::complex<double>> eta3(coeffs.n_comp);
    std::vector<std::complex<double>> phi3(coeffs.n_comp);
    for (int i = 0; i < coeffs.n_comp; ++i) {
      const std::complex<double> amp(coeffs.a3[i], coeffs.b3[i]);
      z3[i] = amp * std::exp(std::complex<double>(0.0, -coeffs.omega3v[i] * inp.t));
      eta3[i] = z3[i] * coeffs.g3[i];
      phi3[i] = z3[i] * std::complex<double>(0.0, coeffs.mu3[i]);
    }
    accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx3, coeffs.ky3, eta3);
    accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx3, coeffs.ky3, phi3);

    const auto accumulate_branch = [&](const std::vector<double>& a_branch, const std::vector<double>& b_branch, const std::vector<double>& omega_branch, const std::vector<double>& kx_branch, const std::vector<double>& ky_branch, const std::vector<double>& g_branch, const std::vector<double>& mu_branch, double amp_scale) {
      const std::size_t len = a_branch.size();
      std::vector<std::complex<double>> z(len);
      std::vector<std::complex<double>> eta(len);
      std::vector<std::complex<double>> phi(len);
      for (std::size_t i = 0; i < len; ++i) {
        const std::complex<double> amp(a_branch[i], b_branch[i]);
        z[i] = amp_scale * amp * std::exp(std::complex<double>(0.0, -omega_branch[i] * inp.t));
        eta[i] = z[i] * g_branch[i];
        phi[i] = z[i] * std::complex<double>(0.0, mu_branch[i]);
      }
      accumulate_spectrum(spec_eta, nx, ny, dkx, dky, kx_branch, ky_branch, eta);
      accumulate_spectrum(spec_phi, nx, ny, dkx, dky, kx_branch, ky_branch, phi);
    };

    accumulate_branch(coeffs.a_np2m, coeffs.b_np2m, coeffs.omega_np2m, coeffs.kx_np2m, coeffs.ky_np2m, coeffs.g_np2m, coeffs.mu_np2m, 1.0);
    accumulate_branch(coeffs.a_2npm, coeffs.b_2npm, coeffs.omega_2npm, coeffs.kx_2npm, coeffs.ky_2npm, coeffs.g_2npm, coeffs.mu_2npm, 1.0);
    accumulate_branch(coeffs.a_npmpp, coeffs.b_npmpp, coeffs.omega_npmpp, coeffs.kx_npmpp, coeffs.ky_npmpp, coeffs.g_npmpp, coeffs.mu_npmpp, 2.0);
  }

  SurfaceFields fields;
  fields.eta = inverse_dft2_unscaled(spec_eta, nx, ny);
  fields.phi = inverse_dft2_unscaled(spec_phi, nx, ny);
  fields.x.rows = 1;
  fields.x.cols = static_cast<std::size_t>(nx);
  fields.x.values.resize(static_cast<std::size_t>(nx));
  fields.y.rows = static_cast<std::size_t>(ny);
  fields.y.cols = 1;
  fields.y.values.resize(static_cast<std::size_t>(ny));

  for (int ix = 0; ix < nx; ++ix) {
    fields.x.values[static_cast<std::size_t>(ix)] = dx * ix;
  }
  for (int iy = 0; iy < ny; ++iy) {
    fields.y.values[static_cast<std::size_t>(iy)] = dy * iy;
  }

  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      fields.phi.values[static_cast<std::size_t>(iy * nx + ix)] += coeffs.Ux * fields.x.values[static_cast<std::size_t>(ix)] + coeffs.Uy * fields.y.values[static_cast<std::size_t>(iy)];
    }
  }
  return fields;
}

}  // namespace

ResultBundle run_case(const LoadedCase& loaded, int repeats, bool warmup) {
  if (repeats < 1) {
    throw std::runtime_error("repeats must be >= 1");
  }

  if (warmup) {
    const auto coeffs = compute_coefficients(loaded);
    (void) reconstruct_surface(coeffs, loaded.manifest.inputs);
  }

  std::vector<double> coeff_times;
  std::vector<double> recon_times;
  SurfaceFields fields;
  for (int i = 0; i < repeats; ++i) {
    const auto t0 = std::chrono::steady_clock::now();
    const auto coeffs = compute_coefficients(loaded);
    const auto t1 = std::chrono::steady_clock::now();
    fields = reconstruct_surface(coeffs, loaded.manifest.inputs);
    const auto t2 = std::chrono::steady_clock::now();
    coeff_times.push_back(std::chrono::duration<double>(t1 - t0).count());
    recon_times.push_back(std::chrono::duration<double>(t2 - t1).count());
  }

  ResultBundle result;
  result.eta = fields.eta;
  result.phi = fields.phi;
  result.x = fields.x;
  result.y = fields.y;
  result.runtime.repeats = repeats;
  result.runtime.warmup = warmup;
  double coeff_sum = 0.0;
  double recon_sum = 0.0;
  result.runtime.best_total_s = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < coeff_times.size(); ++i) {
    coeff_sum += coeff_times[i];
    recon_sum += recon_times[i];
    result.runtime.best_total_s = std::min(result.runtime.best_total_s, coeff_times[i] + recon_times[i]);
  }
  result.runtime.mean_coefficient_s = coeff_sum / static_cast<double>(coeff_times.size());
  result.runtime.mean_reconstruction_s = recon_sum / static_cast<double>(recon_times.size());
  result.runtime.mean_total_s = result.runtime.mean_coefficient_s + result.runtime.mean_reconstruction_s;
  return result;
}

void dump_coefficients_for_case(const LoadedCase& loaded, const std::filesystem::path& output_dir) {
  const auto coeffs = compute_coefficients(loaded);
  const auto& inp = loaded.manifest.inputs;
  std::filesystem::create_directories(output_dir);
  auto save_vec = [&](const std::string& name, const std::vector<double>& values) {
    std::ofstream out(output_dir / (name + ".csv"));
    out.precision(17);
    for (double v : values) {
      out << v << "\n";
    }
  };
  save_vec("omega", coeffs.omega);
  save_vec("omega1", coeffs.omega1);
  save_vec("kappa", coeffs.kappa);
  save_vec("muStar", coeffs.mu_star);
  if (coeffs.order >= 2) {
    save_vec("A_2", coeffs.a2);
    save_vec("B_2", coeffs.b2);
    save_vec("F_2", coeffs.f2);
    save_vec("G_2", coeffs.g2);
    save_vec("mu_2", coeffs.mu2);
    save_vec("A_npm", coeffs.a_npm);
    save_vec("B_npm", coeffs.b_npm);
    save_vec("F_npm", coeffs.f_npm);
    save_vec("G_npm", coeffs.g_npm);
    save_vec("gamma_npm", coeffs.gamma_npm);
    save_vec("mu_npm", coeffs.mu_npm);
    save_vec("omega_npm", coeffs.omega_npm);
    save_vec("kx_npm", coeffs.kx_npm);
    save_vec("ky_npm", coeffs.ky_npm);
    std::vector<double> kappa_npm(coeffs.kx_npm.size(), 0.0);
    for (std::size_t i = 0; i < coeffs.kx_npm.size(); ++i) {
      kappa_npm[i] = std::hypot(coeffs.kx_npm[i], coeffs.ky_npm[i]);
    }
    save_vec("kappa_npm", kappa_npm);
  }
  if (coeffs.order >= 3) {
    save_vec("A_3", coeffs.a3);
    save_vec("B_3", coeffs.b3);
    save_vec("F_3", coeffs.f3);
    save_vec("G_3", coeffs.g3);
    save_vec("mu_3", coeffs.mu3);
    save_vec("A_np2m", coeffs.a_np2m);
    save_vec("B_np2m", coeffs.b_np2m);
    save_vec("F_np2m", coeffs.f_np2m);
    save_vec("G_np2m", coeffs.g_np2m);
    save_vec("mu_np2m", coeffs.mu_np2m);
    save_vec("omega_np2m", coeffs.omega_np2m);
    save_vec("kx_np2m", coeffs.kx_np2m);
    save_vec("ky_np2m", coeffs.ky_np2m);
    save_vec("A_2npm", coeffs.a_2npm);
    save_vec("B_2npm", coeffs.b_2npm);
    save_vec("F_2npm", coeffs.f_2npm);
    save_vec("G_2npm", coeffs.g_2npm);
    save_vec("mu_2npm", coeffs.mu_2npm);
    save_vec("omega_2npm", coeffs.omega_2npm);
    save_vec("kx_2npm", coeffs.kx_2npm);
    save_vec("ky_2npm", coeffs.ky_2npm);
    save_vec("A_npmpp", coeffs.a_npmpp);
    save_vec("B_npmpp", coeffs.b_npmpp);
    save_vec("F_npmpp", coeffs.f_npmpp);
    save_vec("G_npmpp", coeffs.g_npmpp);
    save_vec("mu_npmpp", coeffs.mu_npmpp);
    save_vec("omega_npmpp", coeffs.omega_npmpp);
    save_vec("kx_npmpp", coeffs.kx_npmpp);
    save_vec("ky_npmpp", coeffs.ky_npmpp);
    save_vec("idx_npmpp_nm", coeffs.idx_npmpp_nm);
    save_vec("idx_npmpp_np", coeffs.idx_npmpp_np);
    save_vec("idx_npmpp_mp", coeffs.idx_npmpp_mp);

    const int n_comp = coeffs.n_comp;
    const int num_triplets = n_comp * (n_comp - 1) * (n_comp - 2) / 6;
    std::vector<double> lambda_t0(num_triplets, 0.0);
    std::vector<double> lambda_t1(num_triplets, 0.0);
    std::vector<double> lambda_t2(num_triplets, 0.0);
    std::vector<double> lambda_t3(num_triplets, 0.0);
    std::vector<double> lambda_t4(num_triplets, 0.0);
    std::vector<double> lambda_t5(num_triplets, 0.0);
    std::vector<double> lambda_t6(num_triplets, 0.0);
    std::vector<double> gamma_t0(num_triplets, 0.0);
    std::vector<double> gamma_t1(num_triplets, 0.0);
    std::vector<double> gamma_t2(num_triplets, 0.0);
    std::vector<double> gamma_t3(num_triplets, 0.0);
    std::vector<double> gamma_t4(num_triplets, 0.0);
    std::vector<double> gamma_t5(num_triplets, 0.0);
    std::vector<double> gamma_t6(num_triplets, 0.0);
    std::vector<double> pi_t0(num_triplets, 0.0);
    std::vector<double> pi_t1(num_triplets, 0.0);
    std::vector<double> pi_t2(num_triplets, 0.0);
    std::vector<double> pi_t3(num_triplets, 0.0);

    int c3 = 0;
    for (int n = 0; n < n_comp; ++n) {
      for (int m = n + 1; m < n_comp; ++m) {
        for (int p = m + 1; p < n_comp; ++p) {
          const int idx_nm = static_cast<int>(coeffs.idx_npmpp_nm[c3]);
          const int idx_np = static_cast<int>(coeffs.idx_npmpp_np[c3]);
          const int idx_mp = static_cast<int>(coeffs.idx_npmpp_mp[c3]);
          const double omega_npmpp = coeffs.omega1[n] + coeffs.omega1[m] + coeffs.omega1[p];
          const double kappa_npmpp = std::hypot(coeffs.kx_npmpp[c3], coeffs.ky_npmpp[c3]);
          const double alpha_npmpp = omega_npmpp * std::cosh(coeffs.h * kappa_npmpp);
          const double gamma_npmpp = kappa_npmpp * std::sinh(coeffs.h * kappa_npmpp);
          const double beta_npmpp = omega_npmpp * omega_npmpp * std::cosh(coeffs.h * kappa_npmpp) - coeffs.g * kappa_npmpp * std::sinh(coeffs.h * kappa_npmpp);
          const double knkm = coeffs.kx[n] * coeffs.kx[m] + coeffs.ky[n] * coeffs.ky[m];
          const double knkp = coeffs.kx[n] * coeffs.kx[p] + coeffs.ky[n] * coeffs.ky[p];
          const double kmkp = coeffs.kx[m] * coeffs.kx[p] + coeffs.ky[m] * coeffs.ky[p];
          const double kappanpm = std::hypot(coeffs.kx_npm[idx_nm], coeffs.ky_npm[idx_nm]);
          const double kappanpp = std::hypot(coeffs.kx_npm[idx_np], coeffs.ky_npm[idx_np]);
          const double kappampp = std::hypot(coeffs.kx_npm[idx_mp], coeffs.ky_npm[idx_mp]);
          const double gammanpm = coeffs.gamma_npm[idx_nm];
          const double gammanpp = coeffs.gamma_npm[idx_np];
          const double gammampp = coeffs.gamma_npm[idx_mp];
          const double gnpm = coeffs.g_npm[idx_nm];
          const double gnpp = coeffs.g_npm[idx_np];
          const double gmpp = coeffs.g_npm[idx_mp];
          const double fnpm = coeffs.f_npm[idx_nm];
          const double fnpp = coeffs.f_npm[idx_np];
          const double fmpp = coeffs.f_npm[idx_mp];

          lambda_t0[c3] =
              coeffs.h * coeffs.h / (4.0 * beta_npmpp) *
              (alpha_npmpp *
                   (coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]) +
                    coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]) +
                    coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p])) +
               gamma_npmpp *
                   (coeffs.g / coeffs.omega1[n] * (coeffs.omega1[m] * knkm + coeffs.omega1[p] * knkp - omega_npmpp * coeffs.kappa[n] * coeffs.kappa[n]) +
                    coeffs.g / coeffs.omega1[m] * (coeffs.omega1[n] * knkm + coeffs.omega1[p] * kmkp - omega_npmpp * coeffs.kappa[m] * coeffs.kappa[m]) +
                    coeffs.g / coeffs.omega1[p] * (coeffs.omega1[n] * knkp + coeffs.omega1[m] * kmkp - omega_npmpp * coeffs.kappa[p] * coeffs.kappa[p])));
          lambda_t1[c3] =
              -coeffs.h * fnpm / (2.0 * beta_npmpp) *
              (alpha_npmpp * std::cosh(coeffs.h * kappanpm) * (knkp + kmkp + kappanpm * kappanpm) +
               gamma_npmpp * (coeffs.g / coeffs.omega1[p] * (knkp + kmkp) * std::cosh(coeffs.h * kappanpm) - gammanpm * omega_npmpp));
          lambda_t2[c3] =
              -coeffs.h * fnpp / (2.0 * beta_npmpp) *
              (alpha_npmpp * std::cosh(coeffs.h * kappanpp) * (knkm + kmkp + kappanpp * kappanpp) +
               gamma_npmpp * (coeffs.g / coeffs.omega1[m] * (knkm + kmkp) * std::cosh(coeffs.h * kappanpp) - gammanpp * omega_npmpp));
          lambda_t3[c3] =
              -coeffs.h * fmpp / (2.0 * beta_npmpp) *
              (alpha_npmpp * std::cosh(coeffs.h * kappampp) * (knkm + knkp + kappampp * kappampp) +
               gamma_npmpp * (coeffs.g / coeffs.omega1[n] * (knkm + knkp) * std::cosh(coeffs.h * kappampp) - gammampp * omega_npmpp));
          lambda_t4[c3] = coeffs.h * gnpm / (2.0 * beta_npmpp) * (alpha_npmpp * coeffs.g / coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p]) - gamma_npmpp * coeffs.omega1[p] * coeffs.omega1[p]);
          lambda_t5[c3] = coeffs.h * gnpp / (2.0 * beta_npmpp) * (alpha_npmpp * coeffs.g / coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]) - gamma_npmpp * coeffs.omega1[m] * coeffs.omega1[m]);
          lambda_t6[c3] = coeffs.h * gmpp / (2.0 * beta_npmpp) * (alpha_npmpp * coeffs.g / coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]) - gamma_npmpp * coeffs.omega1[n] * coeffs.omega1[n]);

          gamma_t0[c3] =
              -coeffs.g * coeffs.h * coeffs.h / (4.0 * beta_npmpp) *
              (coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]) +
               coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]) +
               coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p]) +
               omega_npmpp / coeffs.omega1[n] * (coeffs.omega1[m] * knkm + coeffs.omega1[p] * knkp - omega_npmpp * coeffs.kappa[n] * coeffs.kappa[n]) +
               omega_npmpp / coeffs.omega1[m] * (coeffs.omega1[n] * knkm + coeffs.omega1[p] * kmkp - omega_npmpp * coeffs.kappa[m] * coeffs.kappa[m]) +
               omega_npmpp / coeffs.omega1[p] * (coeffs.omega1[n] * knkp + coeffs.omega1[m] * kmkp - omega_npmpp * coeffs.kappa[p] * coeffs.kappa[p]));
          gamma_t1[c3] =
              coeffs.h * fnpm / (2.0 * beta_npmpp) *
              (coeffs.g * std::cosh(coeffs.h * kappanpm) * ((knkp + kmkp + kappanpm * kappanpm) + omega_npmpp / coeffs.omega1[p] * (knkp + kmkp)) -
               gammanpm * omega_npmpp * omega_npmpp);
          gamma_t2[c3] =
              coeffs.h * fnpp / (2.0 * beta_npmpp) *
              (coeffs.g * std::cosh(coeffs.h * kappanpp) * ((knkm + kmkp + kappanpp * kappanpp) + omega_npmpp / coeffs.omega1[m] * (knkm + kmkp)) -
               gammanpp * omega_npmpp * omega_npmpp);
          gamma_t3[c3] =
              coeffs.h * fmpp / (2.0 * beta_npmpp) *
              (coeffs.g * std::cosh(coeffs.h * kappampp) * ((knkm + knkp + kappampp * kappampp) + omega_npmpp / coeffs.omega1[n] * (knkm + knkp)) -
               gammampp * omega_npmpp * omega_npmpp);
          gamma_t4[c3] = coeffs.h * gnpm / (2.0 * beta_npmpp) * (coeffs.omega1[p] * coeffs.omega1[p] * omega_npmpp - coeffs.g * coeffs.g / coeffs.omega1[p] * (knkp + kmkp + coeffs.kappa[p] * coeffs.kappa[p]));
          gamma_t5[c3] = coeffs.h * gnpp / (2.0 * beta_npmpp) * (coeffs.omega1[m] * coeffs.omega1[m] * omega_npmpp - coeffs.g * coeffs.g / coeffs.omega1[m] * (knkm + kmkp + coeffs.kappa[m] * coeffs.kappa[m]));
          gamma_t6[c3] = coeffs.h * gmpp / (2.0 * beta_npmpp) * (coeffs.omega1[n] * coeffs.omega1[n] * omega_npmpp - coeffs.g * coeffs.g / coeffs.omega1[n] * (knkm + knkp + coeffs.kappa[n] * coeffs.kappa[n]));

          pi_t0[c3] = coeffs.f_npmpp[c3] * std::cosh(coeffs.h * kappa_npmpp);
          pi_t1[c3] = -coeffs.g * coeffs.h * coeffs.h / 4.0 * (coeffs.kappa[n] * coeffs.kappa[n] / coeffs.omega1[n] + coeffs.kappa[m] * coeffs.kappa[m] / coeffs.omega1[m] + coeffs.kappa[p] * coeffs.kappa[p] / coeffs.omega1[p]);
          pi_t2[c3] = -coeffs.h / 2.0 * (coeffs.omega1[n] * gmpp + coeffs.omega1[m] * gnpp + coeffs.omega1[p] * gnpm);
          pi_t3[c3] = coeffs.h / 2.0 * (fnpm * gammanpm + fnpp * gammanpp + fmpp * gammampp);
          ++c3;
        }
      }
    }

    save_vec("debug_lambda3_t0_npmpp", lambda_t0);
    save_vec("debug_lambda3_t1_npmpp", lambda_t1);
    save_vec("debug_lambda3_t2_npmpp", lambda_t2);
    save_vec("debug_lambda3_t3_npmpp", lambda_t3);
    save_vec("debug_lambda3_t4_npmpp", lambda_t4);
    save_vec("debug_lambda3_t5_npmpp", lambda_t5);
    save_vec("debug_lambda3_t6_npmpp", lambda_t6);
    save_vec("debug_gamma3_t0_npmpp", gamma_t0);
    save_vec("debug_gamma3_t1_npmpp", gamma_t1);
    save_vec("debug_gamma3_t2_npmpp", gamma_t2);
    save_vec("debug_gamma3_t3_npmpp", gamma_t3);
    save_vec("debug_gamma3_t4_npmpp", gamma_t4);
    save_vec("debug_gamma3_t5_npmpp", gamma_t5);
    save_vec("debug_gamma3_t6_npmpp", gamma_t6);
    save_vec("debug_pi_t0_npmpp", pi_t0);
    save_vec("debug_pi_t1_npmpp", pi_t1);
    save_vec("debug_pi_t2_npmpp", pi_t2);
    save_vec("debug_pi_t3_npmpp", pi_t3);
  }

  const int nx = inp.Nx;
  const int ny = inp.Ny;
  const double dkx = 2.0 * kPi / inp.Lx;
  const double dky = 2.0 * kPi / inp.Ly;
  std::vector<std::complex<double>> spec_eta(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> spec_phi(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
  const auto phase_lin = exp_phase(coeffs.omega, inp.t);
  std::vector<std::complex<double>> z_lin(coeffs.n_comp);
  std::vector<std::complex<double>> phi_lin(coeffs.n_comp);
  for (int i = 0; i < coeffs.n_comp; ++i) {
    const std::complex<double> amp(coeffs.a[i], coeffs.b[i]);
    z_lin[i] = amp * phase_lin[i];
    phi_lin[i] = z_lin[i] * std::complex<double>(0.0, coeffs.mu[i] + coeffs.mu_star[i]);
  }
  accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx, coeffs.ky, z_lin);
  accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx, coeffs.ky, phi_lin);
  if (coeffs.order >= 2) {
    std::vector<std::complex<double>> z2(coeffs.n_comp);
    std::vector<std::complex<double>> eta2(coeffs.n_comp);
    std::vector<std::complex<double>> phi2(coeffs.n_comp);
    for (int i = 0; i < coeffs.n_comp; ++i) {
      const std::complex<double> amp(coeffs.a2[i], coeffs.b2[i]);
      z2[i] = amp * std::exp(std::complex<double>(0.0, -coeffs.omega2[i] * inp.t));
      eta2[i] = z2[i] * coeffs.g2[i];
      phi2[i] = z2[i] * std::complex<double>(0.0, coeffs.mu2[i]);
    }
    accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx2, coeffs.ky2, eta2);
    accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx2, coeffs.ky2, phi2);

    const std::size_t pair_len = coeffs.a_npm.size();
    std::vector<std::complex<double>> znpm(pair_len);
    std::vector<std::complex<double>> etanpm(pair_len);
    std::vector<std::complex<double>> phinpm(pair_len);
    for (std::size_t i = 0; i < pair_len; ++i) {
      const std::complex<double> amp(coeffs.a_npm[i], coeffs.b_npm[i]);
      znpm[i] = amp * std::exp(std::complex<double>(0.0, -coeffs.omega_npm[i] * inp.t));
      etanpm[i] = znpm[i] * coeffs.g_npm[i];
      phinpm[i] = znpm[i] * std::complex<double>(0.0, coeffs.mu_npm[i]);
    }
    accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx_npm, coeffs.ky_npm, etanpm);
    accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx_npm, coeffs.ky_npm, phinpm);
  }
  if (coeffs.order >= 3) {
    std::vector<std::complex<double>> spec_phi_3(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> spec_phi_np2m(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> spec_phi_2npm(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> spec_phi_npmpp(static_cast<std::size_t>(nx * ny), std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> z3(coeffs.n_comp);
    std::vector<std::complex<double>> eta3(coeffs.n_comp);
    std::vector<std::complex<double>> phi3(coeffs.n_comp);
    for (int i = 0; i < coeffs.n_comp; ++i) {
      const std::complex<double> amp(coeffs.a3[i], coeffs.b3[i]);
      z3[i] = amp * std::exp(std::complex<double>(0.0, -coeffs.omega3v[i] * inp.t));
      eta3[i] = z3[i] * coeffs.g3[i];
      phi3[i] = z3[i] * std::complex<double>(0.0, coeffs.mu3[i]);
    }
    accumulate_spectrum(spec_eta, nx, ny, dkx, dky, coeffs.kx3, coeffs.ky3, eta3);
    accumulate_spectrum(spec_phi, nx, ny, dkx, dky, coeffs.kx3, coeffs.ky3, phi3);
    accumulate_spectrum(spec_phi_3, nx, ny, dkx, dky, coeffs.kx3, coeffs.ky3, phi3);

    const auto accumulate_branch = [&](const std::vector<double>& a_branch, const std::vector<double>& b_branch, const std::vector<double>& omega_branch, const std::vector<double>& kx_branch, const std::vector<double>& ky_branch, const std::vector<double>& g_branch, const std::vector<double>& mu_branch, double amp_scale, std::vector<std::complex<double>>& spec_phi_branch) {
      const std::size_t len = a_branch.size();
      std::vector<std::complex<double>> z(len);
      std::vector<std::complex<double>> eta(len);
      std::vector<std::complex<double>> phi(len);
      for (std::size_t i = 0; i < len; ++i) {
        const std::complex<double> amp(a_branch[i], b_branch[i]);
        z[i] = amp_scale * amp * std::exp(std::complex<double>(0.0, -omega_branch[i] * inp.t));
        eta[i] = z[i] * g_branch[i];
        phi[i] = z[i] * std::complex<double>(0.0, mu_branch[i]);
      }
      accumulate_spectrum(spec_eta, nx, ny, dkx, dky, kx_branch, ky_branch, eta);
      accumulate_spectrum(spec_phi, nx, ny, dkx, dky, kx_branch, ky_branch, phi);
      accumulate_spectrum(spec_phi_branch, nx, ny, dkx, dky, kx_branch, ky_branch, phi);
    };
    accumulate_branch(coeffs.a_np2m, coeffs.b_np2m, coeffs.omega_np2m, coeffs.kx_np2m, coeffs.ky_np2m, coeffs.g_np2m, coeffs.mu_np2m, 1.0, spec_phi_np2m);
    accumulate_branch(coeffs.a_2npm, coeffs.b_2npm, coeffs.omega_2npm, coeffs.kx_2npm, coeffs.ky_2npm, coeffs.g_2npm, coeffs.mu_2npm, 1.0, spec_phi_2npm);
    accumulate_branch(coeffs.a_npmpp, coeffs.b_npmpp, coeffs.omega_npmpp, coeffs.kx_npmpp, coeffs.ky_npmpp, coeffs.g_npmpp, coeffs.mu_npmpp, 2.0, spec_phi_npmpp);

    auto save_complex_vec = [&](const std::string& stem, const std::vector<std::complex<double>>& values) {
      std::vector<double> re(values.size(), 0.0);
      std::vector<double> im(values.size(), 0.0);
      for (std::size_t i = 0; i < values.size(); ++i) {
        re[i] = values[i].real();
        im[i] = values[i].imag();
      }
      save_vec(stem + "_real", re);
      save_vec(stem + "_imag", im);
    };
    save_complex_vec("spec_phi_branch_3", spec_phi_3);
    save_complex_vec("spec_phi_branch_np2m", spec_phi_np2m);
    save_complex_vec("spec_phi_branch_2npm", spec_phi_2npm);
    save_complex_vec("spec_phi_branch_npmpp", spec_phi_npmpp);
  }
  std::vector<double> spec_eta_real(spec_eta.size(), 0.0);
  std::vector<double> spec_eta_imag(spec_eta.size(), 0.0);
  std::vector<double> spec_phi_real(spec_phi.size(), 0.0);
  std::vector<double> spec_phi_imag(spec_phi.size(), 0.0);
  for (std::size_t i = 0; i < spec_eta.size(); ++i) {
    spec_eta_real[i] = spec_eta[i].real();
    spec_eta_imag[i] = spec_eta[i].imag();
    spec_phi_real[i] = spec_phi[i].real();
    spec_phi_imag[i] = spec_phi[i].imag();
  }
  save_vec("spec_eta_real", spec_eta_real);
  save_vec("spec_eta_imag", spec_eta_imag);
  save_vec("spec_phi_real", spec_phi_real);
  save_vec("spec_phi_imag", spec_phi_imag);
}

}  // namespace mf12_cpp
