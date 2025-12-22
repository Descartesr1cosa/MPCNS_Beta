#include <algorithm>
#include <cmath>
#include <cstdint>

// Z3_HallMHD
#include "HallMHD_Solver.h"

#if Reconstruction_Scheme == 1
void HallMHDSolver::Reconstruction(double *metric, int32_t direction,
                                   FieldBlock &PV, FieldBlock &U, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k,
                                   double *out_flux)
{
    auto calc_Jac_radius_GCL = [&](double &out, double *uu, double *pv, double *B, double *metric)
    {
        double rho, p, u, v, w;
        double BB2, Bx, By, Bz;
        double K1, K2, K3;
        rho = uu[0];
        u = pv[0];
        v = pv[1];
        w = pv[2];
        p = pv[3];
        Bx = B[0];
        By = B[1];
        Bz = B[2];

        BB2 = Bx * Bx + By * By + Bz * Bz;

        K1 = metric[0];
        K2 = metric[1];
        K3 = metric[2];

        double uvw = K1 * u + K2 * v + K3 * w;
        // double cc1 = sqrt((gamma_ * p / rho) * (K1 * K1 + K2 * K2 + K3 * K3));
        double cc = sqrt((gamma_ * p / rho + BB2 / rho * inver_MA2) * (K1 * K1 + K2 * K2 + K3 * K3));
        // constexpr double C_hall_safe = 1.5;
        // double c_hall = C_hall_safe * ion_inertial_len * sqrt(BB2 / rho * inver_MA2 * (K1 * K1 + K2 * K2 + K3 * K3)); // ≈ Jac * v_A * d_i * |k|
        // out = fabs(uvw) + cc + c_hall;
        out = fabs(uvw) + cc;
        return;
    };

    // 注意这里的B长度为4，最后一个为B_Jac_nabla\xi eta zeta
    auto calc_Jac_Flux_GCL = [&](double *flux, double *uu, double *pv, double *B, double *metric)
    {
        double k1, k2, k3; // GCL 这里为Jac *k1, Jac *k2, Jac *k3
        k1 = metric[0];
        k2 = metric[1];
        k3 = metric[2];

        double rho, p, u, v, w, uvw, rhoe;
        double Bx, By, Bz, B_Jac_nabla, P_B;

        rho = uu[0];
        u = pv[0];
        v = pv[1];
        w = pv[2];
        p = pv[3];

        Bx = B[0];
        By = B[1];
        Bz = B[2];
        B_Jac_nabla = B[3]; // Bx * k1 + By * k2 + Bz * k3;
        uvw = k1 * u + k2 * v + k3 * w;

        P_B = 0.5 * (Bx * Bx + By * By + Bz * Bz) * inver_MA2;
        rhoe = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

        // flux[0] = rho * uvw;
        // flux[1] = (rho * uvw * u + k1 * p); // Euler Flux - Maxwell Tensor (B\otimesB - B^2/2 I)\cdot J\nabla \xi\eta\zeta
        // flux[2] = (rho * uvw * v + k2 * p);
        // flux[3] = (rho * uvw * w + k3 * p);
        // flux[4] = uvw * (rhoe + p); // Euler Flux

        flux[0] = rho * uvw;
        flux[1] = (rho * uvw * u + k1 * (p + P_B) - inver_MA2 * B_Jac_nabla * Bx); // Euler Flux - Maxwell Tensor (B\otimesB - B^2/2 I)\cdot J\nabla \xi\eta\zeta
        flux[2] = (rho * uvw * v + k2 * (p + P_B) - inver_MA2 * B_Jac_nabla * By);
        flux[3] = (rho * uvw * w + k3 * (p + P_B) - inver_MA2 * B_Jac_nabla * Bz);
        flux[4] = uvw * (rhoe + p);                                                          // Euler Flux
        flux[4] += (2.0 * P_B * uvw - inver_MA2 * (Bx * u + By * v + Bz * w) * B_Jac_nabla); // S(Poynting Vector) \cdot J\nabla \xi\eta\zeta
        // u += H[0];
        // v += H[1];
        // w += H[2];
        // uvw = k1 * u + k2 * v + k3 * w; // With Hall Effect

        // flux[4] += inver_MA2 * (2.0 * P_B * uvw - (Bx * u + By * v + Bz * w) * B_Jac_nabla); // S(Poynting Vector) \cdot J\nabla \xi\eta\zeta
        flux[5] = uvw * Bx - B_Jac_nabla * u;
        flux[6] = uvw * By - B_Jac_nabla * v;
        flux[7] = uvw * Bz - B_Jac_nabla * w;
    };

    int i = index_i;
    int j = index_j;
    int k = index_k;

    double UL[8], UR[8], ppvvL[4], ppvvR[4], BL[4], BR[4];
    double radius[2];

    auto fill_state = [&](int ic, int jc, int kc, double *Ucons, double *pv, double *B)
    {
        double rho = U(ic, jc, kc, 0);
        double u = PV(ic, jc, kc, 0);
        double v = PV(ic, jc, kc, 1);
        double w = PV(ic, jc, kc, 2);
        double p = PV(ic, jc, kc, 3);
        double Bx = B_cell(ic, jc, kc, 0); // including B_add
        double By = B_cell(ic, jc, kc, 1); // including B_add
        double Bz = B_cell(ic, jc, kc, 2); // including B_add

        double inner_product = Bx * metric[0] + By * metric[1] + Bz * metric[2];
        // double inner_product_add = B_add_x * metric[0] + B_add_y * metric[1] + B_add_z * metric[2];
        double inver_norm2 = 1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);

        double B_jac_total = B_jac_nabla; // 法向通量仅仅为induced部分
        // 修正Bx By Bz保持法向分量的通量与感应CT值一致
        // double dBx, dBy, dBz;
        // dBx = (inner_product - B_jac_total) * inver_norm2 * metric[0];
        // dBy = (inner_product - B_jac_total) * inver_norm2 * metric[1];
        // dBz = (inner_product - B_jac_total) * inver_norm2 * metric[2];
        // double dB = sqrt(dBx * dBx + dBy * dBy + dBz * dBz);
        // double B_mag = sqrt(Bx * Bx + By * By + Bz * Bz + 1E-20);
        // double factor = 0.2;
        // if (B_mag > 1)
        // {
        //     if (dB / B_mag > factor)
        //     {
        //         double alpha = factor / dB * B_mag;
        //         Bx -= alpha * dBx;
        //         By -= alpha * dBy;
        //         Bz -= alpha * dBz;
        //     }
        //     else
        //     {
        //         Bx -= dBx;
        //         By -= dBy;
        //         Bz -= dBz;
        //     }
        // }
        // Bx -= (inner_product - B_jac_total) * inver_norm2 * metric[0];
        // By -= (inner_product - B_jac_total) * inver_norm2 * metric[1];
        // Bz -= (inner_product - B_jac_total) * inver_norm2 * metric[2];

        pv[0] = u;
        pv[1] = v;
        pv[2] = w;
        pv[3] = p;

        Ucons[0] = rho;
        Ucons[1] = rho * u;
        Ucons[2] = rho * v;
        Ucons[3] = rho * w;
        Ucons[4] = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w) + 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);
        Ucons[5] = Bx;
        Ucons[6] = By;
        Ucons[7] = Bz;

        B[0] = Bx;
        B[1] = By;
        B[2] = Bz;
        B[3] = B_jac_total;
    };

    if (direction == 0)
    {
        int iL = index_i - 1;
        int iR = index_i;
        fill_state(iL, j, k, UL, ppvvL, BL);
        fill_state(iR, j, k, UR, ppvvR, BR);
    }
    else if (direction == 1)
    {
        int jL = index_j - 1;
        int jR = index_j;
        fill_state(i, jL, k, UL, ppvvL, BL);
        fill_state(i, jR, k, UR, ppvvR, BR);
    }
    else
    { // direction == 2
        int kL = index_k - 1;
        int kR = index_k;
        fill_state(i, j, kL, UL, ppvvL, BL);
        fill_state(i, j, kR, UR, ppvvR, BR);
    }

    calc_Jac_radius_GCL(radius[0], UL, ppvvL, BL, metric);
    calc_Jac_radius_GCL(radius[1], UR, ppvvR, BR, metric);

    double radius_max = std::max(radius[0], radius[1]);

    double FL[8], FR[8];
    calc_Jac_Flux_GCL(FL, UL, ppvvL, BL, metric);
    calc_Jac_Flux_GCL(FR, UR, ppvvR, BR, metric);

    for (int m = 0; m < 8; ++m)
        out_flux[m] = 0.5 * (FL[m] + FR[m]) - 0.5 * radius_max * (UR[m] - UL[m]);

    double Elec_flux[3] = {out_flux[5], out_flux[6], out_flux[7]};
    double norm2 = -1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);
    out_flux[5] = norm2 * (metric[1] * Elec_flux[2] - metric[2] * Elec_flux[1]); // Averaged Electric in Face xi eta zeta
    out_flux[6] = norm2 * (metric[2] * Elec_flux[0] - metric[0] * Elec_flux[2]); // Averaged Electric in Face xi eta zeta
    out_flux[7] = norm2 * (metric[0] * Elec_flux[1] - metric[1] * Elec_flux[0]); // Averaged Electric in Face xi eta zeta
}

#endif