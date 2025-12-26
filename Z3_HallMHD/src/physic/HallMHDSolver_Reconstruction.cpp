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
#elif Reconstruction_Scheme == 2
void HallMHDSolver::Reconstruction(double *metric, int32_t direction,
                                   FieldBlock &PV, FieldBlock &U, FieldBlock &B_cell, double B_jac_nabla,
                                   int iblock, int index_i, int index_j, int index_k,
                                   double *out_flux)
{
    auto minmod = [](double a, double b) -> double
    {
        if (a * b <= 0.0)
            return 0.0;
        return (std::abs(a) < std::abs(b)) ? a : b;
    };

    struct Prim8
    {
        double rho, u, v, w, p, Bx, By, Bz;
    };

    auto get_prim = [&](int ic, int jc, int kc) -> Prim8
    {
        Prim8 q;
        q.rho = U(ic, jc, kc, 0);
        q.u = PV(ic, jc, kc, 0);
        q.v = PV(ic, jc, kc, 1);
        q.w = PV(ic, jc, kc, 2);
        q.p = PV(ic, jc, kc, 3);
        q.Bx = B_cell(ic, jc, kc, 0);
        q.By = B_cell(ic, jc, kc, 1);
        q.Bz = B_cell(ic, jc, kc, 2);
        return q;
    };

    // minmod slope in the given direction (in computational index space)
    auto slope_dir = [&](int ic, int jc, int kc) -> Prim8
    {
        int di = (direction == 0) ? 1 : 0;
        int dj = (direction == 1) ? 1 : 0;
        int dk = (direction == 2) ? 1 : 0;

        Prim8 qm = get_prim(ic - di, jc - dj, kc - dk);
        Prim8 q0 = get_prim(ic, jc, kc);
        Prim8 qp = get_prim(ic + di, jc + dj, kc + dk);

        Prim8 s;
        s.rho = minmod(q0.rho - qm.rho, qp.rho - q0.rho);
        s.u = minmod(q0.u - qm.u, qp.u - q0.u);
        s.v = minmod(q0.v - qm.v, qp.v - q0.v);
        s.w = minmod(q0.w - qm.w, qp.w - q0.w);
        s.p = minmod(q0.p - qm.p, qp.p - q0.p);
        s.Bx = minmod(q0.Bx - qm.Bx, qp.Bx - q0.Bx);
        s.By = minmod(q0.By - qm.By, qp.By - q0.By);
        s.Bz = minmod(q0.Bz - qm.Bz, qp.Bz - q0.Bz);
        return s;
    };

    // enforce (B · metric) = B_jac_nabla by projecting along metric
    auto project_B_normal = [&](double &Bx, double &By, double &Bz)
    {
        double n1 = metric[0], n2 = metric[1], n3 = metric[2];
        double dot = Bx * n1 + By * n2 + Bz * n3;
        double invn2 = 1.0 / (n1 * n1 + n2 * n2 + n3 * n3 + 1e-20);
        double delta = (dot - B_jac_nabla) * invn2;
        Bx -= delta * n1;
        By -= delta * n2;
        Bz -= delta * n3;
    };

    auto fill_state_from_prim = [&](const Prim8 &q, double *Ucons, double *pv, double *B)
    {
        double rho = q.rho;
        double u = q.u;
        double v = q.v;
        double w = q.w;
        double p = q.p;
        double Bx = q.Bx;
        double By = q.By;
        double Bz = q.Bz;

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
        B[3] = B_jac_nabla; // face normal (given by CT)
    };

    // ----- your existing lambdas (radius/flux) unchanged -----
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
        double cc = sqrt((gamma_ * p / rho + BB2 / rho * inver_MA2) * (K1 * K1 + K2 * K2 + K3 * K3));
        out = fabs(uvw) + cc;
        return;
    };

    auto calc_Jac_Flux_GCL = [&](double *flux, double *uu, double *pv, double *B, double *metric)
    {
        double k1, k2, k3;
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
        B_Jac_nabla = B[3];
        uvw = k1 * u + k2 * v + k3 * w;

        P_B = 0.5 * (Bx * Bx + By * By + Bz * Bz) * inver_MA2;
        rhoe = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

        flux[0] = rho * uvw;
        flux[1] = (rho * uvw * u + k1 * (p + P_B) - inver_MA2 * B_Jac_nabla * Bx);
        flux[2] = (rho * uvw * v + k2 * (p + P_B) - inver_MA2 * B_Jac_nabla * By);
        flux[3] = (rho * uvw * w + k3 * (p + P_B) - inver_MA2 * B_Jac_nabla * Bz);
        flux[4] = uvw * (rhoe + p);
        flux[4] += (2.0 * P_B * uvw - inver_MA2 * (Bx * u + By * v + Bz * w) * B_Jac_nabla);

        flux[5] = uvw * Bx - B_Jac_nabla * u;
        flux[6] = uvw * By - B_Jac_nabla * v;
        flux[7] = uvw * Bz - B_Jac_nabla * w;
    };

    int i = index_i, j = index_j, k = index_k;

    // Which two cells share this face?
    int di = (direction == 0) ? 1 : 0;
    int dj = (direction == 1) ? 1 : 0;
    int dk = (direction == 2) ? 1 : 0;

    int iL = i - di, jL = j - dj, kL = k - dk; // left cell center
    int iR = i, jR = j, kR = k;                // right cell center

    // Cell-centered primitives
    Prim8 qL0 = get_prim(iL, jL, kL);
    Prim8 qR0 = get_prim(iR, jR, kR);

    // Limited slopes
    Prim8 sL = slope_dir(iL, jL, kL);
    Prim8 sR = slope_dir(iR, jR, kR);

    // Reconstruct to face: q_{L}^{+} and q_{R}^{-}
    Prim8 qLf = qL0;
    Prim8 qRf = qR0;

    qLf.rho += 0.5 * sL.rho;
    qRf.rho -= 0.5 * sR.rho;
    qLf.u += 0.5 * sL.u;
    qRf.u -= 0.5 * sR.u;
    qLf.v += 0.5 * sL.v;
    qRf.v -= 0.5 * sR.v;
    qLf.w += 0.5 * sL.w;
    qRf.w -= 0.5 * sR.w;
    qLf.p += 0.5 * sL.p;
    qRf.p -= 0.5 * sR.p;
    qLf.Bx += 0.5 * sL.Bx;
    qRf.Bx -= 0.5 * sR.Bx;
    qLf.By += 0.5 * sL.By;
    qRf.By -= 0.5 * sR.By;
    qLf.Bz += 0.5 * sL.Bz;
    qRf.Bz -= 0.5 * sR.Bz;

    // Optional but strongly recommended: positivity floors
    // Use your existing p_floor/rho_floor strategy; placeholders here:
    const double rho_floor = 1e-14;
    const double p_floor = 1e-14;
    qLf.rho = std::max(qLf.rho, rho_floor);
    qRf.rho = std::max(qRf.rho, rho_floor);
    qLf.p = std::max(qLf.p, p_floor);
    qRf.p = std::max(qRf.p, p_floor);

    // Enforce face-normal B consistency with CT provided B_jac_nabla
    project_B_normal(qLf.Bx, qLf.By, qLf.Bz);
    project_B_normal(qRf.Bx, qRf.By, qRf.Bz);

    // Build UU/PV/B arrays for flux routines
    double UL[8], UR[8], ppvvL[4], ppvvR[4], BL[4], BR[4];
    fill_state_from_prim(qLf, UL, ppvvL, BL);
    fill_state_from_prim(qRf, UR, ppvvR, BR);

    double radius[2];
    calc_Jac_radius_GCL(radius[0], UL, ppvvL, BL, metric);
    calc_Jac_radius_GCL(radius[1], UR, ppvvR, BR, metric);
    double radius_max = std::max(radius[0], radius[1]);

    double FL[8], FR[8];
    calc_Jac_Flux_GCL(FL, UL, ppvvL, BL, metric);
    calc_Jac_Flux_GCL(FR, UR, ppvvR, BR, metric);

    for (int m = 0; m < 8; ++m)
        out_flux[m] = 0.5 * (FL[m] + FR[m]) - 0.5 * radius_max * (UR[m] - UL[m]);

    // your existing conversion: induction flux -> face-centered electric field
    double Elec_flux[3] = {out_flux[5], out_flux[6], out_flux[7]};
    double norm2 = -1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);
    out_flux[5] = norm2 * (metric[1] * Elec_flux[2] - metric[2] * Elec_flux[1]);
    out_flux[6] = norm2 * (metric[2] * Elec_flux[0] - metric[0] * Elec_flux[2]);
    out_flux[7] = norm2 * (metric[0] * Elec_flux[1] - metric[1] * Elec_flux[0]);
}
#endif