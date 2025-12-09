#include "MHD_Solver.h"

void MHDSolver::inv_fluid()
{
    auto &Jac_ = fld_->field(fid_Jac);
    auto &Xi_ = fld_->field(fid_Xi);
    auto &Eta_ = fld_->field(fid_Eta);
    auto &Zeta_ = fld_->field(fid_Zeta);

    Int3 sub, sup;
    double metric[3];
    double Flux[8];

    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &U = fld_->field(fid_U, iblk);
        auto &PV = fld_->field(fid_PV, iblk);
        auto &Bcell = fld_->field(fid_Bcell, iblk);
        int ncomp = U.descriptor().ncomp;

        // 计算XI面的通量
        {
            auto &flux_xi = fld_->field("F_xi", iblk);
            auto &E_face_xi = fld_->field("E_face_xi", iblk);
            auto &B_face = fld_->field("B_xi", iblk);
            auto &XI = Xi_[iblk];
            sub = flux_xi.inner_lo();
            sup = flux_xi.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        metric[0] = XI(i, j, k, 0);
                        metric[1] = XI(i, j, k, 1);
                        metric[2] = XI(i, j, k, 2);
                        Reconstruction(metric, 0, PV, U, Bcell, B_face(i, j, k, 0), iblk, i, j, k, Flux);
                        for (int m = 0; m < ncomp; m++)
                            flux_xi(i, j, k, m) = Flux[m];

                        E_face_xi(i, j, k, 0) = Flux[5];
                        E_face_xi(i, j, k, 1) = Flux[6];
                        E_face_xi(i, j, k, 2) = Flux[7];
                    }
        }

        // 计算ETA面的通量
        {
            auto &flux_eta = fld_->field("F_eta", iblk);
            auto &E_face_eta = fld_->field("E_face_eta", iblk);
            auto &B_face = fld_->field("B_eta", iblk);
            auto &ET = Eta_[iblk];
            sub = flux_eta.inner_lo();
            sup = flux_eta.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        metric[0] = ET(i, j, k, 0);
                        metric[1] = ET(i, j, k, 1);
                        metric[2] = ET(i, j, k, 2);

                        Reconstruction(metric, 1, PV, U, Bcell, B_face(i, j, k, 0), iblk, i, j, k, Flux);
                        for (int m = 0; m < ncomp; m++)
                            flux_eta(i, j, k, m) = Flux[m];

                        E_face_eta(i, j, k, 0) = Flux[5];
                        E_face_eta(i, j, k, 1) = Flux[6];
                        E_face_eta(i, j, k, 2) = Flux[7];
                    }
        }

        // 计算ZETA面的通量
        {
            auto &flux_zeta = fld_->field("F_zeta", iblk);
            auto &E_face_zeta = fld_->field("E_face_zeta", iblk);
            auto &B_face = fld_->field("B_zeta", iblk);
            auto &ZETA = Zeta_[iblk];
            sub = flux_zeta.inner_lo();
            sup = flux_zeta.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        metric[0] = ZETA(i, j, k, 0);
                        metric[1] = ZETA(i, j, k, 1);
                        metric[2] = ZETA(i, j, k, 2);

                        Reconstruction(metric, 2, PV, U, Bcell, B_face(i, j, k, 0), iblk, i, j, k, Flux);
                        for (int m = 0; m < ncomp; m++)
                            flux_zeta(i, j, k, m) = Flux[m];
                        E_face_zeta(i, j, k, 0) = Flux[5];
                        E_face_zeta(i, j, k, 1) = Flux[6];
                        E_face_zeta(i, j, k, 2) = Flux[7];
                    }
        }
    }
    // 计算RHS
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &RHS = fld_->field("RHS", iblk);
        int ncomp = RHS.descriptor().ncomp;

        auto &Jac = Jac_[iblk];

        auto &flux_xi = fld_->field("F_xi", iblk);
        auto &flux_eta = fld_->field("F_eta", iblk);
        auto &flux_zeta = fld_->field("F_zeta", iblk);
        sub = RHS.inner_lo();
        sup = RHS.inner_hi();

        double inver_jac = 0.0;

        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    inver_jac = 1.0 / Jac(i, j, k, 0);
                    for (int m = 0; m < ncomp; m++)
                    {

                        RHS(i, j, k, m) -= (flux_xi(i + 1, j, k, m) - flux_xi(i, j, k, m)) * inver_jac;
                        RHS(i, j, k, m) -= (flux_eta(i, j + 1, k, m) - flux_eta(i, j, k, m)) * inver_jac;
                        RHS(i, j, k, m) -= (flux_zeta(i, j, k + 1, m) - flux_zeta(i, j, k, m)) * inver_jac;
                    }
                }
    }
}

#if Reconstruction_Scheme == 1
void MHDSolver::Reconstruction(double *metric, int32_t direction,
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
        double inner_product_add = B_add_x * metric[0] + B_add_y * metric[1] + B_add_z * metric[2];
        double inver_norm2 = 1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);

        double B_jac_total = B_jac_nabla + inner_product_add; // 法向通量仅仅为induced部分
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

    // 对最后三个变量，B的通量进行旋转，获得face上的电场
    // E_tangetial  =  - K \times Flux / |K^2|     K=\nabla\xi \eta \zeta
    double Elec_flux[3] = {out_flux[5], out_flux[6], out_flux[7]};
    double norm2 = -1.0 / (metric[0] * metric[0] + metric[1] * metric[1] + metric[2] * metric[2] + 1E-20);
    out_flux[5] = norm2 * (metric[1] * Elec_flux[2] - metric[2] * Elec_flux[1]); // Averaged Electric in Face xi eta zeta
    out_flux[6] = norm2 * (metric[2] * Elec_flux[0] - metric[0] * Elec_flux[2]); // Averaged Electric in Face xi eta zeta
    out_flux[7] = norm2 * (metric[0] * Elec_flux[1] - metric[1] * Elec_flux[0]); // Averaged Electric in Face xi eta zeta
}

#endif
void MHDSolver::inv_induce()
{
    Int3 sub, sup;
    // 插值计算电场E=-u\times B, 存储的E_xi eta zeta均为E\cdot dr的线积分量
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {

        auto &PV = fld_->field(fid_PV, iblk);
        auto &Bcell = fld_->field(fid_Bcell, iblk);
        Double3 vel, B, E, dr;
        double3D &x = fld_->grd->grids(iblk).x;
        double3D &y = fld_->grd->grids(iblk).y;
        double3D &z = fld_->grd->grids(iblk).z;

        {
            auto &Exi = fld_->field("E_xi", iblk);
            auto &E_face_eta = fld_->field("E_face_eta", iblk);
            auto &E_face_zeta = fld_->field("E_face_zeta", iblk);
            sub = Exi.inner_lo();
            sup = Exi.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {

                        E.vector[0] = 0.25 * (E_face_eta(i, j, k, 0) + E_face_eta(i, j, k - 1, 0) + E_face_zeta(i, j, k, 0) + E_face_zeta(i, j - 1, k, 0));
                        E.vector[1] = 0.25 * (E_face_eta(i, j, k, 1) + E_face_eta(i, j, k - 1, 1) + E_face_zeta(i, j, k, 1) + E_face_zeta(i, j - 1, k, 1));
                        E.vector[2] = 0.25 * (E_face_eta(i, j, k, 2) + E_face_eta(i, j, k - 1, 2) + E_face_zeta(i, j, k, 2) + E_face_zeta(i, j - 1, k, 2));

                        dr = {x(i + 1, j, k) - x(i, j, k),
                              y(i + 1, j, k) - y(i, j, k),
                              z(i + 1, j, k) - z(i, j, k)};
                        Exi(i, j, k, 0) = E * dr;
                    }
        }

        {
            auto &Eeta = fld_->field("E_eta", iblk);
            auto &E_face_xi = fld_->field("E_face_xi", iblk);
            auto &E_face_zeta = fld_->field("E_face_zeta", iblk);
            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        E.vector[0] = 0.25 * (E_face_xi(i, j, k, 0) + E_face_xi(i, j, k - 1, 0) + E_face_zeta(i, j, k, 0) + E_face_zeta(i - 1, j, k, 0));
                        E.vector[1] = 0.25 * (E_face_xi(i, j, k, 1) + E_face_xi(i, j, k - 1, 1) + E_face_zeta(i, j, k, 1) + E_face_zeta(i - 1, j, k, 1));
                        E.vector[2] = 0.25 * (E_face_xi(i, j, k, 2) + E_face_xi(i, j, k - 1, 2) + E_face_zeta(i, j, k, 2) + E_face_zeta(i - 1, j, k, 2));

                        dr = {x(i, j + 1, k) - x(i, j, k),
                              y(i, j + 1, k) - y(i, j, k),
                              z(i, j + 1, k) - z(i, j, k)};
                        Eeta(i, j, k, 0) = E * dr;
                    }
        }

        {
            auto &Ezeta = fld_->field("E_zeta", iblk);
            auto &E_face_xi = fld_->field("E_face_xi", iblk);
            auto &E_face_eta = fld_->field("E_face_eta", iblk);
            sub = Ezeta.inner_lo();
            sup = Ezeta.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        E.vector[0] = 0.25 * (E_face_xi(i, j, k, 0) + E_face_xi(i, j - 1, k, 0) + E_face_eta(i, j, k, 0) + E_face_eta(i - 1, j, k, 0));
                        E.vector[1] = 0.25 * (E_face_xi(i, j, k, 1) + E_face_xi(i, j - 1, k, 1) + E_face_eta(i, j, k, 1) + E_face_eta(i - 1, j, k, 1));
                        E.vector[2] = 0.25 * (E_face_xi(i, j, k, 2) + E_face_xi(i, j - 1, k, 2) + E_face_eta(i, j, k, 2) + E_face_eta(i - 1, j, k, 2));

                        dr = {x(i, j, k + 1) - x(i, j, k),
                              y(i, j, k + 1) - y(i, j, k),
                              z(i, j, k + 1) - z(i, j, k)};
                        Ezeta(i, j, k, 0) = E * dr;
                    }
        }
    }

    // for (auto &p : topo_->physical_patches)
    // {
    //     if (p.bc_name != "Solid_Surface")
    //         continue;

    //     int iblk = p.this_block;
    //     int direction = p.direction;
    //     if (abs(direction) == 1)
    //     {
    //         auto &Eeta = fld_->field("E_eta", iblk);
    //         sub = p.this_box_node.lo;
    //         sup = p.this_box_node.hi;
    //         for (int i = sub.i; i < sup.i; i++)
    //             for (int j = sub.j; j < sup.j - 1; j++)
    //                 for (int k = sub.k; k < sup.k; k++)
    //                 {
    //                     Eeta(i, j, k, 0) = 0.0;
    //                 }
    //         auto &Ezeta = fld_->field("E_zeta", iblk);
    //         sub = p.this_box_node.lo;
    //         sup = p.this_box_node.hi;
    //         for (int i = sub.i; i < sup.i; i++)
    //             for (int j = sub.j; j < sup.j; j++)
    //                 for (int k = sub.k; k < sup.k - 1; k++)
    //                 {
    //                     Ezeta(i, j, k, 0) = 0.0;
    //                 }
    //     }
    //     else if (abs(direction) == 2)
    //     {
    //         auto &Exi = fld_->field("E_xi", iblk);
    //         sub = p.this_box_node.lo;
    //         sup = p.this_box_node.hi;
    //         for (int i = sub.i; i < sup.i - 1; i++)
    //             for (int j = sub.j; j < sup.j; j++)
    //                 for (int k = sub.k; k < sup.k; k++)
    //                 {
    //                     Exi(i, j, k, 0) = 0.0;
    //                 }
    //         auto &Ezeta = fld_->field("E_zeta", iblk);
    //         sub = p.this_box_node.lo;
    //         sup = p.this_box_node.hi;
    //         for (int i = sub.i; i < sup.i; i++)
    //             for (int j = sub.j; j < sup.j; j++)
    //                 for (int k = sub.k; k < sup.k - 1; k++)
    //                 {
    //                     Ezeta(i, j, k, 0) = 0.0;
    //                 }
    //     }
    //     else if (abs(direction) == 3)
    //     {
    //         auto &Exi = fld_->field("E_xi", iblk);
    //         sub = p.this_box_node.lo;
    //         sup = p.this_box_node.hi;
    //         for (int i = sub.i; i < sup.i - 1; i++)
    //             for (int j = sub.j; j < sup.j; j++)
    //                 for (int k = sub.k; k < sup.k; k++)
    //                 {
    //                     Exi(i, j, k, 0) = 0.0;
    //                 }
    //         auto &Eeta = fld_->field("E_eta", iblk);
    //         sub = p.this_box_node.lo;
    //         sup = p.this_box_node.hi;
    //         for (int i = sub.i; i < sup.i; i++)
    //             for (int j = sub.j; j < sup.j - 1; j++)
    //                 for (int k = sub.k; k < sup.k; k++)
    //                 {
    //                     Eeta(i, j, k, 0) = 0.0;
    //                 }
    //     }
    // }
    bound_.add_Edge_pole_boundary("E_xi");
    // bound_.add_Edge_pole_boundary("E_eta");
    // bound_.add_Edge_pole_boundary("E_zeta");

    // 计算RHS_xi
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &RHS = fld_->field("RHS_xi", iblk);
        auto &Eeta = fld_->field("E_eta", iblk);
        auto &Ezeta = fld_->field("E_zeta", iblk);
        sub = RHS.inner_lo();
        sup = RHS.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    RHS(i, j, k, 0) -= (Eeta(i, j, k, 0) - Eeta(i, j, k + 1, 0));
                    RHS(i, j, k, 0) -= (Ezeta(i, j + 1, k, 0) - Ezeta(i, j, k, 0));
                }
    }
    // 计算RHS_eta
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &RHS = fld_->field("RHS_eta", iblk);
        auto &Exi = fld_->field("E_xi", iblk);
        auto &Ezeta = fld_->field("E_zeta", iblk);
        sub = RHS.inner_lo();
        sup = RHS.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    RHS(i, j, k, 0) -= (Exi(i, j, k + 1, 0) - Exi(i, j, k, 0));
                    RHS(i, j, k, 0) -= (Ezeta(i, j, k, 0) - Ezeta(i + 1, j, k, 0));
                }
    }
    // 计算RHS_zeta
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &RHS = fld_->field("RHS_zeta", iblk);
        auto &Eeta = fld_->field("E_eta", iblk);
        auto &Exi = fld_->field("E_xi", iblk);
        sub = RHS.inner_lo();
        sup = RHS.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    RHS(i, j, k, 0) -= (Exi(i, j, k, 0) - Exi(i, j + 1, k, 0));
                    RHS(i, j, k, 0) -= (Eeta(i + 1, j, k, 0) - Eeta(i, j, k, 0));
                }
    }
}
