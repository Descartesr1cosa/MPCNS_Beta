#include "Euler_Solver.h"

void EulerSolver::inv_rhs()
{
    auto &Jac_ = fld_->field(fid_Jac);
    auto &Xi_ = fld_->field(fid_Xi);
    auto &Eta_ = fld_->field(fid_Eta);
    auto &Zeta_ = fld_->field(fid_Zeta);

    Int3 sub, sup;
    double metric[3];
    double Flux[5];

    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &U = fld_->field(fid_U, iblk);
        auto &PV = fld_->field(fid_PV, iblk);
        int ncomp = U.descriptor().ncomp;

        // 计算XI面的通量
        {
            auto &flux_xi = fld_->field("F_xi", iblk);
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
                        Reconstruction(metric, 0, PV, U, iblk, i, j, k, Flux);
                        for (int m = 0; m < ncomp; m++)
                            flux_xi(i, j, k, m) = Flux[m];
                    }
        }

        // 计算ETA面的通量
        {
            auto &flux_eta = fld_->field("F_eta", iblk);
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

                        Reconstruction(metric, 1, PV, U, iblk, i, j, k, Flux);
                        for (int m = 0; m < ncomp; m++)
                            flux_eta(i, j, k, m) = Flux[m];
                    }
        }

        // 计算ZETA面的通量
        {
            auto &flux_zeta = fld_->field("F_zeta", iblk);
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

                        Reconstruction(metric, 2, PV, U, iblk, i, j, k, Flux);
                        for (int m = 0; m < ncomp; m++)
                            flux_zeta(i, j, k, m) = Flux[m];
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

void EulerSolver::Reconstruction(double *metric, int32_t direction, FieldBlock &PV, FieldBlock &U, int iblock, int index_i, int index_j, int index_k, double *out_flux)
{
    auto RECON_minmod = [&](double *stencil, double *Left, double *Right)
    {
        double DR, DL, DM, DLeft, DRight;
        constexpr double beta = 1.5;
        // For Left
        DR = beta * (stencil[2] - stencil[1]);
        DL = beta * (stencil[1] - stencil[0]);
        DM = 0.5 * (stencil[2] - stencil[0]);
        if (DR * DL >= 0.0 && DL * DM >= 0.0)
            DLeft = (DR >= 0.0) ? fmin(fmin(DL, DR), DM) : -fmin(fmin(-DL, -DR), -DM);
        else
            DLeft = 0.0;
        // *Left = stencil[1] + 0.5 * DLeft;
        *Left = stencil[1];
        // For Right
        DR = beta * (stencil[3] - stencil[2]);
        DL = beta * (stencil[2] - stencil[1]);
        DM = 0.5 * (stencil[3] - stencil[1]);
        if (DR * DL >= 0.0 && DL * DM >= 0.0)
            DRight = (DR >= 0.0) ? fmin(fmin(DL, DR), DM) : -fmin(fmin(-DL, -DR), -DM);
        else
            DRight = 0.0;
        // *Right = stencil[2] + 0.5 * DRight;
        *Right = stencil[2];
    };

    auto Withdraw = [&](double *stencil, double *UUU, double *ppvv) ///* ,double *Badd_)*/
    {
        double rho, u, v, w, p, Bx, By, Bz, bx, by, bz;
        rho = stencil[0];
        u = stencil[1];
        v = stencil[2];
        w = stencil[3];
        p = stencil[4];
        // bx = stencil[5];
        // by = stencil[6];
        // bz = stencil[7];
        // Bx = stencil[8];
        // By = stencil[9];
        // Bz = stencil[10];
        // Badd_[0] = Bx;
        // Badd_[1] = By;
        // Badd_[2] = Bz;
        ppvv[0] = u;
        ppvv[1] = v;
        ppvv[2] = w;
        ppvv[3] = p;
        UUU[0] = rho;
        UUU[1] = rho * u;
        UUU[2] = rho * v;
        UUU[3] = rho * w;
        // Bx += bx;
        // By += by;
        // Bz += bz;
        UUU[4] = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w); // + 0.5 * inver_MA2 * (Bx * Bx + By * By + Bz * Bz);
        // UUU[5] = bx;
        // UUU[6] = by;
        // UUU[7] = bz;
    };

    auto calc_Jac_radius_GCL = [&](double &out, double *uu, double *pv, double *metric)
    {
        double rho, p, u, v, w;
        // double BB2, Bx, By, Bz;
        double K1, K2, K3;
        rho = uu[0];
        u = pv[0];
        v = pv[1];
        w = pv[2];
        p = pv[3];
        // Bx = uu[5] + B[0];
        // By = uu[6] + B[1];
        // Bz = uu[7] + B[2];

        // BB2 = Bx * Bx + By * By + Bz * Bz;

        K1 = metric[0];
        K2 = metric[1];
        K3 = metric[2];

        double uvw = K1 * u + K2 * v + K3 * w;
        double cc = sqrt((gamma_ * p / rho) * (K1 * K1 + K2 * K2 + K3 * K3));
        // double cc = sqrt((gamma_ * p / rho + BB2 / rho * inver_MA2) * (K1 * K1 + K2 * K2 + K3 * K3));
        // constexpr double C_hall_safe = 1.5;
        // double c_hall = C_hall_safe * ion_inertial_len * sqrt(BB2 / rho * inver_MA2 * (K1 * K1 + K2 * K2 + K3 * K3)); // ≈ Jac * v_A * d_i * |k|
        // out = fabs(uvw) + cc + c_hall;
        out = fabs(uvw) + cc;
        return;
    };

    auto calc_Jac_Flux_GCL = [&](double *flux, double *uu, double *pv, double *metric)
    {
        double k1, k2, k3; // GCL 这里为Jac *k1, Jac *k2, Jac *k3
        k1 = metric[0];
        k2 = metric[1];
        k3 = metric[2];

        double rho, p, u, v, w, uvw, rhoe;
        // double Bx, By, Bz, B_Jac_nabla, P_B;

        rho = uu[0];
        u = pv[0];
        v = pv[1];
        w = pv[2];
        p = pv[3];

        // Bx = B[0] + uu[5];
        // By = B[1] + uu[6];
        // Bz = B[2] + uu[7];
        // B_Jac_nabla = Bx * k1 + By * k2 + Bz * k3;
        uvw = k1 * u + k2 * v + k3 * w;

        // P_B = 0.5 * (Bx * Bx + By * By + Bz * Bz);
        rhoe = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

        flux[0] = rho * uvw;
        flux[1] = (rho * uvw * u + k1 * p); // Euler Flux - Maxwell Tensor (B\otimesB - B^2/2 I)\cdot J\nabla \xi\eta\zeta
        flux[2] = (rho * uvw * v + k2 * p);
        flux[3] = (rho * uvw * w + k3 * p);
        flux[4] = uvw * (rhoe + p); // Euler Flux

        // flux[0] = rho * uvw;
        // flux[1] = (rho * uvw * u + k1 * (p + inver_MA2 * P_B) - inver_MA2 * B_Jac_nabla * Bx); // Euler Flux - Maxwell Tensor (B\otimesB - B^2/2 I)\cdot J\nabla \xi\eta\zeta
        // flux[2] = (rho * uvw * v + k2 * (p + inver_MA2 * P_B) - inver_MA2 * B_Jac_nabla * By);
        // flux[3] = (rho * uvw * w + k3 * (p + inver_MA2 * P_B) - inver_MA2 * B_Jac_nabla * Bz);
        // flux[4] = uvw * (rhoe + p); // Euler Flux

        // u += H[0];
        // v += H[1];
        // w += H[2];
        // uvw = k1 * u + k2 * v + k3 * w; // With Hall Effect

        // flux[4] += inver_MA2 * (2.0 * P_B * uvw - (Bx * u + By * v + Bz * w) * B_Jac_nabla); // S(Poynting Vector) \cdot J\nabla \xi\eta\zeta
        // flux[5] = uvw * Bx - B_Jac_nabla * u;
        // flux[6] = uvw * By - B_Jac_nabla * v;
        // flux[7] = uvw * Bz - B_Jac_nabla * w;
    };

    // double metric[3]; // 度量系数。依次为jac*kx, jac*ky, jac*kz
    int i, j, k;
    i = index_i;
    j = index_j;
    k = index_k;

    double UL[8], UR[8], ppvvL[4], ppvvR[4], BaddL[3], BaddR[3], radius[2], Hall[3];
    double recon[11][4], reconL[11], reconR[11];

    if (direction == 0)
    {
        // 对模板点重构，获得半点左右状态
        i -= 1;
        for (int ii = 0; ii < 4; ii++) // All stencils
        {
            recon[0][ii] = U(i, j, k, 0);  // rho
            recon[1][ii] = PV(i, j, k, 0); // u
            recon[2][ii] = PV(i, j, k, 1); // v
            recon[3][ii] = PV(i, j, k, 2); // w
            recon[4][ii] = PV(i, j, k, 3); // p
            // recon[5][ii] = U(i, j, k, 5);                 // bx
            // recon[6][ii] = U(i, j, k, 6);                 // by
            // recon[7][ii] = U(i, j, k, 7);                 // bz
            // recon[8][ii] = (*B_add)[iblock](i, j, k, 0);  // Baddx
            // recon[9][ii] = (*B_add)[iblock](i, j, k, 1);  // Baddy
            // recon[10][ii] = (*B_add)[iblock](i, j, k, 2); // Baddz

            i++;
        }
        for (int index = 0; index < 11; index++)
            RECON_minmod(recon[index], &(reconL[index]), &(reconR[index]));
        Withdraw(reconL, UL, ppvvL); //, BaddL);
        Withdraw(reconR, UR, ppvvR); //, BaddR);
        calc_Jac_radius_GCL(radius[0], UL, ppvvL, metric);
        calc_Jac_radius_GCL(radius[1], UR, ppvvR, metric);
    }
    else if (direction == 1)
    {
        // 对模板点重构，获得半点左右状态
        j -= 1;
        for (int ii = 0; ii < 4; ii++) // All stencils
        {
            recon[0][ii] = U(i, j, k, 0);  // rho
            recon[1][ii] = PV(i, j, k, 0); // u
            recon[2][ii] = PV(i, j, k, 1); // v
            recon[3][ii] = PV(i, j, k, 2); // w
            recon[4][ii] = PV(i, j, k, 3); // p
            // recon[5][ii] = U(i, j, k, 5);                 // bx
            // recon[6][ii] = U(i, j, k, 6);                 // by
            // recon[7][ii] = U(i, j, k, 7);                 // bz
            // recon[8][ii] = (*B_add)[iblock](i, j, k, 0);  // Baddx
            // recon[9][ii] = (*B_add)[iblock](i, j, k, 1);  // Baddy
            // recon[10][ii] = (*B_add)[iblock](i, j, k, 2); // Baddz

            j++;
        }
        for (int index = 0; index < 5; index++)
            RECON_minmod(recon[index], &(reconL[index]), &(reconR[index]));
        Withdraw(reconL, UL, ppvvL); //, BaddL);
        Withdraw(reconR, UR, ppvvR); //, BaddR);
        calc_Jac_radius_GCL(radius[0], UL, ppvvL, metric);
        calc_Jac_radius_GCL(radius[1], UR, ppvvR, metric);
    }
    else
    {

        // 对模板点重构，获得半点左右状态
        k -= 1;
        for (int ii = 0; ii < 4; ii++) // All stencils
        {
            recon[0][ii] = U(i, j, k, 0);  // rho
            recon[1][ii] = PV(i, j, k, 0); // u
            recon[2][ii] = PV(i, j, k, 1); // v
            recon[3][ii] = PV(i, j, k, 2); // w
            recon[4][ii] = PV(i, j, k, 3); // p
            // recon[5][ii] = U(i, j, k, 5);                 // bx
            // recon[6][ii] = U(i, j, k, 6);                 // by
            // recon[7][ii] = U(i, j, k, 7);                 // bz
            // recon[8][ii] = (*B_add)[iblock](i, j, k, 0);  // Baddx
            // recon[9][ii] = (*B_add)[iblock](i, j, k, 1);  // Baddy
            // recon[10][ii] = (*B_add)[iblock](i, j, k, 2); // Baddz

            k++;
        }
        for (int index = 0; index < 11; index++)
            RECON_minmod(recon[index], &(reconL[index]), &(reconR[index]));
        Withdraw(reconL, UL, ppvvL); //, BaddL);
        Withdraw(reconR, UR, ppvvR); //, BaddR);
        calc_Jac_radius_GCL(radius[0], UL, ppvvL, metric);
        calc_Jac_radius_GCL(radius[1], UR, ppvvR, metric);
    }
    //=============================================================================================
    // Flux Splitting
    double radius_max = 0.0;
    double LeftF[8], RightF[8];
    for (int ii = 0; ii < 2; ii++)
        radius_max = (radius_max < radius[ii]) ? radius[ii] : radius_max;
    calc_Jac_Flux_GCL(LeftF, UL, ppvvL, metric);
    calc_Jac_Flux_GCL(RightF, UR, ppvvR, metric);
    for (int jj = 0; jj < 5; jj++)
        out_flux[jj] = 0.5 * (RightF[jj] + LeftF[jj]) - 0.5 * radius_max * (UR[jj] - UL[jj]);
}
