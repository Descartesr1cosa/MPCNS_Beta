#include "HallMHD_Solver.h"

// 1. 生成 F_xi/F_eta/F_zeta（用于流体方程散度）
// 2. 顺带产出 E_face_xi/eta/zeta（理想 MHD 的 face 电场）
// 3. 用 F_* 组装流体守恒变量 cell 的 RHS (-div F / Jac)
void HallMHDSolver::AssembleRHS_Fluid()
{
    auto &Xi_ = fld_->field(fid_Xi);
    auto &Eta_ = fld_->field(fid_Eta);
    auto &Zeta_ = fld_->field(fid_Zeta);

    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &U = fld_->field(fid_U, iblk);
        auto &PV = fld_->field(fid_PV, iblk);
        auto &Bcell = fld_->field(fid_Bcell, iblk);

        // 1) 三方向：通量 + ideal face EMF
        for (int dir = 0; dir < 3; ++dir)
        {
            auto &flux = fld_->field(fid_F_[dir], iblk);
            auto &E_face = fld_->field(fid_Eface_[dir], iblk);
            auto &B_face = fld_->field(fid_Bface_[dir], iblk);
            auto &metric = fld_->field(fid_metric_[dir], iblk); // Xi_[iblk]/Eta_[iblk]/Zeta_[iblk]

            AssembleOneDirectionFluxAndEMF_(iblk, dir, flux, E_face, B_face, metric, PV, U, Bcell);
        }
    }

    // Calc RHS Terms For Fluid Eqs
    AssembleCellRHSFromFlux_();
}

void HallMHDSolver::AssembleOneDirectionFluxAndEMF_(
    int iblk,
    int dir,                 // 0 xi, 1 eta, 2 zeta
    FieldBlock &flux,        // F_xi / F_eta / F_zeta (ncomp=5)
    FieldBlock &E_face,      // E_face_xi/eta/zeta   (ncomp=3)
    FieldBlock &B_face,      // B_xi/eta/zeta        (ncomp=1)
    FieldBlock &metricField, // Xi_/Eta_/Zeta_       (ncomp=3)
    FieldBlock &PV,
    FieldBlock &U,
    FieldBlock &Bcell)
{
    Int3 sub = flux.inner_lo();
    Int3 sup = flux.inner_hi();

    double metric[3];
    double flux8[8];                        // 注意：Reconstruction 输出 8 个，其中[5..7] 已被旋转成 E_face :contentReference[oaicite:4]{index=4}
    const int ncomp = U.descriptor().ncomp; // 5

    for (int i = sub.i; i < sup.i; ++i)
        for (int j = sub.j; j < sup.j; ++j)
            for (int k = sub.k; k < sup.k; ++k)
            {
                metric[0] = metricField(i, j, k, 0);
                metric[1] = metricField(i, j, k, 1);
                metric[2] = metricField(i, j, k, 2);

                Reconstruction(metric, dir, PV, U, Bcell,
                               B_face(i, j, k, 0), iblk, i, j, k, flux8);

                // 1) 流体通量（只存 0..ncomp-1）
                for (int m = 0; m < ncomp; ++m)
                    flux(i, j, k, m) = flux8[m];

                // 2) 理想 MHD 的 face 电场（Reconstruction 已旋转好）
                E_face(i, j, k, 0) = flux8[5];
                E_face(i, j, k, 1) = flux8[6];
                E_face(i, j, k, 2) = flux8[7];
            }
}

void HallMHDSolver::AssembleCellRHSFromFlux_()
{
    Int3 sub, sup;

    // 计算RHS
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &RHS = fld_->field("RHS", iblk);
        int ncomp = RHS.descriptor().ncomp;

        auto &Jac = fld_->field(fid_Jac, iblk);

        auto &flux_xi = fld_->field(fid_F_[0], iblk);
        auto &flux_eta = fld_->field(fid_F_[1], iblk);
        auto &flux_zeta = fld_->field(fid_F_[2], iblk);

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
