#include "HallMHD_Solver.h"

// 输入：E_face_xi/eta/zeta（已做完物理边界 + halo）
// 输出：RHS_xi/eta/zeta（face 上的 dB/dt 项，CT curl）
void HallMHDSolver::AssembleRHS_Induction()
{
    // Calc RHS Terms For MHD induced Eqs
    AssembleEdgeEMF_FromFaceE_Ideal_();  // 1) face E -> edge (E·dr)
    ApplyBC_EdgeEMF_();                  // 2) edge 的物理边界/极点修补 + halo
    AssembleFaceRHS_FromEdgeEMF_Curl_(); // 3) curl(edge EMF) -> RHS_face
}

// Eelectric face  to  edge
void HallMHDSolver::AssembleEdgeEMF_FromFaceE_Ideal_()
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
}

// Eelectric  edge Boundary
void HallMHDSolver::ApplyBC_EdgeEMF_()
{
    bound_.add_Edge_pole_boundary("E_xi"); // pole边界处理
    // bound_.add_Edge_pole_boundary("E_eta");
    // bound_.add_Edge_pole_boundary("E_zeta");

    // halo_->data_trans("E_xi");
    // halo_->data_trans("E_eta");
    // halo_->data_trans("E_zeta");
}

// curl E
void HallMHDSolver::AssembleFaceRHS_FromEdgeEMF_Curl_()
{
    Int3 sub, sup;

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
