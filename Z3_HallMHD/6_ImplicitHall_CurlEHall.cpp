#include "HallMHD_Solver.h"
#include "CTOperators.h"

// 这里只实现ComputeRHShallFromCurrentBface_： 求解curl E_hall，这将作为函数指针传入ImplicitHall_Sovler中调用
// 从而实现隐式迭代求解
void HallMHDSolver::ComputeRHShallFromCurrentBface_()
{
    // 1) 清零 RHShall_*（不要清 RHS_*）
    ZeroRHShall();

    // 2) 用当前 B_face 计算 J_edge
    ComputeJ_AtEdges_Inner_();

    // 3) 施加J的边界条件
    ApplyBC_EdgeJ_();

    // 3) 用 J_edge 和 B 计算 Ehall_edge
    ComputeHallE_AtEdges_EnergyPreserving_();  // 只填 Ehall_xi/eta/zeta（线积分量）
    bound_.add_Edge_pole_boundary("Ehall_xi"); //   edge 的物理边界/极点修补 + halo

    // 4) curl(Ehall_edge) -> RHShall_face
    AssembleFaceRHSHall_FromEdgeHallEMF_Curl_();
}

void HallMHDSolver::ZeroRHShall()
{
    //-----------------------------------------------------------------------------------------
    // RHS is initialized as 0.0
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        // B_xi
        {
            auto &RHS = fld_->field("RHShall_xi", iblock);
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        RHS(i, j, k, 0) = 0.0;
        }
        // B_eta
        {
            auto &RHS = fld_->field("RHShall_eta", iblock);
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        RHS(i, j, k, 0) = 0.0;
        }
        // B_zeta
        {
            auto &RHS = fld_->field("RHShall_zeta", iblock);
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        RHS(i, j, k, 0) = 0.0;
        }
    }
    //-----------------------------------------------------
}

void HallMHDSolver::AssembleFaceRHSHall_FromEdgeHallEMF_Curl_()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &RHS_xi = fld_->field("RHShall_xi", iblk);
        auto &RHS_eta = fld_->field("RHShall_eta", iblk);
        auto &RHS_zeta = fld_->field("RHShall_zeta", iblk);
        auto &Exi = fld_->field("Ehall_xi", iblk);
        auto &Eeta = fld_->field("Ehall_eta", iblk);
        auto &Ezeta = fld_->field("Ehall_zeta", iblk);
        CTOperators::CurlEdgeToFace(iblk, Exi, Eeta, Ezeta, RHS_xi, RHS_eta, RHS_zeta, /*multiper=*/1.0); // curl E
    }
}