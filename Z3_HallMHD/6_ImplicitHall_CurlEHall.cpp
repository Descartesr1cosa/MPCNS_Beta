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
    ComputeHallE_AtEdges_EnergyPreserving_(); // 只填 Ehall_xi/eta/zeta（线积分量）
    ApplyBC_EdgeHallEMF_();                   //   edge 的物理边界/极点修补 + halo

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

void HallMHDSolver::Modify_TotalEnergy_AfterHall()
{
    // 之前已经计算了U_ PV_,Hall只修改了B
    // 因此需要更新B_cell 然后把U_的能量修改即可，无需（也不能）更新PV
    SyncPrimaryFaceB();  // B_face: BC + halo
    ComputeBcellInner(); // B_cell: inner compute
    SyncDerivedBcell();  // B_cell: derived BC + halo (+corner)
    UpdateTotalEnergy(); // using B_cell and PV --> U
    SyncPrimaryCellU();  // U: BC + halo (+corner)
    calc_divB();         // divB
}

void HallMHDSolver::UpdateTotalEnergy()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &U = fld_->field(fid_U, ib);
        auto &PV = fld_->field(fid_PV, ib);
        auto &Bcell = fld_->field(fid_Bcell, ib);

        Int3 lo = U.get_lo();
        Int3 hi = U.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    double rho = U(i, j, k, 0);
                    double u = PV(i, j, k, 0);
                    double v = PV(i, j, k, 1);
                    double w = PV(i, j, k, 2);
                    double p = PV(i, j, k, 3);

                    double kin = 0.5 * rho * (u * u + v * v + w * w);

                    double Bx = Bcell(i, j, k, 0);
                    double By = Bcell(i, j, k, 1);
                    double Bz = Bcell(i, j, k, 2);
                    double B2 = Bx * Bx + By * By + Bz * Bz;
                    double E_mag = 0.5 * B2 * inver_MA2; //  磁场能量 new

                    U(i, j, k, 4) = kin + p / (gamma_ - 1.0) + E_mag; //  MHD 总能量 new
                }
    }
}

void HallMHDSolver::ApplyBC_EdgeHallEMF_()
{
    // 后续CT只会用到inner的电场，不会使用nghost区域，因此只需对Pole处理即可
    bound_.add_Edge_pole_boundary("Ehall_xi"); // pole边界处理
    bound_.add_Edge_pole_boundary("Ehall_eta");
    // bound_.add_Edge_pole_boundary("Ehall_zeta");

    // halo_->data_trans("Ehall_xi");
    // halo_->data_trans("Ehall_eta");
    // halo_->data_trans("Ehall_zeta");
}