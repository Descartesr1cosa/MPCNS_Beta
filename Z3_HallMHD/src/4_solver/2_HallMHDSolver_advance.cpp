// Core
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"

// Z3_HallMHD
#include "HallMHD_Solver.h"

//=========================================================================
//-------------------------------------------------------------------------
//=========================================================================

void HallMHDSolver::Time_Advance()
{
    ZeroRHS();
    AssembleRHS_Fluid();     // 原 inv_fluid()
    SyncElectricFace();      // E_face_*: BC + halo
    AssembleRHS_Induction(); // 原 inv_induce()，后会在这里插 H
    ApplyTimeUpdate_Euler(); // U += dt*RHS, B += dt*RHS_B
    Update_Physic_Time();    // 记录物理时间
}
//=========================================================================
//-------------------------------------------------------------------------
//=========================================================================

void HallMHDSolver::ZeroRHS()
{
    //---------------------------------------------------------------
    // RHS is initialized as 0.0
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        // Cell
        {
            auto &RHS = fld_->field(fid_.fid_RHS_U, iblock);
            int ncomp = RHS.descriptor().ncomp;
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        for (int m = 0; m < ncomp; m++)
                        {
                            RHS(i, j, k, m) = 0.0;
                        }
        }

        auto zero_face_rhs = [&](int fid)
        {
            auto &R = fld_->field(fid, iblock);
            const auto &lo = R.get_lo();
            const auto &hi = R.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        R(i, j, k, 0) = 0.0;
        };

        zero_face_rhs(fid_.fid_RHS_Bface.xi);
        zero_face_rhs(fid_.fid_RHS_Bface.eta);
        zero_face_rhs(fid_.fid_RHS_Bface.zeta);
    }
    //---------------------------------------------------------------
}

// 对 E_face_xi/eta/zeta 做 BC + halo同步
void HallMHDSolver::SyncElectricFace()
{
    std::vector<std::string> Electric_Face = {"E_face_xi", "E_face_eta", "E_face_zeta"};
    bound_.add_Face_boundary(Electric_Face);
    halo_->data_trans(Electric_Face[0]);
    halo_->data_trans(Electric_Face[1]);
    halo_->data_trans(Electric_Face[2]);
}

// 显式Euler推进
void HallMHDSolver::ApplyTimeUpdate_Euler()
{
    //-----------------------------------------------------
    // RK1 Time Advance
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        // Cell
        {
            auto &RHS = fld_->field(fid_.fid_RHS_U, iblock);
            auto &U = fld_->field(fid_.fid_U, iblock);
            int ncomp = U.descriptor().ncomp;
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        for (int m = 0; m < ncomp; m++)
                        {
                            U(i, j, k, m) += dt * RHS(i, j, k, m);
                        }
        }
        // B_xi
        {
            auto &RHS = fld_->field(fid_.fid_RHS_Bface.xi, iblock);
            auto &U = fld_->field(fid_.fid_Bface.xi, iblock);
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        U(i, j, k, 0) += dt * RHS(i, j, k, 0);
        }
        // B_eta
        {
            auto &RHS = fld_->field(fid_.fid_RHS_Bface.eta, iblock);
            auto &U = fld_->field(fid_.fid_Bface.eta, iblock);
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        U(i, j, k, 0) += dt * RHS(i, j, k, 0);
        }
        // B_zeta
        {
            auto &RHS = fld_->field(fid_.fid_RHS_Bface.zeta, iblock);
            auto &U = fld_->field(fid_.fid_Bface.zeta, iblock);
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        U(i, j, k, 0) += dt * RHS(i, j, k, 0);
        }
    }
}

void HallMHDSolver::Update_Physic_Time()
{
    //-----------------------------------------------------
    // Update time
    control_.Physic_Time += dt;
    par_->AddParam("Physic_Time", control_.Physic_Time);
    //-----------------------------------------------------
}
