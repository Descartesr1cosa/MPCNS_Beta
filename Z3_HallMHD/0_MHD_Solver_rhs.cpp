#include "MHD_Solver.h"

void MHDSolver::time_advance()
{
    //-----------------------------------------------------------------------------------------
    // Step 1
    //-----------------------------------------------------
    // RHS is initialized as 0.0
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        // Cell
        {
            auto &RHS = fld_->field("RHS", iblock);
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
        // B_xi
        {
            auto &RHS = fld_->field("RHS_xi", iblock);
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        RHS(i, j, k, 0) = 0.0;
        }
        // B_eta
        {
            auto &RHS = fld_->field("RHS_eta", iblock);
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        RHS(i, j, k, 0) = 0.0;
        }
        // B_zeta
        {
            auto &RHS = fld_->field("RHS_zeta", iblock);
            const Int3 &sub = RHS.get_lo();
            const Int3 &sup = RHS.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        RHS(i, j, k, 0) = 0.0;
        }
    }
    //-----------------------------------------------------
    // Calc RHS Terms For Fluid Eqs
    inv_fluid();

    std::vector<std::string> Electric_Face = {"E_face_xi", "E_face_eta", "E_face_zeta"};
    bound_.add_Face_boundary(Electric_Face);
    halo_->data_trans(Electric_Face[0]);
    halo_->data_trans(Electric_Face[1]);
    halo_->data_trans(Electric_Face[2]);

    // Calc RHS Terms For MHD induced Eqs
    inv_induce();
    //-----------------------------------------------------
    // RK1 Time Advance
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        // Cell
        {
            auto &RHS = fld_->field("RHS", iblock);
            auto &U = fld_->field("U_", iblock);
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
            auto &RHS = fld_->field("RHS_xi", iblock);
            auto &U = fld_->field("B_xi", iblock);
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        U(i, j, k, 0) += dt * RHS(i, j, k, 0);
        }
        // B_eta
        {
            auto &RHS = fld_->field("RHS_eta", iblock);
            auto &U = fld_->field("B_eta", iblock);
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        U(i, j, k, 0) += dt * RHS(i, j, k, 0);
        }
        // B_zeta
        {
            auto &RHS = fld_->field("RHS_zeta", iblock);
            auto &U = fld_->field("B_zeta", iblock);
            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        U(i, j, k, 0) += dt * RHS(i, j, k, 0);
        }
    }
    //-----------------------------------------------------
    // Update time
    control_.Physic_Time += dt;
    par_->AddParam("Physic_Time", control_.Physic_Time);
    //-----------------------------------------------------
}
