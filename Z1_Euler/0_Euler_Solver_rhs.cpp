#include "Euler_Solver.h"

void EulerSolver::time_advance()
{
    //-----------------------------------------------------------------------------------------
    // Step 1
    //-----------------------------------------------------
    // RHS is initialized as 0.0
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
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
    //-----------------------------------------------------
    // Calc RHS Terms For Fluid Eqs
    inv_rhs();
    // f_vis(RHS, U_temp[0]);
    // f_src(RHS, U_temp[0]);
    //-----------------------------------------------------
    // RK1 Time Advance
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
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
    //-----------------------------------------------------
    // Update time
    control_.Physic_Time += dt;
    par_->AddParam("Physic_Time", control_.Physic_Time);
    //-----------------------------------------------------
}
