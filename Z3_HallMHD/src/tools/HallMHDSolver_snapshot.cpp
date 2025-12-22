#include "HallMHD_Solver.h"

void HallMHDSolver::copy_field()
{
    for (int ib = 0; ib < fld_->num_blocks(); ib++)
    {
        {
            auto &U = fld_->field(fid_.fid_U, ib);
            auto &old_U = fld_->field(fid_.fid_old_U, ib);
            int ncomp = U.descriptor().ncomp;
            const Int3 &sub = old_U.get_lo();
            const Int3 &sup = old_U.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        for (int m = 0; m < ncomp; m++)
                            old_U(i, j, k, m) = U(i, j, k, m);
        }

        {
            auto &U = fld_->field(fid_.fid_Bface.xi, ib);
            auto &old_U = fld_->field(fid_.fid_old_Bface.xi, ib);
            int ncomp = U.descriptor().ncomp;
            const Int3 &sub = old_U.get_lo();
            const Int3 &sup = old_U.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        for (int m = 0; m < ncomp; m++)
                            old_U(i, j, k, m) = U(i, j, k, m);
        }
        {
            auto &U = fld_->field(fid_.fid_Bface.eta, ib);
            auto &old_U = fld_->field(fid_.fid_old_Bface.eta, ib);
            int ncomp = U.descriptor().ncomp;
            const Int3 &sub = old_U.get_lo();
            const Int3 &sup = old_U.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        for (int m = 0; m < ncomp; m++)
                            old_U(i, j, k, m) = U(i, j, k, m);
        }
        {
            auto &U = fld_->field(fid_.fid_Bface.zeta, ib);
            auto &old_U = fld_->field(fid_.fid_old_Bface.zeta, ib);
            int ncomp = U.descriptor().ncomp;
            const Int3 &sub = old_U.get_lo();
            const Int3 &sup = old_U.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                        for (int m = 0; m < ncomp; m++)
                            old_U(i, j, k, m) = U(i, j, k, m);
        }
    }
}