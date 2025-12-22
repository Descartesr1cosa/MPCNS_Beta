#include "HallMHD_Solver.h"

void HallMHDSolver::add_Emag_to_Etotal()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &U = fld_->field(fid_.fid_U, ib);
        auto &PV = fld_->field(fid_.fid_PV, ib);
        auto &Bcell = fld_->field(fid_.fid_Bcell, ib);

        Int3 lo = U.get_lo();
        Int3 hi = U.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // double rho = U(i, j, k, 0);
                    // double u = U(i, j, k, 1) / rho;
                    // double v = U(i, j, k, 2) / rho;
                    // double w = U(i, j, k, 3) / rho;
                    // double kin = 0.5 * rho * (u * u + v * v + w * w);
                    // double E_inner = PV(i, j, k, 3);

                    double Bx = Bcell(i, j, k, 0);
                    double By = Bcell(i, j, k, 1);
                    double Bz = Bcell(i, j, k, 2);
                    double B2 = Bx * Bx + By * By + Bz * Bz;
                    double E_mag = 0.5 * B2 * inver_MA2; //  磁场能量 new

                    U(i, j, k, 4) += E_mag; //  MHD 总能量 new
                }
    }
}