#include "HallMHD_Solver.h"

double HallMHDSolver::ComputeLocalMaxRadius_()
{
    double radius_max = 0.0;

    for (int iblock = 0; iblock < fld_->num_blocks(); ++iblock)
    {
        auto &U = fld_->field(fid_.fid_U, iblock);
        auto &PV = fld_->field(fid_.fid_PV, iblock);
        auto &Bcell = fld_->field(fid_.fid_Bcell, iblock);

        auto &GCL_xi = fld_->field(fid_.fid_metric.xi, iblock);
        auto &GCL_eta = fld_->field(fid_.fid_metric.eta, iblock);
        auto &GCL_zeta = fld_->field(fid_.fid_metric.zeta, iblock);
        auto &jac = fld_->field(fid_.fid_Jac, iblock);

        const Int3 &sub = U.inner_lo();
        const Int3 &sup = U.inner_hi();

        for (int i = sub.i; i < sup.i; ++i)
            for (int j = sub.j; j < sup.j; ++j)
                for (int k = sub.k; k < sup.k; ++k)
                {
                    const double J = jac(i, j, k, 0);

                    double xi[3], eta[3], zeta[3];
                    xi[0] = 0.5 * (GCL_xi(i + 1, j, k, 0) + GCL_xi(i, j, k, 0)) / J;
                    xi[1] = 0.5 * (GCL_xi(i + 1, j, k, 1) + GCL_xi(i, j, k, 1)) / J;
                    xi[2] = 0.5 * (GCL_xi(i + 1, j, k, 2) + GCL_xi(i, j, k, 2)) / J;

                    eta[0] = 0.5 * (GCL_eta(i, j + 1, k, 0) + GCL_eta(i, j, k, 0)) / J;
                    eta[1] = 0.5 * (GCL_eta(i, j + 1, k, 1) + GCL_eta(i, j, k, 1)) / J;
                    eta[2] = 0.5 * (GCL_eta(i, j + 1, k, 2) + GCL_eta(i, j, k, 2)) / J;

                    zeta[0] = 0.5 * (GCL_zeta(i, j, k + 1, 0) + GCL_zeta(i, j, k, 0)) / J;
                    zeta[1] = 0.5 * (GCL_zeta(i, j, k + 1, 1) + GCL_zeta(i, j, k, 1)) / J;
                    zeta[2] = 0.5 * (GCL_zeta(i, j, k + 1, 2) + GCL_zeta(i, j, k, 2)) / J;

                    const double rho = U(i, j, k, 0);

                    const double u = U(i, j, k, 1) / rho;
                    const double v = U(i, j, k, 2) / rho;
                    const double w = U(i, j, k, 3) / rho;

                    const double Bx = Bcell(i, j, k, 0);
                    const double By = Bcell(i, j, k, 1);
                    const double Bz = Bcell(i, j, k, 2);
                    const double B2 = Bx * Bx + By * By + Bz * Bz;

                    const double p = PV(i, j, k, 3);
                    // fast magnetosonic speed（近似）
                    const double c = std::sqrt(gamma_ * p / rho + (B2 / rho) * inver_MA2);

                    auto radius_dir = [&](const double a[3])
                    {
                        const double un = u * a[0] + v * a[1] + w * a[2];
                        const double an = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
                        return std::fabs(un) + c * an;
                    };

                    const double radius = radius_dir(xi) + radius_dir(eta) + radius_dir(zeta);
                    if (radius > radius_max)
                        radius_max = radius;
                }
    }

    return radius_max;
}
