#include "HallMHD_Solver.h"

void HallMHDSolver::Compute_Timestep()
{
    if (control_.Time_step <= 0.0)
    {
        if (control_.CFL * (control_.CFL - 1.0) >= 0.0)
        {
            std::cout << "Input timestep control parameter goes wrong for Explicit Euler\n";
            std::cout << "Input timestep is\t" << control_.Time_step << "\tCFL=\t" << control_.CFL << std::endl;
            exit(-1);
        }
        // 遍历本进程的块，获得最小的稳定时间步长
        double radius, radius_max = 0.0;
        double xi[3], eta[3], zeta[3];
        double u, v, w, rho, p, c, B2, Bx, By, Bz;
        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            auto &U = fld_->field("U_", iblock);
            auto &PV = fld_->field("PV_", iblock);
            auto &Bcell = fld_->field("B_cell", iblock);
            auto &GCL_xi = fld_->field("JDxi", iblock);
            auto &GCL_eta = fld_->field("JDet", iblock);
            auto &GCL_zeta = fld_->field("JDze", iblock);
            auto &jac = fld_->field("Jac", iblock);

            const Int3 &sub = U.inner_lo();
            const Int3 &sup = U.inner_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        xi[0] = 0.5 * (GCL_xi(i + 1, j, k, 0) + GCL_xi(i, j, k, 0)) / jac(i, j, k, 0);
                        xi[1] = 0.5 * (GCL_xi(i + 1, j, k, 1) + GCL_xi(i, j, k, 1)) / jac(i, j, k, 0);
                        xi[2] = 0.5 * (GCL_xi(i + 1, j, k, 2) + GCL_xi(i, j, k, 2)) / jac(i, j, k, 0);

                        eta[0] = 0.5 * (GCL_eta(i, j + 1, k, 0) + GCL_eta(i, j, k, 0)) / jac(i, j, k, 0);
                        eta[1] = 0.5 * (GCL_eta(i, j + 1, k, 1) + GCL_eta(i, j, k, 1)) / jac(i, j, k, 0);
                        eta[2] = 0.5 * (GCL_eta(i, j + 1, k, 2) + GCL_eta(i, j, k, 2)) / jac(i, j, k, 0);

                        zeta[0] = 0.5 * (GCL_zeta(i, j, k + 1, 0) + GCL_zeta(i, j, k, 0)) / jac(i, j, k, 0);
                        zeta[1] = 0.5 * (GCL_zeta(i, j, k + 1, 1) + GCL_zeta(i, j, k, 1)) / jac(i, j, k, 0);
                        zeta[2] = 0.5 * (GCL_zeta(i, j, k + 1, 2) + GCL_zeta(i, j, k, 2)) / jac(i, j, k, 0);

                        rho = U(i, j, k, 0);
                        u = U(i, j, k, 1) / rho;
                        v = U(i, j, k, 2) / rho;
                        w = U(i, j, k, 3) / rho;

                        Bx = Bcell(i, j, k, 0);
                        By = Bcell(i, j, k, 1);
                        Bz = Bcell(i, j, k, 2);
                        B2 = Bx * Bx + By * By + Bz * Bz;

                        p = PV(i, j, k, 3);
                        c = sqrt(gamma_ * p / rho + B2 / rho * inver_MA2); // fast magnetosonic speed

                        // xi
                        radius = fabs(u * xi[0] + v * xi[1] + w * xi[2]) + c * sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
                        // eta
                        radius += fabs(u * eta[0] + v * eta[1] + w * eta[2]) + c * sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
                        // zeta
                        radius += fabs(u * zeta[0] + v * zeta[1] + w * zeta[2]) + c * sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);

                        // find the maximum radius
                        if (radius > radius_max)
                            radius_max = radius;
                    }
        }
        dt = control_.CFL / radius_max;
        // reduce for ALL processes
        double dtemp = dt;
        PARALLEL::mpi_min(&dtemp, &dt, 1);
    }
    else
    {
        dt = control_.Time_step;
    }
    control_.Physic_Time_Step = dt;
    fld_->par->AddParam("Physic_Time_Step", dt);
}