#include "HallMHD_Solver.h"
#include "2_operators/CTOperators.h"
#include "util/Vector_tool.h"

//=========================================================================
//-------------------------------------------------------------------------
//=========================================================================
#if HALL_MODE == 1
void HallMHDSolver::AddExplicitHallToEdgeEMF_()
{
    ComputeJ_AtEdges_Inner_();
    ApplyBC_EdgeJ_();
    ComputeHallE_AtEdges_EnergyPreserving_(); // 只填 Ehall_xi/eta/zeta（线积分量）
    AccumulateHallE_ToTotalEdgeEMF_();        // E_edge += E_hall(B, rho_frozen)
}
#else
void HallMHDSolver::AddExplicitHallToEdgeEMF_()
{
}
#endif
//=========================================================================
//-------------------------------------------------------------------------
//=========================================================================
#if HALL_MODE == 1
// E_* += Ehall_*
void HallMHDSolver::AccumulateHallE_ToTotalEdgeEMF_()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &E_xi = fld_->field(fid_.fid_Eedge.xi, iblk);
        auto &E_eta = fld_->field(fid_.fid_Eedge.eta, iblk);
        auto &E_zeta = fld_->field(fid_.fid_Eedge.zeta, iblk);

        auto &Eh_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
        auto &Eh_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
        auto &Eh_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);

        // EdgeXi
        {
            Int3 lo = E_xi.inner_lo();
            Int3 hi = E_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_xi(i, j, k, 0) += Eh_xi(i, j, k, 0);
        }

        // EdgeEt
        {
            Int3 lo = E_eta.inner_lo();
            Int3 hi = E_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_eta(i, j, k, 0) += Eh_eta(i, j, k, 0);
        }

        // EdgeZe
        {
            Int3 lo = E_zeta.inner_lo();
            Int3 hi = E_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E_zeta(i, j, k, 0) += Eh_zeta(i, j, k, 0);
        }
    }
}
#endif

#if HALL_MODE != 0
void HallMHDSolver::ComputeJ_AtEdges_Inner_()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Bxi = fld_->field(fid_.fid_Bface.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_Bface.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_Bface.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_Jedge.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_Jedge.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_Jedge.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       Jxi, Jeta, Jzeta,
                                       /*multiper=*/1.0);
    }
}

void HallMHDSolver::ApplyBC_EdgeJ_()
{
    std::string Electric_Current[3] = {"J_xi", "J_eta", "J_zeta"};
    halo_->data_trans(Electric_Current[0]);
    halo_->data_trans(Electric_Current[1]);
    halo_->data_trans(Electric_Current[2]);

    bound_.add_Edge_pole_boundary(Electric_Current[0]);
    bound_.add_Edge_copy_boundary(Electric_Current[0]);
    bound_.add_Edge_copy_boundary(Electric_Current[1]);
    bound_.add_Edge_copy_boundary(Electric_Current[2]);
}

void HallMHDSolver::ComputeHallE_AtEdges_EnergyPreserving_()
{
    // ------------------------------------------------------------
    // Hall coefficient:
    //   E_hall = alpha * (J x B)
    //   alpha typically = +1/(n e)  (sign can be absorbed here)
    // ------------------------------------------------------------
    const double hall_coeff = hall_coef;
    const double rho_floor = 1e-300; // 防止除零

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &U = fld_->field(fid_.fid_U, iblk);

        auto &Bxi = fld_->field(fid_.fid_Bface.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_Bface.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_Bface.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_Jedge.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_Jedge.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_Jedge.zeta, iblk);

        auto &Ehall_xi = fld_->field(fid_.fid_Ehall.xi, iblk);
        auto &Ehall_eta = fld_->field(fid_.fid_Ehall.eta, iblk);
        auto &Ehall_zeta = fld_->field(fid_.fid_Ehall.zeta, iblk);

        // edge cache: 9 comps (row-major 3x3)
        auto &pinvGT_xi = fld_->field(fid_.fid_pinvGT.xi, iblk);
        auto &pinvGT_eta = fld_->field(fid_.fid_pinvGT.eta, iblk);
        auto &pinvGT_zeta = fld_->field(fid_.fid_pinvGT.zeta, iblk);

        auto &pinvAT_xi = fld_->field(fid_.fid_pinvAT.xi, iblk);
        auto &pinvAT_eta = fld_->field(fid_.fid_pinvAT.eta, iblk);
        auto &pinvAT_zeta = fld_->field(fid_.fid_pinvAT.zeta, iblk);

        auto &x = grd_->grids(iblk).x;
        auto &y = grd_->grids(iblk).y;
        auto &z = grd_->grids(iblk).z;

        // small helpers
        auto matvec3 = [&](FieldBlock &M9, int i, int j, int k,
                           double c0, double c1, double c2) -> Vec3
        {
            Vec3 v;
            v.vec[0] = M9(i, j, k, 0) * c0 + M9(i, j, k, 1) * c1 + M9(i, j, k, 2) * c2;
            v.vec[1] = M9(i, j, k, 3) * c0 + M9(i, j, k, 4) * c1 + M9(i, j, k, 5) * c2;
            v.vec[2] = M9(i, j, k, 6) * c0 + M9(i, j, k, 7) * c1 + M9(i, j, k, 8) * c2;
            return v;
        };

        auto cache_is_valid = [&](FieldBlock &M9, int i, int j, int k) -> bool
        {
            double s = 0.0;
            for (int m = 0; m < 9; ++m)
                s += std::abs(M9(i, j, k, m));
            return s > 0.0;
        };

        // ============================================================
        // 1) EdgeXi : Ehall_xi(i,j,k) = (alpha * (J x B)) · dr_xi
        // ============================================================
        {
            Int3 lo = Ehall_xi.inner_lo();
            Int3 hi = Ehall_xi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_xi, i, j, k) || !cache_is_valid(pinvAT_xi, i, j, k))
                        {
                            Ehall_xi(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // rho at xi-edge: average 4 surrounding cells (j- and k- directions)
                        double rho = 0.25 * (U(i, j, k, 0) +
                                             U(i, j - 1, k, 0) +
                                             U(i, j, k - 1, 0) +
                                             U(i, j - 1, k - 1, 0));
                        double alpha = hall_coeff / (rho + rho_floor);

                        // Phi (2-form) co-located at xi-edge center
                        double Phi_xi = 0.0;
                        for (int di : {0, 1})
                            for (int dj : {0, -1})
                                for (int dk : {0, -1})
                                    Phi_xi += Bxi(i + di, j + dj, k + dk, 0);
                        Phi_xi *= 0.125;

                        double Phi_eta = 0.5 * (Beta(i, j, k, 0) + Beta(i, j, k - 1, 0));
                        double Phi_zeta = 0.5 * (Bzeta(i, j, k, 0) + Bzeta(i, j - 1, k, 0));

                        // j (1-form) co-located at same xi-edge center
                        double j_xi = Jxi(i, j, k, 0);
                        double j_eta = 0.25 * (Jeta(i, j, k, 0) + Jeta(i + 1, j, k, 0) + Jeta(i, j - 1, k, 0) + Jeta(i + 1, j - 1, k, 0));
                        double j_zeta = 0.25 * (Jzeta(i, j, k, 0) + Jzeta(i + 1, j, k, 0) + Jzeta(i, j, k - 1, 0) + Jzeta(i + 1, j, k - 1, 0));

                        // map to physical vectors
                        Vec3 Jvec = matvec3(pinvGT_xi, i, j, k, j_xi, j_eta, j_zeta);
                        Vec3 Bvec = matvec3(pinvAT_xi, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

                        Vec3 Evec = (Jvec ^ Bvec);
                        Evec *= alpha;

                        Vec3 dr;
                        dr.vec[0] = x(i + 1, j, k) - x(i, j, k);
                        dr.vec[1] = y(i + 1, j, k) - y(i, j, k);
                        dr.vec[2] = z(i + 1, j, k) - z(i, j, k);

                        Ehall_xi(i, j, k, 0) = Evec * dr; // line integral
                    }
        }

        // ============================================================
        // 2) EdgeEt
        // ============================================================
        {
            Int3 lo = Ehall_eta.inner_lo();
            Int3 hi = Ehall_eta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_eta, i, j, k) || !cache_is_valid(pinvAT_eta, i, j, k))
                        {
                            Ehall_eta(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // rho at eta-edge: average 4 surrounding cells (i- and k- directions)
                        double rho = 0.25 * (U(i, j, k, 0) +
                                             U(i - 1, j, k, 0) +
                                             U(i, j, k - 1, 0) +
                                             U(i - 1, j, k - 1, 0));
                        double alpha = hall_coeff / (rho + rho_floor);

                        // Phi co-located at eta-edge center
                        double Phi_eta = 0.0;
                        for (int di : {0, -1})
                            for (int dj : {0, 1})
                                for (int dk : {0, -1})
                                    Phi_eta += Beta(i + di, j + dj, k + dk, 0);
                        Phi_eta *= 0.125;

                        double Phi_xi = 0.5 * (Bxi(i, j, k, 0) + Bxi(i, j, k - 1, 0));
                        double Phi_zeta = 0.5 * (Bzeta(i, j, k, 0) + Bzeta(i - 1, j, k, 0));

                        // j co-located at eta-edge center
                        double j_eta = Jeta(i, j, k, 0);
                        double j_xi = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j + 1, k, 0) + Jxi(i - 1, j, k, 0) + Jxi(i - 1, j + 1, k, 0));
                        double j_zeta = 0.25 * (Jzeta(i, j, k, 0) + Jzeta(i, j + 1, k, 0) + Jzeta(i, j, k - 1, 0) + Jzeta(i, j + 1, k - 1, 0));

                        Vec3 Jvec = matvec3(pinvGT_eta, i, j, k, j_xi, j_eta, j_zeta);
                        Vec3 Bvec = matvec3(pinvAT_eta, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

                        Vec3 Evec = (Jvec ^ Bvec);
                        Evec *= alpha;

                        Vec3 dr;
                        dr.vec[0] = x(i, j + 1, k) - x(i, j, k);
                        dr.vec[1] = y(i, j + 1, k) - y(i, j, k);
                        dr.vec[2] = z(i, j + 1, k) - z(i, j, k);

                        Ehall_eta(i, j, k, 0) = Evec * dr;
                    }
        }

        // ============================================================
        // 3) EdgeZe
        // ============================================================
        {
            Int3 lo = Ehall_zeta.inner_lo();
            Int3 hi = Ehall_zeta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        if (!cache_is_valid(pinvGT_zeta, i, j, k) || !cache_is_valid(pinvAT_zeta, i, j, k))
                        {
                            Ehall_zeta(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // rho at zeta-edge: average 4 surrounding cells (i- and j- directions)
                        double rho = 0.25 * (U(i, j, k, 0) +
                                             U(i - 1, j, k, 0) +
                                             U(i, j - 1, k, 0) +
                                             U(i - 1, j - 1, k, 0));
                        double alpha = hall_coeff / (rho + rho_floor);

                        // Phi co-located at zeta-edge center
                        double Phi_zeta = 0.0;
                        for (int di : {0, -1})
                            for (int dj : {0, -1})
                                for (int dk : {0, 1})
                                    Phi_zeta += Bzeta(i + di, j + dj, k + dk, 0);
                        Phi_zeta *= 0.125;

                        double Phi_xi = 0.5 * (Bxi(i, j, k, 0) + Bxi(i, j - 1, k, 0));
                        double Phi_eta = 0.5 * (Beta(i, j, k, 0) + Beta(i - 1, j, k, 0));

                        // j co-located at zeta-edge center
                        double j_zeta = Jzeta(i, j, k, 0);
                        double j_xi = 0.25 * (Jxi(i, j, k, 0) + Jxi(i, j, k + 1, 0) + Jxi(i - 1, j, k, 0) + Jxi(i - 1, j, k + 1, 0));
                        double j_eta = 0.25 * (Jeta(i, j, k, 0) + Jeta(i, j, k + 1, 0) + Jeta(i, j - 1, k, 0) + Jeta(i, j - 1, k + 1, 0));

                        Vec3 Jvec = matvec3(pinvGT_zeta, i, j, k, j_xi, j_eta, j_zeta);
                        Vec3 Bvec = matvec3(pinvAT_zeta, i, j, k, Phi_xi, Phi_eta, Phi_zeta);

                        Vec3 Evec = (Jvec ^ Bvec);
                        Evec *= alpha;

                        Vec3 dr;
                        dr.vec[0] = x(i, j, k + 1) - x(i, j, k);
                        dr.vec[1] = y(i, j, k + 1) - y(i, j, k);
                        dr.vec[2] = z(i, j, k + 1) - z(i, j, k);

                        Ehall_zeta(i, j, k, 0) = Evec * dr;
                    }
        }
    }
}
#endif