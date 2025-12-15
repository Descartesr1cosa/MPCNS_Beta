#pragma once
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"
#include "Control.h"
#include "Output.h"
#include "Boundary.h"
#include "Initial.h"

#include <array>

class MHDSolver
{
public:
    MHDSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par, std::vector<std::string> Solver_Name)
    {
        grd_ = grd;
        fld_ = fld;
        halo_ = halo;
        par_ = par;
        topo_ = topo;

        Solver_Name_ = Solver_Name;

        initial_.Initialization(fld_);

        fid_U = fld_->field_id("U_");
        fid_PV = fld_->field_id("PV_");
        fid_Jac = fld_->field_id("Jac");
        fid_Xi = fld_->field_id("JDxi");
        fid_Eta = fld_->field_id("JDet");
        fid_Zeta = fld_->field_id("JDze");
        fid_Bcell = fld_->field_id("B_cell");
        fid_Bxi = fld_->field_id("B_xi");
        fid_Beta = fld_->field_id("B_eta");
        fid_Bzeta = fld_->field_id("B_zeta");

        gamma_ = par_->GetDou_List("constant").data["gamma"];

        B_add_x = par->GetDou("B_add_x");
        B_add_y = par->GetDou("B_add_y");
        B_add_z = par->GetDou("B_add_z");

        inver_MA2 = par->GetDou("inver_MA2");

        control_.SetUp(par_, 8);
        output_.SetUp(par_, fld_);
        bound_.SetUp(grd_, fld_, topo_, par_);

        fld_->register_field(FieldDescriptor{"old_U_", StaggerLocation::Cell, 5, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"divB", StaggerLocation::Cell, 1, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"old_B_xi", StaggerLocation::FaceXi, 1, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"old_B_eta", StaggerLocation::FaceEt, 1, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"old_B_zeta", StaggerLocation::FaceZe, 1, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"RHS", StaggerLocation::Cell, 5, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"RHS_xi", StaggerLocation::FaceXi, 1, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"RHS_eta", StaggerLocation::FaceEt, 1, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"RHS_zeta", StaggerLocation::FaceZe, 1, par->GetInt("ngg")});
    }

    void Advance()
    {
        // 0. Preparation
        PrepareStep();

        while (true)
        {

            // 1. 计算时间步长
            compute_timestep();

            // 2. 时间推进（e.g. 最简单的 forward Euler）:计算 RHS: dU/dt = -div F(U)
            time_advance();

            // 3. 计算残差
            calc_Residual();

            // 4. Preparation
            PrepareStep();

            // 5. 更新控制
            control_.Update();
            if (control_.if_outfile)
                output_.output_field();
            // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
            // output_.output_plt_field(); // output_.output_plt_cell_field(output_.var_defaut_plt_name); // output_field();
            if (control_.if_stop)
            {
                // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
                output_.output_field();
                break;
            }
        }
    }

private:
    Grid *grd_;
    TOPO::Topology *topo_;
    Field *fld_;
    Halo *halo_;
    Param *par_;

    MHD_Control control_;
    MHD_Output output_;
    MHD_Boundary bound_;
    MHD_Initial initial_;

    int fid_Jac, fid_Xi, fid_Eta, fid_Zeta;

    int fid_U;
    int fid_PV;
    int fid_Bcell, fid_Bxi, fid_Beta, fid_Bzeta;
    std::vector<std::string> Solver_Name_;
    double gamma_, B_add_x, B_add_y, B_add_z, inver_MA2;

    double dt;

    void compute_timestep();
    void time_advance();
    void calc_Residual();

    void inv_fluid();
    void inv_induce();

    void Reconstruction(double *metric, int32_t direction, FieldBlock &PV, FieldBlock &U, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k, double *out_flux);

    void PrepareStep()
    {
        // A) Face 主量的物理边界 + halo
        std::vector<std::string> face_name = {"B_xi", "B_eta", "B_zeta"};
        bound_.add_Face_boundary(face_name);
        halo_->data_trans(face_name[0]);
        halo_->data_trans(face_name[1]);
        halo_->data_trans(face_name[2]);

        // B) 只计算 inner （计算域）网格的派生量 B_cell
        calc_Bcell();

        // C) B_cell 的“派生边界策略”先用最简单 copy + halo (暂时不用角区通信)
        std::string Bcell_string = "B_cell";
        bound_.add_derived_Cell_boundary(Bcell_string);
        halo_->data_trans(Bcell_string);
        halo_->data_trans_2DCorner(Bcell_string);
        halo_->data_trans_3DCorner(Bcell_string);

        // D) U_ 的物理边界（此时可安全用 B_cell） + halo (暂时不用角区通信)
        std::string U_string = "U_";
        bound_.add_Cell_boundary(U_string);
        halo_->data_trans(U_string);
        halo_->data_trans_2DCorner(U_string);
        halo_->data_trans_3DCorner(U_string);

        // E) 原始量 + divB
        calc_PV();
        calc_divB();

        // F) 复制物理场，用于残差计算
        copy_field();
    }

    struct Double3
    {
        double vector[3];
        Double3() : vector{0.0, 0.0, 0.0} {}
        Double3(double x, double y, double z) : vector{x, y, z} {}

        Double3 &operator+=(const Double3 &add)
        {
            vector[0] += add.vector[0];
            vector[1] += add.vector[1];
            vector[2] += add.vector[2];
            return *this;
        }
        Double3 &operator-=(const Double3 &add)
        {
            vector[0] -= add.vector[0];
            vector[1] -= add.vector[1];
            vector[2] -= add.vector[2];
            return *this;
        }
        Double3 &operator*=(double add)
        {
            vector[0] *= add;
            vector[1] *= add;
            vector[2] *= add;
            return *this;
        }
        // 点积
        friend double operator*(const Double3 &a, const Double3 &b)
        {
            return a.vector[0] * b.vector[0] + a.vector[1] * b.vector[1] + a.vector[2] * b.vector[2];
        }

        // 叉积，用 ^
        friend Double3 operator^(const Double3 &a, const Double3 &b)
        {
            return Double3(
                a.vector[1] * b.vector[2] - a.vector[2] * b.vector[1],
                a.vector[2] * b.vector[0] - a.vector[0] * b.vector[2],
                a.vector[0] * b.vector[1] - a.vector[1] * b.vector[0]);
        }
        void set(double a, double b, double c)
        {
            vector[0] = a;
            vector[1] = b;
            vector[2] = c;
        }
    };

    void calc_PV()
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
                        double u = U(i, j, k, 1) / rho;
                        double v = U(i, j, k, 2) / rho;
                        double w = U(i, j, k, 3) / rho;

                        PV(i, j, k, 0) = u;
                        PV(i, j, k, 1) = v;
                        PV(i, j, k, 2) = w;

                        double kin = 0.5 * rho * (u * u + v * v + w * w);

                        double Bx = Bcell(i, j, k, 0);
                        double By = Bcell(i, j, k, 1);
                        double Bz = Bcell(i, j, k, 2);
                        double B2 = Bx * Bx + By * By + Bz * Bz;

                        double E = U(i, j, k, 4);            //  MHD 总能量
                        double E_mag = 0.5 * B2 * inver_MA2; //  磁场能量

                        double p = (gamma_ - 1.0) * (E - kin - 0.5 * B2 * inver_MA2);
                        const double p_floor = 1e-8;

                        if (p < p_floor)
                        {
                            // 希望：保持 E、B 不变，把一部分动能变成内能，使 p >= p_floor

                            // 目标内能对应的 e_floor
                            double e_floor = p_floor / (gamma_ - 1.0);

                            // 对应的目标动能：E = K' + E_mag + e_floor
                            double K_target = E - E_mag - e_floor;

                            if (K_target < 0.0)
                            {
                                // 能量不够支撑这么大的内能，只能把所有非磁能都变成内能
                                K_target = 0.0;
                                e_floor = E - E_mag;
                                p = (gamma_ - 1.0) * e_floor;
                            }
                            else
                            {
                                // 正常情况：达到指定的 floor
                                p = p_floor;
                            }

                            if (kin > 0.0)
                            {
                                // 按 K' = beta^2 * K 缩放速度，从动能里挪出一部分作为内能
                                double beta = sqrt(K_target / kin);

                                u *= beta;
                                v *= beta;
                                w *= beta;

                                PV(i, j, k, 0) = u;
                                PV(i, j, k, 1) = v;
                                PV(i, j, k, 2) = w;

                                U(i, j, k, 1) = rho * u;
                                U(i, j, k, 2) = rho * v;
                                U(i, j, k, 3) = rho * w;
                            }

                            // 注意：E 不变，保持总能量守恒
                            // U(i, j, k, 4) = E;  // 不需要改
                        }
                        PV(i, j, k, 3) = p;
                        // PV(i, j, k, 3) = (E - kin - 0.5 * B2 * inver_MA2) * (gamma_ - 1.0);
                    }
        }
    }

    void calc_Bcell()
    {
        const int nblock = fld_->num_blocks();

        const double eps = 1e-300;
        const double delta = 1e-15; // same spirit as your reference

        auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        };
        auto norm = [&](const std::array<double, 3> &a)
        {
            return std::sqrt(dot(a, a));
        };

        for (int ib = 0; ib < nblock; ++ib)
        {
            auto &Bcell = fld_->field(fid_Bcell, ib);
            auto &U = fld_->field(fid_U, ib);

            auto &Bxi = fld_->field(fid_Bxi, ib);
            auto &Beta = fld_->field(fid_Beta, ib);
            auto &Bzeta = fld_->field(fid_Bzeta, ib);

            auto &Jac = fld_->field(fid_Jac, ib);
            auto &A_xi = fld_->field(fid_Xi, ib);   // JDxi
            auto &A_eta = fld_->field(fid_Eta, ib); // JDet
            auto &A_ze = fld_->field(fid_Zeta, ib); // JDze

            auto &x = grd_->grids(ib).x;
            auto &y = grd_->grids(ib).y;
            auto &z = grd_->grids(ib).z;

            auto &cx = grd_->grids(ib).dual_x;
            auto &cy = grd_->grids(ib).dual_y;
            auto &cz = grd_->grids(ib).dual_z;

            // face centers (consistent with your Jac construction style)
            auto xfc_xi = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {
                    0.25 * (x(i, j, k) + x(i, j + 1, k) + x(i, j, k + 1) + x(i, j + 1, k + 1)),
                    0.25 * (y(i, j, k) + y(i, j + 1, k) + y(i, j, k + 1) + y(i, j + 1, k + 1)),
                    0.25 * (z(i, j, k) + z(i, j + 1, k) + z(i, j, k + 1) + z(i, j + 1, k + 1))};
            };
            auto xfc_eta = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {
                    0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j, k + 1) + x(i + 1, j, k + 1)),
                    0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j, k + 1) + y(i + 1, j, k + 1)),
                    0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j, k + 1) + z(i + 1, j, k + 1))};
            };
            auto xfc_zeta = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {
                    0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k)),
                    0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k)),
                    0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k))};
            };

            Int3 lo = Jac.inner_lo();
            Int3 hi = Jac.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // cell center
                        std::array<double, 3> Xc = {
                            cx(i + 1, j + 1, k + 1),
                            cy(i + 1, j + 1, k + 1),
                            cz(i + 1, j + 1, k + 1)};

                        struct FaceEq
                        {
                            std::array<double, 3> n; // normalized S
                            double phi;              // normalized Phi
                            double w;                // weight
                        };
                        FaceEq eqs[6];
                        int K = 0;

                        auto push = [&](const std::array<double, 3> &S,
                                        double Phi,
                                        const std::array<double, 3> &Xf)
                        {
                            double s_norm = norm(S) + eps;
                            std::array<double, 3> nvec = {S[0] / s_norm, S[1] / s_norm, S[2] / s_norm};
                            double phi_hat = Phi / s_norm;

                            double dx = Xf[0] - Xc[0];
                            double dy = Xf[1] - Xc[1];
                            double dz = Xf[2] - Xc[2];
                            double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                            double w = 1.0 / std::sqrt(dist * dist + delta * delta);

                            eqs[K++] = {nvec, phi_hat, w};
                        };

                        // ---------- xi- face (at i) ----------
                        std::array<double, 3> S_xm = {
                            -A_xi(i, j, k, 0),
                            -A_xi(i, j, k, 1),
                            -A_xi(i, j, k, 2)};
                        double Phi_xm = -Bxi(i, j, k, 0);
                        push(S_xm, Phi_xm, xfc_xi(i, j, k));

                        // ---------- xi+ face (at i+1) ----------
                        std::array<double, 3> S_xp = {
                            A_xi(i + 1, j, k, 0),
                            A_xi(i + 1, j, k, 1),
                            A_xi(i + 1, j, k, 2)};
                        double Phi_xp = Bxi(i + 1, j, k, 0);
                        push(S_xp, Phi_xp, xfc_xi(i + 1, j, k));

                        // ---------- eta- face (at j) ----------
                        std::array<double, 3> S_em = {
                            -A_eta(i, j, k, 0),
                            -A_eta(i, j, k, 1),
                            -A_eta(i, j, k, 2)};
                        double Phi_em = -Beta(i, j, k, 0);
                        push(S_em, Phi_em, xfc_eta(i, j, k));

                        // ---------- eta+ face (at j+1) ----------
                        std::array<double, 3> S_ep = {
                            A_eta(i, j + 1, k, 0),
                            A_eta(i, j + 1, k, 1),
                            A_eta(i, j + 1, k, 2)};
                        double Phi_ep = Beta(i, j + 1, k, 0);
                        push(S_ep, Phi_ep, xfc_eta(i, j + 1, k));

                        // ---------- zeta- face (at k) ----------
                        std::array<double, 3> S_zm = {
                            -A_ze(i, j, k, 0),
                            -A_ze(i, j, k, 1),
                            -A_ze(i, j, k, 2)};
                        double Phi_zm = -Bzeta(i, j, k, 0);
                        push(S_zm, Phi_zm, xfc_zeta(i, j, k));

                        // ---------- zeta+ face (at k+1) ----------
                        std::array<double, 3> S_zp = {
                            A_ze(i, j, k + 1, 0),
                            A_ze(i, j, k + 1, 1),
                            A_ze(i, j, k + 1, 2)};
                        double Phi_zp = Bzeta(i, j, k + 1, 0);
                        push(S_zp, Phi_zp, xfc_zeta(i, j, k + 1));

                        // ---------- build normal equations N = A^T W A, r = A^T W phi ----------
                        double N00 = 0, N01 = 0, N02 = 0, N11 = 0, N12 = 0, N22 = 0;
                        double rx = 0, ry = 0, rz = 0;

                        for (int t = 0; t < K; ++t)
                        {
                            double w = eqs[t].w;
                            const auto &n = eqs[t].n;
                            double phi = eqs[t].phi;

                            N00 += w * n[0] * n[0];
                            N01 += w * n[0] * n[1];
                            N02 += w * n[0] * n[2];
                            N11 += w * n[1] * n[1];
                            N12 += w * n[1] * n[2];
                            N22 += w * n[2] * n[2];

                            rx += w * phi * n[0];
                            ry += w * phi * n[1];
                            rz += w * phi * n[2];
                        }

                        auto det3 = [&](double a, double b, double c, double d, double e, double f)
                        {
                            // | a b c |
                            // | b d e |
                            // | c e f |
                            return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
                        };

                        double det = det3(N00, N01, N02, N11, N12, N22);
                        double reg = 1e-14 * (N00 + N11 + N22);

                        if (std::abs(det) < reg)
                        {
                            N00 += reg;
                            N11 += reg;
                            N22 += reg;
                            det = det3(N00, N01, N02, N11, N12, N22);
                        }

                        // cofactors of symmetric matrix
                        double C00 = (N11 * N22 - N12 * N12);
                        double C01 = (N02 * N12 - N01 * N22);
                        double C02 = (N01 * N12 - N02 * N11);
                        double C11 = (N00 * N22 - N02 * N02);
                        double C12 = (N01 * N02 - N00 * N12);
                        double C22 = (N00 * N11 - N01 * N01);

                        double inv = 1.0 / det;

                        double Bx_ind = inv * (C00 * rx + C01 * ry + C02 * rz);
                        double By_ind = inv * (C01 * rx + C11 * ry + C12 * rz);
                        double Bz_ind = inv * (C02 * rx + C12 * ry + C22 * rz);

                        // background field if you want it
                        const double Bx_tot = Bx_ind + B_add_x;
                        const double By_tot = By_ind + B_add_y;
                        const double Bz_tot = Bz_ind + B_add_z;

                        // energy consistency (keep your existing style)
                        const double Bx_old = Bcell(i, j, k, 0);
                        const double By_old = Bcell(i, j, k, 1);
                        const double Bz_old = Bcell(i, j, k, 2);

                        const double Delta_Eb =
                            0.5 * inver_MA2 * (Bx_tot * Bx_tot + By_tot * By_tot + Bz_tot * Bz_tot) -
                            0.5 * inver_MA2 * (Bx_old * Bx_old + By_old * By_old + Bz_old * Bz_old);

                        // U(i, j, k, 4) += Delta_Eb;

                        Bcell(i, j, k, 0) = Bx_tot;
                        Bcell(i, j, k, 1) = By_tot;
                        Bcell(i, j, k, 2) = Bz_tot;
                    }
        }
    }

    void calc_divB()
    {
        const int nblock = fld_->num_blocks();

        for (int ib = 0; ib < nblock; ++ib)
        {
            auto &divB = fld_->field("divB", ib);
            auto &Bxi = fld_->field("B_xi", ib);
            auto &Beta = fld_->field("B_eta", ib);
            auto &Bzeta = fld_->field("B_zeta", ib);
            auto &Jac = fld_->field("Jac", ib);

            Int3 lo = divB.inner_lo();
            Int3 hi = divB.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        divB(i, j, k, 0) = (Bxi(i + 1, j, k, 0) - Bxi(i, j, k, 0) +
                                            Beta(i, j + 1, k, 0) - Beta(i, j, k, 0) +
                                            Bzeta(i, j, k + 1, 0) - Bzeta(i, j, k, 0)) /
                                           Jac(i, j, k, 0);
                    }
        }
    }

    void copy_field()
    {
        for (int ib = 0; ib < fld_->num_blocks(); ib++)
        {
            {
                auto &U = fld_->field("U_", ib);
                auto &old_U = fld_->field("old_U_", ib);
                int ncomp = U.descriptor().ncomp;
                const Int3 &sub = U.get_lo();
                const Int3 &sup = U.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            for (int m = 0; m < ncomp; m++)
                                old_U(i, j, k, m) = U(i, j, k, m);
            }

            {
                auto &U = fld_->field("B_xi", ib);
                auto &old_U = fld_->field("old_B_xi", ib);
                int ncomp = U.descriptor().ncomp;
                const Int3 &sub = U.get_lo();
                const Int3 &sup = U.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            for (int m = 0; m < ncomp; m++)
                                old_U(i, j, k, m) = U(i, j, k, m);
            }
            {
                auto &U = fld_->field("B_eta", ib);
                auto &old_U = fld_->field("old_B_eta", ib);
                int ncomp = U.descriptor().ncomp;
                const Int3 &sub = U.get_lo();
                const Int3 &sup = U.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            for (int m = 0; m < ncomp; m++)
                                old_U(i, j, k, m) = U(i, j, k, m);
            }
            {
                auto &U = fld_->field("B_zeta", ib);
                auto &old_U = fld_->field("old_B_zeta", ib);
                int ncomp = U.descriptor().ncomp;
                const Int3 &sub = U.get_lo();
                const Int3 &sup = U.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            for (int m = 0; m < ncomp; m++)
                                old_U(i, j, k, m) = U(i, j, k, m);
            }
        }
    }
};