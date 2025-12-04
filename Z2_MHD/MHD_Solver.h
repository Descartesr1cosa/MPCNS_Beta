#pragma once
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"
#include "Control.h"
#include "Output.h"
#include "Boundary.h"
#include "Initial.h"

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
        std::string U_string = "U_";
        std::string Bcell_string = "B_cell";
        // 0. halo + 物理边界 + 计算Bcell + PV_
        {
            // 计算更新B_cell
            calc_Bcell();

            // 传递更新虚网格的主控变量 U_ B_xi eta zeta
            bound_.add_boundary(Solver_Name_);
            for (auto &fld_name : Solver_Name_)
            {
                halo_->data_trans(fld_name);
            }
            halo_->data_trans_2DCorner(U_string);
            halo_->data_trans_3DCorner(U_string);

            // 计算原始变量
            calc_PV();

            copy_field();
        }

        while (true)
        {

            // 1. 计算时间步长
            compute_timestep();

            // 2. 时间推进（e.g. 最简单的 forward Euler）:计算 RHS: dU/dt = -div F(U)
            time_advance();

            // 3. 计算残差
            calc_Residual();

            // 4. halo + 物理边界+ 计算Bcell + PV_
            {
                // 计算更新B_cell
                calc_Bcell();

                // 传递更新虚网格的主控变量 U_ B_xi eta zeta
                bound_.add_boundary(Solver_Name_);
                for (auto &fld_name : Solver_Name_)
                {
                    halo_->data_trans(fld_name);
                }
                halo_->data_trans_2DCorner(U_string);
                halo_->data_trans_3DCorner(U_string);

                // 计算原始变量
                calc_PV();
                calc_divB();
                copy_field();
            }

            // 5. 更新控制
            control_.Update();
            if (control_.if_outfile)
                output_.output_plt_cell_field(output_.var_defaut_plt_name);
            // output_.output_plt_field(); // output_.output_plt_cell_field(output_.var_defaut_plt_name); // output_field();
            if (control_.if_stop)
            {
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

                        double E = U(i, j, k, 4); //  MHD 总能量
                        PV(i, j, k, 3) = (E - kin - 0.5 * B2 * inver_MA2) * (gamma_ - 1.0);
                    }
        }
    }

    void calc_Bcell()
    {
        const int nblock = fld_->num_blocks();

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

            // 这里用几何体积的 inner_lo/inner_hi 作为 cell 中心循环范围
            Int3 lo = Jac.inner_lo();
            Int3 hi = Jac.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double J = Jac(i, j, k, 0);

                        // ---- 2.1 cell 中心的对偶分量 B^ξ, B^η, B^ζ （face 值平均）----

                        const double Bxi_m = Bxi(i, j, k, 0);
                        const double Bxi_p = Bxi(i + 1, j, k, 0);
                        const double Beta_m = Beta(i, j, k, 0);
                        const double Beta_p = Beta(i, j + 1, k, 0);
                        const double Bze_m = Bzeta(i, j, k, 0);
                        const double Bze_p = Bzeta(i, j, k + 1, 0);

                        const double Bxi_c = 0.5 * (Bxi_m + Bxi_p);
                        const double Beta_c = 0.5 * (Beta_m + Beta_p);
                        const double Bzeta_c = 0.5 * (Bze_m + Bze_p);

                        // ---- 2.2 cell 中心的面积向量 Aξ, Aη, Aζ （JDxi/JDet/JDze 平均）----

                        // Aξ 在 i 和 i+1 两个 xi-face 之间平均
                        double Axi_c[3] = {
                            0.5 * (A_xi(i, j, k, 0) + A_xi(i + 1, j, k, 0)),
                            0.5 * (A_xi(i, j, k, 1) + A_xi(i + 1, j, k, 1)),
                            0.5 * (A_xi(i, j, k, 2) + A_xi(i + 1, j, k, 2))};

                        // Aη 在 j 和 j+1 两个 eta-face 之间平均
                        double Aeta_c[3] = {
                            0.5 * (A_eta(i, j, k, 0) + A_eta(i, j + 1, k, 0)),
                            0.5 * (A_eta(i, j, k, 1) + A_eta(i, j + 1, k, 1)),
                            0.5 * (A_eta(i, j, k, 2) + A_eta(i, j + 1, k, 2))};

                        // Aζ 在 k 和 k+1 两个 zeta-face 之间平均
                        double Aze_c[3] = {
                            0.5 * (A_ze(i, j, k, 0) + A_ze(i, j, k + 1, 0)),
                            0.5 * (A_ze(i, j, k, 1) + A_ze(i, j, k + 1, 1)),
                            0.5 * (A_ze(i, j, k, 2) + A_ze(i, j, k + 1, 2))};

                        // ---- 2.3 按 B = (1/J) [ Bξ Aξ + Bη Aη + Bζ Aζ ] 组装物理 B_ind ----
                        double Bx_ind =
                            (Bxi_c * Axi_c[0] + Beta_c * Aeta_c[0] + Bzeta_c * Aze_c[0]) / J;

                        double By_ind =
                            (Bxi_c * Axi_c[1] + Beta_c * Aeta_c[1] + Bzeta_c * Aze_c[1]) / J;

                        double Bz_ind =
                            (Bxi_c * Axi_c[2] + Beta_c * Aeta_c[2] + Bzeta_c * Aze_c[2]) / J;

                        // ---- 2.4 加上背景场----
                        const double Bx_tot = Bx_ind + B_add_x;
                        const double By_tot = By_ind + B_add_y;
                        const double Bz_tot = Bz_ind + B_add_z;

                        //----- 2.5 修改磁场和磁能----
                        const double Bx_tot_old = Bcell(i, j, k, 0);
                        const double By_tot_old = Bcell(i, j, k, 1);
                        const double Bz_tot_old = Bcell(i, j, k, 2);
                        const double Delta_Eb = 0.5 * inver_MA2 * (Bx_tot * Bx_tot + By_tot * By_tot + Bz_tot * Bz_tot) - 0.5 * inver_MA2 * (Bx_tot_old * Bx_tot_old + By_tot_old * By_tot_old + Bz_tot_old * Bz_tot_old);

                        U(i, j, k, 4) += Delta_Eb;
                        Bcell(i, j, k, 0) = Bx_tot;
                        Bcell(i, j, k, 1) = By_tot;
                        Bcell(i, j, k, 2) = Bz_tot;
                    }
        }
        std::string Bcell_string = "B_cell";
        bound_.cell_copy_boundary(Bcell_string);
        halo_->data_trans(Bcell_string);
        halo_->data_trans_2DCorner(Bcell_string);
        halo_->data_trans_3DCorner(Bcell_string);
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