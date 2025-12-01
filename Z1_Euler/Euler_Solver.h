#pragma once
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"
#include "Control.h"
#include "Output.h"
#include "Boundary.h"
#include "Initial.h"

class EulerSolver
{
public:
    EulerSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par, std::string Solver_Name)
    {
        grd_ = grd;
        fld_ = fld;
        halo_ = halo;
        par_ = par;
        topo_ = topo;

        Solver_Name_ = Solver_Name;

        fid_U = fld_->field_id(Solver_Name_);
        fid_PV = fld_->field_id("PV_");

        fid_Jac = fld_->field_id("Jac");
        fid_Xi = fld_->field_id("JDxi");
        fid_Eta = fld_->field_id("JDet");
        fid_Zeta = fld_->field_id("JDze");

        gamma_ = par_->GetDou_List("constant").data["gamma"];

        control_.SetUp(par_, fld_->field(Solver_Name_)[0].descriptor().ncomp);
        output_.SetUp(par_, fld_);
        bound_.SetUp(grd_, fld_, topo_, par_);
        initial_.Initialization(fld_);

        fld_->register_field(FieldDescriptor{"old_U_", StaggerLocation::Cell, 5, par->GetInt("ngg")});
        fld_->register_field(FieldDescriptor{"RHS", StaggerLocation::Cell, 5, par->GetInt("ngg")});
    }

    void Advance()
    {
        // 0. halo + 物理边界 + 计算PV_
        calc_PV();
        halo_->exchange_inner(Solver_Name_);
        halo_->exchange_parallel(Solver_Name_);
        bound_.add_boundary(Solver_Name_);
        calc_PV();

        while (true)
        {

            // // 1. 计算时间步长
            compute_timestep();

            // // 2. 时间推进（e.g. 最简单的 forward Euler）:计算 RHS: dU/dt = -div F(U)
            time_advance();

            // 3. 计算残差
            calc_Residual();

            // 4. halo + 物理边界+ 计算PV_
            halo_->exchange_inner(Solver_Name_);
            halo_->exchange_parallel(Solver_Name_);
            bound_.add_boundary(Solver_Name_);
            calc_PV();
            copy_field();

            // 5. 更新控制
            control_.Update();
            if (control_.if_outfile)
                output_.output_plt_field(); // output_field();
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

    Euler_Control control_;
    Euler_Output output_;
    Euler_Boundary bound_;
    Euler_Initial initial_;

    int fid_Jac, fid_Xi, fid_Eta, fid_Zeta;

    int fid_U;
    int fid_PV;
    std::string Solver_Name_;
    double gamma_;

    double dt;

    void compute_timestep();
    void time_advance();
    void calc_Residual();

    void inv_rhs();

    struct Double3
    {
        double vector[3];
        Double3 &operator+(Double3 add)
        {
            vector[0] += add.vector[0];
            vector[1] += add.vector[1];
            vector[2] += add.vector[2];
            return *this;
        }
        Double3 &operator-(Double3 add)
        {
            vector[0] -= add.vector[0];
            vector[1] -= add.vector[1];
            vector[2] -= add.vector[2];
            return *this;
        }
        Double3 &operator*(double add)
        {
            vector[0] *= add;
            vector[1] *= add;
            vector[2] *= add;
            return *this;
        }
        double operator*(Double3 add)
        {
            return vector[0] * add.vector[0] + vector[1] * add.vector[1] + vector[2] * add.vector[2];
        }
        Double3 &operator()(double a, double b, double c)
        {
            vector[0] = a;
            vector[1] = b;
            vector[2] = c;
            return *this;
        }
    };

    void Reconstruction(double *metric, int32_t direction, FieldBlock &PV, FieldBlock &U, int iblock, int index_i, int index_j, int index_k, double *out_flux);

    void calc_PV()
    {
        for (int ib = 0; ib < fld_->num_blocks(); ib++)
        {
            auto &U = fld_->field(fid_U, ib);
            auto &PV = fld_->field(fid_PV, ib);
            const Int3 &sub = U.get_lo();
            const Int3 &sup = U.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        double rho = U(i, j, k, 0);
                        double u = U(i, j, k, 1) / rho;
                        double v = U(i, j, k, 2) / rho;
                        double w = U(i, j, k, 3) / rho;
                        PV(i, j, k, 0) = u;
                        PV(i, j, k, 1) = v;
                        PV(i, j, k, 2) = w;
                        PV(i, j, k, 3) = (U(i, j, k, 4) - 0.5 * rho * (u * u + v * v + w * w)) * (gamma_ - 1.0);
                    }
        }
    }

    void copy_field()
    {
        for (int ib = 0; ib < fld_->num_blocks(); ib++)
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
    }
};