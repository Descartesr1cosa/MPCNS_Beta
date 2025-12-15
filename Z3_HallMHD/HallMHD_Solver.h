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

class HallMHDSolver
{
public:
    HallMHDSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par, std::vector<std::string> Solver_Name)
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

        fid_F_[0] = fld_->field_id("F_xi");
        fid_F_[1] = fld_->field_id("F_eta");
        fid_F_[2] = fld_->field_id("F_zeta");

        fid_Eface_[0] = fld_->field_id("E_face_xi");
        fid_Eface_[1] = fld_->field_id("E_face_eta");
        fid_Eface_[2] = fld_->field_id("E_face_zeta");

        fid_Bface_[0] = fid_Bxi;
        fid_Bface_[1] = fid_Beta;
        fid_Bface_[2] = fid_Bzeta;

        fid_metric_[0] = fid_Xi;
        fid_metric_[1] = fid_Eta;
        fid_metric_[2] = fid_Zeta;

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
            if (StepOnce())
                break;
        }
    }

private:
    // --- core pointers ---
    Grid *grd_;
    TOPO::Topology *topo_;
    Field *fld_;
    Halo *halo_;
    Param *par_;

    // --- components ---
    MHD_Control control_;
    MHD_Output output_;
    MHD_Boundary bound_;
    MHD_Initial initial_;

    // --- cached field ids / constants ---
    int fid_Jac, fid_Xi, fid_Eta, fid_Zeta;
    int fid_U;
    int fid_PV;
    int fid_Bcell, fid_Bxi, fid_Beta, fid_Bzeta;

    int fid_F_[3], fid_Eface_[3], fid_Bface_[3], fid_metric_[3];

    std::vector<std::string> Solver_Name_;
    double gamma_, B_add_x, B_add_y, B_add_z, inver_MA2;
    double dt;

    // ========== Step driver ==========
    bool StepOnce()
    {
        Compute_Timestep();
        Time_Advance(); // 内部只编排：ZeroRHS/Assemble/Update
        Calc_Residual();
        PrepareStep(); // 用于下一步/输出前字段一致
        return UpdateControlAndOutput();
    };

    void Compute_Timestep();
    void Calc_Residual();
    bool UpdateControlAndOutput();

    // ========== Prepare / sync pipeline ==========
    void PrepareStep()
    {
        SyncPrimaryFaceB();       // B_face: BC + halo
        ComputeBcellInner();      // B_cell: inner compute
        SyncDerivedBcell();       // B_cell: derived BC + halo (+corner)
        SyncPrimaryCellU();       // U: BC + halo (+corner)
        UpdateDerivedPVandDivB(); // PV + divB
        SnapshotOldFields();      // old_* for residual/monitor
    };

    void SyncPrimaryFaceB();       // B_face_*（ghost）
    void ComputeBcellInner();      // 在 inner 域从 B_face_* 重建 cell 磁场。
    void SyncDerivedBcell();       // 更新 B_cell 的 ghost
    void SyncPrimaryCellU();       // 守恒变量U添加边界条件
    void UpdateDerivedPVandDivB(); // 计算原始变量等被动量
    void SnapshotOldFields();      // 拷贝保存当前场，用于残差计算

    // ========== RHS assembly ==========
    void Time_Advance()
    {
        ZeroRHS();
        AssembleRHS_Fluid();     // 原 inv_fluid()
        SyncElectricFace();      // E_face_*: BC + halo
        AssembleRHS_Induction(); // 原 inv_induce()，后会在这里插 Hall
        ApplyTimeUpdate_Euler(); // U += dt*RHS, B += dt*RHS_B
        Update_Physic_Time();    // 记录物理时间
    };

    void ZeroRHS();
    void AssembleRHS_Fluid();     // inv_fluid()
    void SyncElectricFace();      // E_face_* BC+halo
    void AssembleRHS_Induction(); // inv_induce()
    void ApplyTimeUpdate_Euler(); // += dt*RHS
    void Update_Physic_Time();    // record physical time

    //------------------------------------
    // helper for AssembleRHS_Fluid()
    void AssembleOneDirectionFluxAndEMF_(int iblk,
                                         int dir,                 // 0 xi, 1 eta, 2 zeta
                                         FieldBlock &flux,        // F_xi / F_eta / F_zeta (ncomp=5)
                                         FieldBlock &E_face,      // E_face_xi/eta/zeta   (ncomp=3)
                                         FieldBlock &B_face,      // B_xi/eta/zeta        (ncomp=1)
                                         FieldBlock &metricField, // Xi_/Eta_/Zeta_       (ncomp=3)
                                         FieldBlock &PV,
                                         FieldBlock &U,
                                         FieldBlock &Bcell);
    void AssembleCellRHSFromFlux_();
    // Reconstruction / flux
    void Reconstruction(double *metric, int32_t direction, FieldBlock &PV, FieldBlock &U, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k, double *out_flux);

    //------------------------------------
    // helper for AssembleRHS_Induction()
    void AssembleEdgeEMF_FromFaceE_Ideal_();
    void ApplyBC_EdgeEMF_();
    void AssembleFaceRHS_FromEdgeEMF_Curl_();

    // ========== Derived / diagnostics ==========
    void calc_PV();
    void calc_Bcell();
    void calc_divB();
    void copy_field();

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
};