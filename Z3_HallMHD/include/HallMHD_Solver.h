#pragma once

#include <cstdint>

// Z3_HallMHD
#include "0_state/SolverFields.h"
#include "1_boundary/Boundary.h"
#include "3_io/Output.h"
#include "3_io/Initial.h"
#include "4_solver/Control.h"
#include "0_state/HallConfig.h"

// ---- forward declarations (avoid heavy includes in header) ----
class Grid;
namespace TOPO
{
    class Topology;
}
class Field;
class Halo;
class Param;
class FieldBlock;
class ImplicitHall_Solver;

class HallMHDSolver
{
public:
    HallMHDSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par, ImplicitHall_Solver *hall_imp);

    void Advance();

private:
    // --- core pointers ---
    Grid *grd_{nullptr};
    TOPO::Topology *topo_{nullptr};
    Field *fld_{nullptr};
    Halo *halo_{nullptr};
    Param *par_{nullptr};
    ImplicitHall_Solver *hall_{nullptr};

    // --- components ---
    MHD_Control control_;
    MHD_Output output_;
    MHD_Boundary bound_;
    MHD_Initial initial_;

    // --- cached field ids  ---
    SolverFields fid_;

    // --- constants ---
    double gamma_{0.0};
    double inver_MA2{0.0};
    double hall_coef{0.0};
    double dt{0.0};

private:
    // =============================== Step driver ============================
    bool StepOnce();
    //-------------------------------------------------------------------------
    void Compute_Timestep();
    void Calc_Residual();
    bool UpdateControlAndOutput();
    //=========================================================================

    // ======================== Prepare / Sync pipeline =======================
    void PrepareStep();
    void PrepareSubstep_NoSnapshot();
    //-------------------------------------------------------------------------
    void SyncPrimaryFaceB();       // B_face_*（ghost）
    void ComputeBcellInner();      // 在 inner 域从 B_face_* 重建 cell 磁场。
    void SyncDerivedBcell();       // 更新 B_cell 的 ghost
    void SyncPrimaryCellU();       // 守恒变量U添加边界条件
    void UpdateDerivedPVandDivB(); // 计算原始变量等被动量
    void SnapshotOldFields();      // 拷贝保存当前场，用于残差计算
    //=========================================================================

    // ============================= RHS assembly =============================
    void Time_Advance();
    //-------------------------------------------------------------------------
    void ZeroRHS();
    void AssembleRHS_Fluid();     // inv_fluid()
    void SyncElectricFace();      // E_face_* BC+halo
    void AssembleRHS_Induction(); // inv_induce()
    void ApplyTimeUpdate_Euler(); // += dt*RHS
    void Update_Physic_Time();    // record physical time
    //=========================================================================

    // ========================== helper for Fluid ============================
    void AssembleOneDirectionFluxAndEMF_(int iblk,
                                         int dir,                 // 0 xi, 1 eta, 2 zeta
                                         FieldBlock &flux,        // F_xi / F_eta / F_zeta (ncomp=5)
                                         FieldBlock &E_face,      // E_face_xi/eta/zeta   (ncomp=3)
                                         FieldBlock &B_face,      // B_xi/eta/zeta        (ncomp=1)
                                         FieldBlock &B_face_add,  // B_xi/eta/zeta add        (ncomp=1)
                                         FieldBlock &metricField, // Xi_/Eta_/Zeta_       (ncomp=3)
                                         FieldBlock &PV,
                                         FieldBlock &U,
                                         FieldBlock &Bcell);
    void AssembleCellRHSFromFlux_();
    // Reconstruction / flux
    void Reconstruction(double *metric, int32_t direction, FieldBlock &PV, FieldBlock &U, FieldBlock &B_cell, double B_jac_nabla, int iblock, int index_i, int index_j, int index_k, double *out_flux);
    //=========================================================================

    // ========================= helper for Induction =========================
    void AssembleEdgeEMF_FromFaceE_Ideal_();
    void AddExplicitHallToEdgeEMF_(); // 只在 hall_explicit.cpp 实现由宏控制, 对于Ideal Implicit均为空
    void ApplyBC_EdgeEMF_();
    void AssembleFaceRHS_FromEdgeEMF_Curl_();

    // for Hall MHD
#if HALL_MODE != 0
    void ComputeJ_AtEdges_Inner_();
    void ApplyBC_EdgeJ_();
    void ComputeHallE_AtEdges_EnergyPreserving_(); // 生成 Ehall_xi/Ehall_eta/Ehall_zeta（标量线积分）
#endif

    // for Hall MHD Explicit
#if HALL_MODE == 1
    void AccumulateHallE_ToTotalEdgeEMF_(); // E_* += Ehall_*
#endif

    // for Hall MHD Implicit
#if HALL_MODE == 2
    void ComputeRHShallFromCurrentBface_();  // 真正干活的成员函数：计算 RHShall_*（写入 face fields）
    static void RHShallCallback_(void *ctx); // 静态 wrapper：签名匹配 void(*)(void*)
    void Modify_TotalEnergy_AfterHall();     // Hall只修改了B, 更新B_cell 然后把U_的能量修改即可
    //------------------------------------
    void ZeroRHShall();
    void ApplyBC_EdgeHallEMF_(); // Eelectric edge Boundary
    void AssembleFaceRHSHall_FromEdgeHallEMF_Curl_();
    void UpdateTotalEnergy();
#endif
    //=========================================================================

    // ================================== TOOLS ==============================
    void calc_physical_constant(Param *par);
    void calc_PV();
    void calc_Bcell();
    void calc_divB();
    void copy_field();
    void PrintMinMaxDiagnostics_();
    void add_Emag_to_Etotal();
    double ComputeLocalMaxRadius_();
    //=========================================================================
};