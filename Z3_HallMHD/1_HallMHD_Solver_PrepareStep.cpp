#include "HallMHD_Solver.h"

// 目的：让所有用到B_face的地方都不再思考halo
// 输入：B_face_*（本地可能已更新）
// 输出：B_face_*（ghost 最新）
void HallMHDSolver::SyncPrimaryFaceB()
{
    // A) Face 主量的物理边界 + halo
    std::vector<std::string> face_name = {"B_xi", "B_eta", "B_zeta"};
    bound_.add_Face_boundary(face_name);
    halo_->data_trans(face_name[0]);
    halo_->data_trans(face_name[1]);
    halo_->data_trans(face_name[2]);
}

// 目的：只在 inner 计算域从 face 重建 cell 磁场。
// 前置：B_face_* ghost 已有效（来自 SyncPrimaryFaceB）
// 输出：B_cell（只写 inner）
void HallMHDSolver::ComputeBcellInner()
{
    // B) 只计算 inner （计算域）网格的派生量 B_cell
    calc_Bcell();
}

// 目的：让 B_cell 对后续 calc_PV/divB 可用。
// 前置：B_cell(inner) 已写好
// 操作：boundary + halo + corner halo
// 输出：B_cell ghost 最新
void HallMHDSolver::SyncDerivedBcell()
{
    // C) B_cell 的“派生边界策略”先用最简单 copy + halo (及角区通信)
    std::string Bcell_string = "B_cell";
    bound_.add_derived_Cell_boundary(Bcell_string);
    halo_->data_trans(Bcell_string);
    halo_->data_trans_2DCorner(Bcell_string);
    halo_->data_trans_3DCorner(Bcell_string);
}

// 操作：U_ 物理边界 + halo(+corner)
// 输出：U_ ghost 最新
void HallMHDSolver::SyncPrimaryCellU()
{
    // D) U_ 的物理边界（此时可安全用 B_cell） + halo (及角区通信)
    std::string U_string = "U_";
    bound_.add_Cell_boundary(U_string);
    halo_->data_trans(U_string);
    halo_->data_trans_2DCorner(U_string);
    halo_->data_trans_3DCorner(U_string);
}

// 前置：U_ ghost 最新、B_cell ghost 最新
// 操作：calc_PV()、calc_divB()
// 输出：PV_、divB
void HallMHDSolver::UpdateDerivedPVandDivB()
{
    // E) 原始量 + divB
    calc_PV();
    calc_divB();
}

// 只做 copy_field()（保存 old_*）
void HallMHDSolver::SnapshotOldFields()
{
    // F) 复制物理场，用于残差计算
    copy_field();
}
