#include "HallMHD_Solver.h"

// 这里只实现ComputeRHShallFromCurrentBface_： 求解curl E_hall，这将作为函数指针传入ImplicitHall_Sovler中调用
// 从而实现隐式迭代求解
void HallMHDSolver::ComputeRHShallFromCurrentBface_()
{
    // 1) 清零 RHShall_*（不要清 RHS_*）
    // 2) 用当前 B_face 计算 J_edge
    // 3) 用 J_edge 和 B 计算 Ehall_edge
    // 4) curl(Ehall_edge) -> RHShall_face
}
