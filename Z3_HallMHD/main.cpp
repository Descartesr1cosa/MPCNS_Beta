//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"
#include "HallMHD_Solver.h"
//==============================================================================

//==============================================================================

//==============================================================================

int main(int arg, char **argv)
{
    //=============================================================================================
    // MPI initialization
    int myid;
    PARALLEL::mpi_initial(arg, argv);
    PARALLEL::mpi_rank(&myid);
    //=============================================================================================

    //=============================================================================================
    //--------------------------------------------------------------------------
    // 读入控制参数
    Param *par = new Param;
    par->ReadParam(myid);
    //--------------------------------------------------------------------------
    // 读入网格并作预处理
    Grid *grd = new Grid;
    grd->Grid_Preprocess(par);
    //--------------------------------------------------------------------------
    // 建立topology
    TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
    //--------------------------------------------------------------------------
    // 建立Field
    Field *fld = new Field(grd, par);
    //-------------------------------------
    // 加入求解物理场
    int ngg = par->GetInt("ngg");
    // 守恒变量、独立变量，用于构建CT方法
    fld->register_field({"U_", StaggerLocation::Cell, 5, ngg});
    fld->register_field({"B_xi", StaggerLocation::FaceXi, 1, ngg});
    fld->register_field({"B_eta", StaggerLocation::FaceEt, 1, ngg});
    fld->register_field({"B_zeta", StaggerLocation::FaceZe, 1, ngg});

    // 辅助电场变量
    fld->register_field({"E_xi", StaggerLocation::EdgeXi, 1, ngg});
    fld->register_field({"E_eta", StaggerLocation::EdgeEt, 1, ngg});
    fld->register_field({"E_zeta", StaggerLocation::EdgeZe, 1, ngg});
#if HALL_MODE != 0
    // 辅助Hall电场变量
    fld->register_field({"Ehall_xi", StaggerLocation::EdgeXi, 1, ngg});
    fld->register_field({"Ehall_eta", StaggerLocation::EdgeEt, 1, ngg});
    fld->register_field({"Ehall_zeta", StaggerLocation::EdgeZe, 1, ngg});
    fld->register_field({"J_xi", StaggerLocation::EdgeXi, 1, ngg});
    fld->register_field({"J_eta", StaggerLocation::EdgeEt, 1, ngg});
    fld->register_field({"J_zeta", StaggerLocation::EdgeZe, 1, ngg});
#endif

    // 辅助磁场
    fld->register_field({"B_cell", StaggerLocation::Cell, 3, ngg});
    // 流体方程通量，辅助
    fld->register_field(
        FieldDescriptor{"F_xi", StaggerLocation::FaceXi, 5, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"F_eta", StaggerLocation::FaceEt, 5, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"F_zeta", StaggerLocation::FaceZe, 5, par->GetInt("ngg")});
    // 原始变量
    fld->register_field(
        FieldDescriptor{"PV_", StaggerLocation::Cell, 4, par->GetInt("ngg")});
    // 网格界面(Face)的电场向量，用于CT方法
    fld->register_field(
        FieldDescriptor{"E_face_xi", StaggerLocation::FaceXi, 3, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"E_face_eta", StaggerLocation::FaceEt, 3, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"E_face_zeta", StaggerLocation::FaceZe, 3, par->GetInt("ngg")});
    //--------------------------------------------------------------------------
    // 建立Halo通信
    Halo *hal = new Halo(fld, &topology);
    //=============================================================================================

    //=============================================================================================
    HallMHDSolver solver(grd, &topology, fld, hal, par, {"U_", "B_xi", "B_eta", "B_zeta"});
    solver.Advance();
    //=============================================================================================

    //=============================================================================================
    //--------------------------------------------------------------------------
    // MPI终止
    PARALLEL::mpi_finalize();
    //--------------------------------------------------------------------------
    // 释放所分配的空间，建议按照创建顺序逆序释放
    delete hal;
    delete fld;
    delete par;
    delete grd;
    //=============================================================================================
    return 0;
}