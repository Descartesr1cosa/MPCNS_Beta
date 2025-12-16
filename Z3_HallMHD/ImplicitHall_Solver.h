// ImplicitHall_Solver.h
#pragma once
#include <iostream>
#include <cstdlib>

#include "petscksp.h"
#include "petscsnes.h"

#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"
#include "Boundary.h"

class ImplicitHall_Solver
{
public:
    // 构造函数：在这里做 PetscInitialize
    // main 里是： HallImplicitSolver* hall_imp = new HallImplicitSolver(arg, argv);
    ImplicitHall_Solver(int arg, char **argv)
    {
        PetscErrorCode ierr = PetscInitialize(&arg, &argv, nullptr, nullptr);
        if (ierr != 0)
        {
            // 简单粗暴一点，可以直接抛异常或 abort
            std::cout << "Fatal Error! !  PetscInitialize failed in ImplicitHall_Solver ! ! !\n";
            exit(-999);
        }
        initialized_ = true;
    }
    ~ImplicitHall_Solver()
    {
        // 这里先什么都不做，真正的 PetscFinalize 由 Finalization() 负责，这里防止异常退出
        Finalization();
    };

    // 显式的收尾接口，在 main 里调用：
    // hall_imp->Finalization();
    void Finalization()
    {
        if (!initialized_)
            return;

        destroy_petsc_objects_();

        PetscBool finalized = PETSC_FALSE;
        PetscFinalized(&finalized);
        if (!finalized)
        {
            destroy_petsc_objects_();
            PetscFinalize();
        }
        initialized_ = false;
    }

    // 之后再实现：传入 grd/fld，建立 dof、创建 SNES/Vec 等
    void setup(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par, MHD_Boundary *bound)
    {
        grd_ = grd;
        topo_ = topo;
        fld_ = fld;
        halo_ = halo;
        par_ = par;
        bound_ = bound;

        //-------------------------------------------------------------------------
        // --- cache field ids (minimum) ---
        fid_Bxi = fld_->field_id("B_xi");
        fid_Beta = fld_->field_id("B_eta");
        fid_Bzeta = fld_->field_id("B_zeta");
        fid_RHShall_xi = fld_->field_id("RHShall_xi");
        fid_RHShall_eta = fld_->field_id("RHShall_eta");
        fid_RHShall_zeta = fld_->field_id("RHShall_zeta");
        PARALLEL::mpi_rank(&myid);
        PARALLEL::mpi_size(&rank_num);

        B_face_name.resize(3);
        B_face_name[0] = "B_xi";
        B_face_name[1] = "B_eta";
        B_face_name[2] = "B_zeta";
        //-------------------------------------------------------------------------

        //-------------------------------------------------------------------------
        // --- build owner and map ---
        build_owner_mask_faces_(); // 创建owner，初始化为-1 0 1分别代表被动量，非DOF和DOF

        assign_global_ids_to_owners_(); // 对owner = 1的进行统计，获得全局唯一的DOF编号
        build_local_dof_map_faces_owner_only_();

        fill_passiveowner_inner_();    // 对owner = -1的修改为真实DOF的负值
        fill_passiveowner_parallel_(); // 对owner = -1的修改为真实DOF的负值
        //-------------------------------------------------------------------------

        //-------------------------------------------------------------------------
        // --- create petsc ---
        create_petsc_objects_();
        //-------------------------------------------------------------------------
    };

    void solve_implicit_hall(double dt);

    void vec2B(Vec X)
    {
        unpack_vec_to_Bface_owner_only_(X);
        InterfaceExchange_Bface();
        bound_->add_Face_boundary(B_face_name);
        halo_->data_trans(B_face_name[0]);
        halo_->data_trans(B_face_name[1]);
        halo_->data_trans(B_face_name[2]);
    }

    void B2vec(Vec X)
    {
        pack_Bface_owner_only_to_vec_(X);
    }

private:
    bool initialized_ = false;
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Halo *halo_ = nullptr;
    Param *par_ = nullptr;
    MHD_Boundary *bound_ = nullptr;

    int fid_Bxi, fid_Beta, fid_Bzeta;
    int fid_RHShall_xi = -1, fid_RHShall_eta = -1, fid_RHShall_zeta = -1;
    std::vector<std::string> B_face_name;

    int32_t myid, rank_num;
    double dt_;

    //-------------------------------------------------------------------------
    // ownership array
    using OwnerFlag = FieldScalar<PetscInt>;
    std::vector<OwnerFlag> owner_bxi_, owner_beta_, owner_bzeta_;
    //=========================================================================

    //=========================================================================
    // build owner and be initialized as -1(passive dof) 0 (not dof) 1(dof)
    void build_owner_mask_faces_()
    {
        alloc_owner_by_B_inner_and_fill1_();
        apply_physical_bc_mask_();
        apply_interface_nonowner_mask_();
    }
    //-------------------------------------------------------------------------
    // 开辟owner空间，初始化为1
    void alloc_owner_by_B_inner_and_fill1_();

    // 将物理边界的face都设置为0
    void apply_physical_bc_mask_();

    // 将inner parallel边界的face按照优先级(rank block i j k)取较小的保留为1，否则为被动量并设置为-1
    void apply_interface_nonowner_mask_();

    //-------------------------------------------------------------------------
    // tool find out which face the Box3 is
    int face_axis_from_node_box_(const Box3 &b);

    // tool mark the whole face (range of which is Box3) as val
    void mark_face_patch_(int ib, const Box3 &box, int axis, PetscInt val);

    // tool 用于比较优先级，确定ownership
    struct PatchKey
    {
        int rank, block;
        int axis;       // 该 patch 是 Xi/Eta/Zeta 哪类面
        int fixed;      // 该面固定坐标（Xi 面用 lo.i，Eta 用 lo.j，Zeta 用 lo.k）
        int lo_u, hi_u; // 按照ijk顺序，除去固定坐标之后的高优先级坐标
        int lo_v, hi_v; // 按照ijk顺序，除去固定坐标之后的低优先级坐标

        bool operator<(const PatchKey &o) const
        {
            if (rank != o.rank)
                return rank < o.rank;
            if (block != o.block)
                return block < o.block;
            if (axis != o.axis)
                return axis < o.axis;
            if (fixed != o.fixed)
                return fixed < o.fixed;
            if (lo_u != o.lo_u)
                return lo_u < o.lo_u;
            if (hi_u != o.hi_u)
                return hi_u < o.hi_u;
            if (lo_v != o.lo_v)
                return lo_v < o.lo_v;
            return hi_v < o.hi_v;
        }
    };

    // tool 构造一个patch的PatchKey
    PatchKey make_key_(int rank, int block, const Box3 &b);

    // tool 根据PatchKey判断自己是否是owner
    bool this_side_is_owner_(const TOPO::InterfacePatch &p, bool is_parallel);
    //=========================================================================

    //=========================================================================
    // 给所有1的owner进行编号，获得全局唯一的DOF号码，从1开始
    void assign_global_ids_to_owners_();
    // 统计本进程的总DOF数量
    PetscInt count_owned_local_();
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // 对owner = -1的修改为真实DOF的负值
    void fill_passiveowner_inner_();
    void fill_passiveowner_parallel_();
    //-------------------------------------------------------------------------
    // tool ijk的坐标变换
    Int3 apply_trans_node_(const TOPO::IndexTransform &T, int i, int j, int k);
    // tool 计box中face的数量
    int face_count_from_node_box_(const Box3 &b, int axis);
    // 线性下标按“接收侧 natural loop 顺序”定义：
    // axis=0: j outer, k inner
    // axis=1: i outer, k inner
    // axis=2: i outer, j inner
    int face_linear_index_in_nodebox_(const Box3 &b, int axis, int i, int j, int k);
    // 取 owner 数组的引用（按 comp=0/1/2）
    FieldScalar<PetscInt> &owner_arr_(int comp, int ib)
    {
        // 取 owner 数组的引用（按 comp=0/1/2）
        if (comp == 0)
            return owner_bxi_[ib];
        if (comp == 1)
            return owner_beta_[ib];
        return owner_bzeta_[ib];
    }
    //=========================================================================
    void InterfaceExchange_Bface()
    {
        InterfaceExchange_Bface_Inner();
        InterfaceExchange_Bface_Parallel();
    }
    void InterfaceExchange_Bface_Inner();
    void InterfaceExchange_Bface_Parallel();
    // 访问 B_face 的 helper
    FieldBlock &B_face_(int axis, int ib)
    {
        if (axis == 0)
            return fld_->field(fid_Bxi, ib);
        if (axis == 1)
            return fld_->field(fid_Beta, ib);
        return fld_->field(fid_Bzeta, ib);
    }
    //=========================================================================
    struct FaceDofRef
    {
        int axis; // 0 xi, 1 eta, 2 zeta
        int ib;   // block id
        int i, j, k;
    };
    PetscInt ndof_local_owned_ = 0;      // 本 rank owned DOF 数
    PetscInt ndof_global_ = 0;           // 全局总DOF
    PetscInt start_gid_ = 1;             // 本 rank owned gid 起点（1-based）
    std::vector<FaceDofRef> local_dofs_; // size = ndof_local_owned_
    void build_local_dof_map_faces_owner_only_();
    //=========================================================================
    void pack_Bface_owner_only_to_vec_(Vec X)
    {
        PetscScalar *x = nullptr;
        VecGetArray(X, &x);

        for (PetscInt lid = 0; lid < ndof_local_owned_; ++lid)
        {
            const auto &d = local_dofs_[(size_t)lid];
            FieldBlock &B = (d.axis == 0   ? fld_->field(fid_Bxi, d.ib)
                             : d.axis == 1 ? fld_->field(fid_Beta, d.ib)
                                           : fld_->field(fid_Bzeta, d.ib));

            x[lid] = (PetscScalar)B(d.i, d.j, d.k, 0);
        }

        VecRestoreArray(X, &x);
    }
    void unpack_vec_to_Bface_owner_only_(Vec X)
    {
        const PetscScalar *x = nullptr;
        VecGetArrayRead(X, &x);

        for (PetscInt lid = 0; lid < ndof_local_owned_; ++lid)
        {
            const auto &d = local_dofs_[(size_t)lid];
            FieldBlock &B = (d.axis == 0   ? fld_->field(fid_Bxi, d.ib)
                             : d.axis == 1 ? fld_->field(fid_Beta, d.ib)
                                           : fld_->field(fid_Bzeta, d.ib));

            B(d.i, d.j, d.k, 0) = (double)x[lid];
        }

        VecRestoreArrayRead(X, &x);
    }
    //=========================================================================

    //=========================================================================
    // PETSc objects
    SNES snes_ = nullptr;
    Vec X_ = nullptr;
    Vec F_ = nullptr;
    Vec Xref_ = nullptr; // 框架自检用（以后可替换为 B*）
    bool petsc_objs_built_ = false;
    //-------------------------------------------------------------------------
    void create_petsc_objects_();
    void destroy_petsc_objects_();
    PetscErrorCode FormFunction_(Vec X, Vec F);
    static PetscErrorCode FormFunction_PETSc(SNES, Vec X, Vec F, void *ctx)
    {
        return static_cast<ImplicitHall_Solver *>(ctx)->FormFunction_(X, F);
    }
    //=========================================================================

    //=========================================================================
    // 一个回调函数类型：输入一个 void *上下文，不返回值
    using RHShallCallback = void (*)(void *);

    RHShallCallback rhshall_cb_ = nullptr;
    void *rhshall_ctx_ = nullptr;
    // 注册回调：cb 是函数地址；ctx 通常传入 HallMHDSolver*（也就是 this）
public:
    void set_rhshall_callback(RHShallCallback cb, void *ctx)
    {
        rhshall_cb_ = cb;
        rhshall_ctx_ = ctx;
    }
    //=========================================================================
};