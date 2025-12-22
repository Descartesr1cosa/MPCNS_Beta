#include "4_solver/ImplicitHall_Solver.h"

//=========================================================================
//-------------------------------------------------------------------------
// 开辟owner空间，初始化为1
void ImplicitHall_Solver::alloc_owner_by_B_inner_and_fill1_()
{
    owner_bxi_.resize(grd_->nblock);
    owner_beta_.resize(grd_->nblock);
    owner_bzeta_.resize(grd_->nblock);

    for (int ib = 0; ib < grd_->nblock; ++ib)
    {
        auto &Bxi = fld_->field(fid_Bxi, ib);
        auto &Beta = fld_->field(fid_Beta, ib);
        auto &Bze = fld_->field(fid_Bzeta, ib);

        auto lo = Bxi.inner_lo();
        auto hi = Bxi.inner_hi();
        owner_bxi_[ib].SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, 0);
        owner_bxi_[ib] = (PetscInt)1;

        lo = Beta.inner_lo();
        hi = Beta.inner_hi();
        owner_beta_[ib].SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, 0);
        owner_beta_[ib] = (PetscInt)1;

        lo = Bze.inner_lo();
        hi = Bze.inner_hi();
        owner_bzeta_[ib].SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, 0);
        owner_bzeta_[ib] = (PetscInt)1;
    }
}

// 将物理边界的face都设置为0
void ImplicitHall_Solver::apply_physical_bc_mask_()
{
    // 第一版最稳策略：物理边界上的 B_face 都不作为未知量（=0），每次残差评估由 BC 覆盖它们
    for (auto &p : topo_->physical_patches)
    {
        const int axis = face_axis_from_node_box_(p.this_box_node);
        if (axis < 0)
            continue;

        // topo_里面的范围是基于node的，这里会自动处理为对应的范围
        mark_face_patch_(p.this_block, p.this_box_node, axis, (PetscInt)1);
    }
}

// 将inner parallel边界的face按照优先级(rank block i j k)取较小的保留为1，否则为被动量并设置为-1
void ImplicitHall_Solver::apply_interface_nonowner_mask_()
{
    // Inner interface: non-owner => -1
    for (auto &p : topo_->inner_patches)
    {
        const int axis = face_axis_from_node_box_(p.this_box_node);
        if (axis < 0)
            continue;
        if (!this_side_is_owner_(p, /*is_parallel=*/false))
            mark_face_patch_(p.this_block, p.this_box_node, axis, (PetscInt)-1);
    }

    // Parallel interface: non-owner => -1
    for (auto &p : topo_->parallel_patches)
    {
        const int axis = face_axis_from_node_box_(p.this_box_node);
        if (axis < 0)
            continue;
        if (!this_side_is_owner_(p, /*is_parallel=*/true))
            mark_face_patch_(p.this_block, p.this_box_node, axis, (PetscInt)-1);
    }
}

//-------------------------------------------------------------------------
// tool find out which face the Box3 is
int ImplicitHall_Solver::face_axis_from_node_box_(const Box3 &b)
{
    const int di = b.hi.i - b.lo.i;
    const int dj = b.hi.j - b.lo.j;
    const int dk = b.hi.k - b.lo.k;

    if (di == 1 && dj > 1 && dk > 1)
        return 0; // Xi-normal face
    if (dj == 1 && di > 1 && dk > 1)
        return 1; // Eta-normal face
    if (dk == 1 && di > 1 && dj > 1)
        return 2; // Zeta-normal face
    return -1;
}

// tool mark the whole face (range of which is Box3) as val
void ImplicitHall_Solver::mark_face_patch_(int ib, const Box3 &box, int axis, PetscInt val)
{
    if (axis == 0)
    { // Xi-face: i fixed, (j,k) run over faces => node-range minus 1
        const int i_face = box.lo.i;
        for (int j = box.lo.j; j < box.hi.j - 1; ++j)
            for (int k = box.lo.k; k < box.hi.k - 1; ++k)
                if (owner_bxi_[ib](i_face, j, k) != 0)
                    owner_bxi_[ib](i_face, j, k) = val;
    }
    else if (axis == 1)
    { // Eta-face: j fixed
        const int j_face = box.lo.j;
        for (int i = box.lo.i; i < box.hi.i - 1; ++i)
            for (int k = box.lo.k; k < box.hi.k - 1; ++k)
                if (owner_beta_[ib](i, j_face, k) != 0)
                    owner_beta_[ib](i, j_face, k) = val;
    }
    else if (axis == 2)
    { // Zeta-face: k fixed
        const int k_face = box.lo.k;
        for (int i = box.lo.i; i < box.hi.i - 1; ++i)
            for (int j = box.lo.j; j < box.hi.j - 1; ++j)
                if (owner_bzeta_[ib](i, j, k_face) != 0)
                    owner_bzeta_[ib](i, j, k_face) = val;
    }
}

// tool 构造一个patch的PatchKey
ImplicitHall_Solver::PatchKey ImplicitHall_Solver::make_key_(int rank, int block, const Box3 &b)
{
    const int axis = face_axis_from_node_box_(b);

    PatchKey k{};
    k.rank = rank;
    k.block = block;
    k.axis = axis;

    if (axis == 0)
    { // Xi face: i fixed, (j,k) vary
        k.fixed = b.lo.i;
        k.lo_u = b.lo.j;
        k.hi_u = b.hi.j;
        k.lo_v = b.lo.k;
        k.hi_v = b.hi.k;
    }
    else if (axis == 1)
    { // Eta face: j fixed
        k.fixed = b.lo.j;
        k.lo_u = b.lo.i;
        k.hi_u = b.hi.i;
        k.lo_v = b.lo.k;
        k.hi_v = b.hi.k;
    }
    else
    { // axis == 2, Zeta face: k fixed
        k.fixed = b.lo.k;
        k.lo_u = b.lo.i;
        k.hi_u = b.hi.i;
        k.lo_v = b.lo.j;
        k.hi_v = b.hi.j;
    }
    return k;
}

// tool 根据PatchKey判断自己是否是owner
bool ImplicitHall_Solver::this_side_is_owner_(const TOPO::InterfacePatch &p, bool is_parallel)
{
    // 1) parallel: rank tie-break
    if (is_parallel)
        return (p.this_rank < p.nb_rank);

    // 2) same rank: block tie-break
    if (p.this_block != p.nb_block)
        return (p.this_block < p.nb_block);

    // 3) same rank + same block: patch-key tie-break (periodic/self-connect/pole fold)
    const PatchKey kt = make_key_(p.this_rank, p.this_block, p.this_box_node);
    const PatchKey kn = make_key_(p.nb_rank, p.nb_block, p.nb_box_node);

    // 如果 key 完全相同（极少见：数据重复/退化），就让 this 作为 owner，避免两边都标 non-owner
    if (!(kt < kn) && !(kn < kt))
        return true;

    return (kt < kn);
}
//=========================================================================

//=========================================================================
//-------------------------------------------------------------------------
// 给所有1的owner进行编号，获得全局唯一的DOF号码，从1开始
void ImplicitHall_Solver::assign_global_ids_to_owners_()
{
    PetscInt nloc = count_owned_local_();

    ndof_local_owned_ = nloc; // 记录本 rank owned DOF 数

    PetscInt prefix = 0;
    MPI_Exscan(&nloc, &prefix, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);
    if (myid == 0)
        prefix = 0;

    PetscInt gid = prefix + 1;

    start_gid_ = gid; // 记录本 rank 的 owned gid 起点（1-based）

    for (int ib = 0; ib < grd_->nblock; ++ib)
    {
        auto lo = fld_->field(fid_Bxi, ib).inner_lo();
        auto hi = fld_->field(fid_Bxi, ib).inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (owner_bxi_[ib](i, j, k) == 1)
                        owner_bxi_[ib](i, j, k) = gid++;

        lo = fld_->field(fid_Beta, ib).inner_lo();
        hi = fld_->field(fid_Beta, ib).inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (owner_beta_[ib](i, j, k) == 1)
                        owner_beta_[ib](i, j, k) = gid++;

        lo = fld_->field(fid_Bzeta, ib).inner_lo();
        hi = fld_->field(fid_Bzeta, ib).inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (owner_bzeta_[ib](i, j, k) == 1)
                        owner_bzeta_[ib](i, j, k) = gid++;
    }

    // 保存全局 dof 数
    PetscInt nglob = 0;
    MPI_Allreduce(&nloc, &nglob, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);
    ndof_global_ = nglob; // 记录全局的总DOF数量
}

//-------------------------------------------------------------------------
// 统计本进程的总DOF数量
PetscInt ImplicitHall_Solver::count_owned_local_()
{
    PetscInt c = 0;
    for (int ib = 0; ib < grd_->nblock; ++ib)
    {
        auto lo = fld_->field(fid_Bxi, ib).inner_lo();
        auto hi = fld_->field(fid_Bxi, ib).inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (owner_bxi_[ib](i, j, k) == 1)
                        ++c;

        lo = fld_->field(fid_Beta, ib).inner_lo();
        hi = fld_->field(fid_Beta, ib).inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (owner_beta_[ib](i, j, k) == 1)
                        ++c;

        lo = fld_->field(fid_Bzeta, ib).inner_lo();
        hi = fld_->field(fid_Bzeta, ib).inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (owner_bzeta_[ib](i, j, k) == 1)
                        ++c;
    }
    return c;
}
//=========================================================================

//=========================================================================
//-------------------------------------------------------------------------
// 对owner = -1的修改为真实DOF的负值
void ImplicitHall_Solver::fill_passiveowner_inner_()
{
    for (const auto &p : topo_->inner_patches)
    {
        const int axis_this = face_axis_from_node_box_(p.this_box_node);
        if (axis_this < 0)
            continue;

        const int axis_nb = p.trans.perm[axis_this]; // 关键：this 的法向轴映射到 nb 的法向轴
        const int ibA = p.this_block;
        const int ibB = p.nb_block;

        auto &A = owner_arr_(axis_this, ibA);
        auto &B = owner_arr_(axis_nb, ibB);

        const Box3 &boxA = p.this_box_node;

        if (axis_this == 0)
        {
            const int i_face = boxA.lo.i;
            for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                {
                    // this face dof anchor node
                    const Int3 nbnode = apply_trans_node_(p.trans, i_face, j, k);

                    int iB = nbnode.i, jB = nbnode.j, kB = nbnode.k;
                    // nb face indices depend on axis_nb; 但用 (iB,jB,kB) 作为 face 的索引三元组是成立的（固定轴那一维就是常数）
                    PetscInt va = A(i_face, j, k);
                    PetscInt vb = B(iB, jB, kB);

                    if (va > 0 && vb == -1)
                        B(iB, jB, kB) = -va;
                    else if (vb > 0 && va == -1)
                        A(i_face, j, k) = -vb;
                    else if (va > 0 && vb < 0)
                    { /* 已填过：可选检查 vb==-va */
                    }
                    else if (vb > 0 && va < 0)
                    { /* 已填过 */
                    }
                }
        }
        else if (axis_this == 1)
        {
            const int j_face = boxA.lo.j;
            for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                {
                    const Int3 nbnode = apply_trans_node_(p.trans, i, j_face, k);
                    int iB = nbnode.i, jB = nbnode.j, kB = nbnode.k;

                    PetscInt va = A(i, j_face, k);
                    PetscInt vb = B(iB, jB, kB);

                    if (va > 0 && vb == -1)
                        B(iB, jB, kB) = -va;
                    else if (vb > 0 && va == -1)
                        A(i, j_face, k) = -vb;
                }
        }
        else // axis_this == 2
        {
            const int k_face = boxA.lo.k;
            for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                {
                    const Int3 nbnode = apply_trans_node_(p.trans, i, j, k_face);
                    int iB = nbnode.i, jB = nbnode.j, kB = nbnode.k;

                    PetscInt va = A(i, j, k_face);
                    PetscInt vb = B(iB, jB, kB);

                    if (va > 0 && vb == -1)
                        B(iB, jB, kB) = -va;
                    else if (vb > 0 && va == -1)
                        A(i, j, k_face) = -vb;
                }
        }
    }
}

// 对owner = -1的修改为真实DOF的负值
void ImplicitHall_Solver::fill_passiveowner_parallel_()
{
    struct RecvTask
    {
        const TOPO::InterfacePatch *p;
        int axis_this;
        int ib;
        std::vector<PetscInt> buf;
        MPI_Request req;
    };
    struct SendTask
    {
        std::vector<PetscInt> buf;
        MPI_Request req;
    };

    std::vector<RecvTask> recvs;
    std::vector<SendTask> sends;

    // 1) 先把所有 Irecv post 掉（避免死锁）
    for (const auto &p : topo_->parallel_patches)
    {
        const int axis_this = face_axis_from_node_box_(p.this_box_node);
        if (axis_this < 0)
            continue;

        auto &A = owner_arr_(axis_this, p.this_block);

        // 用 patch 上第一个 face dof 判定角色
        PetscInt v0 = 0;
        if (axis_this == 0)
            v0 = A(p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k);
        if (axis_this == 1)
            v0 = A(p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k);
        if (axis_this == 2)
            v0 = A(p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k);

        if (v0 == -1)
        {
            RecvTask t;
            t.p = &p;
            t.axis_this = axis_this;
            t.ib = p.this_block;

            const int nface = face_count_from_node_box_(p.this_box_node, axis_this);
            t.buf.resize(nface, (PetscInt)0);

            MPI_Irecv(t.buf.data(), nface, MPIU_INT, p.nb_rank, p.recv_flag, PETSC_COMM_WORLD, &t.req);
            recvs.push_back(std::move(t));
        }
    }

    // 2) 再做所有 Isend
    for (const auto &p : topo_->parallel_patches)
    {
        const int axis_this = face_axis_from_node_box_(p.this_box_node);
        if (axis_this < 0)
            continue;

        auto &A = owner_arr_(axis_this, p.this_block);

        PetscInt v0 = 0;
        if (axis_this == 0)
            v0 = A(p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k);
        if (axis_this == 1)
            v0 = A(p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k);
        if (axis_this == 2)
            v0 = A(p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k);

        if (v0 > 0)
        {
            const int axis_nb = p.trans.perm[axis_this];

            Box3 nb_box = p.nb_box_node;

            const int nface = face_count_from_node_box_(nb_box, axis_nb);

            SendTask s;
            s.buf.resize(nface, (PetscInt)0);

            const Box3 &boxA = p.this_box_node;

            if (axis_this == 0)
            {
                const int i_face = boxA.lo.i;
                for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                    for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                    {
                        const PetscInt gid = A(i_face, j, k);
                        const Int3 nbnode = apply_trans_node_(p.trans, i_face, j, k);
                        const int idx = face_linear_index_in_nodebox_(nb_box, axis_nb, nbnode.i, nbnode.j, nbnode.k);
                        s.buf[idx] = gid;
                    }
            }
            else if (axis_this == 1)
            {
                const int j_face = boxA.lo.j;
                for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                    for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                    {
                        const PetscInt gid = A(i, j_face, k);
                        const Int3 nbnode = apply_trans_node_(p.trans, i, j_face, k);
                        const int idx = face_linear_index_in_nodebox_(nb_box, axis_nb, nbnode.i, nbnode.j, nbnode.k);
                        s.buf[idx] = gid;
                    }
            }
            else
            {
                const int k_face = boxA.lo.k;
                for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                    for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                    {
                        const PetscInt gid = A(i, j, k_face);
                        const Int3 nbnode = apply_trans_node_(p.trans, i, j, k_face);
                        const int idx = face_linear_index_in_nodebox_(nb_box, axis_nb, nbnode.i, nbnode.j, nbnode.k);
                        s.buf[idx] = gid;
                    }
            }

            MPI_Isend(s.buf.data(), nface, MPIU_INT, p.nb_rank, p.send_flag, PETSC_COMM_WORLD, &s.req);
            sends.push_back(std::move(s));
        }
    }

    // 3) 等接收完成并写回（写成负值）
    for (auto &t : recvs)
    {
        MPI_Wait(&t.req, MPI_STATUS_IGNORE);

        const auto &p = *t.p;
        auto &A = owner_arr_(t.axis_this, t.ib);

        int idx = 0;
        const Box3 &box = p.this_box_node;

        if (t.axis_this == 0)
        {
            const int i_face = box.lo.i;
            for (int j = box.lo.j; j < box.hi.j - 1; ++j)
                for (int k = box.lo.k; k < box.hi.k - 1; ++k, ++idx)
                    if (A(i_face, j, k) == -1)
                        A(i_face, j, k) = -t.buf[idx];
        }
        else if (t.axis_this == 1)
        {
            const int j_face = box.lo.j;
            for (int i = box.lo.i; i < box.hi.i - 1; ++i)
                for (int k = box.lo.k; k < box.hi.k - 1; ++k, ++idx)
                    if (A(i, j_face, k) == -1)
                        A(i, j_face, k) = -t.buf[idx];
        }
        else
        {
            const int k_face = box.lo.k;
            for (int i = box.lo.i; i < box.hi.i - 1; ++i)
                for (int j = box.lo.j; j < box.hi.j - 1; ++j, ++idx)
                    if (A(i, j, k_face) == -1)
                        A(i, j, k_face) = -t.buf[idx];
        }
    }

    // 4) 等发送完成（保证 buffer 生命周期）
    for (auto &s : sends)
        MPI_Wait(&s.req, MPI_STATUS_IGNORE);
}

//-------------------------------------------------------------------------
Int3 ImplicitHall_Solver::apply_trans_node_(const TOPO::IndexTransform &T, int i, int j, int k)
{
    const int loc[3] = {i, j, k};
    const int off[3] = {T.offset.i, T.offset.j, T.offset.k};

    int nb[3] = {0, 0, 0};
    for (int a = 0; a < 3; ++a)
    {
        nb[T.perm[a]] = T.sign[a] * loc[a] + off[a];
    }
    return Int3{nb[0], nb[1], nb[2]};
}

int ImplicitHall_Solver::face_count_from_node_box_(const Box3 &b, int axis)
{
    if (axis == 0)
    {
        const int nj = b.hi.j - b.lo.j;
        const int nk = b.hi.k - b.lo.k;
        return (nj - 1) * (nk - 1);
    }
    else if (axis == 1)
    {
        const int ni = b.hi.i - b.lo.i;
        const int nk = b.hi.k - b.lo.k;
        return (ni - 1) * (nk - 1);
    }
    else
    {
        const int ni = b.hi.i - b.lo.i;
        const int nj = b.hi.j - b.lo.j;
        return (ni - 1) * (nj - 1);
    }
}

int ImplicitHall_Solver::face_linear_index_in_nodebox_(const Box3 &b, int axis, int i, int j, int k)
{
    // 线性下标按“接收侧 natural loop 顺序”定义：
    // axis=0: j outer, k inner
    // axis=1: i outer, k inner
    // axis=2: i outer, j inner
    if (axis == 0)
    {
        const int njf = (b.hi.j - b.lo.j - 1);
        const int nkf = (b.hi.k - b.lo.k - 1);
        return (j - b.lo.j) * nkf + (k - b.lo.k);
    }
    else if (axis == 1)
    {
        const int nif = (b.hi.i - b.lo.i - 1);
        const int nkf = (b.hi.k - b.lo.k - 1);
        return (i - b.lo.i) * nkf + (k - b.lo.k);
    }
    else
    {
        const int nif = (b.hi.i - b.lo.i - 1);
        const int njf = (b.hi.j - b.lo.j - 1);
        return (i - b.lo.i) * njf + (j - b.lo.j);
    }
}
//=========================================================================

void ImplicitHall_Solver::InterfaceExchange_Bface_Inner()
{
    for (const auto &p : topo_->inner_patches)
    {
        const int axis_this = face_axis_from_node_box_(p.this_box_node);
        if (axis_this < 0)
            continue;

        const int axis_nb = p.trans.perm[axis_this];
        const int ibA = p.this_block;
        const int ibB = p.nb_block;

        auto &ownA = owner_arr_(axis_this, ibA);
        auto &ownB = owner_arr_(axis_nb, ibB);

        auto &BA = B_face_(axis_this, ibA);
        auto &BB = B_face_(axis_nb, ibB);

        const Box3 &boxA = p.this_box_node;

        // 用一个代表点判定方向：本 patch 必须是一边 owner 一边 non-owner
        PetscInt sA = 0;
        if (axis_this == 0)
            sA = ownA(boxA.lo.i, boxA.lo.j, boxA.lo.k);
        if (axis_this == 1)
            sA = ownA(boxA.lo.i, boxA.lo.j, boxA.lo.k);
        if (axis_this == 2)
            sA = ownA(boxA.lo.i, boxA.lo.j, boxA.lo.k);
        if (sA == 0)
            continue;

        // 遍历 this-side face dof，用 trans 找到 nb-side 对应 face dof
        if (axis_this == 0)
        {
            const int i_face = boxA.lo.i;
            for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                {
                    const Int3 nb = apply_trans_node_(p.trans, i_face, j, k);

                    const PetscInt oa = ownA(i_face, j, k);
                    const PetscInt ob = ownB(nb.i, nb.j, nb.k);

                    if (oa > 0 && ob < 0)
                        BB(nb.i, nb.j, nb.k, 0) = BA(i_face, j, k, 0);
                    else if (ob > 0 && oa < 0)
                        BA(i_face, j, k, 0) = BB(nb.i, nb.j, nb.k, 0);
                }
        }
        else if (axis_this == 1)
        {
            const int j_face = boxA.lo.j;
            for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                {
                    const Int3 nb = apply_trans_node_(p.trans, i, j_face, k);

                    const PetscInt oa = ownA(i, j_face, k);
                    const PetscInt ob = ownB(nb.i, nb.j, nb.k);

                    if (oa > 0 && ob < 0)
                        BB(nb.i, nb.j, nb.k, 0) = BA(i, j_face, k, 0);
                    else if (ob > 0 && oa < 0)
                        BA(i, j_face, k, 0) = BB(nb.i, nb.j, nb.k, 0);
                }
        }
        else
        {
            const int k_face = boxA.lo.k;
            for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                {
                    const Int3 nb = apply_trans_node_(p.trans, i, j, k_face);

                    const PetscInt oa = ownA(i, j, k_face);
                    const PetscInt ob = ownB(nb.i, nb.j, nb.k);

                    if (oa > 0 && ob < 0)
                        BB(nb.i, nb.j, nb.k, 0) = BA(i, j, k_face, 0);
                    else if (ob > 0 && oa < 0)
                        BA(i, j, k_face, 0) = BB(nb.i, nb.j, nb.k, 0);
                }
        }
    }
}

void ImplicitHall_Solver::InterfaceExchange_Bface_Parallel()
{
    struct RecvTask
    {
        const TOPO::InterfacePatch *p;
        int axis_this;
        int axis_nb; // 仅用于一致性
        int ib;
        Box3 box_this;
        std::vector<double> buf;
        MPI_Request req;
    };
    struct SendTask
    {
        std::vector<double> buf;
        MPI_Request req;
    };

    std::vector<RecvTask> recvs;
    std::vector<SendTask> sends;

    // 先 post Irecv
    for (const auto &p : topo_->parallel_patches)
    {
        const int axis_this = face_axis_from_node_box_(p.this_box_node);
        if (axis_this < 0)
            continue;

        auto &own = owner_arr_(axis_this, p.this_block);
        const Box3 &box = p.this_box_node;

        PetscInt v0 = own(box.lo.i, box.lo.j, box.lo.k);
        if (v0 == 0)
            continue;

        if (v0 < 0)
        {
            RecvTask t;
            t.p = &p;
            t.axis_this = axis_this;
            t.axis_nb = p.trans.perm[axis_this];
            t.ib = p.this_block;
            t.box_this = box;

            const int nface = face_count_from_node_box_(box, axis_this);
            t.buf.resize(nface, 0.0);

            MPI_Irecv(t.buf.data(), nface, MPI_DOUBLE, p.nb_rank, p.recv_flag, PETSC_COMM_WORLD, &t.req);
            recvs.push_back(std::move(t));
        }
    }

    // 再做 Isend
    for (const auto &p : topo_->parallel_patches)
    {
        const int axis_this = face_axis_from_node_box_(p.this_box_node);
        if (axis_this < 0)
            continue;

        auto &ownA = owner_arr_(axis_this, p.this_block);
        const Box3 &boxA = p.this_box_node;

        PetscInt v0 = ownA(boxA.lo.i, boxA.lo.j, boxA.lo.k);
        if (v0 == 0)
            continue;

        if (v0 > 0)
        {
            const int axis_nb = p.trans.perm[axis_this];
            const Box3 nb_box = p.nb_box_node; // transform_box_node_(p.this_box_node, p.trans);

            const int nface = face_count_from_node_box_(nb_box, axis_nb);

            SendTask s;
            s.buf.assign(nface, 0.0);

            auto &BA = B_face_(axis_this, p.this_block);

            if (axis_this == 0)
            {
                const int i_face = boxA.lo.i;
                for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                    for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                    {
                        const Int3 nb = apply_trans_node_(p.trans, i_face, j, k);
                        const int idx = face_linear_index_in_nodebox_(nb_box, axis_nb, nb.i, nb.j, nb.k);
                        s.buf[idx] = BA(i_face, j, k, 0);
                    }
            }
            else if (axis_this == 1)
            {
                const int j_face = boxA.lo.j;
                for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                    for (int k = boxA.lo.k; k < boxA.hi.k - 1; ++k)
                    {
                        const Int3 nb = apply_trans_node_(p.trans, i, j_face, k);
                        const int idx = face_linear_index_in_nodebox_(nb_box, axis_nb, nb.i, nb.j, nb.k);
                        s.buf[idx] = BA(i, j_face, k, 0);
                    }
            }
            else
            {
                const int k_face = boxA.lo.k;
                for (int i = boxA.lo.i; i < boxA.hi.i - 1; ++i)
                    for (int j = boxA.lo.j; j < boxA.hi.j - 1; ++j)
                    {
                        const Int3 nb = apply_trans_node_(p.trans, i, j, k_face);
                        const int idx = face_linear_index_in_nodebox_(nb_box, axis_nb, nb.i, nb.j, nb.k);
                        s.buf[idx] = BA(i, j, k_face, 0);
                    }
            }

            MPI_Isend(s.buf.data(), nface, MPI_DOUBLE, p.nb_rank, p.send_flag, PETSC_COMM_WORLD, &s.req);
            sends.push_back(std::move(s));
        }
    }

    // 收到后写回 non-owner
    for (auto &t : recvs)
    {
        MPI_Wait(&t.req, MPI_STATUS_IGNORE);

        auto &own = owner_arr_(t.axis_this, t.ib);
        auto &B = B_face_(t.axis_this, t.ib);

        int idx = 0;
        const Box3 &box = t.box_this;

        if (t.axis_this == 0)
        {
            const int i_face = box.lo.i;
            for (int j = box.lo.j; j < box.hi.j - 1; ++j)
                for (int k = box.lo.k; k < box.hi.k - 1; ++k, ++idx)
                    if (own(i_face, j, k) < 0)
                        B(i_face, j, k, 0) = t.buf[idx];
        }
        else if (t.axis_this == 1)
        {
            const int j_face = box.lo.j;
            for (int i = box.lo.i; i < box.hi.i - 1; ++i)
                for (int k = box.lo.k; k < box.hi.k - 1; ++k, ++idx)
                    if (own(i, j_face, k) < 0)
                        B(i, j_face, k, 0) = t.buf[idx];
        }
        else
        {
            const int k_face = box.lo.k;
            for (int i = box.lo.i; i < box.hi.i - 1; ++i)
                for (int j = box.lo.j; j < box.hi.j - 1; ++j, ++idx)
                    if (own(i, j, k_face) < 0)
                        B(i, j, k_face, 0) = t.buf[idx];
        }
    }

    for (auto &s : sends)
        MPI_Wait(&s.req, MPI_STATUS_IGNORE);
}

//=========================================================================

void ImplicitHall_Solver::build_local_dof_map_faces_owner_only_()
{

    // 你应当已在 assign_global_ids_to_owners_() 里算过这些
    // ndof_local_owned_, ndof_global_, start_gid_
    // 这里假设它们已赋值。如果没赋值，就在此函数里再 Exscan/Allreduce 一次也行。

    local_dofs_.assign((size_t)ndof_local_owned_, FaceDofRef{-1, -1, -1, -1, -1});

    for (int ib = 0; ib < grd_->nblock; ++ib)
    {
        // ---- B_xi ----
        {
            auto lo = fld_->field(fid_Bxi, ib).inner_lo();
            auto hi = fld_->field(fid_Bxi, ib).inner_hi();
            auto &own = owner_bxi_[ib];

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const PetscInt gid = own(i, j, k);
                        if (gid > 0)
                        {
                            const PetscInt lid = gid - start_gid_;
                            if (lid < 0 || lid >= ndof_local_owned_)
                            {
                                std::cout << "[ImplicitHall] Fatal: gid out of local range. gid=" << gid
                                          << " start_gid=" << start_gid_
                                          << " ndof_local=" << ndof_local_owned_ << "\n";
                                MPI_Abort(PETSC_COMM_WORLD, -911);
                            }
                            local_dofs_[(size_t)lid] = FaceDofRef{0, ib, i, j, k};
                        }
                    }
        }

        // ---- B_eta ----
        {
            auto lo = fld_->field(fid_Beta, ib).inner_lo();
            auto hi = fld_->field(fid_Beta, ib).inner_hi();
            auto &own = owner_beta_[ib];

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const PetscInt gid = own(i, j, k);
                        if (gid > 0)
                        {
                            const PetscInt lid = gid - start_gid_;
                            if (lid < 0 || lid >= ndof_local_owned_)
                            {
                                MPI_Abort(PETSC_COMM_WORLD, -912);
                            }
                            local_dofs_[(size_t)lid] = FaceDofRef{1, ib, i, j, k};
                        }
                    }
        }

        // ---- B_zeta ----
        {
            auto lo = fld_->field(fid_Bzeta, ib).inner_lo();
            auto hi = fld_->field(fid_Bzeta, ib).inner_hi();
            auto &own = owner_bzeta_[ib];

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const PetscInt gid = own(i, j, k);
                        if (gid > 0)
                        {
                            const PetscInt lid = gid - start_gid_;
                            if (lid < 0 || lid >= ndof_local_owned_)
                            {
                                MPI_Abort(PETSC_COMM_WORLD, -913);
                            }
                            local_dofs_[(size_t)lid] = FaceDofRef{2, ib, i, j, k};
                        }
                    }
        }
    }

    // 简单一致性检查：看有没有未填项
    for (PetscInt lid = 0; lid < ndof_local_owned_; ++lid)
    {
        if (local_dofs_[(size_t)lid].axis < 0)
        {
            std::cout << "[ImplicitHall] Fatal: local_dofs_ has hole at lid=" << lid << "\n";
            MPI_Abort(PETSC_COMM_WORLD, -914);
        }
    }
}

//=========================================================================
