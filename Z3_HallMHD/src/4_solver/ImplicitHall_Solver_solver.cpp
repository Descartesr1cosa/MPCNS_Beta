#include "4_solver/ImplicitHall_Solver.h"

void ImplicitHall_Solver::create_petsc_objects_()
{
    if (petsc_objs_built_)
        return;

    // Vec sizes: local=ndof_local_owned_, global=ndof_global_
    PetscErrorCode ierr;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, ndof_local_owned_, ndof_global_, &X_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecDuplicate(X_, &F_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecDuplicate(X_, &Xref_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    ierr = SNESCreate(PETSC_COMM_WORLD, &snes_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = SNESSetFunction(snes_, F_, &ImplicitHall_Solver::FormFunction_PETSc, this);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // 先用不需要 Jacobian 的类型把壳跑通；后面再切 newtonls + mf + ksp
    ierr = SNESSetType(snes_, SNESQN);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // 允许命令行覆盖（例如后面你想切 -snes_type newtonls 等）
    ierr = SNESSetFromOptions(snes_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    petsc_objs_built_ = true;
}

void ImplicitHall_Solver::destroy_petsc_objects_()
{
    if (!petsc_objs_built_)
        return;

    if (snes_)
    {
        SNESDestroy(&snes_);
        snes_ = nullptr;
    }
    if (X_)
    {
        VecDestroy(&X_);
        X_ = nullptr;
    }
    if (F_)
    {
        VecDestroy(&F_);
        F_ = nullptr;
    }
    if (Xref_)
    {
        VecDestroy(&Xref_);
        Xref_ = nullptr;
    }

    petsc_objs_built_ = false;
}

PetscErrorCode ImplicitHall_Solver::FormFunction_(Vec X, Vec F)
{
    // // 这一步会反复被 SNES 调用：用你的现成链路准备 B_face 状态
    // // （owner-only unpack + interface exchange + BC + halo）
    // vec2B(X);
    // // dummy residual: F = X - Xref_
    // PetscErrorCode ierr;
    // ierr = VecCopy(X, F);
    // CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // ierr = VecAXPY(F, -1.0, Xref_);
    // CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // return 0;

    // 1) 把候选解 X 写回 B_face，并做 interface/BC/halo 一致化
    vec2B(X); // owner-only unpack + InterfaceExchange + BC + halo :contentReference[oaicite:4]{index=4}

    // 2) 计算 RHShall_* = RHS_Hall(B_face)   (这里你需要自己填：清零 + CT 链路)
    // RHS_Hall为curl E_Hall
    // if (!rhshall_cb_)
    //     SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "RHShall callback not set.");
    rhshall_cb_(rhshall_ctx_);

    // 3) 组装 VecF：F_lid = X_lid - Xref_lid + dt * RHShall(dof)
    const PetscScalar *x = nullptr, *xref = nullptr;
    PetscScalar *f = nullptr;
    VecGetArrayRead(X, &x);
    VecGetArrayRead(Xref_, &xref);
    VecGetArray(F, &f);

    for (PetscInt lid = 0; lid < ndof_local_owned_; ++lid)
    {
        const auto &d = local_dofs_[(size_t)lid];

        FieldBlock &rhs = (d.axis == 0   ? fld_->field(fid_RHShall_xi, d.ib)
                           : d.axis == 1 ? fld_->field(fid_RHShall_eta, d.ib)
                                         : fld_->field(fid_RHShall_zeta, d.ib));

        f[lid] = x[lid] - xref[lid] + (PetscScalar)(dt_ * rhs(d.i, d.j, d.k, 0));
    }

    VecRestoreArray(F, &f);
    VecRestoreArrayRead(Xref_, &xref);
    VecRestoreArrayRead(X, &x);
    return 0;
}

void ImplicitHall_Solver::solve_implicit_hall(double dt)
{
    dt_ = dt;
    // 初值：从当前 B_face pack 到 X_
    pack_Bface_owner_only_to_vec_(X_);

    // 参考解：令 Xref_ = X_，这样 dummy residual 立刻为 0，SNES 会正常退出
    VecCopy(X_, Xref_);

    // 运行 SNES
    SNESSolve(snes_, nullptr, X_);

    // 写回（也可不写回；这里写回方便你验证整套链）
    vec2B(X_);
}