#include <cmath>
#include <cstdio>

// Z3_HallMHD
#include "HallMHD_Solver.h"
#include "4_solver/ImplicitHall_Solver.h"
//=========================================================================
//-------------------------------------------------------------------------
//=========================================================================
void HallMHDSolver::Advance()
{
    // 0. Preparation
    PrepareStep();

    while (true)
    {
        if (StepOnce())
            break;
    }
}
//=========================================================================
//-------------------------------------------------------------------------
//=========================================================================
bool HallMHDSolver::StepOnce()
{
    Compute_Timestep();
    Time_Advance(); // 内部只编排：ZeroRHS/Assemble/Update

#if HALL_MODE == 2
    PrepareSubstep_NoSnapshot(); // 把 * 步状态同步到一致（但不要 SnapshotOldFields）
    hall_->solve_implicit_hall(dt);
    Modify_TotalEnergy_AfterHall();
    Calc_Residual();
    SnapshotOldFields();
#else
    Calc_Residual();
    PrepareStep(); // 用于下一步/输出前字段一致
#endif

    return UpdateControlAndOutput();
}

void HallMHDSolver::Compute_Timestep()
{
    if (control_.Time_step > 0.0)
    {
        dt = control_.Time_step;
    }
    else
    {
        // CFL sanity check
        if (!(control_.CFL > 0.0 && control_.CFL < 1.0))
        {
            std::fprintf(stderr,
                         "[Timestep] CFL invalid for Explicit Euler: CFL=%g Time_step=%g\n",
                         control_.CFL, control_.Time_step);
            std::abort();
        }

        double radius_max = ComputeLocalMaxRadius_();

        // 全局取 max 再算 dt
        double radius_global = radius_max;
        PARALLEL::mpi_max(&radius_max, &radius_global, 1);

        if (radius_global <= 0.0 || std::isnan(radius_global))
        {
            std::fprintf(stderr, "[Timestep] radius_global <= 0 or is NAN, cannot compute dt.\n");
            std::abort();
        }

        dt = control_.CFL / radius_global;
    }

    control_.Physic_Time_Step = dt;
    fld_->par->AddParam("Physic_Time_Step", dt);
}

bool HallMHDSolver::UpdateControlAndOutput()
{
    // 更新控制
    control_.Update();

    if (control_.Nstep % control_.Resouput_step == 0)
        PrintMinMaxDiagnostics_();

    if (control_.if_outfile)
        output_.output_field();

    if (control_.if_stop)
    {
        // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
        output_.output_field();
        return true;
    }
    return false;
}

void HallMHDSolver::Calc_Residual()
{
    int ncomp = fld_->field(fid_.fid_U, 0).descriptor().ncomp;
    double temp;
    for (int32_t i = 0; i < control_.residual.Getsize1(); i++)
        control_.residual(i) = 0.0;
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        auto &U_calc = fld_->field(fid_.fid_U, iblock);
        auto &field_U = fld_->field(fid_.fid_old_U, iblock);

        const Int3 &sub = U_calc.inner_lo();
        const Int3 &sup = U_calc.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    if (U_calc(i, j, k, 0) < 0.0)
                    {
                        std::cout << "Negative Density Appears ! ! !   Myid=\t"
                                  << par_->GetInt("myid") << "\tn, i, j, k=\t" << iblock << "\t" << i << "\t" << j << "\t" << k << "\n";
                        exit(-1);
                    }

                    for (int32_t l = 0; l < ncomp; l++)
                    {
                        if (std::isnan(U_calc(i, j, k, l)))
                        {
                            std::cout << "NAN Appear ! ! !   Myid=\t"
                                      << par_->GetInt("myid") << "\tn, i, j, k=\t" << iblock << "\t" << i << "\t" << j << "\t" << k << "For componet " << l << "\n";
                            exit(-1);
                        }
                        temp = U_calc(i, j, k, l) - field_U(i, j, k, l);
                        control_.residual(l) += (temp * temp);
                    }
                }
    }

    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        auto &U_calc = fld_->field(fid_.fid_Bface.xi, iblock);
        auto &field_U = fld_->field(fid_.fid_old_Bface.xi, iblock);

        const Int3 &sub = U_calc.inner_lo();
        const Int3 &sup = U_calc.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    if (std::isnan(U_calc(i, j, k, 0)))
                    {
                        std::cout << "NAN Appear ! ! !   Myid=\t"
                                  << par_->GetInt("myid") << "\tn, i, j, k=\t" << iblock << "\t" << i << "\t" << j << "\t" << k << "For componet B_xi" << "\n";
                        exit(-1);
                    }
                    temp = U_calc(i, j, k, 0) - field_U(i, j, k, 0);
                    control_.residual(5) += (temp * temp);
                }
    }
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        auto &U_calc = fld_->field(fid_.fid_Bface.eta, iblock);
        auto &field_U = fld_->field(fid_.fid_old_Bface.eta, iblock);

        const Int3 &sub = U_calc.inner_lo();
        const Int3 &sup = U_calc.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    if (std::isnan(U_calc(i, j, k, 0)))
                    {
                        std::cout << "NAN Appear ! ! !   Myid=\t"
                                  << par_->GetInt("myid") << "\tn, i, j, k=\t" << iblock << "\t" << i << "\t" << j << "\t" << k << "For componet B_eta" << "\n";
                        exit(-1);
                    }
                    temp = U_calc(i, j, k, 0) - field_U(i, j, k, 0);
                    control_.residual(6) += (temp * temp);
                }
    }
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        auto &U_calc = fld_->field(fid_.fid_Bface.zeta, iblock);
        auto &field_U = fld_->field(fid_.fid_old_Bface.zeta, iblock);

        const Int3 &sub = U_calc.inner_lo();
        const Int3 &sup = U_calc.inner_hi();
        for (int i = sub.i; i < sup.i; i++)
            for (int j = sub.j; j < sup.j; j++)
                for (int k = sub.k; k < sup.k; k++)
                {
                    if (std::isnan(U_calc(i, j, k, 0)))
                    {
                        std::cout << "NAN Appear ! ! !   Myid=\t"
                                  << par_->GetInt("myid") << "\tn, i, j, k=\t" << iblock << "\t" << i << "\t" << j << "\t" << k << "For componet B_zeta" << "\n";
                        exit(-1);
                    }
                    temp = U_calc(i, j, k, 0) - field_U(i, j, k, 0);
                    control_.residual(7) += (temp * temp);
                }
    }

    // 所有进程残差规约操作
    double *send = new double[control_.residual.Getsize1()];
    double *recv = new double[control_.residual.Getsize1()];
    for (int32_t l = 0; l < control_.residual.Getsize1(); l++)
        send[l] = control_.residual(l);
    PARALLEL::mpi_sum(send, recv, control_.residual.Getsize1());
    for (int32_t l = 0; l < control_.residual.Getsize1(); l++)
        control_.residual(l) = sqrt(recv[l]);
    delete[] send;
    delete[] recv;
    // par->AddParam("continue_calc_index", true);
    if (control_.Nstep == 0 || par_->GetBoo("continue_calc") && par_->GetBoo("continue_calc_index"))
    {
        par_->AddParam("continue_calc_index", false);
        for (int32_t l = 0; l < control_.residual.Getsize1(); l++)
        {
            control_.residual_reference(l) = control_.residual(l);
            if (control_.residual_reference(l) < 1e-20)
            {
                control_.residual_reference(l) = 1e-20;
            }
        }
        par_->AddParam("Res_Ref_0", control_.residual_reference(0));
        par_->AddParam("Res_Ref_1", control_.residual_reference(1));
        par_->AddParam("Res_Ref_2", control_.residual_reference(2));
        par_->AddParam("Res_Ref_3", control_.residual_reference(3));
        par_->AddParam("Res_Ref_4", control_.residual_reference(4));
        par_->AddParam("Res_Ref_5", control_.residual_reference(5));
        par_->AddParam("Res_Ref_6", control_.residual_reference(6));
        par_->AddParam("Res_Ref_7", control_.residual_reference(7));
    }
    for (int32_t l = 0; l < control_.residual.Getsize1(); l++)
        control_.residual(l) /= control_.residual_reference(l);
}