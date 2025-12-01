#include "Euler_Solver.h"
void EulerSolver::calc_Residual()
{
    int ncomp = fld_->field("U_", 0).descriptor().ncomp;
    double temp;
    for (int32_t i = 0; i < control_.residual.Getsize1(); i++)
        control_.residual(i) = 0.0;
    for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
    {
        auto &U_calc = fld_->field("U_", iblock);
        auto &field_U = fld_->field("old_U_", iblock);

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
                                      << par_->GetInt("myid") << "\tn, i, j, k=\t" << iblock << "\t" << i << "\t" << j << "\t" << k << "\n";
                            exit(-1);
                        }
                        temp = U_calc(i, j, k, l) - field_U(i, j, k, l);
                        control_.residual(l) += (temp * temp);
                    }
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
    }
    for (int32_t l = 0; l < control_.residual.Getsize1(); l++)
        control_.residual(l) /= control_.residual_reference(l);
}