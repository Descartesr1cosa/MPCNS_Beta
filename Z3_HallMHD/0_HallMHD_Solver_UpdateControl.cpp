#include "HallMHD_Solver.h"

bool HallMHDSolver::UpdateControlAndOutput()
{
    // 更新控制
    control_.Update();
    if (control_.if_outfile)
        output_.output_field();
    // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
    // output_.output_plt_field(); // output_.output_plt_cell_field(output_.var_defaut_plt_name); // output_field();
    if (control_.if_stop)
    {
        // output_.output_plt_cell_field(output_.var_defaut_plt_name); //For debug
        output_.output_field();
        return true;
    }
    return false;
}