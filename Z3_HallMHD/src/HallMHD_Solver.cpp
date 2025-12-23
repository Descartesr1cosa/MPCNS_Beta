
// Core
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"

// Z3_HallMHD
#include "HallMHD_Solver.h"
#include "4_solver/ImplicitHall_Solver.h"

namespace
{
    constexpr HallMode CompiledHallMode()
    {
#if HALL_MODE == 0
        return HallMode::Ideal;
#elif HALL_MODE == 1
        return HallMode::Explicit;
#else
        return HallMode::Implicit;
#endif
    }
}

HallMHDSolver::HallMHDSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo,
                             Param *par, ImplicitHall_Solver *hall_imp)
    : grd_(grd),
      topo_(topo),
      fld_(fld),
      halo_(halo),
      par_(par),
      hall_(hall_imp)
{
    // ---- Cache field ids (compiled mode) ----
    fid_.Init(fld_, CompiledHallMode());

    // ---- Calc Constants ----
    calc_physical_constant(par);

    // ---- components ----
    control_.SetUp(par_, 8);
    output_.SetUp(par_, fld_);
    bound_.SetUp(grd_, fld_, topo_, par_);

    // ---- Initialization ----
    initial_.Initialization(fld_);
    ComputeBcellInner();
    SyncDerivedBcell();
    add_Emag_to_Etotal();

    // ---- hall initial ----
#if HALL_MODE == 2
    if (hall_)
    {
        hall_->setup(grd_, topo_, fld_, halo_, par_, &bound_);
        hall_->set_rhshall_callback(&HallMHDSolver::RHShallCallback_, this);
    }
#endif
}
