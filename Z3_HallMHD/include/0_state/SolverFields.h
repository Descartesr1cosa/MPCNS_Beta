// 0_state/SolverFields.h
#pragma once

#include <array>
#include <cstdio>
#include <cstdlib>

// Core
#include "3_field/2_MPCNS_Field.h"

// Z3_HallMHD
#include "0_state/StateTypes.h" // 已有 Triplet(FieldT*)；这里我们再定义一个 IdTriplet
#include "0_state/HallConfig.h"

struct SolverFields
{
    Field *fld = nullptr;
    HallMode hall_mode = HallMode::Ideal;

    // ---- geometry ----
    int fid_Jac = -1;
    IdTriplet fid_metric; // (xi,eta,zeta) <- (JDxi,JDet,JDze)
    IdTriplet fid_pinvGT; // (pinvGT_xi, pinvGT_eta, pinvGT_zeta)  each ncomp=9
    IdTriplet fid_pinvAT; // (pinvAT_xi, pinvAT_eta, pinvAT_zeta)  each ncomp=9

    // ---- field ids ----
    int fid_U = -1;
    IdTriplet fid_Bface; // (B_xi, B_eta, B_zeta)

    // ---- auxiliary ----
    int fid_PV = -1;
    int fid_Bcell = -1;
    IdTriplet fid_Badd; // (Badd_xi,Badd_eta,Badd_zeta)

    // ---- fluid flux and face E ----
    std::array<int, 3> fid_F{{-1, -1, -1}};     // F_xi/F_eta/F_zeta
    std::array<int, 3> fid_Eface{{-1, -1, -1}}; // E_face_xi/E_face_eta/E_face_zeta

    // ---- edge fields for CT / Hall ----
    IdTriplet fid_Eedge; // (E_xi,E_eta,E_zeta)
    IdTriplet fid_Jedge; // (J_xi,J_eta,J_zeta)
    IdTriplet fid_Ehall; // (Ehall_xi,Ehall_eta,Ehall_zeta)

    // -------- buffers / time advance ----------
    int fid_old_U = -1;
    int fid_divB = -1;
    IdTriplet fid_old_Bface; // old_B_xi/eta/zeta

    int fid_RHS_U = -1;          // RHS (cell,5)
    IdTriplet fid_RHS_Bface;     // RHS_xi/eta/zeta (face,1)
    IdTriplet fid_RHShall_Bface; // RHShall_xi/eta/zeta

    void Init(Field *fld_in, HallMode mode)
    {
        fld = fld_in;
        hall_mode = mode;

        if (!fld)
        {
            std::fprintf(stderr, "[SolverFields] fld is null\n");
            std::abort();
        }

        // ---- geometry ----
        fid_Jac = fld->field_id("Jac");
        fid_metric.xi = fld->field_id("JDxi");
        fid_metric.eta = fld->field_id("JDet");
        fid_metric.zeta = fld->field_id("JDze");
        fid_pinvGT.xi = fld->field_id("pinvGT_xi");
        fid_pinvGT.eta = fld->field_id("pinvGT_eta");
        fid_pinvGT.zeta = fld->field_id("pinvGT_zeta");
        fid_pinvAT.xi = fld->field_id("pinvAT_xi");
        fid_pinvAT.eta = fld->field_id("pinvAT_eta");
        fid_pinvAT.zeta = fld->field_id("pinvAT_zeta");

        // ---- field ids ----
        fid_U = fld->field_id("U_");
        fid_Bface.xi = fld->field_id("B_xi");
        fid_Bface.eta = fld->field_id("B_eta");
        fid_Bface.zeta = fld->field_id("B_zeta");

        // ---- auxiliary ----
        fid_PV = fld->field_id("PV_");
        fid_Bcell = fld->field_id("B_cell");
        fid_Badd.xi = fld->field_id("Badd_xi");
        fid_Badd.eta = fld->field_id("Badd_eta");
        fid_Badd.zeta = fld->field_id("Badd_zeta");

        // ---- fluid flux and face E ----
        fid_F[0] = fld->field_id("F_xi");
        fid_F[1] = fld->field_id("F_eta");
        fid_F[2] = fld->field_id("F_zeta");
        fid_Eface[0] = fld->field_id("E_face_xi");
        fid_Eface[1] = fld->field_id("E_face_eta");
        fid_Eface[2] = fld->field_id("E_face_zeta");

        // ---- edge fields for CT / Hall ----
        fid_Eedge.xi = fld->field_id("E_xi");
        fid_Eedge.eta = fld->field_id("E_eta");
        fid_Eedge.zeta = fld->field_id("E_zeta");

#if HALL_MODE != 0
        if (hall_mode != HallMode::Ideal)
        {
            fid_Jedge.xi = fld->field_id("J_xi");
            fid_Jedge.eta = fld->field_id("J_eta");
            fid_Jedge.zeta = fld->field_id("J_zeta");

            fid_Ehall.xi = fld->field_id("Ehall_xi");
            fid_Ehall.eta = fld->field_id("Ehall_eta");
            fid_Ehall.zeta = fld->field_id("Ehall_zeta");
        }
#endif

        // -------- buffers / time advance ----------
        fid_old_U = fld->field_id("old_U_");
        fid_divB = fld->field_id("divB");
        fid_old_Bface.xi = fld->field_id("old_B_xi");
        fid_old_Bface.eta = fld->field_id("old_B_eta");
        fid_old_Bface.zeta = fld->field_id("old_B_zeta");

        fid_RHS_U = fld->field_id("RHS");
        fid_RHS_Bface.xi = fld->field_id("RHS_xi");
        fid_RHS_Bface.eta = fld->field_id("RHS_eta");
        fid_RHS_Bface.zeta = fld->field_id("RHS_zeta");

        // ---- RHS hall only for implicit ----
#if HALL_MODE == 2
        fid_RHShall_Bface.xi = fld->field_id("RHShall_xi");
        fid_RHShall_Bface.eta = fld->field_id("RHShall_eta");
        fid_RHShall_Bface.zeta = fld->field_id("RHShall_zeta");
#endif
        // -------- Check --------
        Validate();
    }

    void Validate() const
    {
        if (!fld)
        {
            std::fprintf(stderr, "[SolverFields] fld is null\n");
            std::abort();
        }

        auto require_id = [](int id, const char *name)
        {
            if (id < 0)
            {
                std::fprintf(stderr, "[SolverFields] missing field id: %s\n", name ? name : "(null)");
                std::abort();
            }
        };

        // ---- geometry ----
        require_id(fid_Jac, "Jac");
        fid_metric.require_all("metric(JDxi/JDet/JDze)");
        fid_pinvGT.require_all("pinvGT(edge)");
        fid_pinvAT.require_all("pinvAT(edge)");

        // ---- primary / auxiliary ----
        require_id(fid_U, "U_");
        fid_Bface.require_all("B_face(B_xi/B_eta/B_zeta)");

        require_id(fid_PV, "PV_");
        require_id(fid_Bcell, "B_cell");
        fid_Badd.require_all("Badd(Badd_xi/Badd_eta/Badd_zeta)");

        // ---- fluid flux and face E ----
        for (int d = 0; d < 3; ++d)
        {
            if (fid_F[d] < 0)
            {
                std::fprintf(stderr, "[SolverFields] missing F_dir field id at d=%d\n", d);
                std::abort();
            }
            if (fid_Eface[d] < 0)
            {
                std::fprintf(stderr, "[SolverFields] missing E_face_dir field id at d=%d\n", d);
                std::abort();
            }
        }

        // ---- edge fields (CT / Hall) ----
        fid_Eedge.require_all("E_edge(E_xi/E_eta/E_zeta)");

#if HALL_MODE != 0
        if (hall_mode != HallMode::Ideal)
        {
            fid_Jedge.require_all("J_edge(J_xi/J_eta/J_zeta)");
            fid_Ehall.require_all("Ehall_edge(Ehall_xi/Ehall_eta/Ehall_zeta)");
        }
#else
        // 编译时不含 Hall，但运行时选了 Hall：必须报错
        if (hall_mode != HallMode::Ideal)
        {
            std::fprintf(stderr,
                         "[SolverFields] Hall requested at runtime (mode=%d) but compiled with HALL_MODE==0\n",
                         static_cast<int>(hall_mode));
            std::abort();
        }
#endif

        // ---- buffers / time advance ----
        require_id(fid_old_U, "old_U_");
        require_id(fid_divB, "divB");
        fid_old_Bface.require_all("old_B_face(old_B_xi/eta/zeta)");

        require_id(fid_RHS_U, "RHS");
        fid_RHS_Bface.require_all("RHS_B_face(RHS_xi/eta/zeta)");

        // ---- Implicit Hall only: compile-time gated by your registration ----
#if HALL_MODE == 2
        if (hall_mode == HallMode::Implicit)
        {
            fid_RHShall_Bface.require_all("RHShall_B_face(RHShall_xi/eta/zeta)");
        }
#else
        if (hall_mode == HallMode::Implicit)
        {
            std::fprintf(stderr,
                         "[SolverFields] Implicit Hall requested at runtime but compiled with HALL_MODE!=2\n");
            std::abort();
        }
#endif
    }
};