#pragma once
#include "3_field/2_MPCNS_Field.h"

class MHD_Boundary
{
public:
    void add_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 未来可以用于分辨Cell Face Edge

        if (desc.location == StaggerLocation::Cell)
            // 遍历所有物理边界 patch
            for (const auto &patch : topo_->physical_patches)
            {
                apply_cell_patch_U(patch, field_id);
            }
        else if (desc.location == StaggerLocation::FaceXi ||
                 desc.location == StaggerLocation::FaceEt ||
                 desc.location == StaggerLocation::FaceZe)
        {
            for (const auto &patch : topo_->physical_patches)
                apply_face_patch_B(patch, field_id, desc.location);
        }
    };

    void add_boundary(std::vector<std::string> field_name_set)
    {
        for (auto fieldname : field_name_set)
            add_boundary(fieldname);
    }

    void cell_copy_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 未来可以用于分辨Cell Face Edge

        if (desc.location == StaggerLocation::Cell)
            // 遍历所有物理边界 patch
            for (const auto &patch : topo_->physical_patches)
            {
                const int ib = patch.this_block;
                FieldBlock &U = fld_->field(field_id, ib); // 该块上的 U
                apply_cell_copy(U, patch);
            }
    }

    // 初始化，把需要用到的指针和 U_ 的 fid 记起来
    void SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par)
    {
        grd_ = grd;
        fld_ = fld;
        topo_ = topo;
        par_ = par;

        calc_farfield_data();
    }

private:
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Param *par_ = nullptr;

    // 一个 patch 代表某块上的一个物理边界片段
    void apply_cell_patch_U(const TOPO::PhysicalPatch &patch, int32_t field_id);

    // 一个 patch 代表某块上的一个物理边界片段
    void apply_face_patch_B(const TOPO::PhysicalPatch &patch, int32_t field_id, StaggerLocation location);

    // 针对不同边界类型做分发（根据 bc_name）
    void apply_cell_wall(FieldBlock &U, FieldBlock &Bcell, const TOPO::PhysicalPatch &patch);
    // void apply_bc_outflow(FieldBlock &U, const TOPO::PhysicalPatch &patch);
    void apply_cell_farfield(FieldBlock &U, const TOPO::PhysicalPatch &patch);
    void apply_cell_copy(FieldBlock &U, const TOPO::PhysicalPatch &patch);
    void apply_cell_pole(FieldBlock &U, const TOPO::PhysicalPatch &patch);

    void apply_face_copy(FieldBlock &U, const TOPO::PhysicalPatch &patch, int loc);
    void apply_face_farfield(FieldBlock &U, const TOPO::PhysicalPatch &patch, int loc);

    void calc_farfield_data()
    {
        List<double> temp;
        temp = par_->GetDou_List("INITIAL");

        gamma_ = par_->GetDou_List("constant").data["gamma"];
        double NA = par_->GetDou_List("constant").data["NA"];
        double R_uni = par_->GetDou_List("constant").data["R_uni"];

        double c_y = temp.data["c_y"];
        double c_z = temp.data["c_z"];
        double c_x = sqrt(1.0 - c_y * c_y - c_z * c_z);

        double Bx = temp.data["Bx"];
        double By = temp.data["By"];
        double Bz = temp.data["Bz"];
        double B_ref = par_->GetDou("B_ref");
        Bx /= B_ref;
        By /= B_ref;
        Bz /= B_ref;

        numdensity_ref = temp.data["n"];
        T_ref = temp.data["T"];
        double Molecular_mass = temp.data["Molecular_mass"];

        double k_Boltz = R_uni / NA;
        rho_ref = Molecular_mass / NA * numdensity_ref;
        p_ref = numdensity_ref * k_Boltz * T_ref;

        double mu_mag = par_->GetDou_List("constant").data["mu_mag"];
        c_A = B_ref / sqrt(mu_mag * rho_ref);
        Velocity_ref = par_->GetDou_List("INITIAL").data["U"];
        M_A = Velocity_ref / c_A;
        inver_MA2 = 1.0 / (M_A * M_A);

        //=========================================================================
        farfield_rho = par_->GetDou_List("Boundary_Farfield").data["n"] / numdensity_ref;
        double vvv = par_->GetDou_List("Boundary_Farfield").data["U"] / Velocity_ref;
        double cy = par_->GetDou_List("Boundary_Farfield").data["c_y"];
        double cz = par_->GetDou_List("Boundary_Farfield").data["c_z"];
        double temperature = par_->GetDou_List("Boundary_Farfield").data["T"];
        double cx = sqrt(1 - cy * cy - cz * cz);
        double c_sound = sqrt(gamma_ * p_ref / rho_ref);
        double Ma = Velocity_ref / c_sound;
        farfield_u = vvv * cx;
        farfield_v = vvv * cy;
        farfield_w = vvv * cz;
        farfield_T = temperature / T_ref;
        farfield_p = farfield_rho * farfield_T / (gamma_ * Ma * Ma);
        farfield_bx = Bx;
        farfield_by = By;
        farfield_bz = Bz;
    };

    double Velocity_ref, numdensity_ref, T_ref, rho_ref, p_ref, M_A, c_A, gamma_, inver_MA2;
    double farfield_u, farfield_T, farfield_rho, farfield_v, farfield_w, farfield_p, farfield_bx, farfield_by, farfield_bz;

public:
    MHD_Boundary() {};
    ~MHD_Boundary() = default;
};