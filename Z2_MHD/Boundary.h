#pragma once
#include "3_field/2_MPCNS_Field.h"

class MHD_Boundary
{
private:
    struct PhysicalRegion
    {
        //--------------------------------------------------------------------------
        // Block_ID
        int this_block; // 我是哪个 block
        //--------------------------------------------------------------------------
        // IJK_range
        Box3 box_bound; // 真正的网格边界,一般为inner计算区域，不包含ghost 区
        Int3 cycle;     // 往虚网格的方向（朝block外）
        //--------------------------------------------------------------------------

        std::string this_block_name;
        int bc_id;                              // boundary_num
        std::string bc_name;                    // boundary_name
        const Physical_Boundary *raw = nullptr; // 回指原始结构
    };

    struct PhysicalPattern
    {
        StaggerLocation location; // Cell / FaceXi / ...
        std::vector<PhysicalRegion> regions;
    };

    std::map<StaggerLocation, PhysicalPattern> phy_patterns_;

    void build_boundary_pattern()
    {
        if (par_->GetInt("myid") == 0)
            std::cout << "---->Start building boundary pattern...\n";
        // Cell
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::Cell;
            for (auto &p : topo_->physical_patches)
            {
                PhysicalRegion phy;
                phy.this_block = p.this_block;
                phy.this_block_name = p.this_block_name;
                phy.bc_id = p.bc_id;
                phy.bc_name = p.bc_name;
                phy.raw = p.raw;

                phy.cycle.i = p.raw->cycle[0];
                phy.cycle.j = p.raw->cycle[1];
                phy.cycle.k = p.raw->cycle[2];

                int sub[3] = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                int sup[3] = {p.raw->sup[0], p.raw->sup[1], p.raw->sup[2]}; // 开区间，正好形成cell范围

                sub[abs(p.direction) - 1] += (p.direction > 0) ? -1 : 0;
                sup[abs(p.direction) - 1] += (p.direction > 0) ? 0 : 1;

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0]; // 开区间
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1]; // 开区间
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2]; // 开区间

                tmp.regions.push_back(phy);
            }
            phy_patterns_[StaggerLocation::Cell] = std::move(tmp);
        }

        // FaceXi
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::FaceXi;
            for (auto &p : topo_->physical_patches)
            {
                PhysicalRegion phy;
                phy.this_block = p.this_block;
                phy.this_block_name = p.this_block_name;
                phy.bc_id = p.bc_id;
                phy.bc_name = p.bc_name;
                phy.raw = p.raw;
                phy.cycle.i = p.raw->cycle[0];
                phy.cycle.j = p.raw->cycle[1];
                phy.cycle.k = p.raw->cycle[2];

                std::vector<int> sub, sup;
                sub.resize(3);
                sup.resize(3);
                sub = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                sup = {p.raw->sup[0] + 1, p.raw->sup[1], p.raw->sup[2]}; // 开区间
                if (abs(p.direction) == 1)
                {
                }
                else if (abs(p.direction) == 2)
                {
                    if (p.direction < 0)
                        sup[1] += 1;
                    else
                        sub[1] -= 1;
                }
                else if (abs(p.direction) == 3)
                {
                    if (p.direction < 0)
                        sup[2] += 1;
                    else
                        sub[2] -= 1;
                }

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0]; // 开区间
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1]; // 开区间
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2]; // 开区间
                tmp.regions.push_back(phy);
            }
            phy_patterns_[StaggerLocation::FaceXi] = std::move(tmp);
        }

        // FaceEt
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::FaceEt;
            for (auto &p : topo_->physical_patches)
            {
                PhysicalRegion phy;
                phy.this_block = p.this_block;
                phy.this_block_name = p.this_block_name;
                phy.bc_id = p.bc_id;
                phy.bc_name = p.bc_name;
                phy.raw = p.raw;
                phy.cycle.i = p.raw->cycle[0];
                phy.cycle.j = p.raw->cycle[1];
                phy.cycle.k = p.raw->cycle[2];

                std::vector<int> sub, sup;
                sub.resize(3);
                sup.resize(3);
                sub = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                sup = {p.raw->sup[0], p.raw->sup[1] + 1, p.raw->sup[2]}; // 开区间
                if (abs(p.direction) == 2)
                {
                }
                else if (abs(p.direction) == 1)
                {
                    if (p.direction < 0)
                        sup[0] += 1;
                    else
                        sub[0] -= 1;
                }
                else if (abs(p.direction) == 3)
                {
                    if (p.direction < 0)
                        sup[2] += 1;
                    else
                        sub[2] -= 1;
                }

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0]; // 开区间
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1]; // 开区间
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2]; // 开区间

                tmp.regions.push_back(phy);
            }
            phy_patterns_[StaggerLocation::FaceEt] = std::move(tmp);
        }

        // FaceZe
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::FaceZe;
            for (auto &p : topo_->physical_patches)
            {
                PhysicalRegion phy;
                phy.this_block = p.this_block;
                phy.this_block_name = p.this_block_name;
                phy.bc_id = p.bc_id;
                phy.bc_name = p.bc_name;
                phy.raw = p.raw;
                phy.cycle.i = p.raw->cycle[0];
                phy.cycle.j = p.raw->cycle[1];
                phy.cycle.k = p.raw->cycle[2];

                std::vector<int> sub, sup;
                sub.resize(3);
                sup.resize(3);
                sub = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                sup = {p.raw->sup[0], p.raw->sup[1], p.raw->sup[2] + 1}; // 开区间
                if (abs(p.direction) == 3)
                {
                }
                else if (abs(p.direction) == 1)
                {
                    if (p.direction < 0)
                        sup[0] += 1;
                    else
                        sub[0] -= 1;
                }
                else if (abs(p.direction) == 2)
                {
                    if (p.direction < 0)
                        sup[1] += 1;
                    else
                        sub[1] -= 1;
                }

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0]; // 开区间
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1]; // 开区间
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2]; // 开区间

                tmp.regions.push_back(phy);
            }
            phy_patterns_[StaggerLocation::FaceZe] = std::move(tmp);
        }
        if (par_->GetInt("myid") == 0)
            std::cout << "***********Finish the Physical Pattern building Process! !************\n\n"
                      << std::flush;
    }

public:
    void add_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 未来可以用于分辨Cell Face Edge

        if (desc.location == StaggerLocation::Cell)
            // 遍历所有物理边界 patch
            for (auto &patch : phy_patterns_[desc.location].regions)
            {
                apply_cell_patch_U(patch, field_id);
            }
        else if (desc.location == StaggerLocation::FaceXi ||
                 desc.location == StaggerLocation::FaceEt ||
                 desc.location == StaggerLocation::FaceZe)
        {
            for (auto &patch : phy_patterns_[desc.location].regions)
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
            for (auto &patch : phy_patterns_[desc.location].regions)
            {
                int ib = patch.this_block;
                FieldBlock &U = fld_->field(field_id, ib); // 该块上的 U
                if (patch.bc_name == "Pole")
                    apply_cell_pole(U, patch);
                // else if (patch.bc_name == "Solid_Surface")
                //     apply_cell_wall_B(U, patch);
                else
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

        build_boundary_pattern();
    }

private:
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Param *par_ = nullptr;

    // 一个 patch 代表某块上的一个物理边界片段
    void apply_cell_patch_U(PhysicalRegion &patch, int32_t field_id);

    // 一个 patch 代表某块上的一个物理边界片段
    void apply_face_patch_B(PhysicalRegion &patch, int32_t field_id, StaggerLocation location);

    // 针对不同边界类型做分发（根据 bc_name）
    void apply_cell_wall(FieldBlock &U, FieldBlock &Bcell, PhysicalRegion &patch);
    void apply_cell_farfield(FieldBlock &U, PhysicalRegion &patch);
    void apply_cell_copy(FieldBlock &U, PhysicalRegion &patch);
    void apply_cell_pole(FieldBlock &U, PhysicalRegion &patch);

    void apply_cell_wall_B(FieldBlock &U, PhysicalRegion &patch);

    void apply_face_copy(FieldBlock &U, PhysicalRegion &patch);

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

        ngg = fld_->grd->ngg;
    };

    int ngg;
    double Velocity_ref, numdensity_ref, T_ref, rho_ref, p_ref, M_A, c_A, gamma_, inver_MA2;
    double farfield_u, farfield_T, farfield_rho, farfield_v, farfield_w, farfield_p, farfield_bx, farfield_by, farfield_bz;

public:
    MHD_Boundary() {};
    ~MHD_Boundary() = default;
};