#pragma once
#include "3_field/2_MPCNS_Field.h"

class MHD_Boundary
{
public:
    //===================================================================================
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

    // 仅仅对Cell施加边界条件
    void add_Cell_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 用于分辨Cell Face Edge

        // 用于辅助施加边界条件，这里不做修改
        int32_t field_id_Bcell = fld_->field_id("B_cell");

        if (desc.location == StaggerLocation::Cell)
        {
            // 遍历所有物理边界 patch
            for (auto &patch : phy_patterns_[desc.location].regions)
            {
                const int ib = patch.this_block;
                FieldBlock &U = fld_->field(field_id, ib);            // 该块上的 U
                FieldBlock &B_cell = fld_->field(field_id_Bcell, ib); // 该块上的 B_cell

                // 这里你可以根据 patch.bc_name 判断是什么类型的边界
                if (patch.bc_name == "Solid_Surface")
                {
                    // apply_cell_wall(U, B_cell, patch);
                    apply_cell_wall_lunar(U, B_cell, patch);
                }
                else if (patch.bc_name == "Outflow")
                {
                    apply_cell_copy(U, patch);
                }
                else if (patch.bc_name == "Farfield")
                {
                    apply_cell_farfield(U, patch);
                }
                else if (patch.bc_name == "Pole")
                {
                    // apply_cell_pole(U, patch);
                    apply_cell_copy(U, patch);
                }
                else
                {
                    // 默认给一个简单的拷贝边界
                    apply_cell_copy(U, patch);
                }
            }
        }
        else
        {
            std::cout << "Fatal Error! ! ! add_Cell_boundary only for Cell !\n"
                      << std::flush;
            exit(-1);
        }
    };

    // 仅仅对Face施加边界条件
    void add_Face_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 未来可以用于分辨Cell Face Edge

        if (desc.location == StaggerLocation::FaceXi ||
            desc.location == StaggerLocation::FaceEt ||
            desc.location == StaggerLocation::FaceZe)
        {
            for (auto &patch : phy_patterns_[desc.location].regions)
            {
                const int ib = patch.this_block;
                FieldBlock &U = fld_->field(field_id, ib); // 该face上的 U

                // 默认给Face使用简单的拷贝边界
                apply_face_copy(U, patch);
            }
        }
        else
        {
            std::cout << "Fatal Error! ! ! add_Face_boundary only for Face xi eta zeta !\n"
                      << std::flush;
            exit(-1);
        }
    };

    // 针对“派生量”的轻量策略，主要用于B_cell
    void add_derived_Cell_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 可用于分辨Cell Face Edge

        if (desc.location == StaggerLocation::Cell && field_name == "B_cell")
        {
            // 遍历所有物理边界 patch
            for (auto &patch : phy_patterns_[desc.location].regions)
            {
                int ib = patch.this_block;
                FieldBlock &U = fld_->field(field_id, ib); // 该块上的 U
                if (patch.bc_name == "Solid_Surface")
                    apply_cell_copy(U, patch);
                // apply_derived_cell_wall(U, patch);
                else if (patch.bc_name == "Pole")
                    // apply_cell_pole(U, patch);
                    apply_cell_copy(U, patch);
                else
                    apply_cell_copy(U, patch);
            }
        }
        else
        {
            std::cout << "Error! ! ! add_derived_Cell_boundary is only for Cell and only for B_cell now! !\n"
                      << std::flush;
            exit(-1);
        }
    };

    void add_Cell_boundary(std::vector<std::string> field_name_set)
    {
        for (auto fieldname : field_name_set)
            add_Cell_boundary(fieldname);
    }

    void add_Face_boundary(std::vector<std::string> field_name_set)
    {
        for (auto fieldname : field_name_set)
            add_Face_boundary(fieldname);
    }
    //===================================================================================

    void add_Edge_pole_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        auto patern = fld_->descriptor(field_id);
        for (auto &top : topo_->physical_patches)
        {
            if (top.bc_name == "Pole")
            {
                const int ib = top.this_block;
                FieldBlock &U = fld_->field(field_id, ib); // 该face上的 U
                if (patern.location == StaggerLocation::EdgeZe)
                {
                    continue; // 一定为0
                }
                else if (patern.location == StaggerLocation::EdgeXi && top.this_box_node.lo.i == top.this_box_node.hi.i - 1)
                {
                    continue;
                }
                else if (patern.location == StaggerLocation::EdgeXi && top.this_box_node.lo.j == top.this_box_node.hi.j - 1)
                {
                    double temp, count;
                    for (int i = top.this_box_node.lo.i; i < top.this_box_node.hi.i - 1; i++)
                        for (int j = top.this_box_node.lo.j; j < top.this_box_node.hi.j; j++)
                        {
                            temp = 0.0;
                            count = 0.0;

                            // Sum
                            for (int k = top.this_box_node.lo.k; k < top.this_box_node.hi.k - 1; k++)
                            {
                                temp += U(i, j, k, 0);
                                count += 1.0;
                            }

                            // Average
                            temp /= count;

                            // Unified
                            for (int k = top.this_box_node.lo.k; k < top.this_box_node.hi.k; k++)
                                U(i, j, k, 0) = temp;
                        }
                }
                else if (patern.location == StaggerLocation::EdgeEt && top.this_box_node.lo.i == top.this_box_node.hi.i - 1)
                {
                    double temp, count;
                    for (int i = top.this_box_node.lo.i; i < top.this_box_node.hi.i; i++)
                        for (int j = top.this_box_node.lo.j; j < top.this_box_node.hi.j - 1; j++)
                        {
                            temp = 0.0;
                            count = 0.0;

                            // Sum
                            for (int k = top.this_box_node.lo.k; k < top.this_box_node.hi.k - 1; k++)
                            {
                                temp += U(i, j, k, 0);
                                count += 1.0;
                            }

                            // Average
                            temp /= count;

                            // Unified
                            for (int k = top.this_box_node.lo.k; k < top.this_box_node.hi.k; k++)
                                U(i, j, k, 0) = temp;
                        }
                }
                else if (patern.location == StaggerLocation::EdgeEt && top.this_box_node.lo.j == top.this_box_node.hi.j - 1)
                {
                    continue;
                }
            }
        }
    }

    void add_Edge_copy_boundary(std::string field_name)
    {
        int32_t field_id = fld_->field_id(field_name);
        const FieldDescriptor &desc = fld_->descriptor(field_id);
        // 未来可以用于分辨Cell Face Edge
        for (auto &patch : phy_patterns_[desc.location].regions)
        {
            const int ib = patch.this_block;
            FieldBlock &U = fld_->field(field_id, ib); // 该face上的 U

            // 默认给Face使用简单的拷贝边界
            apply_edge_copy(U, patch);
        }
    }

private:
    //===================================================================================
    // Physical Pattern Data
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

        // EdgeXi
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::EdgeXi;

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

                std::vector<int> sub(3), sup(3);
                sub = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                // EdgeXi: i cell-based => no +1; j,k node-based => +1
                sup = {p.raw->sup[0], p.raw->sup[1] + 1, p.raw->sup[2] + 1}; // 开区间

                // 只有法向为 xi（abs=1）时，i 为 cell-based，需要扩展一层
                if (abs(p.direction) == 1)
                {
                    if (p.direction < 0)
                        sup[0] += 1; // 低端边界：向外扩展到 ghost
                    else
                        sub[0] -= 1; // 高端边界：向外扩展到 ghost
                }

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0];
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1];
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2];

                tmp.regions.push_back(phy);
            }

            phy_patterns_[StaggerLocation::EdgeXi] = std::move(tmp);
        }

        // EdgeEt
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::EdgeEt;

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

                std::vector<int> sub(3), sup(3);
                sub = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                // EdgeEt: j cell-based => no +1; i,k node-based => +1
                sup = {p.raw->sup[0] + 1, p.raw->sup[1], p.raw->sup[2] + 1}; // 开区间

                // 只有法向为 eta（abs=2）时，j 为 cell-based，需要扩展一层
                if (abs(p.direction) == 2)
                {
                    if (p.direction < 0)
                        sup[1] += 1;
                    else
                        sub[1] -= 1;
                }

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0];
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1];
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2];

                tmp.regions.push_back(phy);
            }

            phy_patterns_[StaggerLocation::EdgeEt] = std::move(tmp);
        }

        // EdgeZe
        {
            PhysicalPattern tmp;
            tmp.location = StaggerLocation::EdgeZe;

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

                std::vector<int> sub(3), sup(3);
                sub = {p.raw->sub[0], p.raw->sub[1], p.raw->sub[2]};
                // EdgeZe: k cell-based => no +1; i,j node-based => +1
                sup = {p.raw->sup[0] + 1, p.raw->sup[1] + 1, p.raw->sup[2]}; // 开区间

                // 只有法向为 zeta（abs=3）时，k 为 cell-based，需要扩展一层
                if (abs(p.direction) == 3)
                {
                    if (p.direction < 0)
                        sup[2] += 1;
                    else
                        sub[2] -= 1;
                }

                phy.box_bound.lo.i = sub[0];
                phy.box_bound.hi.i = sup[0];
                phy.box_bound.lo.j = sub[1];
                phy.box_bound.hi.j = sup[1];
                phy.box_bound.lo.k = sub[2];
                phy.box_bound.hi.k = sup[2];

                tmp.regions.push_back(phy);
            }

            phy_patterns_[StaggerLocation::EdgeZe] = std::move(tmp);
        }

        if (par_->GetInt("myid") == 0)
            std::cout << "***********Finish the Physical Pattern building Process! !************\n\n"
                      << std::flush;
    }
    //===================================================================================

    //===================================================================================
    // Data Tools
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Param *par_ = nullptr;
    //===================================================================================

    //===================================================================================
    // Actual Boundary Conditions
    //-------------------------------------------------------------------------
    // Cell
    //-------------------------------------------------------------------------
    // 针对不同边界类型做分发（根据 bc_name）
    void apply_cell_wall(FieldBlock &U, FieldBlock &Bcell, PhysicalRegion &patch);
    void apply_cell_wall_lunar(FieldBlock &U, FieldBlock &Bcell, PhysicalRegion &patch);
    void apply_cell_farfield(FieldBlock &U, PhysicalRegion &patch);
    void apply_cell_copy(FieldBlock &U, PhysicalRegion &patch);
    void apply_cell_pole(FieldBlock &U, PhysicalRegion &patch);

    // B_cell是导出量derived，壁面需要特殊处理
    void apply_derived_cell_wall(FieldBlock &U, PhysicalRegion &patch);

    //-------------------------------------------------------------------------
    // Face
    //-------------------------------------------------------------------------
    // 针对Face暂时采用拷贝边界
    void apply_face_copy(FieldBlock &U, PhysicalRegion &patch);

    //-------------------------------------------------------------------------
    // Edge
    //-------------------------------------------------------------------------
    // 针对Edge暂时采用拷贝边界
    void apply_edge_copy(FieldBlock &U, PhysicalRegion &patch);
    //===================================================================================

    //===================================================================================
    // 数据初始化 通用无量纲数据
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

        // ngg = fld_->grd->ngg;
    };

    // int ngg;
    double Velocity_ref, numdensity_ref, T_ref, rho_ref, p_ref, M_A, c_A, gamma_, inver_MA2;
    double farfield_u, farfield_T, farfield_rho, farfield_v, farfield_w, farfield_p, farfield_bx, farfield_by, farfield_bz;
    //===================================================================================
public:
    MHD_Boundary() {};
    ~MHD_Boundary() = default;
};