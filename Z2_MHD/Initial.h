#pragma once
#include "3_field/2_MPCNS_Field.h"

class MHD_Initial
{
public:
    void Initialization(Field *fld)
    {
        Param *par = fld->par;
        par->AddParam("Nstep", 0);

        if (par->GetBoo("continue_calc"))
        {
            Read_bin_Initial(fld);
        }
        else
        {
            Common_Initial(fld);
        }

        if (par->GetInt("myid") == 0)
        {
            std::cout << "********************Finish the Initial Process! !*********************\n\n\n";
            std::cout << "====================================="
                      << "===============================================================\n";
            std::cout << "\t \t Multi-Physics Coupling Numerical Simulation Solver Begins NOW!\n";
            std::cout << "====================================="
                      << "===============================================================\n\n\n";
        }
    }

    /**
     * @brief   从par的LIST中提取用户输入的INITIAL LIST信息，存入初始化数组U_initial中
     */
    void Initial_Transfer(std::map<std::string, double1D> &U_initial_General, Param *par)
    {
        List<double> temp;
        temp = par->GetDou_List("INITIAL");
        double1D &U_initial = U_initial_General["U_"];

        // ---- 常数 ----
        double gamma = par->GetDou_List("constant").data["gamma"];
        double NA = par->GetDou_List("constant").data["NA"];
        double R_uni = par->GetDou_List("constant").data["R_uni"];

        // ---- 流向单位向量 ----
        double c_y = temp.data["c_y"];
        double c_z = temp.data["c_z"];
        double c_x = std::sqrt(1.0 - c_y * c_y - c_z * c_z);

        // ---- 背景磁场（物理单位）----
        double Bx_phy = temp.data["Bx"];
        double By_phy = temp.data["By"];
        double Bz_phy = temp.data["Bz"];

        double B_ref = par->GetDou("B_ref");
        // 无量纲背景磁场
        double Bx = Bx_phy / B_ref;
        double By = By_phy / B_ref;
        double Bz = Bz_phy / B_ref;

        // 把无量纲 B_add 存进 Param，后面求解器/通量都可以直接用
        par->AddParam("B_add_x", Bx);
        par->AddParam("B_add_y", By);
        par->AddParam("B_add_z", Bz);

        // ---- 流体参考量 ----
        double U_ref = temp.data["U"];

        double n_ref = temp.data["n"];
        double T_ref = temp.data["T"];
        double Molecular_mass = temp.data["Molecular_mass"];

        double k_Boltz = R_uni / NA;
        double rho_ref = Molecular_mass / NA * n_ref;
        double p_ref = n_ref * k_Boltz * T_ref;

        double mu_mag = par->GetDou_List("constant").data["mu_mag"];
        double c_A = B_ref / std::sqrt(mu_mag * rho_ref);
        double Velocity_ref = U_ref;
        double M_A = Velocity_ref / c_A; // Alfven Mach 数

        // 把无量纲 M_A^(-2)加入par便于调用
        par->AddParam("inver_MA2", 1.0 / (M_A * M_A));

        // ---- 无量纲原始量 ----
        double rho0 = 1.0;
        double u0 = c_x;
        double v0 = c_y;
        double w0 = c_z;
        double v2 = u0 * u0 + v0 * v0 + w0 * w0; // 一般应为 1

        // 无量纲气体压力
        double p0 = p_ref / (rho_ref * U_ref * U_ref);

        // 外加磁场对应的无量纲磁能密度
        double B2 = Bx * Bx + By * By + Bz * Bz;
        double Emag0 = 0.5 * B2 / (M_A * M_A);

        // ---- 填充 U_initial ----
        U_initial(0) = rho0;
        U_initial(1) = rho0 * u0;
        U_initial(2) = rho0 * v0;
        U_initial(3) = rho0 * w0;
        U_initial(4) = 0.5 * rho0 * v2      // 动能
                       + p0 / (gamma - 1.0) // 内能
                       + Emag0;             // 总磁场（这里只是 B_add）的能量

        // 感应磁场一律从 0 开始：
        //   B_ind = 0  ->  B_total = B_add
        //   具体 Face 上的 B_ind 下面在 Common_Initial 里设。
    };

public:
    MHD_Initial() {};
    ~MHD_Initial() = default;

private:
    void Common_Initial(Field *fld)
    {
        Param *par = fld->par;
        if (par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        int ngg = par->GetInt("ngg");
        std::map<std::string, double1D> U_initial;

        const FieldDescriptor &desc = fld->descriptor(fld->field_id("U_"));
        double1D tmp;
        U_initial.insert(std::pair<std::string, double1D>{desc.name, tmp});
        U_initial[desc.name].SetSize(desc.ncomp);

        // 1) 构造 U_initial（包含 B_add 的能量）
        Initial_Transfer(U_initial, par);

        // 无量纲 B_add
        double Bx = par->GetDou("B_add_x");
        double By = par->GetDou("B_add_y");
        double Bz = par->GetDou("B_add_z");

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            double1D &uinitial = U_initial["U_"];

            FieldBlock &Ublk = fld->field("U_", iblock);
            FieldBlock &B_cell = fld->field("B_cell", iblock);
            const Int3 &sub = Ublk.get_lo();
            const Int3 &sup = Ublk.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < desc.ncomp; ll++)
                            Ublk(i, j, k, ll) = uinitial(ll);
                        B_cell(i, j, k, 0) = Bx;
                        B_cell(i, j, k, 1) = By;
                        B_cell(i, j, k, 2) = Bz;
                    }
        }

        par->AddParam("Physic_Time", 0.0);

        // MHD Part
        int fid_Bxi = fld->field_id("B_xi");
        int fid_Beta = fld->field_id("B_eta");
        int fid_Bze = fld->field_id("B_zeta");
        int fid_Ex = fld->field_id("E_xi");
        int fid_Ey = fld->field_id("E_eta");
        int fid_Ez = fld->field_id("E_zeta");

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            // --- FaceXi 上的感应 B^ξ ---
            {
                FieldBlock &Bxi_blk = fld->field(fid_Bxi, iblock);
                const Int3 &lo = Bxi_blk.get_lo();
                const Int3 &hi = Bxi_blk.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Bxi_blk(i, j, k, 0) = 0.0; // 感应场从 0 开始
            }
            // --- FaceEt 上的感应 B^η ---
            {
                FieldBlock &Beta_blk = fld->field(fid_Beta, iblock);
                const Int3 &lo = Beta_blk.get_lo();
                const Int3 &hi = Beta_blk.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Beta_blk(i, j, k, 0) = 0.0;
            }
            // --- FaceZe 上的感应 B^ζ ---
            {
                FieldBlock &Bze_blk = fld->field(fid_Bze, iblock);
                const Int3 &lo = Bze_blk.get_lo();
                const Int3 &hi = Bze_blk.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Bze_blk(i, j, k, 0) = 0.0;
            }

            // --- FaceXi 上的感应 B^ξ ---
            {
                FieldBlock &Bxi_blk = fld->field(fid_Bxi, iblock);
                FieldBlock &xi = fld->field("JDxi", iblock);
                const Int3 &lo = Bxi_blk.inner_lo();
                const Int3 &hi = Bxi_blk.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Bxi_blk(i, j, k, 0) = Bx * xi(i, j, k, 0) + By * xi(i, j, k, 1) + Bz * xi(i, j, k, 2); // 感应场从 0 开始
            }
            // --- FaceEt 上的感应 B^η ---
            {
                FieldBlock &Beta_blk = fld->field(fid_Beta, iblock);
                FieldBlock &xi = fld->field("JDet", iblock);
                const Int3 &lo = Beta_blk.inner_lo();
                const Int3 &hi = Beta_blk.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Beta_blk(i, j, k, 0) = Bx * xi(i, j, k, 0) + By * xi(i, j, k, 1) + Bz * xi(i, j, k, 2);
            }
            // --- FaceZe 上的感应 B^ζ ---
            {
                FieldBlock &Bze_blk = fld->field(fid_Bze, iblock);
                FieldBlock &xi = fld->field("JDze", iblock);
                const Int3 &lo = Bze_blk.inner_lo();
                const Int3 &hi = Bze_blk.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Bze_blk(i, j, k, 0) = Bx * xi(i, j, k, 0) + By * xi(i, j, k, 1) + Bz * xi(i, j, k, 2);
            }

            // --- Edge 上的 EMF 也设零（可选，但推荐） ---
            {
                FieldBlock &Ex_blk = fld->field(fid_Ex, iblock);
                const Int3 &lo = Ex_blk.get_lo();
                const Int3 &hi = Ex_blk.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Ex_blk(i, j, k, 0) = 0.0;
            }

            {
                FieldBlock &Ey_blk = fld->field(fid_Ey, iblock);
                const Int3 &lo = Ey_blk.get_lo();
                const Int3 &hi = Ey_blk.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Ey_blk(i, j, k, 0) = 0.0;
            }

            {
                FieldBlock &Ez_blk = fld->field(fid_Ez, iblock);
                const Int3 &lo = Ez_blk.get_lo();
                const Int3 &hi = Ez_blk.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Ez_blk(i, j, k, 0) = 0.0;
            }
        }
    }

    void Read_bin_Initial(Field *fld)
    {
        Param *par = fld->par;
        if (par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        //----------------------------------------------------------------------
        // 打开二进制文件  read the "flow_field   *.bin" binary files
        std::string _my_id_s;
        int my_id = par->GetInt("myid");
        if (my_id < 10)
        {
            _my_id_s = "   " + std::to_string(my_id);
        }
        else if (my_id < 100)
        {
            _my_id_s = "  " + std::to_string(my_id);
        }
        else if (my_id < 1000)
        {
            _my_id_s = " " + std::to_string(my_id);
        }
        else // 这说明并行进程数不得超过9999
        {
            _my_id_s = std::to_string(my_id);
        }
        std::ifstream inFile_bin("./DATA/flow_field" + _my_id_s + ".bin",
                                 std::ios_base::in | std::ios::binary);
        if (!inFile_bin)
        {
            std::cerr << "Error: cannot open ./DATA/flow_field"
                      << _my_id_s << ".bin" << std::endl;
            return;
        }

        //----------------------------------------------------------------------
        // 一些小工具：读名字（带长度）
        auto bin_read_name = [](std::ifstream &file) -> std::string
        {
            int32_t temp_i;
            //  1、读入名称长度
            file.read((char *)&temp_i, sizeof(int32_t));
            //  2、分配空间并读入名称，存入string后释放
            char *temp_char = new char[temp_i];
            file.read(temp_char, temp_i);
            std::string physical_name;
            physical_name.reserve(temp_i);
            physical_name = temp_char;
            delete[] temp_char;
            return physical_name;
        };

        // 读 int 参数： (length + name + value)
        auto bin_read_int_param = [&](std::ifstream &file)
        {
            std::string pname = bin_read_name(file);
            int value = 0;
            file.read(reinterpret_cast<char *>(&value), sizeof(int));
            // 根据你的 Param 接口写回去，这里假定有 SetInt / SetDou
            par->AddParam(pname, value);
        };

        // 读 double 参数： (length + name + value)
        auto bin_read_double_param = [&](std::ifstream &file)
        {
            std::string pname = bin_read_name(file);
            double value = 0.0;
            file.read(reinterpret_cast<char *>(&value), sizeof(double));
            par->AddParam(pname, value);
        };

        // kind(int) -> StaggerLocation 的反向映射
        auto map_int2loc = [](int kind) -> StaggerLocation
        {
            switch (kind)
            {
            case 0:
                return StaggerLocation::Cell;
            case 1:
                return StaggerLocation::FaceXi;
            case 2:
                return StaggerLocation::FaceEt;
            case 3:
                return StaggerLocation::FaceZe;
            case 4:
                return StaggerLocation::EdgeXi;
            case 5:
                return StaggerLocation::EdgeEt;
            case 6:
                return StaggerLocation::EdgeZe;
            case 7:
                return StaggerLocation::Node;
            default:
                return StaggerLocation::Cell; // fallback
            }
        };

        //----------------------------------------------------------------------
        // 1. 读入 int 类型记录变量
        int32_t num_int = 0;
        inFile_bin.read(reinterpret_cast<char *>(&num_int), sizeof(int32_t));
        for (int32_t n = 0; n < num_int; ++n)
        {
            bin_read_int_param(inFile_bin);
        }

        //----------------------------------------------------------------------
        // 2. 读入 double 类型记录变量
        int32_t num_double = 0;
        inFile_bin.read(reinterpret_cast<char *>(&num_double), sizeof(int32_t));
        for (int32_t n = 0; n < num_double; ++n)
        {
            bin_read_double_param(inFile_bin);
        }

        //----------------------------------------------------------------------
        // 3. 读入 data 信息
        // 3.1 网格块数
        int32_t nblock_file = 0;
        inFile_bin.read(reinterpret_cast<char *>(&nblock_file), sizeof(int32_t));

        int32_t nblock = fld->num_blocks();
        if (nblock_file != nblock)
        {
            std::cerr << "Error: nblock in file (" << nblock_file
                      << ") != fld->num_blocks() (" << nblock << ")\n";
            exit(-1);
        }

        // 3.2 物理场个数（文件中实际输出的）
        int32_t n_field_file = 0;
        inFile_bin.read(reinterpret_cast<char *>(&n_field_file), sizeof(int32_t));

        // 注意：文件里只包含你写出时没有被 skip 的物理场，
        // 我们按名字匹配到当前 Field 中的 fid 再回填数据。
        for (int32_t ifld = 0; ifld < n_field_file; ++ifld)
        {
            // 3.2.1 读 FieldDescriptor：name + kind + ncomp + nghost
            std::string fname = bin_read_name(inFile_bin);

            int32_t kind = 0;
            int32_t ncomp_file = 0;
            int32_t nghost_file = 0;
            inFile_bin.read(reinterpret_cast<char *>(&kind), sizeof(int32_t));
            inFile_bin.read(reinterpret_cast<char *>(&ncomp_file), sizeof(int32_t));
            inFile_bin.read(reinterpret_cast<char *>(&nghost_file), sizeof(int32_t));

            // 在当前 Field 中查这个名字的 field ID
            int fid = -1;
            try
            {
                fid = fld->field_id(fname); // Field::field_id(name):contentReference[oaicite:1]{index=1}
            }
            catch (const std::out_of_range &)
            {
                std::cerr << "Error: field name \"" << fname
                          << "\" not found in current Field. Skip this field.\n";
                // 如果找不到，只能把这块数据从文件中 “读掉丢弃”
                // 以保证后面能对齐。
                // 这里简单实现一个丢弃过程：
                for (int ib = 0; ib < nblock_file; ++ib)
                {
                    int32_t lo[3], hi[3];
                    inFile_bin.read(reinterpret_cast<char *>(&lo[0]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&lo[1]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&lo[2]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&hi[0]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&hi[1]), sizeof(int32_t));
                    inFile_bin.read(reinterpret_cast<char *>(&hi[2]), sizeof(int32_t));

                    const int Ni = hi[0] - lo[0];
                    const int Nj = hi[1] - lo[1];
                    const int Nk = hi[2] - lo[2];
                    const std::size_t nval =
                        static_cast<std::size_t>(Ni) *
                        static_cast<std::size_t>(Nj) *
                        static_cast<std::size_t>(Nk) *
                        static_cast<std::size_t>(ncomp_file);

                    // 直接跳过这些 double
                    inFile_bin.seekg(static_cast<std::streamoff>(
                                         nval * sizeof(double)),
                                     std::ios_base::cur);
                }
                continue; // 进入下一个 field
            }

            const FieldDescriptor &desc = fld->descriptor(fid); //: contentReference[oaicite:2]{index=2}

            // 一些健壮性检查（可选）
            if (desc.ncomp != ncomp_file)
            {
                std::cerr << "Error: field \"" << fname << "\" ncomp mismatch: file="
                          << ncomp_file << ", current=" << desc.ncomp << std::endl;
                exit(-1);
            }
            if (desc.nghost != nghost_file)
            {
                std::cerr << "Error: field \"" << fname << "\" nghost mismatch: file="
                          << nghost_file << ", current=" << desc.nghost << std::endl;
                exit(-1);
            }
            // 如果你希望严格一致，也可以这里直接 return 或 throw

            // 3.2.2 对每个 block 读 lo/hi 和数据
            for (int ib = 0; ib < nblock_file; ++ib)
            {
                int32_t lo[3], hi[3];
                inFile_bin.read(reinterpret_cast<char *>(&lo[0]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&lo[1]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&lo[2]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&hi[0]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&hi[1]), sizeof(int32_t));
                inFile_bin.read(reinterpret_cast<char *>(&hi[2]), sizeof(int32_t));

                FieldBlock &fb = fld->field(fid, ib); //: contentReference[oaicite:3]{index=3}

                // 可选：检查 lo/hi 是否和当前 fb 一致
                // const Int3 &cur_lo = fb.get_lo();
                // const Int3 &cur_hi = fb.get_hi();

                // 实际数据读取：完全按写出时的循环顺序读回来
                for (int i = lo[0]; i < hi[0]; ++i)
                    for (int j = lo[1]; j < hi[1]; ++j)
                        for (int k = lo[2]; k < hi[2]; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                            {
                                double val = 0.0;
                                inFile_bin.read(reinterpret_cast<char *>(&val), sizeof(double));
                                fb(i, j, k, m) = val;
                            }
            }
        }

        inFile_bin.close();
    }
};
