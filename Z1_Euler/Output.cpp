#include "Output.h"

//======================================================================
// 主函数：xyz 为节点，物理量为 cell-centered 的 Tecplot PLT 输出
//======================================================================
void Euler_Output::output_plt_field()
{
    Grid *grd = fld->grd;
    Param *par = fld->par;
    const double gamma = par->GetDou_List("constant").data["gamma"]; // 理想气体比热比

    //============================== 1. 准备文件 ==============================
    std::string path = "./DATA/flow_field" + filename + ".plt";
    const char *filenameall = path.c_str();

    // 确保 DATA 目录存在
    {
        std::ofstream test(path, std::ios::out);
        if (!test)
        {
            std::string cmd = "mkdir -p DATA";
            std::system(cmd.c_str());
        }
    }

    FILE *fp = std::fopen(filenameall, "wb");
    if (!fp)
    {
        std::printf("Failed to open file %s\n", filenameall);
        std::exit(-1);
    }

    //============================== 2. 一些 lambda 工具 ==============================
    auto write_int32 = [&](int32_t v)
    { std::fwrite(&v, 4, 1, fp); };
    auto write_float = [&](float v)
    { std::fwrite(&v, 4, 1, fp); };
    auto write_double = [&](double v)
    { std::fwrite(&v, 8, 1, fp); };

    // Tecplot 的字符串写法：每个字符写成 int，一直写到 '\0' 的 int
    auto write_string = [&](const std::string &s)
    {
        int value = 0;

        for (char c : s)
        {
            value = static_cast<int>(c);
            std::fwrite(&value, sizeof(int), 1, fp);
        }
        value = static_cast<int>('\0');
        std::fwrite(&value, sizeof(int), 1, fp);
    };

    auto write_file_header = [&]()
    {
        // magic number "#!TDV112"
        const char magic[8] = {'#', '!', 'T', 'D', 'V', '1', '1', '2'};
        std::fwrite(magic, 8, 1, fp);

        // byte order: 1 = little endian
        int32_t byte_order = 1;
        write_int32(byte_order);

        // file type: 0 = FULL
        int32_t file_type = 0;
        write_int32(file_type);

        // title
        std::string title = "flow_field";
        write_string(title);

        // 变量名：x,y,z,rho,u,v,w,p
        std::vector<std::string> var_names = {
            "x", "y", "z",
            "rho", "u", "v", "w", "p"};

        int32_t n_var = static_cast<int32_t>(var_names.size());
        write_int32(n_var);

        for (const auto &name : var_names)
            write_string(name);

        return n_var; // 返回变量个数，后面还要用
    };

    auto write_zone_header = [&](int iblock, int num_of_var, int *var_location)
    {
        const Block &blk = grd->grids(iblock);

        // Tecplot zone marker
        write_float(299.0f);

        std::string zone_name = "Zone" + std::to_string(iblock);
        write_string(zone_name);

        // parent zone id
        // parent zone id：0 表示没有父区（更保守）
        int32_t parentZone = -1;
        write_int32(parentZone);
        // strand id
        int32_t strandId = -2;
        write_int32(strandId);
        // solution time
        // solution time：写成真实物理时间
        double solTime = par->GetDou("Physic_Time");
        write_double(solTime);
        // zone color
        write_int32(-1);
        // zone type: 0 = ORDERED
        write_int32(0);
        // 是否指定每个变量的位置：1=是
        int32_t specify_var_location = 1;
        write_int32(specify_var_location);
        // 写每个变量的 location（0=nodal, 1=cell-centered）
        for (int v = 0; v < num_of_var; ++v)
            write_int32(var_location[v]);

        // raw face neighbor, misc face（这里不用）
        // Are raw local 1-to-1 face neighbors supplied? 0=FALSE（对 ORDERED 必须是 0）
        write_int32(0);
        // Number of miscellaneous user-defined face neighbor connections: 0
        write_int32(0);

        // Ordered zone 的 IMax,JMax,KMax = 节点数 = mx+1, my+1, mz+1
        int32_t IMax = blk.mx + 1;
        int32_t JMax = blk.my + 1;
        int32_t KMax = blk.mz + 1;
        write_int32(IMax);
        write_int32(JMax);
        write_int32(KMax);

        // aux data count = 0
        // AuxData 结束标志：0 = 没有 AuxData
        write_int32(0);
    };

    auto write_zone_data_header = [&](int num_of_var)
    {
        // zone marker
        write_float(299.0f);

        // 每个变量的 data type: 1 = float 2=Double, 3=LongInt,  4=ShortInt, 5=Byte,   6=Bit
        int32_t data_type = 1;
        for (int v = 0; v < num_of_var; ++v)
            write_int32(data_type);

        // passive variables: 0 = 无
        write_int32(0);
        // shared variables: 0 = 无
        write_int32(0);
        // shared connectivity: -1 = 无
        write_int32(-1);
    };
    //============================== 3. 文件 Header（全局） ==============================

    const int n_var = write_file_header();

    // 每个变量的位置：0= nodal, 1 = cell-centered
    // 我们约定：前 3 个 (x,y,z) 为节点变量，其余为 cell-centered
    std::vector<int32_t> var_location(n_var, 1);
    if (n_var >= 3)
    {
        var_location[0] = 0; // x
        var_location[1] = 0; // y
        var_location[2] = 0; // z
    }

    //============================== 4. 找到 U 场 ==============================
    const int fid_U = fld->field_id("U_");
    const int fid_PV = fld->field_id("PV_");
    const FieldDescriptor &descU = fld->descriptor(fid_U);
    if (descU.location != StaggerLocation::Cell || descU.ncomp < 5)
    {
        std::cout << "Fatal Error!!! Field \"U\" must be Cell-centered with at least 5 components (rho, rho*u, rho*v, rho*w, E).\n";
        std::exit(-1);
    }

    //============================== 5. Zone Header（每个 block 一个） ==============================

    for (int ib = 0; ib < grd->nblock; ++ib)
    {
        const Block &blk = grd->grids(ib);
        if (blk.block_name != "Fluid")
            continue;

        write_zone_header(ib, n_var, var_location.data());
    }

    // 结束 Header 的 marker
    write_float(357.0f);

    //============================== 6. Data Section（每个 zone 的数据） ==============================

    // 变量在 var_names 中的索引约定
    enum VarIndex
    {
        VID_X = 0,
        VID_Y = 1,
        VID_Z = 2,
        VID_RHO = 3,
        VID_U = 4,
        VID_V = 5,
        VID_W = 6,
        VID_P = 7
    };

    //----------------- 逐 block 写数据 -----------------
    for (int iblock = 0; iblock < grd->nblock; ++iblock)
    {
        Block &blk = grd->grids(iblock);
        if (blk.block_name != "Fluid")
            continue;

        const int Ni = blk.mx; // cell 数
        const int Nj = blk.my;
        const int Nk = blk.mz;
        const int Ni_node = Ni + 1;
        const int Nj_node = Nj + 1;
        const int Nk_node = Nk + 1;

        FieldBlock &Ublk = fld->field(fid_U, iblock);
        FieldBlock &PVblk = fld->field(fid_PV, iblock);

        // 写这个 zone 的数据 header
        write_zone_data_header(n_var);

        // 对每个变量逐个写：
        for (int v = 0; v < n_var; ++v)
        {
            const bool is_nodal = (var_location[v] == 0);

            double minv = 0.0;
            double maxv = 0.0;
            bool first = true;

            //================== 6.1 先算 min/max ==================
            if (is_nodal)
            {
                // ---- 节点变量：x,y,z ----
                for (int k = 0; k < Nk_node; ++k)
                {
                    for (int j = 0; j < Nj_node; ++j)
                    {
                        for (int i = 0; i < Ni_node; ++i)
                        {
                            double val = 0.0;
                            if (v == VID_X)
                                val = blk.x(i, j, k);
                            else if (v == VID_Y)
                                val = blk.y(i, j, k);
                            else if (v == VID_Z)
                                val = blk.z(i, j, k);
                            else
                                val = 0.0; // 理论上不会来

                            if (first)
                            {
                                minv = maxv = val;
                                first = false;
                            }
                            else
                            {
                                if (val < minv)
                                    minv = val;
                                if (val > maxv)
                                    maxv = val;
                            }
                        }
                    }
                }
            }
            else
            {
                // ---- cell-centered 变量：rho,u,v,w,p ----
                for (int kc = 0; kc < Nk; ++kc)
                {
                    for (int jc = 0; jc < Nj; ++jc)
                    {
                        for (int ic = 0; ic < Ni; ++ic)
                        {
                            double val = 0.0;
                            if (v == VID_RHO)
                                val = Ublk(ic, jc, kc, 0);
                            else if (v == VID_U)
                                val = PVblk(ic, jc, kc, 0);
                            else if (v == VID_V)
                                val = PVblk(ic, jc, kc, 1);
                            else if (v == VID_W)
                                val = PVblk(ic, jc, kc, 2);
                            else if (v == VID_P)
                                val = PVblk(ic, jc, kc, 3);
                            else
                                val = 0.0;

                            if (first)
                            {
                                minv = maxv = val;
                                first = false;
                            }
                            else
                            {
                                if (val < minv)
                                    minv = val;
                                if (val > maxv)
                                    maxv = val;
                            }
                        }
                    }
                }
            }

            // 写 min/max
            write_double(minv);
            write_double(maxv);
        }

        //================== 6.2 再写具体数据 ==================
        // 对每个变量逐个写：
        for (int v = 0; v < n_var; ++v)
        {
            const bool is_nodal = (var_location[v] == 0);

            if (is_nodal)
            {
                // ---- 节点变量：遍历所有节点 (i,j,k)，节点个数 = (mx+1)*(my+1)*(mz+1) ----
                for (int k = 0; k < Nk_node; ++k)
                {
                    for (int j = 0; j < Nj_node; ++j)
                    {
                        for (int i = 0; i < Ni_node; ++i)
                        {
                            double val = 0.0;
                            if (v == VID_X)
                                val = blk.x(i, j, k);
                            else if (v == VID_Y)
                                val = blk.y(i, j, k);
                            else if (v == VID_Z)
                                val = blk.z(i, j, k);
                            else
                                val = 0.0;

                            const float fval = static_cast<float>(val);
                            write_float(fval);
                        }
                    }
                }
            }
            else
            {
                // ---- cell-centered 变量：遍历所有 cell (ic,jc,kc)，个数 = mx*my*mz ----
                // Tecplot 存储 cell-centered 变量时，为了高效读取，实际写的是 IMax * JMax * (KMax - 1) 个数，
                // IMax, JMax, KMax 是点数（节点数）。 所以会在最快变化的索引方向（I 方向）写入一些 0 作为 ghost 值。
                int IMax = Ni + 1;
                int JMax = Nj + 1;
                int KMax = Nk + 1;

                // 对每个 cell-centered 变量 v：
                for (int k = 0; k < KMax - 1; ++k)
                { // 0 .. Nk-1
                    for (int j = 0; j < JMax; ++j)
                    { // 0 .. Nj
                        for (int i = 0; i < IMax; ++i)
                        { // 0 .. Ni
                            double val = 0.0;
                            if (i < Ni && j < Nj)
                            {
                                // (i,j,k) 在真正的 cell 范围内
                                int ic = i;
                                int jc = j;
                                int kc = k;

                                if (v == VID_RHO)
                                    val = Ublk(ic, jc, kc, 0);
                                else if (v == VID_U)
                                    val = PVblk(ic, jc, kc, 0);
                                else if (v == VID_V)
                                    val = PVblk(ic, jc, kc, 1);
                                else if (v == VID_W)
                                    val = PVblk(ic, jc, kc, 2);
                                else if (v == VID_P)
                                    val = PVblk(ic, jc, kc, 3);
                            }
                            else
                            {
                                // ghost，保持 0.0
                            }
                            write_float((float)val);
                        }
                    }
                }
            }
        } // end for(v)
    } // end for(iblock)

    std::fclose(fp);
}

void Euler_Output::output_bin_field()
{
    Param *par = fld->par;
    //--------------------------------------------------------------------------
    // 输出二进制文件
    std::ofstream outFile_bin("./DATA/flow_field" + filename + ".bin", std::ios_base::out | std::ios::binary);
    if (!outFile_bin)
    {
        std::string command;
        command = "mkdir DATA";
        std::system(command.c_str());
        outFile_bin.close();
        outFile_bin.open("./DATA/flow_field" + filename + ".bin", std::ios_base::out | std::ios::binary);
    }

    int32_t length;          // 用来计量string类型的长度
    std::string temp_string; // 用来临时存储string类型

    auto bin_write_int = [](std::ofstream &file, const std::string &parameter_name, const int &parameter_value)
    {
        int32_t length = parameter_name.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(parameter_name.c_str(), length);
        file.write((char *)&parameter_value, sizeof(int));
    };
    auto bin_write_double = [](std::ofstream &file, const std::string &parameter_name, const double &parameter_value)
    {
        int32_t length = parameter_name.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(parameter_name.c_str(), length);
        file.write((char *)&parameter_value, sizeof(double));
    };
    auto bin_write_string = [](std::ofstream &file, const std::string &parameter_name)
    {
        int32_t length = parameter_name.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(parameter_name.c_str(), length);
    };
    //--------------------------------------------------------------------------
    // 1、输出int类型的记录变量
    int32_t num_int = 1;
    outFile_bin.write((char *)&num_int, sizeof(int32_t));

    bin_write_int(outFile_bin, "Nstep", fld->par->GetInt("Nstep")); // 输出Nstep
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // 2、输出double类型的记录变量
    int32_t num_double = 6;
    outFile_bin.write((char *)&num_double, sizeof(int32_t));
    // 输出物理时间
    bin_write_double(outFile_bin, "Physic_Time", fld->par->GetDou("Physic_Time")); //  输出Physic_Time
    bin_write_double(outFile_bin, "Res_Ref_0", fld->par->GetDou("Res_Ref_0"));     //  输出Res_Ref_0
    bin_write_double(outFile_bin, "Res_Ref_1", fld->par->GetDou("Res_Ref_1"));     //  输出Res_Ref_1
    bin_write_double(outFile_bin, "Res_Ref_2", fld->par->GetDou("Res_Ref_2"));     //  输出Res_Ref_2
    bin_write_double(outFile_bin, "Res_Ref_3", fld->par->GetDou("Res_Ref_3"));     //  输出Res_Ref_3
    bin_write_double(outFile_bin, "Res_Ref_4", fld->par->GetDou("Res_Ref_4"));     //  输出Res_Ref_4
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // 3、输出data信息
    int32_t nblock = fld->num_blocks();
    // 一共有多少个网格块
    outFile_bin.write((char *)&nblock, sizeof(int32_t));

    int32_t multi_phy_num = fld->num_fields();
    int num = 0;
    for (int fid = 0; fid < multi_phy_num; fid++)
    {
        auto desc = fld->descriptor(fid);
        if (desc.name == "PV_" || desc.name == "Jac" || desc.name == "JDxi" ||
            desc.name == "JDet" || desc.name == "JDze")
            continue;
        num++;
    }
    // 一共有多少个物理场
    outFile_bin.write((char *)&num, sizeof(int32_t));

    // 正式输出data数据
    auto map_loc2int = [&](StaggerLocation &input_loc) -> int
    {
        if (input_loc == StaggerLocation::Cell)
            return 0;
        else if (input_loc == StaggerLocation::FaceXi)
            return 1;
        else if (input_loc == StaggerLocation::FaceEt)
            return 2;
        else if (input_loc == StaggerLocation::FaceZe)
            return 3;
        else if (input_loc == StaggerLocation::EdgeXi)
            return 4;
        else if (input_loc == StaggerLocation::EdgeEt)
            return 5;
        else if (input_loc == StaggerLocation::EdgeZe)
            return 6;
        else if (input_loc == StaggerLocation::Node)
            return 7;
        else
            return -1;
    };
    for (int fid = 0; fid < multi_phy_num; fid++)
    {
        auto desc = fld->descriptor(fid);
        if (desc.name == "PV_" || desc.name == "Jac" || desc.name == "JDxi" ||
            desc.name == "JDet" || desc.name == "JDze")
            continue;
        // 开始输出 FieldDescriptor
        bin_write_string(outFile_bin, desc.name);
        int kind = map_loc2int(desc.location);
        outFile_bin.write((char *)&kind, sizeof(int32_t));
        outFile_bin.write((char *)&desc.ncomp, sizeof(int32_t));
        outFile_bin.write((char *)&desc.nghost, sizeof(int32_t));

        for (int ib = 0; ib < nblock; ib++)
        {
            auto &field = fld->field(fid)[ib];
            int lo[3] = {field.get_lo().i,
                         field.get_lo().j,
                         field.get_lo().k};
            int hi[3] = {field.get_hi().i,
                         field.get_hi().j,
                         field.get_hi().k};
            outFile_bin.write((char *)&lo[0], sizeof(int32_t));
            outFile_bin.write((char *)&lo[1], sizeof(int32_t));
            outFile_bin.write((char *)&lo[2], sizeof(int32_t));
            outFile_bin.write((char *)&hi[0], sizeof(int32_t));
            outFile_bin.write((char *)&hi[1], sizeof(int32_t));
            outFile_bin.write((char *)&hi[2], sizeof(int32_t));
            for (int i = lo[0]; i < hi[0]; i++)
                for (int j = lo[1]; j < hi[1]; j++)
                    for (int k = lo[2]; k < hi[2]; k++)
                        for (int32_t ll = 0; ll < desc.ncomp; ll++)
                            outFile_bin.write((char *)&(field(i, j, k, ll)), sizeof(double));
        }
    }
    outFile_bin.close();
    //--------------------------------------------------------------------------
}

void Euler_Output::get_filename(Param *par)
{
    int32_t my_id = par->GetInt("myid");
    std::string _my_id_s;
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
    filename = _my_id_s;
}

void Euler_Output::plt_write_Headersection(const std::string &fieldname, std::vector<std::string> &var_name, FILE *file)
{
    constexpr auto magic_number{"#!TDV112"};
    fwrite(magic_number, 8, 1, file);
    constexpr int32_t byte_order{1};
    fwrite(&byte_order, 4, 1, file);
    constexpr int32_t file_type{0};
    fwrite(&file_type, 4, 1, file);

    const char *field_name = fieldname.c_str();
    plt_write_str(field_name, file);

    int32_t n_var = var_name.size();
    // std::cout << "  nvar   " << n_var << std::endl;
    fwrite(&n_var, 4, 1, file);

    for (const std::string &str : var_name)
    {
        const char *varname = str.c_str();
        plt_write_str(varname, file);
    }
};

void Euler_Output::plt_write_str(const char *str, FILE *file)
{
    int value = 0;
    while (*str != '\0')
    {
        value = static_cast<int>(*str);
        fwrite(&value, sizeof(int), 1, file);
        ++str;
    }
    constexpr char null_char = '\0';
    value = static_cast<int>(null_char);
    fwrite(&value, sizeof(int), 1, file);
};

void Euler_Output::plt_write_ZoneHeadersection(const std::string &Zone_name, int *mxyz, FILE *file)
{
    constexpr float zone_marker{299.0f};
    fwrite(&zone_marker, 4, 1, file);

    const char *Zonename = Zone_name.c_str();

    plt_write_str(Zonename, file);
    constexpr int32_t parent_zone{-1};
    fwrite(&parent_zone, 4, 1, file);
    constexpr int32_t strand_id{-2};
    fwrite(&strand_id, 4, 1, file);
    constexpr double solution_time{0.0};
    fwrite(&solution_time, 8, 1, file);
    constexpr int32_t zone_color{-1};
    fwrite(&zone_color, 4, 1, file);
    constexpr int32_t zone_type{0};
    fwrite(&zone_type, 4, 1, file);
    constexpr int32_t var_location{0};
    fwrite(&var_location, 4, 1, file);
    constexpr int32_t raw_face_neighbor{0};
    fwrite(&raw_face_neighbor, 4, 1, file);
    constexpr int32_t miscellaneous_face{0};
    fwrite(&miscellaneous_face, 4, 1, file);

    // For ordered zone, specify IMax, JMax, KMax
    fwrite(&mxyz[0], 4, 1, file);
    fwrite(&mxyz[1], 4, 1, file);
    fwrite(&mxyz[2], 4, 1, file);

    constexpr int32_t auxi_data{0};
    fwrite(&auxi_data, 4, 1, file);
};

void Euler_Output::plt_write_DataSection(int &n_var, FILE *file)
{
    constexpr float zone_marker{299.0f};
    fwrite(&zone_marker, 4, 1, file);
    constexpr int32_t data_format{1};

    for (int l = 0; l < n_var; ++l)
        fwrite(&data_format, 4, 1, file);
    constexpr int32_t passive_var{0};
    fwrite(&passive_var, 4, 1, file);
    constexpr int32_t shared_var{0};
    fwrite(&shared_var, 4, 1, file);
    constexpr int32_t shared_connect{-1};
    fwrite(&shared_connect, 4, 1, file);
};
