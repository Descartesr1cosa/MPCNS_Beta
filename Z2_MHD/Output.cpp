#include "Output.h"

void MHD_Output::output_bin_field(const std::vector<std::string> &field_names)
{
    Param *par = fld->par;
    //----------------------------------------------------------------------
    // 打开/准备文件
    std::ofstream outFile_bin("./DATA/flow_field" + filename + ".bin",
                              std::ios_base::out | std::ios::binary);
    if (!outFile_bin)
    {
        std::string command = "mkdir DATA";
        std::system(command.c_str());
        outFile_bin.close();
        outFile_bin.open("./DATA/flow_field" + filename + ".bin",
                         std::ios_base::out | std::ios::binary);
    }

    auto bin_write_int = [](std::ofstream &file,
                            const std::string &parameter_name,
                            const int &parameter_value)
    {
        int32_t length = parameter_name.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(parameter_name.c_str(), length);
        file.write((char *)&parameter_value, sizeof(int));
    };
    auto bin_write_double = [](std::ofstream &file,
                               const std::string &parameter_name,
                               const double &parameter_value)
    {
        int32_t length = parameter_name.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(parameter_name.c_str(), length);
        file.write((char *)&parameter_value, sizeof(double));
    };
    auto bin_write_string = [](std::ofstream &file,
                               const std::string &parameter_name)
    {
        int32_t length = parameter_name.size() + 1;
        file.write((char *)&length, sizeof(int32_t));
        file.write(parameter_name.c_str(), length);
    };

    //----------------------------------------------------------------------
    // 1. 输出 int 记录变量
    int32_t num_int = 1;
    outFile_bin.write((char *)&num_int, sizeof(int32_t));

    bin_write_int(outFile_bin, "Nstep", fld->par->GetInt("Nstep"));
    //----------------------------------------------------------------------
    // 2. 输出 double 记录变量（保持原格式）
    int32_t num_double = 9;
    outFile_bin.write((char *)&num_double, sizeof(int32_t));
    bin_write_double(outFile_bin, "Physic_Time", fld->par->GetDou("Physic_Time"));
    bin_write_double(outFile_bin, "Res_Ref_0", fld->par->GetDou("Res_Ref_0"));
    bin_write_double(outFile_bin, "Res_Ref_1", fld->par->GetDou("Res_Ref_1"));
    bin_write_double(outFile_bin, "Res_Ref_2", fld->par->GetDou("Res_Ref_2"));
    bin_write_double(outFile_bin, "Res_Ref_3", fld->par->GetDou("Res_Ref_3"));
    bin_write_double(outFile_bin, "Res_Ref_4", fld->par->GetDou("Res_Ref_4"));
    bin_write_double(outFile_bin, "Res_Ref_5", fld->par->GetDou("Res_Ref_5"));
    bin_write_double(outFile_bin, "Res_Ref_6", fld->par->GetDou("Res_Ref_6"));
    bin_write_double(outFile_bin, "Res_Ref_7", fld->par->GetDou("Res_Ref_7"));
    //----------------------------------------------------------------------
    // 3. 输出 data 信息

    // 网格块数量
    int32_t nblock = fld->num_blocks();
    outFile_bin.write((char *)&nblock, sizeof(int32_t));

    // 先把 field_names 转成 fid 列表（有可能有名字找不到，顺便跳过）
    std::vector<int> fids;
    fids.reserve(field_names.size());

    for (auto &name : field_names)
    {
        try
        {
            int fid = fld->field_id(name);
            fids.push_back(fid);
        }
        catch (const std::out_of_range &)
        {
            // 找不到这个 field，就跳过。可以在 myid==0 时打印一下 warning
            if (par->GetInt("myid") == 0)
            {
                std::cout << "[MHD_Output::output_bin_field] Warning: field \""
                          << name << "\" not found, skip.\n";
            }
        }
    }

    int32_t num = static_cast<int32_t>(fids.size());
    // 一共有多少个物理场
    outFile_bin.write((char *)&num, sizeof(int32_t));

    // 位置映射
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

    // 正式输出每个选中的 field
    for (int idx = 0; idx < num; ++idx)
    {
        int fid = fids[idx];
        auto desc = fld->descriptor(fid);

        // 1) 输出 FieldDescriptor 信息
        bin_write_string(outFile_bin, desc.name);
        int kind = map_loc2int(desc.location);
        outFile_bin.write((char *)&kind, sizeof(int32_t));
        outFile_bin.write((char *)&desc.ncomp, sizeof(int32_t));
        outFile_bin.write((char *)&desc.nghost, sizeof(int32_t));

        // 2) 输出每个 block 的数据
        for (int ib = 0; ib < nblock; ++ib)
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
                        {
                            double val = field(i, j, k, ll);
                            outFile_bin.write((char *)&val, sizeof(double));
                        }
        }
    }

    outFile_bin.close();
    //----------------------------------------------------------------------
}

void MHD_Output::output_plt_field(const std::vector<std::string> &var_list)
{
    Grid *grd = fld->grd;
    Param *par = fld->par;

    //============================== 1. 准备文件 ==============================
    std::string path = "./DATA/flow_field" + filename + ".plt";
    const char *filenameall = path.c_str();

    {
        std::ofstream test(path, std::ios::out);
        if (!test)
        {
            std::string cmd = "mkdir DATA";
            std::system(cmd.c_str());
        }
    }

    FILE *fp = std::fopen(filenameall, "wb");
    if (!fp)
    {
        std::printf("Failed to open file %s\n", filenameall);
        std::exit(-1);
    }

    auto write_int32 = [&](int32_t v)
    { std::fwrite(&v, 4, 1, fp); };
    auto write_float = [&](float v)
    { std::fwrite(&v, 4, 1, fp); };
    auto write_double = [&](double v)
    { std::fwrite(&v, 8, 1, fp); };

    //============================== 2. 变量名列表 ==============================
    std::vector<std::string> var_names;
    var_names.reserve(3 + var_list.size());
    var_names.push_back("x");
    var_names.push_back("y");
    var_names.push_back("z");
    for (auto &s : var_list)
        var_names.push_back(s);

    const int n_var = static_cast<int>(var_names.size());

    // 写 Header
    plt_write_Headersection("flow_field", var_names, fp);

    //============================== 3. Zone Header ==============================
    for (int ib = 0; ib < grd->nblock; ++ib)
    {
        const Block &blk = grd->grids(ib);
        if (blk.block_name != "Fluid")
            continue;

        int mxyz[3] = {blk.mx + 1, blk.my + 1, blk.mz + 1};
        std::string zone_name = "Zone" + std::to_string(ib);
        plt_write_ZoneHeadersection(zone_name, mxyz, fp);
    }

    // Header结束marker
    write_float(357.0f);

    //============================== 4. 找到需要的 Field ==============================
    const int fid_U = fld->field_id("U_");
    const int fid_PV = fld->field_id("PV_");

    const FieldDescriptor &descU = fld->descriptor(fid_U);
    if (descU.location != StaggerLocation::Cell || descU.ncomp < 5)
    {
        std::cout << "Fatal Error: \"U_\" 必须是 Cell-centered 且至少 5 分量（rho, rho*u, rho*v, rho*w, E）。\n";
        std::exit(-1);
    }

    // 看看有没有请求 Bx/By/Bz，如果有就需要 B_cell
    bool need_B = false;
    for (auto &s : var_list)
        if (s == "Bx" || s == "By" || s == "Bz")
            need_B = true;

    int fid_Bc = -1;
    if (need_B)
    {
        fid_Bc = fld->field_id("B_cell");
        const FieldDescriptor &descBc = fld->descriptor(fid_Bc);
        if (descBc.location != StaggerLocation::Cell || descBc.ncomp < 3)
        {
            std::cout << "Fatal Error: Output Bx/By/Bz, but \"B_cell\" is not Cell-centered 3 component \n";
            std::cout << descBc.ncomp << "\tfid_Bc=\t" << fid_Bc << std::endl;
            std::exit(-1);
        }
    }

    //============================== 5. Data Section ==============================
    for (int iblock = 0; iblock < grd->nblock; ++iblock)
    {
        Block &blk = grd->grids(iblock);
        if (blk.block_name != "Fluid")
            continue;

        const int Ni = blk.mx;
        const int Nj = blk.my;
        const int Nk = blk.mz;
        const int Ni_node = Ni + 1;
        const int Nj_node = Nj + 1;
        const int Nk_node = Nk + 1;

        FieldBlock &Ublk = fld->field(fid_U, iblock);
        FieldBlock &PVblk = fld->field(fid_PV, iblock);
        FieldBlock *Bcblk_ptr = nullptr;
        if (need_B)
            Bcblk_ptr = &fld->field(fid_Bc, iblock);

        // zone data header
        plt_write_DataSection(const_cast<int &>(n_var), fp);

        // ---------- 5.1 先写每个变量的 min/max ----------
        for (int v = 0; v < n_var; ++v)
        {
            bool is_nodal = (v < 3);

            double minv = 0.0;
            double maxv = 0.0;
            bool first = true;

            if (is_nodal)
            {
                for (int k = 0; k < Nk_node; ++k)
                    for (int j = 0; j < Nj_node; ++j)
                        for (int i = 0; i < Ni_node; ++i)
                        {
                            double val = 0.0;
                            if (v == 0)
                                val = blk.x(i, j, k);
                            else if (v == 1)
                                val = blk.y(i, j, k);
                            else if (v == 2)
                                val = blk.z(i, j, k);

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
            else
            {
                const std::string &name = var_names[v];
                for (int kc = 0; kc < Nk; ++kc)
                    for (int jc = 0; jc < Nj; ++jc)
                        for (int ic = 0; ic < Ni; ++ic)
                        {
                            double val = 0.0;

                            // 根据名字映射到底层 field
                            if (name == "rho")
                            {
                                val = Ublk(ic, jc, kc, 0);
                            }
                            else if (name == "u")
                            {
                                val = PVblk(ic, jc, kc, 0);
                            }
                            else if (name == "v")
                            {
                                val = PVblk(ic, jc, kc, 1);
                            }
                            else if (name == "w")
                            {
                                val = PVblk(ic, jc, kc, 2);
                            }
                            else if (name == "p")
                            {
                                val = PVblk(ic, jc, kc, 3);
                            }
                            else if (name == "Bx" && need_B)
                            {
                                val = (*Bcblk_ptr)(ic, jc, kc, 0);
                            }
                            else if (name == "By" && need_B)
                            {
                                val = (*Bcblk_ptr)(ic, jc, kc, 1);
                            }
                            else if (name == "Bz" && need_B)
                            {
                                val = (*Bcblk_ptr)(ic, jc, kc, 2);
                            }
                            else
                            {
                                // 未支持的名字，全部输出 0
                                val = 0.0;
                            }

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

            write_double(minv);
            write_double(maxv);
        }

        // ---------- 5.2 再写具体数据 ----------
        for (int v = 0; v < n_var; ++v)
        {
            bool is_nodal = (v < 3);

            if (is_nodal)
            {
                for (int k = 0; k < Nk_node; ++k)
                    for (int j = 0; j < Nj_node; ++j)
                        for (int i = 0; i < Ni_node; ++i)
                        {
                            double val = 0.0;
                            if (v == 0)
                                val = blk.x(i, j, k);
                            else if (v == 1)
                                val = blk.y(i, j, k);
                            else if (v == 2)
                                val = blk.z(i, j, k);

                            write_float(static_cast<float>(val));
                        }
            }
            else
            {
                int IMax = Ni + 1;
                int JMax = Nj + 1;
                int KMax = Nk + 1;
                const std::string &name = var_names[v];

                for (int k = 0; k < KMax - 1; ++k)
                {
                    for (int j = 0; j < JMax; ++j)
                    {
                        for (int i = 0; i < IMax; ++i)
                        {
                            double val = 0.0;
                            if (i < Ni && j < Nj)
                            {
                                int ic = i;
                                int jc = j;
                                int kc = k;

                                if (name == "rho")
                                {
                                    val = Ublk(ic, jc, kc, 0);
                                }
                                else if (name == "u")
                                {
                                    val = PVblk(ic, jc, kc, 0);
                                }
                                else if (name == "v")
                                {
                                    val = PVblk(ic, jc, kc, 1);
                                }
                                else if (name == "w")
                                {
                                    val = PVblk(ic, jc, kc, 2);
                                }
                                else if (name == "p")
                                {
                                    val = PVblk(ic, jc, kc, 3);
                                }
                                else if (name == "Bx" && need_B)
                                {
                                    val = (*Bcblk_ptr)(ic, jc, kc, 0);
                                }
                                else if (name == "By" && need_B)
                                {
                                    val = (*Bcblk_ptr)(ic, jc, kc, 1);
                                }
                                else if (name == "Bz" && need_B)
                                {
                                    val = (*Bcblk_ptr)(ic, jc, kc, 2);
                                }
                                else
                                {
                                    val = 0.0;
                                }
                            }
                            write_float(static_cast<float>(val));
                        }
                    }
                }
            }
        }
    }

    std::fclose(fp);
}

void MHD_Output::output_plt_cell_field(const std::vector<std::string> &var_list)
{
    Grid *grd = fld->grd;
    Param *par = fld->par;

    //============================== 1. 准备文件 ==============================
    std::string path = "./DATA/flow_field" + filename + ".plt";
    const char *filenameall = path.c_str();

    {
        std::ofstream test(path, std::ios::out);
        if (!test)
        {
            std::string cmd = "mkdir DATA";
            std::system(cmd.c_str());
        }
    }

    FILE *fp = std::fopen(filenameall, "wb");
    if (!fp)
    {
        std::printf("Failed to open file %s\n", filenameall);
        std::exit(-1);
    }

    auto write_int32 = [&](int32_t v)
    { std::fwrite(&v, 4, 1, fp); };
    auto write_float = [&](float v)
    { std::fwrite(&v, 4, 1, fp); };
    auto write_double = [&](double v)
    { std::fwrite(&v, 8, 1, fp); };

    auto my_plt_write_ZoneHeadersection = [](const std::string &Zone_name,
                                             int *mxyz,
                                             FILE *file, int n_var_plt)
    {
        auto write_int32 = [&](int32_t v)
        { std::fwrite(&v, 4, 1, file); };
        auto write_float = [&](float v)
        { std::fwrite(&v, 4, 1, file); };
        auto write_double = [&](double v)
        { std::fwrite(&v, 8, 1, file); };
        auto write_string = [&](const char *str, FILE *file)
        {
            int value = 0;
            const char *p = str;
            while (*p != '\0')
            {
                value = static_cast<int>(*p);
                std::fwrite(&value, sizeof(int), 1, file);
                ++p;
            }
            value = static_cast<int>('\0');
            std::fwrite(&value, sizeof(int), 1, file);
        };
        // Tecplot zone marker
        write_float(299.0f);

        // Zone 名字
        write_string(Zone_name.c_str(), file);

        // parent zone id
        int32_t parentZone = -1;
        write_int32(parentZone);
        // strand id
        int32_t strandId = -2;
        write_int32(strandId);
        // solution time
        double solTime = 0.0;
        write_double(solTime);
        // zone color
        write_int32(-1);
        // zone type: 0 = ORDERED
        write_int32(0);

        // 是否指定每个变量的位置：1=是
        int32_t specify_var_location = 1;
        write_int32(specify_var_location);

        // 每个变量的位置：0= nodal, 1 = cell-centered
        // 约定：全部 nodal
        for (int v = 0; v < n_var_plt; ++v)
        {
            // int32_t loc = (v < 3) ? 0 : 1;
            int32_t loc = 0;
            write_int32(loc);
        }

        // Are raw local 1-to-1 face neighbors supplied? 0=FALSE
        write_int32(0);
        // Number of miscellaneous user-defined face neighbor connections: 0
        write_int32(0);

        // IMax, JMax, KMax
        write_int32(mxyz[0]);
        write_int32(mxyz[1]);
        write_int32(mxyz[2]);

        // AuxData count = 0
        write_int32(0);
    };

    //============================== 2. 变量名列表 ==============================
    std::vector<std::string> var_names;
    var_names.reserve(3 + var_list.size());
    var_names.push_back("x");
    var_names.push_back("y");
    var_names.push_back("z");
    for (auto &s : var_list)
        var_names.push_back(s);

    const int n_var = static_cast<int>(var_names.size());

    // 写 Header
    plt_write_Headersection("flow_field", var_names, fp);
#if if_Debug_output_ngg == 1
    int temp_ngg = 2;
#endif
    //============================== 3. Zone Header ==============================
    for (int ib = 0; ib < grd->nblock; ++ib)
    {
        const Block &blk = grd->grids(ib);
        if (blk.block_name != "Fluid")
            continue;
#if if_Debug_output_ngg == 0
        int mxyz[3] = {blk.mx, blk.my, blk.mz};
#elif if_Debug_output_ngg == 1
        int mxyz[3] = {blk.mx + 2 * temp_ngg, blk.my + 2 * temp_ngg, blk.mz + 2 * temp_ngg};
#endif
        std::string zone_name = "Zone" + std::to_string(ib);
        // plt_write_ZoneHeadersection(zone_name, mxyz, fp);
        my_plt_write_ZoneHeadersection(zone_name, mxyz, fp, n_var_plt);
    }

    // Header结束marker
    write_float(357.0f);

    //============================== 4. 找到需要的 Field ==============================
    const int fid_U = fld->field_id("U_");
    const int fid_PV = fld->field_id("PV_");

    const FieldDescriptor &descU = fld->descriptor(fid_U);
    if (descU.location != StaggerLocation::Cell || descU.ncomp < 5)
    {
        std::cout << "Fatal Error: \"U_\" 必须是 Cell-centered 且至少 5 分量（rho, rho*u, rho*v, rho*w, E）。\n";
        std::exit(-1);
    }

    // 看看有没有请求 Bx/By/Bz，如果有就需要 B_cell
    bool need_B = false;
    for (auto &s : var_list)
        if (s == "Bx" || s == "By" || s == "Bz")
            need_B = true;

    int fid_Bc = -1;
    if (need_B)
    {
        fid_Bc = fld->field_id("B_cell");
        const FieldDescriptor &descBc = fld->descriptor(fid_Bc);
        if (descBc.location != StaggerLocation::Cell || descBc.ncomp < 3)
        {
            std::cout << "Fatal Error: Output Bx/By/Bz, but \"B_cell\" is not Cell-centered 3 component \n";
            std::cout << descBc.ncomp << "\tfid_Bc=\t" << fid_Bc << std::endl;
            std::exit(-1);
        }
    }

    //============================== 5. Data Section ==============================
    for (int iblock = 0; iblock < grd->nblock; ++iblock)
    {
        Block &blk = grd->grids(iblock);
        if (blk.block_name != "Fluid")
            continue;

        const int Ni = blk.mx;
        const int Nj = blk.my;
        const int Nk = blk.mz;

        FieldBlock &Ublk = fld->field(fid_U, iblock);
        FieldBlock &PVblk = fld->field(fid_PV, iblock);
        FieldBlock *Bcblk_ptr = nullptr;
        if (need_B)
            Bcblk_ptr = &fld->field(fid_Bc, iblock);

        // zone data header
        plt_write_DataSection(const_cast<int &>(n_var), fp);

        // ---------- 5.1 先写每个变量的 min/max ----------
        for (int v = 0; v < n_var; ++v)
        {
            double minv = 0.0;
            double maxv = 0.0;
            bool first = true;

            const std::string &name = var_names[v];
#if if_Debug_output_ngg == 0
            for (int kc = 0; kc < Nk; ++kc)
                for (int jc = 0; jc < Nj; ++jc)
                    for (int ic = 0; ic < Ni; ++ic)
#elif if_Debug_output_ngg == 1
            for (int kc = -temp_ngg; kc < Nk + temp_ngg; ++kc)
                for (int jc = -temp_ngg; jc < Nj + temp_ngg; ++jc)
                    for (int ic = -temp_ngg; ic < Ni + temp_ngg; ++ic)
#endif
                    {
                        double val = 0.0;

                        // 根据名字映射到底层 field
                        if (name == "x")
                            val = blk.dual_x(ic + 1, jc + 1, kc + 1);
                        else if (name == "y")
                            val = blk.dual_y(ic + 1, jc + 1, kc + 1);
                        else if (name == "z")
                            val = blk.dual_z(ic + 1, jc + 1, kc + 1);
                        else if (name == "rho")
                        {
                            val = Ublk(ic, jc, kc, 0);
                        }
                        else if (name == "u")
                        {
                            val = PVblk(ic, jc, kc, 0);
                        }
                        else if (name == "v")
                        {
                            val = PVblk(ic, jc, kc, 1);
                        }
                        else if (name == "w")
                        {
                            val = PVblk(ic, jc, kc, 2);
                        }
                        else if (name == "p")
                        {
                            val = PVblk(ic, jc, kc, 3);
                        }
                        else if (name == "Bx" && need_B)
                        {
                            val = (*Bcblk_ptr)(ic, jc, kc, 0);
                        }
                        else if (name == "By" && need_B)
                        {
                            val = (*Bcblk_ptr)(ic, jc, kc, 1);
                        }
                        else if (name == "Bz" && need_B)
                        {
                            val = (*Bcblk_ptr)(ic, jc, kc, 2);
                        }
                        else
                        {
                            // 未支持的名字，全部输出 0
                            val = 0.0;
                        }

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

            write_double(minv);
            write_double(maxv);
        }

        // ---------- 5.2 再写具体数据 ----------
        for (int v = 0; v < n_var; ++v)
        {

            const std::string &name = var_names[v];
#if if_Debug_output_ngg == 0
            for (int kc = 0; kc < Nk; ++kc)
                for (int jc = 0; jc < Nj; ++jc)
                    for (int ic = 0; ic < Ni; ++ic)
#elif if_Debug_output_ngg == 1
            for (int kc = -temp_ngg; kc < Nk + temp_ngg; ++kc)
                for (int jc = -temp_ngg; jc < Nj + temp_ngg; ++jc)
                    for (int ic = -temp_ngg; ic < Ni + temp_ngg; ++ic)
#endif
                    {
                        double val = 0.0;

                        // 根据名字映射到底层 field
                        if (name == "x")
                            val = blk.dual_x(ic + 1, jc + 1, kc + 1);
                        else if (name == "y")
                            val = blk.dual_y(ic + 1, jc + 1, kc + 1);
                        else if (name == "z")
                            val = blk.dual_z(ic + 1, jc + 1, kc + 1);
                        else if (name == "rho")
                        {
                            val = Ublk(ic, jc, kc, 0);
                        }
                        else if (name == "u")
                        {
                            val = PVblk(ic, jc, kc, 0);
                        }
                        else if (name == "v")
                        {
                            val = PVblk(ic, jc, kc, 1);
                        }
                        else if (name == "w")
                        {
                            val = PVblk(ic, jc, kc, 2);
                        }
                        else if (name == "p")
                        {
                            val = PVblk(ic, jc, kc, 3);
                        }
                        else if (name == "Bx" && need_B)
                        {
                            val = (*Bcblk_ptr)(ic, jc, kc, 0);
                        }
                        else if (name == "By" && need_B)
                        {
                            val = (*Bcblk_ptr)(ic, jc, kc, 1);
                        }
                        else if (name == "Bz" && need_B)
                        {
                            val = (*Bcblk_ptr)(ic, jc, kc, 2);
                        }
                        else
                        {
                            // 未支持的名字，全部输出 0
                            val = 0.0;
                        }

                        write_float(static_cast<float>(val));
                    }
        }
    }

    std::fclose(fp);
}

void MHD_Output::get_filename(Param *par)
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

// 写 Tecplot string：字符逐个写成 int，最后写 '\0'
void MHD_Output::plt_write_str(const char *str, FILE *file)
{
    int value = 0;
    const char *p = str;
    while (*p != '\0')
    {
        value = static_cast<int>(*p);
        std::fwrite(&value, sizeof(int), 1, file);
        ++p;
    }
    value = static_cast<int>('\0');
    std::fwrite(&value, sizeof(int), 1, file);
}

// Header 区：magic + byte order + file type + title + 变量名
void MHD_Output::plt_write_Headersection(const std::string &fieldname,
                                         std::vector<std::string> &var_name,
                                         FILE *file)
{
    auto write_int32 = [&](int32_t v)
    { std::fwrite(&v, 4, 1, file); };

    // magic number "#!TDV112"
    const char magic[8] = {'#', '!', 'T', 'D', 'V', '1', '1', '2'};
    std::fwrite(magic, 8, 1, file);

    // byte order: 1 = little endian
    int32_t byte_order = 1;
    write_int32(byte_order);

    // file type: 0 = FULL
    int32_t file_type = 0;
    write_int32(file_type);

    // title：用 fieldname 做 title
    plt_write_str(fieldname.c_str(), file);

    // 变量名
    int32_t n_var = static_cast<int32_t>(var_name.size());
    n_var_plt = n_var; // 记住变量个数，ZoneHeader会用
    write_int32(n_var);

    for (auto &name : var_name)
        plt_write_str(name.c_str(), file);
}

// Zone header：每个 block 一个
// mxyz = {IMax, JMax, KMax} = {mx+1, my+1, mz+1}
void MHD_Output::plt_write_ZoneHeadersection(const std::string &Zone_name,
                                             int *mxyz,
                                             FILE *file)
{
    auto write_int32 = [&](int32_t v)
    { std::fwrite(&v, 4, 1, file); };
    auto write_float = [&](float v)
    { std::fwrite(&v, 4, 1, file); };
    auto write_double = [&](double v)
    { std::fwrite(&v, 8, 1, file); };

    Param *par = fld->par;

    // Tecplot zone marker
    write_float(299.0f);

    // Zone 名字
    plt_write_str(Zone_name.c_str(), file);

    // parent zone id
    int32_t parentZone = -1;
    write_int32(parentZone);
    // strand id
    int32_t strandId = -2;
    write_int32(strandId);
    // solution time
    double solTime = par->GetDou("Physic_Time");
    write_double(solTime);
    // zone color
    write_int32(-1);
    // zone type: 0 = ORDERED
    write_int32(0);

    // 是否指定每个变量的位置：1=是
    int32_t specify_var_location = 1;
    write_int32(specify_var_location);

    // 每个变量的位置：0= nodal, 1 = cell-centered
    // 约定：前3个 (x,y,z) 节点，其余全部 cell-centered
    for (int v = 0; v < n_var_plt; ++v)
    {
        int32_t loc = (v < 3) ? 0 : 1;
        write_int32(loc);
    }

    // Are raw local 1-to-1 face neighbors supplied? 0=FALSE
    write_int32(0);
    // Number of miscellaneous user-defined face neighbor connections: 0
    write_int32(0);

    // IMax, JMax, KMax
    write_int32(mxyz[0]);
    write_int32(mxyz[1]);
    write_int32(mxyz[2]);

    // AuxData count = 0
    write_int32(0);
}

// Zone 数据头（DataSection 的前半部分）
void MHD_Output::plt_write_DataSection(int &n_var, FILE *file)
{
    auto write_int32 = [&](int32_t v)
    { std::fwrite(&v, 4, 1, file); };
    auto write_float = [&](float v)
    { std::fwrite(&v, 4, 1, file); };

    // zone marker
    write_float(299.0f);

    // 每个变量的 data type: 1 = float
    int32_t data_type = 1;
    for (int v = 0; v < n_var; ++v)
        write_int32(data_type);

    // passive variables: 0 = 无
    write_int32(0);
    // shared variables: 0 = 无
    write_int32(0);
    // shared connectivity: -1 = 无
    write_int32(-1);
}