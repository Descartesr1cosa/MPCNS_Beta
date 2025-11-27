#include "2_topology/2_MPCNS_Topology.h"
#include "0_basic/MPI_WRAPPER.h"

namespace TOPO
{
    Topology build_topology(Grid &grid, int my_rank, int dimension)
    {
        Topology topo;

        //=======================================================================
        // Physic
        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &blk = grid.grids(ib);

            for (const auto &phy : blk.physical_bc)
            {
                PhysicalPatch p;
                p.this_rank = my_rank;
                p.this_block = phy.this_block_num;
                p.this_block_name = phy.this_block_name;
                p.bc_id = phy.boundary_num;
                p.bc_name = phy.boundary_name;
                p.direction = phy.direction;
                p.raw = &phy;

                node_box_from_subsup(phy.sub, phy.sup, p.this_box_node);

                topo.physical_patches.push_back(p);
            }
        }
        //=======================================================================

        //=======================================================================
        // Inner
        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &blk = grid.grids(ib);

            for (const auto &inner : blk.inner_bc)
            {
                InterfacePatch interface;

                interface.kind = PatchKind::Inner;

                interface.this_rank = my_rank;
                interface.nb_rank = my_rank;

                interface.this_block = inner.this_block_num;
                interface.nb_block = inner.tar_block_num;

                interface.this_block_name = inner.this_block_name;
                interface.nb_block_name = inner.target_block_name;

                node_box_from_subsup(inner.sub, inner.sup, interface.this_box_node);
                node_box_from_subsup(inner.tar_sub, inner.tar_sup, interface.nb_box_node);

                //=======================================================================
                // 获取IndexTransform
                //-------------------------------------------------------------
                // 根据映射Transform -> mid state -> inver_Transform -> tar
                // 获取perm
                int inver_Transform[3];
                for (int d = 0; d < 3; ++d)
                    inver_Transform[inner.tar_Transform[d]] = d;
                for (int d = 0; d < 3; ++d)
                    interface.trans.perm[d] = inver_Transform[inner.Transform[d]];
                //-------------------------------------------------------------
                // 获取sign[3]
                for (int d = 0; d < 3; ++d)
                {
                    int my_sub, my_sup, tar_sub, tar_sup;
                    my_sub = abs(inner.sub[d]);
                    my_sup = abs(inner.sup[d]);
                    tar_sub = abs(inner.tar_sub[interface.trans.perm[d]]);
                    tar_sup = abs(inner.tar_sup[interface.trans.perm[d]]);

                    if (my_sub == my_sup && d < dimension)
                    {
                        if (tar_sub != tar_sup)
                        {
                            // For Check Security
                            std::cout << "\tDirection is not right when building inner interface\n";
                            exit(-1);
                        }

                        interface.trans.sign[d] = (inner.direction * inner.tar_direction > 0) ? -1 : 1;
                        continue;
                    }

                    if (dimension == 2 && d == 2)
                    {
                        interface.trans.sign[d] = 0;
                        interface.trans.offset.k = 0;
                        continue;
                    }

                    if (abs(my_sup - my_sub) != abs(tar_sup - tar_sub))
                    {
                        // For Check Security
                        std::cout << "\tDirection is not right when building inner interface\t" << d << std::endl;
                        exit(-1);
                    }
                    else
                        interface.trans.sign[d] = ((my_sup - my_sub) * (tar_sup - tar_sub) > 0) ? 1 : -1;
                }
                //-------------------------------------------------------------
                // 获取offset
                int offset[3] = {0, 0, 0};
                for (int d = 0; d < 3; ++d)
                {
                    int my_sub, tar_sub;
                    my_sub = abs(inner.sub[d]);
                    tar_sub = abs(inner.tar_sub[interface.trans.perm[d]]);

                    offset[d] = -interface.trans.sign[d] * my_sub + tar_sub;
                }
                interface.trans.offset = {offset[0], offset[1], offset[2]};
                //-------------------------------------------------------------
                topo.inner_patches.push_back(interface);
            }
        }
        //=======================================================================

        //=======================================================================
        //  Parallel
        std::vector<std::vector<Parallel_Boundary>> parallel_bc_all;
        // Read All para ** .txt

        {
            int rank_num;
            PARALLEL::mpi_size(&rank_num);
            // 清空并按 rank_num 分配
            parallel_bc_all.clear();
            parallel_bc_all.resize(rank_num);

            for (int r = 0; r < rank_num; ++r)
            {
                std::string _my_id_s;
                if (r < 10)
                {
                    _my_id_s = "   " + std::to_string(r);
                }
                else if (r < 100)
                {
                    _my_id_s = "  " + std::to_string(r);
                }
                else if (r < 1000)
                {
                    _my_id_s = " " + std::to_string(r);
                }
                else // 这说明并行进程数不得超过9999
                {
                    _my_id_s = std::to_string(r);
                }

                std::string filename = "./CASE/geometry/boundary_condition/parallel" + _my_id_s + ".txt";

                std::ifstream grdfile(filename, std::ios_base::in);
                if (!grdfile.is_open())
                {
                    std::cout << "ERROR: cannot open " << filename << " when Read_All_Para(), myid = "
                              << my_rank << std::endl;
                    exit(-1);
                }

                int blk_num_file = 0;
                grdfile >> blk_num_file; // 该 rank 的块数

                // 注意：blk_num_file 不一定等于本 rank 的 nblock（因为这是 rank r 的文件）
                // 我们这里只是“搬运” parallel 信息，不需要和本地 nblock 比较

                int num_parallel_face = 0;
                grdfile >> num_parallel_face; // 该 rank 总共有多少个通信面（暂时用不到）

                // 每个块的 parallel 面数
                std::vector<int> nface(blk_num_file);
                for (int izone = 0; izone < blk_num_file; ++izone)
                    grdfile >> nface[izone];

                std::string dummy;
                std::getline(grdfile, dummy); // 读掉一行尾巴
                std::getline(grdfile, dummy); // 注释行

                // 依次读每个 block 的 parallel 面
                for (int izone = 0; izone < blk_num_file; ++izone)
                {
                    int iface_num = nface[izone];

                    for (int iiface = 0; iiface < iface_num; ++iiface)
                    {
                        Parallel_Boundary pbc;
                        int pointst[3], pointed[3];
                        int srid, sflag, rflag;
                        std::string nameed;

                        grdfile >> pointst[0] >> pointed[0] >> pointst[1] >> pointed[1] >> pointst[2] >> pointed[2] >> srid >> sflag >> rflag;
                        grdfile >> nameed;

                        // 这是文件中写的：本块在 rank r 上，目标是 rank srid
                        pbc.this_myid = r;
                        pbc.tar_myid = srid;
                        pbc.this_block_num = izone;

                        for (int i = 0; i < 3; ++i)
                        {
                            pbc.sub[i] = pointst[i];
                            pbc.sup[i] = pointed[i];
                        }

                        pbc.send_flag = sflag;
                        pbc.rece_flag = rflag;

                        // 块名暂时只能从文件读到 target_block_name，
                        // this_block_name 在原来的 Grid 里是从 Block 取的，这里并不知道，
                        // 可以先留空或用一个占位符
                        pbc.this_block_name = ""; // 如果需要可以后续通过 block_name 映射补
                        pbc.target_block_name = nameed;

                        // 和 Grid::Read_Parallel_Boundary 一样，做预处理
                        pbc.Pre_process(dimension);

                        // 存入 parallel_bc_all[r]
                        parallel_bc_all[r].push_back(pbc);
                    }
                }

                grdfile.close();
            }

            if (my_rank == 0)
                std::cout << "\t--> Building TOPO Para Interface: Read all parallel_*.txt successfully!\n";
        }

        // Building parallel_patches
        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &blk = grid.grids(ib);

            for (const auto &para : blk.parallel_bc)
            {
                //-------------------------------------------------------------
                // Find tar_para
                Parallel_Boundary *tar_para_pointer = nullptr;
                bool if_find = false;
                for (auto &tar_para_temp : parallel_bc_all[para.tar_myid])
                {
                    if (tar_para_temp.rece_flag == para.send_flag && tar_para_temp.send_flag == para.rece_flag)
                    {
                        tar_para_pointer = &tar_para_temp;
                        if_find = true;
                        break;
                    }
                }
                if (!if_find)
                {
                    std::cout << "FATAL Error, Can not find corresponding PARA info!\n";
                    exit(-1);
                }
                Parallel_Boundary &tar_para = *tar_para_pointer;
                tar_para_pointer = nullptr;
                //-------------------------------------------------------------

                InterfacePatch interface;

                interface.kind = PatchKind::Parallel;

                interface.this_rank = my_rank;
                interface.nb_rank = para.tar_myid;

                interface.this_block = para.this_block_num;
                interface.nb_block = tar_para.this_block_num;

                interface.this_block_name = para.this_block_name;
                interface.nb_block_name = para.target_block_name;

                node_box_from_subsup(para.sub, para.sup, interface.this_box_node);
                node_box_from_subsup(tar_para.sub, tar_para.sup, interface.nb_box_node);

                //=======================================================================
                // 获取IndexTransform
                //-------------------------------------------------------------
                // 根据映射Transform -> mid state -> inver_Transform -> tar
                // 获取perm
                int inver_Transform[3];
                for (int d = 0; d < 3; ++d)
                    inver_Transform[tar_para.Transform[d]] = d;
                for (int d = 0; d < 3; ++d)
                    interface.trans.perm[d] = inver_Transform[para.Transform[d]];
                //-------------------------------------------------------------
                // 获取sign[3]
                for (int d = 0; d < 3; ++d)
                {
                    int my_sub, my_sup, tar_sub, tar_sup;
                    my_sub = abs(para.sub[d]);
                    my_sup = abs(para.sup[d]);
                    tar_sub = abs(tar_para.sub[interface.trans.perm[d]]);
                    tar_sup = abs(tar_para.sup[interface.trans.perm[d]]);

                    if (my_sub == my_sup && d < dimension)
                    {
                        if (tar_sub != tar_sup)
                        {
                            // For Check Security
                            std::cout << "\tDirection is not right when building para interface\n";
                            exit(-1);
                        }

                        interface.trans.sign[d] = (para.direction * tar_para.direction > 0) ? -1 : 1;
                        continue;
                    }

                    if (dimension == 2 && d == 2)
                    {
                        interface.trans.sign[d] = 0;
                        interface.trans.offset.k = 0;
                        continue;
                    }

                    if (abs(my_sup - my_sub) != abs(tar_sup - tar_sub))
                    {
                        // For Check Security
                        std::cout << "\tDirection is not right when building para interface\t" << d << std::endl;
                        exit(-1);
                    }
                    else
                        interface.trans.sign[d] = ((my_sup - my_sub) * (tar_sup - tar_sub) > 0) ? 1 : -1;
                }
                //-------------------------------------------------------------
                // 获取offset
                int offset[3] = {0, 0, 0};
                for (int d = 0; d < 3; ++d)
                {
                    int my_sub, tar_sub;
                    my_sub = abs(para.sub[d]);
                    tar_sub = abs(tar_para.sub[interface.trans.perm[d]]);

                    offset[d] = -interface.trans.sign[d] * my_sub + tar_sub;
                }
                interface.trans.offset = {offset[0], offset[1], offset[2]};
                //-------------------------------------------------------------
                topo.parallel_patches.push_back(interface);
            }
        }
        //=======================================================================

        if (my_rank == 0)
        {
            std::cout << "*************Finish the Topology Manipulating Process! !**************\n\n";
        }
        return topo;
    }
}