// 我和谁连，在什么地方，怎么对齐（拓扑/连接）
// 把 1_grid 里“块之间的几何连接关系”抽出来，变成一组统一的接口描述（patch），供以后任何 field / halo 使用。
#pragma once
#include <vector>
#include <string>

#include "0_basic/TYPES.h"        // Int3, Box3
#include "1_grid/1_MPCNS_Grid.h"  // Grid, Block
#include "1_grid/Grid_Boundary.h" // Inner/Parallel/Physical_Boundary

namespace TOPO
{

    // 接口类型：内部同 rank、跨 rank、物理外边界
    enum class PatchKind
    {
        Inner,
        Parallel
    };

    // point_local -> point_nb 的索引变换：
    // nb[ perm[a] ] = sign[a] * local[ a ] + offset[a]
    struct IndexTransform
    {
        int perm[3]; // {0,1,2} 的排列，表示本地哪个坐标对应到邻居的 i/j/k
        int sign[3]; // +1 或 -1（二维时可以把第三维设成 0）
        Int3 offset; // 整数偏移（大多数情况为 0）
    };

    // 块-块接口（Inner + Parallel 都用这个）
    struct InterfacePatch
    {
        PatchKind kind;

        int this_rank; // this myid（0-based）
        int nb_rank;   // neighbor myid（Inner: 同 rank；Parallel: tar_myid）

        int this_block; // 0-based: Block index in this rank
        int nb_block;   // Inner: 对方 block index；Parallel: 如不知道可先设为 -1
        std::string this_block_name;
        std::string nb_block_name;

        // 在“节点 index 空间里的 box”，半开区间 [lo, hi)
        Box3 this_box_node; // 本块 interface 上的 node 区域
        Box3 nb_box_node;   // 对方块对应区域（Inner 可填，Parallel 暂时可留空或等于 this_box）

        // 从本块 (i,j,k) 到邻居块 (i',j',k') 的映射
        IndexTransform trans;
    };

    // 物理边界 patch：只在本块/本 rank 的一侧
    struct PhysicalPatch
    {
        int this_rank;
        int this_block;
        std::string this_block_name;

        int bc_id;           // boundary_num
        std::string bc_name; // boundary_name

        Box3 this_box_node; // 外边界在 node 空间中的区域 [lo,hi)
        int direction;      // ±1,±2,±3（和 Physical_Boundary 相同）

        const Physical_Boundary *raw = nullptr; // 回指原始结构（可选）
    };

    // 汇总：以后 3_field 只拿 Topology 这一个对象
    struct Topology
    {
        std::vector<InterfacePatch> inner_patches;
        std::vector<InterfacePatch> parallel_patches;
        std::vector<PhysicalPatch> physical_patches;
    };

    // 从 Grid 构造 topology（每个 rank 本地调用一次）
    Topology build_topology(Grid &grid, int my_rank, int dimension);

    void node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box);
} // namespace TOPO
