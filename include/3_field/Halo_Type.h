#pragma once
#include "0_basic/types.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/Field_Type.h"
#include <vector>

struct HaloRegion
{
    //--------------------------------------------------------------------------
    // 对于 parallel（跨 rank）情况，还可以加：
    int this_rank;
    int neighbor_rank;
    //--------------------------------------------------------------------------
    // Block_ID
    int this_block;     // 我是哪个 block（接收方）
    int neighbor_block; // 对方 block（发送方）
    //--------------------------------------------------------------------------
    // IJK_range
    Box3 recv_box; // 我这边要填的 ghost 区 (在 this_block 的 FieldBlock 索引空间里)
    Box3 send_box; // 对方的 interior 区 (在 neighbor_block 的 FieldBlock 索引空间里)
    //--------------------------------------------------------------------------
    // Mapping patern
    TOPO::IndexTransform trans; // 从 this_block 逻辑坐标 -> neighbor_block 逻辑坐标

    int send_flag, recv_flag;
};

struct HaloPattern
{
    StaggerLocation location; // Cell / FaceXi / ...
    int nghost;               // ghost 层数

    std::vector<HaloRegion> regions;
};

enum class Direction
{
    XMinus,
    XPlus,
    YMinus,
    YPlus,
    ZMinus,
    ZPlus
};