#pragma once
#include "3_field/Halo_Type.h"
#include "3_field/2_MPCNS_Field.h"
#include "0_basic/MPI_WRAPPER.h"

class Halo
{
public:
    Halo(Field *field, TOPO::Topology *topo)
    {
        fld_ = field;
        topo_ = topo;
        build_pattern();
    };

    void data_trans(std::string &field_name)
    {
        // 顺序上 inner/parallel 没强依赖，一般都可以：
        exchange_inner(field_name);
        exchange_parallel(field_name);
    }

    // 对所有场做 inner halo
    void exchange_inner_all();

    void exchange_inner(std::string field_name);

    // 对所有场做 para halo（MPI）
    void exchange_parallel_all();

    void exchange_parallel(std::string field_name);

private:
    Field *fld_;
    TOPO::Topology *topo_;

    // key = (StaggerLocation, nghost)，value = HaloPattern
    using PatternKey = std::pair<StaggerLocation, int>;
    std::map<PatternKey, HaloPattern> inner_patterns_;
    std::map<PatternKey, HaloPattern> parallel_patterns_;

    // 复用的 send / recv 缓冲区（MPI 并行用）
    std::vector<std::vector<double>> send_buf;
    std::vector<std::vector<double>> recv_buf;
    std::vector<MPI_Request> req_send, req_recv;
    std::vector<MPI_Status> stat_send, stat_recv;
    std::vector<int32_t> length;

    //=========================================================================

    void build_pattern();

    void build_inner_pattern(StaggerLocation loc, int nghost);
    void build_parallel_pattern(StaggerLocation loc, int nghost);

    void switch_range_ghost(Box3 &box, Direction &dir, int nghost);
    void switch_range_inner(Box3 &box, Direction &dir, int nghost);

    Box3 make_cell_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_cell_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facex_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facex_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facey_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facey_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facez_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facez_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgex_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgex_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgey_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgey_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgez_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgez_inner_box(Box3 &face, Direction &dir, int nghost);
};
