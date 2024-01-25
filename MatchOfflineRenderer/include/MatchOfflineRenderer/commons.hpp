#pragma once

namespace MatchOfflineRenderer {
    #define no_copy_construction(cls_name) cls_name(const cls_name &) = delete;
    #define no_move_construction(cls_name) cls_name(cls_name &&) = delete;
    #define default_construction(cls_name) public: \
        cls_name() = default;
    #define no_copy_move_construction(cls_name) no_copy_construction(cls_name) no_move_construction(cls_name)
    #define default_no_copy_move_construction(cls_name) default_construction(cls_name) no_copy_construction(cls_name) no_move_construction(cls_name)
}

#include <MatchOfflineRenderer/core/logger.hpp>
