#pragma once

#include <cassert>

namespace MatchOfflineRenderer {
    #define no_copy_construction(cls_name) cls_name(const cls_name &) = delete;
    #define no_move_construction(cls_name) cls_name(cls_name &&) = delete;
    #define default_construction(cls_name) public: \
        cls_name() = default;
    #define no_copy_move_construction(cls_name) no_copy_construction(cls_name) no_move_construction(cls_name)
    #define default_no_copy_move_construction(cls_name) default_construction(cls_name) no_copy_construction(cls_name) no_move_construction(cls_name)

    #define MCH_ASSERT(expr) (assert(expr));
    #if defined (MATCH_DEBUG)
    #define MCH_DASSERT(expr) MCH_ASSERT(expr)
    #else
    #define MCH_DASSERT(expr)
    #endif
}

#include <MatchOfflineRenderer/core/logger.hpp>
