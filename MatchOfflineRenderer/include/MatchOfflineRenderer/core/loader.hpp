#pragma once

namespace MatchOfflineRenderer {
    struct InitInfo {
        bool generate_rgb_to_spectrum_table { false };
        int gauss_newton_iteration_step { 15 };
    };

    void Initialize(const InitInfo &init_info = {});
    void Destroy();
}
