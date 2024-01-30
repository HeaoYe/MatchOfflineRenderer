#include <MatchOfflineRenderer/core/loader.hpp>
#include <MatchOfflineRenderer/core/logger.hpp>
#include <MatchOfflineRenderer/spectra/loader.hpp>

namespace MatchOfflineRenderer {
    void Initialize(const InitInfo &init_info) {
        g_logger.initialize();
        spectra::Initialize(init_info);
    }

    void Destroy() {
        spectra::Destroy();
        g_logger.destroy();
    }
}
