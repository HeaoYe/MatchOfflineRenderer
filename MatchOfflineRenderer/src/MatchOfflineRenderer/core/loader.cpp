#include <MatchOfflineRenderer/core/loader.hpp>
#include <MatchOfflineRenderer/core/logger.hpp>

namespace MatchOfflineRenderer {
    void Initialize() {
        g_logger.initialize();
    }

    void Destroy() {
        g_logger.destroy();
    }
}
