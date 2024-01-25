#include <MatchOfflineRenderer/core/logger.hpp>

namespace MatchOfflineRenderer {
    Logger g_logger;

    void Logger::initialize() {
        if (is_initialized) {
            return;
        }
        is_initialized = true;
        spd_logger = spdlog::stdout_color_mt("MatchOfflineRenderer Core");
        spd_logger->set_level(level);
        spd_logger->info("Create Logger");
    }

    void Logger::destroy() {
        if (!is_initialized) {
            return;
        }
        is_initialized = false;
        spd_logger->info("Destroy Logger");
        // spdlog::shutdown();
    }

    void Logger::set_level(spdlog::level::level_enum level) {
        this->level = level;
        if (is_initialized) {
            spd_logger->set_level(level);
        }
    }

    void set_log_level(spdlog::level::level_enum level) {
        g_logger.set_level(level);
    }
}
