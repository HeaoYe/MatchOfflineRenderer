project(MatchOfflineRenderer)

file(GLOB_RECURSE source CONFIGURE_DEPENDS src/MatchOfflineRenderer/*.cpp)
add_library(MatchOfflineRenderer STATIC ${source})
target_include_directories(MatchOfflineRenderer PUBLIC include)

# Link spdlog
target_link_libraries(MatchOfflineRenderer PUBLIC spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)

if (MATCH_DEBUG_MODE)
    target_compile_definitions(MatchOfflineRenderer PUBLIC MATCH_DEBUG)
endif()
