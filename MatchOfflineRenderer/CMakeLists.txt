project(MatchOfflineRenderer)

file(GLOB_RECURSE source CONFIGURE_DEPENDS src/MatchOfflineRenderer/*.cpp)

add_executable(MatchOfflineRenderer ${source})
target_include_directories(MatchOfflineRenderer PUBLIC include)

# Link spdlog
target_link_libraries(MatchOfflineRenderer PUBLIC spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)
