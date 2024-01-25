#include <MatchOfflineRenderer/MatchOfflineRenderer.hpp>

int main() {
    MatchOfflineRenderer::Initialize();
    
    MCH_INFO("Hello World!")

    MatchOfflineRenderer::Destroy();
    return 0;
}
