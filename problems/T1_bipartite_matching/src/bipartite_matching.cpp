#include "bipartite_matching.hpp"

#include <algorithm>
#include <functional>
#include <queue>
#include <utility>
using namespace std;
struct BipartiteMatcher::Impl {
    uint32_t left_count, right_count;
    vector<Edge> edges;
    Impl(uint32_t l, uint32_t r, const vector<Edge>& e)
        : left_count(l), right_count(r), edges(e) {}
};

BipartiteMatcher::BipartiteMatcher(uint32_t l, uint32_t r,
                                 const vector<Edge>& edges)
    : impl_(make_unique<Impl>(l, r, edges)) {}
BipartiteMatcher::~BipartiteMatcher() = default;

vector<int32_t> BipartiteMatcher::maximumWeightMatching() const {
    using Wide = __int128;
    const int l = static_cast<int>(impl_->left_count);
    const int r = static_cast<int>(impl_->right_count);
    vector<int32_t> answer(l, -1);
    if (impl_->edges.empty()) return answer;
    const auto best = max_element(impl_->edges.begin(), impl_->edges.end(),
        [](const Edge& a, const Edge& b) { return a.weight < b.weight; });
    if (best->weight <= 0) {
        answer[best->left] = static_cast<int32_t>(best->right);
        return answer;
    }

    struct Arc { int to, next; Wide cost; bool available; };
    const int source = l + r, sink = source + 1, n = sink + 1;
    vector<int> head(n, -1);
    vector<Arc> arcs;
    arcs.reserve(2 * (impl_->edges.size() + l + r));
    auto add = [&](int u, int v, Wide cost) {
        int id = static_cast<int>(arcs.size());
        arcs.push_back({v, head[u], cost, true}); head[u] = id;
        arcs.push_back({u, head[v], -cost, false}); head[v] = id + 1;
    };
    vector<Wide> potential(n, 0);
    for (int i = 0; i < l; ++i) add(source, i, 0);
    for (int j = 0; j < r; ++j) add(l + j, sink, 0);
    for (const auto& e : impl_->edges) {
        if (e.weight <= 0) continue;
        const Wide cost = -static_cast<Wide>(e.weight);
        add(static_cast<int>(e.left), l + static_cast<int>(e.right), cost);
        potential[l + e.right] = min(potential[l + e.right], cost);
    }
 
    for (int j = 0; j < r; ++j)
        potential[sink] = min(potential[sink], potential[l + j]);
    const Wide inf = static_cast<Wide>(1) << 120;
    vector<Wide> distance(n);
    vector<int> previous(n);
    using Item = pair<Wide, int>;
    vector<Item> storage;
    storage.reserve(arcs.size());
    priority_queue<Item, vector<Item>, greater<Item>>
        queue(greater<Item>{}, move(storage));
    while (true) {
        fill(distance.begin(), distance.end(), inf);
        fill(previous.begin(), previous.end(), -1);
        distance[source] = 0;
        queue.emplace(0, source);
        while (!queue.empty()) {
            const auto [d, u] = queue.top(); queue.pop();
            if (d != distance[u]) continue;
            if (u == sink) break;
            for (int id = head[u]; id != -1; id = arcs[id].next) {
                const auto& a = arcs[id];
                if (!a.available) continue;
                const Wide nd = d + a.cost + potential[u] - potential[a.to];
                if (nd < distance[a.to]) {
                    distance[a.to] = nd;
                    previous[a.to] = id;
                    queue.emplace(nd, a.to);
                }
            }
        }
        while (!queue.empty()) queue.pop();
        if (distance[sink] == inf ||
            distance[sink] + potential[sink] - potential[source] >= 0) break;
 
        for (int v = 0; v < n; ++v)
            potential[v] += min(distance[v], distance[sink]);
        for (int v = sink; v != source;) {
            const int id = previous[v];
            arcs[id].available = false;
            arcs[id ^ 1].available = true;
            v = arcs[id ^ 1].to;
        }
    }
    for (int i = 0; i < l; ++i)
        for (int id = head[i]; id != -1; id = arcs[id].next)
            if (arcs[id].to >= l && arcs[id].to < l + r &&
                !arcs[id].available) answer[i] = arcs[id].to - l;
    return answer;
}
