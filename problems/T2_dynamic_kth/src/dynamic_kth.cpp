#include "dynamic_kth.hpp"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <random>
#include <utility>

using namespace std;

struct DynamicKth::Impl {
    struct Node {
        int value;
        int left = -1;
        int right = -1;
        int size = 1;
        uint64_t priority;
    };
    vector<Node> nodes;
    vector<int> free_nodes;
    vector<pair<int, bool>> path;
    vector<int> traversal;
    vector<int> values;
    mt19937_64 rng{712367821};
    int root = -1;

    int Size(int u) const { return u == -1 ? 0 : nodes[u].size; }
    void Pull(int u) {
        nodes[u].size = 1 + Size(nodes[u].left) + Size(nodes[u].right);
    }
    bool Higher(int a, int b) const {
        return nodes[a].priority != nodes[b].priority
                   ? nodes[a].priority > nodes[b].priority : a > b;
    }
    int NewNode(int value) {
        Node node{value, -1, -1, 1, rng()};
        if (free_nodes.empty()) {
            nodes.push_back(node);
            return static_cast<int>(nodes.size()) - 1;
        }
        int u = free_nodes.back();
        free_nodes.pop_back();
        nodes[u] = node;
        return u;
    }
    explicit Impl(const vector<int>& initial) {
        nodes.reserve(initial.size());
        traversal.reserve(initial.size());
        for (int value : initial) {
            int u = NewNode(value);
            int last = -1;
            while (!traversal.empty() && Higher(u, traversal.back())) {
                last = traversal.back();
                traversal.pop_back();
                Pull(last);
            }
            nodes[u].left = last;
            if (!traversal.empty()) nodes[traversal.back()].right = u;
            traversal.push_back(u);
        }
        if (!traversal.empty()) root = traversal.front();
        for (auto it = traversal.rbegin(); it != traversal.rend(); ++it) Pull(*it);
        traversal.clear();
    }
    
    void Prepare() {
        auto required = static_cast<size_t>(Size(root)) + 1;
        if (path.capacity() < required)
            path.reserve(max(required, path.capacity() * 2));
    }
    pair<int, int> Split(int u, size_t count) {
        path.clear();
        while (u != -1) {
            auto left_size = static_cast<size_t>(Size(nodes[u].left));
            if (count <= left_size) {
                path.emplace_back(u, false);
                u = nodes[u].left;
            } else {
                count -= left_size + 1;
                path.emplace_back(u, true);
                u = nodes[u].right;
            }
        }
        int a = -1, b = -1;
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            int v = it->first;
            if (it->second) { nodes[v].right = a; Pull(v); a = v; }
            else { nodes[v].left = b; Pull(v); b = v; }
        }
        return {a, b};
    }
    int Merge(int a, int b) {
        path.clear();
        while (a != -1 && b != -1) {
            if (Higher(a, b)) {
                path.emplace_back(a, true);
                a = nodes[a].right;
            } else {
                path.emplace_back(b, false);
                b = nodes[b].left;
            }
        }
        int u = a == -1 ? b : a;
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            int v = it->first;
            if (it->second) nodes[v].right = u;
            else nodes[v].left = u;
            Pull(v);
            u = v;
        }
        return u;
    }
    void Collect(int u) {
        traversal.clear();
        values.clear();
        while (u != -1 || !traversal.empty()) {
            while (u != -1) {
                traversal.push_back(u);
                u = nodes[u].left;
            }
            u = traversal.back();
            traversal.pop_back();
            values.push_back(nodes[u].value);
            u = nodes[u].right;
        }
    }
};

DynamicKth::DynamicKth(const vector<int>& initial)
    : impl_(make_unique<Impl>(initial)) {}
DynamicKth::~DynamicKth() = default;

void DynamicKth::insert(size_t index, int value) {
    auto& t = *impl_;
    t.Prepare();
    int u = t.NewNode(value);
    auto parts = t.Split(t.root, index);
    t.root = t.Merge(t.Merge(parts.first, u), parts.second);
}
void DynamicKth::erase(size_t index) {
    auto& t = *impl_;
    t.Prepare();
    if (t.free_nodes.size() == t.free_nodes.capacity())
        t.free_nodes.reserve(max<size_t>(16, t.free_nodes.capacity() * 2));
    auto a = t.Split(t.root, index);
    auto b = t.Split(a.second, 1);
    t.root = t.Merge(a.first, b.second);
    t.free_nodes.push_back(b.first);
}
void DynamicKth::update(size_t index, int value) {
    auto& t = *impl_;
    int u = t.root;
    while (u != -1) {
        auto left_size = static_cast<size_t>(t.Size(t.nodes[u].left));
        if (index == left_size) { t.nodes[u].value = value; return; }
        if (index < left_size) u = t.nodes[u].left;
        else { index -= left_size + 1; u = t.nodes[u].right; }
    }
}
int DynamicKth::kthLargest(size_t left, size_t right,
                           size_t k) const {
    auto& t = *impl_;
    const auto length = right - left + 1;
    t.Prepare();
    t.values.reserve(length);
    t.traversal.reserve(length);
    auto a = t.Split(t.root, left);
    auto b = t.Split(a.second, length);
    t.Collect(b.first);
 
    t.root = t.Merge(a.first, t.Merge(b.first, b.second));
    nth_element(t.values.begin(), t.values.begin() + (k - 1),
                     t.values.end(), greater<int>());
    return t.values[k - 1];
}


