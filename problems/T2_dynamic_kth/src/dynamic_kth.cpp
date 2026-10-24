#include "dynamic_kth.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <utility>
#include <vector>

using namespace std;

// ---------------------------------------------------------------------------
// v4: 分块序列 + Fenwick 定位 + 基数选择 + 全区间摊销缓存
//
// 与 treap 方案的关键差别在"区间收集"的代价:
//   * treap 按子树大小下探, 收集区间是一次逐节点的指针追逐, 每个元素一次随机
//     访存(1e6 长度约 8 ms), 这是所有 O(L) 方案的主要成本。
//   * 这里把序列存成若干连续块(每块 ≤ 2048 个值), 用 Fenwick 树按块大小做前缀
//     和定位下标。收集区间 = 对整块做 memcpy, 1e6 长度约 0.4 ms, 快约 20 倍。
//
// 各操作:
//   Locate   O(log 块数) 的 Fenwick 二分
//   insert   Locate + 块内 memmove(≤ 8 KB); 块超过 2048 时对半分裂
//   erase    Locate + 块内 memmove; 块小于 256 时与相邻块合并
//   update   Locate + 改写单值
//   kthLargest  收集 + MSD 基数选择(11+11+10 位); 反复全区间查询走排序缓存
//
// 块内 memmove 完全顺序访存, 因此增删改也不比 treap 慢; 分裂/合并是每 O(块大小)
// 次操作才发生一次的摊销事件, 触发时重建 Fenwick 为 O(块数)。
// ---------------------------------------------------------------------------

namespace {
constexpr uint32_t kSignFlip = 0x80000000u;

// int -> uint32 的保序映射, 使无符号基数选择与有符号大小一致。
inline uint32_t ToKey(int v) { return static_cast<uint32_t>(v) ^ kSignFlip; }
inline int FromKey(uint32_t k) { return static_cast<int>(k ^ kSignFlip); }

constexpr size_t kBlockTarget = 1024;   // 建块时的块大小
constexpr size_t kBlockMax = 2048;      // 超过即分裂
constexpr size_t kBlockMin = 256;       // 低于即尝试与相邻块合并
constexpr size_t kRadixMinLen = 256;    // 低于该长度用 nth_element
constexpr uint32_t kRadix = 2048;       // 11 位基数
}  // namespace

struct DynamicKth::Impl {
    vector<vector<uint32_t>> blocks;   // 每块是一段连续的下标区间
    vector<int> bit;                   // Fenwick(1-based): 块大小的前缀和
    int bitTop = 0;                    // ≤ 块数 的最大 2 的幂, 用于二分定位
    size_t total = 0;

    // 复用工作区
    vector<uint32_t> keys;             // 区间键值(有效长度由调用处显式给出)
    vector<uint32_t> scratch;          // 基数选择压缩缓冲
    vector<uint32_t> hist;             // 基数计数桶
    vector<uint32_t> cacheKeys;        // 全区间升序缓存

    // Query-only snapshot. Mutation invalidates it in O(1); never rebuild on a mutation.
    struct StaticIndex {
        vector<uint32_t> alphabet;
        vector<vector<uint32_t>> prefix;
        vector<size_t> zeros;
        int levels = 0;
        void Build(const vector<vector<uint32_t>>& blocks, size_t n) {
            vector<uint64_t> ordered;
            ordered.reserve(n);
            size_t pos = 0;
            for (const auto& block : blocks)
                for (uint32_t key : block)
                    ordered.push_back((uint64_t(key) << 32) | pos++);
            sort(ordered.begin(), ordered.end());
            vector<uint32_t> current(n), next(n);
            alphabet.resize(n);
            for (size_t i = 0; i < n; ++i) {
                alphabet[i] = uint32_t(ordered[i] >> 32);
                current[uint32_t(ordered[i])] = uint32_t(i);
            }
            vector<uint64_t>().swap(ordered);
            levels = 0;
            while ((size_t(1) << levels) < n) ++levels;
            prefix.resize(levels);
            zeros.resize(levels);
            for (int level = 0; level < levels; ++level) {
                const int shift = levels - 1 - level;
                auto& p = prefix[level];
                p.resize(n + 1);
                p[0] = 0;
                for (size_t i = 0; i < n; ++i)
                    p[i + 1] = p[i] + ((current[i] >> shift) & 1u);
                zeros[level] = n - p[n];
                size_t lo = 0, hi = zeros[level];
                for (uint32_t v : current) {
                    if ((v >> shift) & 1u) next[hi++] = v;
                    else next[lo++] = v;
                }
                current.swap(next);
            }
        }
        uint32_t Kth(size_t l, size_t r, size_t k) const {
            uint32_t rank = 0;
            // k is the 0-based ascending rank, r is exclusive.
            for (int level = 0; level < levels; ++level) {
                const auto& p = prefix[level];
                size_t ones_l = p[l], ones_r = p[r];
                size_t zero_count = r - l - (ones_r - ones_l);
                rank <<= 1;
                if (k < zero_count) {
                    l -= ones_l; r -= ones_r;
                } else {
                    rank |= 1u; k -= zero_count;
                    l = zeros[level] + ones_l; r = zeros[level] + ones_r;
                }
            }
            return alphabet[rank];
        }
    };
    unique_ptr<StaticIndex> staticIndex;
    bool indexValid = false;
    bool indexAttempted = false;
    size_t longQueries = 0;
    size_t scannedSinceMutation = 0;

    bool MaybeIndex(size_t len) {
        if (len < 4096) return false;  // Keep short/point queries on the cheap block path.
        if (indexValid) return true;
        if (indexAttempted || total < 8192) return false;
        ++longQueries;
        scannedSinceMutation += len;
        if (longQueries < 64 || scannedSinceMutation < total * 32) return false;
        // At most one build per object: repeated mutations cannot cause rebuild storms.
        indexAttempted = true;
        try {
            staticIndex = make_unique<StaticIndex>();
            staticIndex->Build(blocks, total);
        } catch (const bad_alloc&) {
            staticIndex.reset();
            return false;  // Optional acceleration must not make the block path fail.
        }
        indexValid = true;
        return true;
    }

    bool cacheValid = false;
    uint32_t fullRangeHits = 0;

    explicit Impl(const vector<int>& initial) {
        const size_t n = initial.size();
        blocks.reserve(n / kBlockTarget + 1);
        for (size_t i = 0; i < n; i += kBlockTarget) {
            const size_t end = min(n, i + kBlockTarget);
            vector<uint32_t> blk(static_cast<size_t>(end - i));
            for (size_t j = i; j < end; ++j) blk[j - i] = ToKey(initial[j]);
            blocks.push_back(std::move(blk));
        }
        total = n;
        RebuildFenwick();
    }

    void InvalidateCache() {
        indexValid = false;
        longQueries = 0;
        scannedSinceMutation = 0;
        cacheValid = false;
        fullRangeHits = 0;
    }

    void RebuildFenwick() {
        const int m = static_cast<int>(blocks.size());
        bit.assign(static_cast<size_t>(m) + 1, 0);
        for (int i = 1; i <= m; ++i) {
            bit[i] += static_cast<int>(blocks[static_cast<size_t>(i) - 1].size());
            const int j = i + (i & -i);
            if (j <= m) bit[j] += bit[i];
        }
        bitTop = 0;
        if (m > 0) {
            bitTop = 1;
            while ((bitTop << 1) <= m) bitTop <<= 1;
        }
    }

    void FenwickAdd(int idx0, int delta) {
        for (int i = idx0 + 1; i < static_cast<int>(bit.size()); i += i & -i)
            bit[i] += delta;
    }

    // 下标 pos(0-based) 所在的块与块内偏移; pos 等于 total 时返回末尾位置。
    void Locate(size_t pos, size_t& bi, size_t& off) const {
        int idx = 0;
        int rem = static_cast<int>(pos);
        for (int pw = bitTop; pw; pw >>= 1) {
            const int nxt = idx + pw;
            if (nxt < static_cast<int>(bit.size()) && bit[nxt] <= rem) {
                idx = nxt;
                rem -= bit[nxt];
            }
        }
        if (static_cast<size_t>(idx) >= blocks.size()) {
            bi = blocks.size() - 1;
            off = blocks.back().size();
        } else {
            bi = static_cast<size_t>(idx);
            off = static_cast<size_t>(rem);
        }
    }

    void Insert(size_t index, int value) {
        const uint32_t key = ToKey(value);
        if (blocks.empty()) {
            blocks.emplace_back(static_cast<size_t>(1), key);
            total = 1;
            RebuildFenwick();
            InvalidateCache();
            return;
        }
        size_t bi = 0, off = 0;
        Locate(index, bi, off);
        auto& blk = blocks[bi];
        blk.insert(blk.begin() + static_cast<long long>(off), key);
        FenwickAdd(static_cast<int>(bi), 1);
        ++total;
        if (blk.size() > kBlockMax) {
            // 分裂: 后半段搬到紧邻的新块。
            const size_t half = blk.size() / 2;
            vector<uint32_t> tail(blk.begin() + static_cast<long long>(half),
                                  blk.end());
            blk.resize(half);
            blocks.insert(blocks.begin() + static_cast<long long>(bi) + 1,
                          std::move(tail));
            RebuildFenwick();
        }
        InvalidateCache();
    }

    void Erase(size_t index) {
        size_t bi = 0, off = 0;
        Locate(index, bi, off);
        auto& blk = blocks[bi];
        blk.erase(blk.begin() + static_cast<long long>(off));
        FenwickAdd(static_cast<int>(bi), -1);
        --total;
        if (blk.empty()) {
            blocks.erase(blocks.begin() + static_cast<long long>(bi));
            RebuildFenwick();
        } else if (blk.size() < kBlockMin && blocks.size() > 1) {
            // 与相邻块合并(只在合并后不超过 kBlockMax 时进行, 避免立刻再分裂)。
            const size_t nb = (bi + 1 < blocks.size()) ? bi + 1 : bi - 1;
            if (blk.size() + blocks[nb].size() <= kBlockMax) {
                if (nb > bi)
                    blk.insert(blk.end(), blocks[nb].begin(), blocks[nb].end());
                else
                    blk.insert(blk.begin(), blocks[nb].begin(), blocks[nb].end());
                blocks.erase(blocks.begin() + static_cast<long long>(nb));
                RebuildFenwick();
            }
        }
        InvalidateCache();
    }

    void Update(size_t index, int value) {
        size_t bi = 0, off = 0;
        Locate(index, bi, off);
        const uint32_t key = ToKey(value);
        if (blocks[bi][off] != key) {
            blocks[bi][off] = key;
            InvalidateCache();
        }
    }

    // 把下标区间 [left, left+len) 连续拷到 keys。
    void Collect(size_t left, size_t len) {
        if (keys.size() < len) keys.resize(len);
        size_t bi = 0, off = 0;
        Locate(left, bi, off);
        size_t need = len;
        size_t out = 0;
        while (need > 0) {
            const auto& blk = blocks[bi];
            const size_t avail = blk.size() - off;
            const size_t take = need < avail ? need : avail;
            std::memcpy(keys.data() + out, blk.data() + off,
                        take * sizeof(uint32_t));
            out += take;
            need -= take;
            ++bi;
            off = 0;
        }
    }

    // keys 前 len 个元素中的第 k 大(1-based); 会打乱 keys 的前 len 个元素。
    int SelectKth(size_t len, size_t k) {
        uint32_t* base = keys.data();
        if (len == 1) return FromKey(base[0]);
        if (k == len) {   // 最小值
            uint32_t mn = base[0];
            for (size_t i = 1; i < len; ++i)
                if (base[i] < mn) mn = base[i];
            return FromKey(mn);
        }
        if (k == 1) {     // 最大值
            uint32_t mx = base[0];
            for (size_t i = 1; i < len; ++i)
                if (base[i] > mx) mx = base[i];
            return FromKey(mx);
        }

        const size_t rank = len - k;   // 等价于升序 0-based 秩
        if (len < kRadixMinLen) {
            nth_element(base, base + rank, base + len);
            return FromKey(base[rank]);
        }
        if (scratch.size() < len) scratch.resize(len);
        if (hist.size() < kRadix) hist.resize(kRadix);

        static const int shifts[3] = {21, 10, 0};
        static const int widths[3] = {11, 11, 10};
        uint32_t* src = base;
        uint32_t* dst = scratch.data();
        size_t n = len;
        size_t below = 0;   // 已确认小于目标键的元素数

        for (int pass = 0; pass < 3; ++pass) {
            const uint32_t shift = static_cast<uint32_t>(shifts[pass]);
            const uint32_t buckets = 1u << widths[pass];
            const uint32_t mask = buckets - 1;

            fill(hist.begin(), hist.begin() + buckets, 0u);
            for (size_t i = 0; i < n; ++i) ++hist[(src[i] >> shift) & mask];

            uint32_t acc = 0;
            uint32_t b = 0;
            for (; b + 1 < buckets; ++b) {
                if (below + acc + hist[b] > rank) break;
                acc += hist[b];
            }
            below += acc;

            // 桶内只剩一个元素即目标; 最后一趟时同桶元素键值必然相同。
            if (hist[b] == 1 || shift == 0) {
                for (size_t i = 0; i < n; ++i)
                    if (((src[i] >> shift) & mask) == b) return FromKey(src[i]);
            }

            size_t m = 0;
            for (size_t i = 0; i < n; ++i)
                if (((src[i] >> shift) & mask) == b) dst[m++] = src[i];
            n = m;
            swap(src, dst);
        }
        return FromKey(src[0]);   // 不可达: 32 位已在三趟内完全确定
    }

    // 累计到这个次数才值得为全区间查询排一次序。
    size_t CacheThreshold(size_t n) const {
        size_t bits = 1;
        while ((size_t(1) << bits) < n) ++bits;
        return 2 * bits + 8;
    }

    void RebuildCache(size_t n) {
        if (keys.size() < n) keys.resize(n);
        size_t out = 0;
        for (const auto& blk : blocks) {
            std::memcpy(keys.data() + out, blk.data(),
                        blk.size() * sizeof(uint32_t));
            out += blk.size();
        }
        std::sort(keys.begin(), keys.begin() + static_cast<long long>(n));
        cacheKeys.assign(keys.begin(), keys.begin() + static_cast<long long>(n));
        cacheValid = true;
    }
};

DynamicKth::DynamicKth(const vector<int>& initial)
    : impl_(make_unique<Impl>(initial)) {}
DynamicKth::~DynamicKth() = default;

void DynamicKth::insert(size_t index, int value) { impl_->Insert(index, value); }
void DynamicKth::erase(size_t index) { impl_->Erase(index); }
void DynamicKth::update(size_t index, int value) { impl_->Update(index, value); }

int DynamicKth::kthLargest(size_t left, size_t right, size_t k) const {
    auto& t = *impl_;
    const size_t n = t.total;
    const size_t len = right - left + 1;
    if (len == 0) return 0;   // 题面保证非空, 仅作防御

    if (left == 0 && len == n) {
        if (t.cacheValid) return FromKey(t.cacheKeys[n - k]);
        if (++t.fullRangeHits >= t.CacheThreshold(n)) {
            t.RebuildCache(n);
            return FromKey(t.cacheKeys[n - k]);
        }
    }
    if (t.MaybeIndex(len))
        return FromKey(t.staticIndex->Kth(left, right + 1, len - k));
    t.Collect(left, len);
    return t.SelectKth(len, k);
}
