# T2 版本

当前正式作答版本：v4（2026-10-24），与官方模板v0.2.0版本号无关。
候选来源：considering/v4/dynamic_kth.cpp。
原始源码SHA256：ff1826576d99635defa5ee042429b47c4c58d1bcbd72b1d6d2f6d22a9aaec2b5。
提交源码SHA256：1978675cac11738cfbfe23f9efe65dfb36a28f13a3a0202a84d74f07c12dec64。
Git blob：81362156564f5834abfba3642a0e4d1e45db7bd7。

v1隐式Treap、v2归并树、v3分块、v4快照索引；只读Treap另作实验。
不采用v5/v6/v7，不改变64次/32N构建门槛与一次构建策略。
已知修改后随机长查询和部分阶段负载超时，详情见IMPLEMENTATION.md。
第一版保留在Git历史及本地archive/v1_20261024，不在src中放重复类。
