# protstruct-core

protstruct 的批量计算核心，包含：

- `align`：BLOSUM62 + 仿射空位 Needleman–Wunsch 两两全局比对；
- `disulfide`：最多 16 个 Cys 的位掩码动态规划，更多 Cys 使用贪心回退。

Python 通过标准输入一次提交多条记录，并通过 `PROTSTRUCT_CORE_BIN` 指定二进制。
未构建此核心时，Python 自动使用功能一致的纯 Python 实现。

## 构建

```powershell
cd software/rust/protstruct_core
cargo +stable-x86_64-pc-windows-gnu build --release
```
