# ABACUS PEXSI 能量记账修复 + 精度验证（本地复现）

日期：2026-10-08 ｜ 补丁：`pexsi_fix.patch`（3 文件，+19/−2）
基线源码：ABACUS **v3.10.1**（克隆 `<abacus-pexsi-checkout>`）
PEXSI：conda-forge **1.2.0** + SuperLU_DIST ｜ 环境 `<pexsi-conda-env>`
验证体系：金刚石 C，Γ-only LCAO，ecutwfc 100 Ry，ONCV-PBE 赝势 + 7au 2s2p1d 轨道，np=4

> 注意：技术栈与 yeesuan（ABACUS 3.9.0.25 + PEXSI 2.0.0）**完全不同**，因此结论可互相印证。

## 1. 补丁内容

| 文件 | 改动 | 目的 |
|---|---|---|
| `module_elecstate/elecstate_lcao.cpp:169` | `kin_r` → `kin_r[0]` | 3.10 上 PEXSI 分支编译错误（`double**` 传给要 `const double*` 的构造函数） |
| `module_hamilt_lcao/hamilt_lcaodft/edm.cpp` | `#ifdef __PEXSI` 内补 `#include "module_elecstate/elecstate_lcao.h"` | 3.10 上 PEXSI 分支编译错误（`ElecStateLCAO` 未声明） |
| `module_hsolver/hsolver_lcao.cpp` | ① 把 `totalEnergyH/S/FreeEnergy` **广播到全部 rank**；② `f_en.eband = totalEnergyH`（真带能 Tr[H·DM]）而不是 `totalFreeEnergy`；③ 新增 `f_en.demet = totalFreeEnergy − totalEnergyH` | 修 rank 间总能量不一致 + 能量项错记（自由能塞进带能、熵项恒为 0） |

## 2. 验证结果

### 2.1 rank 一致性（`out_alllog 1`，2 原子 Γ，method 2）

| | rank0 | rank1 | rank2 | rank3 |
|---|---|---|---|---|
| **修复前** | −279.6757230584205 | −358.9418162859744 | −358.9418162859744 | −358.9418162859744 |
| **修复后** | **−280.0089654573225** | **−280.0089654573225** | **−280.0089654573225** | **−280.0089654573225** |

修复前 rank 间差值恰为一个 E_band（79.2661 eV）；修复后 E_band/E_entropy 在
各 rank 上也完全一致（78.9328606658 / −0.0000274553）。

### 2.2 能量项分解（2 原子 Γ）

| | E_band (eV) | E_entropy(−TS) (eV) | 总能量 (eV) |
|---|---|---|---|
| 稠密 scalapack_gvx（参考） | 78.9328695827 | −0.0000000000 | −280.0089789108313 |
| PEXSI method 1（修复前） | 79.2660932276 | 0.0000000000 | −279.6757230584356 |
| PEXSI method 1（修复后） | **78.9324204981** | **+0.3336727294** | −279.6757230584856 |

修复后 E_band 与稠密解一致到 **0.4 meV**；但**总能量不变**（78.9324 + 0.3337 = 79.2661，
代数上守恒）——所以**单靠记账修复不能让 PEXSI 与 genelpa 一致**。

### 2.3 决定性发现：`pexsi_method 2`

| 设置（2 原子 Γ） | 总能量 (eV) | 与稠密参考之差 |
|---|---|---|
| `pexsi_method 1` + npole 40（**ABACUS 默认**） | −279.6757230584546 | **+333.3 meV** |
| **`pexsi_method 2` + npole 40** | **−280.0089654571346** | **+0.013 meV** ✅ |
| `pexsi_method 2` + npole 60 | −279.9945595706079 | +14.4 meV |
| `pexsi_method 2` + npole 60 + 收紧 μ 区间 | −279.9945593349852 | +14.4 meV（μ 区间在 2 原子上无影响） |

16 原子 Γ（稠密求解器在该 Γ 点不收敛，仅作方法间对照）：

| 设置 | 总能量 (eV) | E_band (eV) | E_entropy (eV) |
|---|---|---|---|
| method 1 | −2399.6864 | 23.7333 | 6.0299 |
| method 2 | −2402.3519 | 23.7336 | 5.8336 |

## 3. 结论

1. **两个真 bug 已修**：PEXSI 路径的总能量在 MPI 各 rank 上不一致（差值恰为 E_band）；
   `eband` 被赋成自由能、`demet`（−TS）恒为 0。修复后各 rank 一致、分解与稠密求解器同口径。
2. **记账修复不改变总能量**（band+entropy 之和守恒），因此**精度问题的根因不在这里**。
3. **精度根因是 PEXSI 的极点方法**：ABACUS 默认 `pexsi_method 1`（轮廓积分）在 40 极点下
   使总能量偏 **333 meV**；换成 `pexsi_method 2`（Moussa）后与稠密求解器吻合到 **0.013 meV**。
   这与 yeesuan 上 `pexsi_bench.c` 的 Laplacian 基准（method 1 比 method 2/3 差 3–6 个数量级）
   完全一致，且这里是在 **PEXSI 1.2.0 + ABACUS 3.10.1** 上独立复现的。
4. **yeesuan 建议**：`lcao_1q_pexsi` 预设里固化 `pexsi_method 2`，再叠加你们已验证的
   `pexsi_mu_lower/upper` 收紧（54 原子必需）；同时把本补丁的能量记账改动移植到 3.9.0.25 构建。
   然后用你们的 G2 验收（genelpa 收敛密度 → PEXSI 5 步，判据：第 1 步 drho < 1e-4、
   E 与 −730768.58 eV 一致）复测 54 原子。

## 4. 复现命令

```bash
# 环境
E=<pexsi-conda-env>
# 二进制（已含补丁）
AB=<abacus-pexsi-checkout>/build/abacus
# 2 原子 Γ：method 1 vs method 2
cd <bench-dir>/mu_test
env LD_LIBRARY_PATH=$E/lib $E/bin/mpirun -np 4 $AB   # INPUT 选 INPUT_m1 / INPUT_m2
```
