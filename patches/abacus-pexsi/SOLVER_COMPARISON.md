# ABACUS ks_solver 对比：pexsi vs genelpa（以碳为例）

生成时间：2026-10-08 ｜ 工作目录：`<bench-dir>`

## 1. 体系与参数

| 项 | 值 |
|---|---|
| 赝势 | `<pseudo-dir>/C_ONCV_PBE-1.0.upf`（ONCV, NC, PBE, ZVAL=4） |
| 轨道 | `<orbital-dir>/C_gga_7au_100Ry_2s2p1d.orb`（13 轨道/原子） |
| ecutwfc | 100 Ry |
| 体系 | ① 2 原子金刚石原胞（a=6.74 Bohr）；② 16 原子超胞（2×2×2 常规胞, a=13.48 Bohr） |
| k 点 | Γ 点（PEXSI 强制）或 2×2×2 MP 网格（仅基准对照） |
| SCF | scf_thr 1e-7, scf_nmax 100；混合/展宽用默认（除非另注） |

## 2. PEXSI 构建（本地原本没有任何带 PEXSI 的 ABACUS）

- 环境：`conda create -p <pexsi-conda-env> -c conda-forge pexsi superlu_dist parmetis metis scalapack elpa fftw cereal cxx-compiler fortran-compiler`
  （PEXSI 1.2.0；ABACUS 用 `c_pexsi_interface.h` 即 PEXSI 1.x C 接口）
- 源码：`<abacus-pexsi-checkout>`（tag **v3.10.1** 的本地克隆）
- 配置：`-DENABLE_PEXSI=ON -DUSE_ELPA=OFF -DPEXSI_DIR=$ENV -DGIT_SUBMODULE=OFF`
- **必须打 2 个补丁才能编译通过**（develop HEAD 与 v3.10.1 同样报错，说明 PEXSI 路径未被维护）：
  1. `source/module_hamilt_lcao/hamilt_lcaodft/edm.cpp`：`#ifdef __PEXSI` 块用到 `elecstate::ElecStateLCAO`，需补 `#include "module_elecstate/elecstate_lcao.h"`
  2. `source/module_elecstate/elecstate_lcao.cpp:169`：`charge->kin_r` 是 `double**`，构造函数要 `const double*` → 改 `kin_r[0]`

## 3. 功能限制（实测）

- **PEXSI 只支持 Γ 点**：多 k（复数）分支在 `source/module_hsolver/diago_pexsi.cpp:92` 直接
  `WARNING_QUIT("DiagoPEXSI", "PEXSI is not completed for multi-k case")`。
- `pexsi_mu_thr` 实测**无效果**（0.05 → 1e-5 结果按位相同）。

## 4. 结果汇总

### 4.1 能量（eV）

| 体系 | genelpa | scalapack_gvx | pexsi | 收敛 |
|---|---|---|---|---|
| 2 原子 Γ | -280.0089789108313 | -280.0089789108313 | **-279.6757230584356** | 全部 yes |
| 16 原子 Γ | -2397.677678589701 | -2397.703069328551 | -2399.686377613799 | genelpa/scalapack **NO**，pexsi yes |
| 2 原子 2×2×2k | -309.7735070153198 | -309.7735070153198 | 不支持 | yes |
| 16 原子 2×2×2k | -2433.750343679361 (np1) / -2433.750343679248 (np4) | -2433.750343679362 (np1) / -2433.750343679272 (np4) | 不支持 | yes |

**精度结论**
- `genelpa` 与 `scalapack_gvx` 在每一个收敛算例上**按位一致**（差 ≤1e-10 eV，跨二进制亦然）。
- `pexsi` 在唯一收敛可比的小体系（2 原子 Γ）上偏高 **+0.333 eV**（≈42 meV/电子）。
- PEXSI 参数扫描（16 原子 Γ）：

| 设置 | ETOT/eV | 相对参考(-2397.703) |
|---|---|---|
| npole 40, mu_thr 0.05, temp 0.015（默认） | -2399.6864 | -1.98 |
| npole 120 | -2402.3520 | -4.65 |
| npole 120 + mu_thr 1e-5 | -2402.3520 | -4.65（mu_thr 无效） |
| npole 120 + temp 0.001 | -2396.7409 | +0.96 |
| temp 0.001 + 固定 pexsi_mu=1.1499 | -2396.7409 | +0.96 |

→ 能量强烈依赖 PEXSI 参数，且**不向稠密求解器的值收敛**；差异出现在密度相关项
（16 原子：E_xc 差 4.5 eV、E_Hartree 差 9.7 eV，E_Ewald 完全一致），即 **PEXSI 路径给出的密度本身不同**。

### 4.2 时间（墙钟，np=4 除非另注）

| 体系 | genelpa | scalapack_gvx | pexsi |
|---|---|---|---|
| 2 原子 Γ | 1.9 s | 1.3 s | 1.8 s |
| 16 原子 Γ（SCF 未收敛，仅供逐迭代比较） | 259.8 s | 49.3 s | **5.8 s**（np=1: 17 s） |
| 2 原子 2×2×2k | 12.4 s | 13.8 s | — |
| 16 原子 2×2×2k（已收敛） | 87.2 s (np1) / 54.5 s (np4) | 20.3 s (np1) / 15.3 s (np4) | — |

- 小体系三者相当（对角化不是瓶颈）。
- 16 原子：`pexsi` 最快，`scalapack_gvx` 比 `genelpa` 快约 4 倍（本机单节点、单进程/4 进程条件下 ELPA 并未体现优势）。

## 5. 已知问题与注意事项

1. 16 原子超胞在 **Γ 点**下 SCF 无法收敛（折叠带在 Γ 处带隙很小 → 近金属性电荷振荡；已试 broyden(0.1~0.2, ndim 20) + gaussian smearing 0.02~0.05 Ry + 200 步，仍不收敛）。同一体系用 2×2×2 k 网格则收敛。因此 16 原子 Γ 的**能量不可用**，仅时间可比。
2. PEXSI 只支持 Γ 点，故与 genelpa/scalapack 的严格对比必须在 Γ 点下、且选用 SCF 稳健的体系（绝缘体/分子）。
3. 结论：ABACUS v3.10.1 的 PEXSI 路径**速度优势明显但精度与可维护性存疑**（需打补丁才能编译、多 k 未实现、mu_thr 失效、能量与稠密求解器不一致）。

## 6. 产物

- 基准输入/输出：本目录（`gamma_bench/`、`big/`、`pexsi_study/`、`conv_test/`）
- PEXSI 版源码（含补丁）：`<abacus-pexsi-checkout>`
- PEXSI 环境：`<pexsi-conda-env>`
