# UAS-SMC Documentation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 生成三份与当前 MATLAB 实现、论文公式和已保存仿真结果一致的中文 Markdown 文档。

**Architecture:** 三份文档按“如何运行”“控制器如何得到”“在什么条件下稳定”分层。MATLAB 源码和 MAT 文件决定当前实现事实，论文与倒立摆 PDF 提供理论来源，稳定性文档严格区分连续未饱和证明与采样、饱和条件下的数值证据。

**Tech Stack:** MATLAB R2025a、GitHub Markdown、LaTeX 数学公式、PowerShell、Poppler/PyPDF 文本与页面核对。

---

当前目录不是 Git 仓库，因此计划中的成果不能提交 commit；所有验证均在工作区直接执行。

## 文件结构

**Create:**

- `docs/01_仿真使用说明.md`：文件职责、结果目录、运行与配置方法。
- `docs/02_UAS-SMC理论推导.md`：模型、虚拟状态、切换面、辅助面、趋近律和控制律推导。
- `docs/03_UAS-SMC稳定性证明.md`：理想闭环证明、内部动态证明及工程边界。

**Read and cross-check:**

- `uas_smc_controller.m`
- `run_uas_smc_sim.m`
- `optimize_uas_smc.m`
- `test_uas_smc.m`
- `uas-smc.pdf`
- `倒立摆状态模型-完整版.pdf`
- `推导手稿/page1.jpg`
- `推导手稿/page2.jpg`
- `results/uas_smc_best_params.mat`
- `results/uas_smc_results.mat`
- `results/optimization_history.mat`

### Task 1: 建立可核验的事实基线

**Files:**

- Read: `run_uas_smc_sim.m`
- Read: `uas_smc_controller.m`
- Read: `optimize_uas_smc.m`
- Read: `results/uas_smc_results.mat`
- Read: `results/uas_smc_best_params.mat`

- [ ] **Step 1: 运行目标文件缺失检查**

Run:

```powershell
$targets = @(
  'docs/01_仿真使用说明.md',
  'docs/02_UAS-SMC理论推导.md',
  'docs/03_UAS-SMC稳定性证明.md'
)
$existing = $targets | Where-Object { Test-Path $_ }
if ($existing) { throw "Target docs unexpectedly exist: $existing" }
'EXPECTED_FAIL_BASELINE: three target documents are absent'
```

Expected: 输出 `EXPECTED_FAIL_BASELINE`，证明验收目标尚未实现。

- [ ] **Step 2: 从 MAT 文件打印当前参数、矩阵和性能数据**

Run:

```powershell
& 'D:\matlab\matlabr2025a\bin\matlab.exe' -batch "d=load('results/uas_smc_results.mat'); p=load('results/uas_smc_best_params.mat'); disp(d.sim_results.model.A); disp(d.sim_results.model.B); disp(d.sim_results.selected.design.S); disp(p.best_ctrl); disp(d.sim_results.selected.metrics); fprintf('validation=%.6f\n',p.optimization_report.validation.passRate);"
```

Expected: 打印数值矩阵 `A`、`B`、`S`，优化控制参数、名义性能指标和验证通过率；文档只使用这次输出中的数值。

- [ ] **Step 3: 核对资料中的公式页面**

Inspect:

- `倒立摆状态模型-完整版.pdf` 第 3-7 页和第 11 页：力学、电机模型、状态空间矩阵及参数。
- `uas-smc.pdf` PDF 第 3-5 页：切换面、子空间和辅助面。
- `uas-smc.pdf` PDF 第 8-9 页：稳定性定理及正不变集结论。
- `uas-smc.pdf` PDF 第 13-18 页：连续控制、对称简化条件和式 (2.57) 趋近律。
- 两页推导手稿：倒立摆模型到标量虚拟状态的适配过程。

Expected: 确认文档符号映射为论文标量状态 `x_i` 对应本项目 `sigma`，论文积分项对应本项目 `eta`。

### Task 2: 编写仿真使用说明

**Files:**

- Create: `docs/01_仿真使用说明.md`
- Read: `run_uas_smc_sim.m`
- Read: `optimize_uas_smc.m`
- Read: `test_uas_smc.m`
- Read: `results/*`

- [ ] **Step 1: 写入文档头和项目快速开始**

文档必须依次包含以下一级或二级标题：

```markdown
# UAS-SMC 倒立摆 MATLAB 仿真使用说明
## 1. 项目目标与适用范围
## 2. 目录与文件职责
## 3. 快速运行
## 4. 仿真配置
## 5. 参数优化
## 6. 输出结果说明
## 7. 自动化测试
## 8. 当前标准结果
## 9. 常见问题与限制
```

快速运行命令必须是：

```matlab
cd('D:\inverted pendulum\UAS-SMC');
run('run_uas_smc_sim.m');
```

- [ ] **Step 2: 写入四个 MATLAB 文件的职责表**

表格至少包含“文件”“职责”“主要输入”“主要输出”“是否通常直接运行”五列，并准确说明：

- `uas_smc_controller.m` 是单步控制器函数。
- `run_uas_smc_sim.m` 是标准仿真入口。
- `optimize_uas_smc.m` 是无工具箱 PSO 加局部搜索入口。
- `test_uas_smc.m` 是九项回归与验收测试。

- [ ] **Step 3: 写入可复制的配置示例**

至少提供以下完整示例：优化控制器、基线控制器、自定义初始角、自定义控制参数、重新优化、运行测试。所有调用都使用脚本支持的 `simulation_options` 或 `optimizer_options` 字段。

- [ ] **Step 4: 逐项解释结果文件**

覆盖当前 `results` 中全部 9 个文件，并说明两张新增图：

- `uas_smc_phase_plane.png`
- `uas_smc_reaching_law.png`

同时给出 `uas_smc_results.mat` 中 `sim_results.model`、`options`、`baseline`、`selected` 和诊断字段的读取示例。

- [ ] **Step 5: 验证文档一**

Run:

```powershell
$p='docs/01_仿真使用说明.md'
if (-not (Test-Path $p)) { throw 'Missing simulation guide' }
$t=Get-Content -Raw -Encoding utf8 $p
@('200 Hz','0.005','run_uas_smc_sim.m','optimize_uas_smc.m','test_uas_smc.m','uas_smc_phase_plane.png','uas_smc_reaching_law.png') | ForEach-Object {
  if (-not $t.Contains($_)) { throw "Missing guide item: $_" }
}
'GUIDE_VERIFY: pass'
```

Expected: `GUIDE_VERIFY: pass`。

### Task 3: 编写 UAS-SMC 理论推导

**Files:**

- Create: `docs/02_UAS-SMC理论推导.md`
- Read: `uas_smc_controller.m`
- Read: `run_uas_smc_sim.m`
- Read: `optimize_uas_smc.m`
- Read: `uas-smc.pdf`
- Read: `倒立摆状态模型-完整版.pdf`
- Read: `推导手稿/page1.jpg`
- Read: `推导手稿/page2.jpg`

- [ ] **Step 1: 写入模型和符号约定**

文档使用以下结构：

```markdown
# 倒立摆 UAS-SMC 控制理论推导
## 1. 推导目标与资料来源
## 2. 倒立摆小角度线性模型
## 3. 单输入系统的标量虚拟状态
## 4. 切换面的选择
## 5. 辅助面的选择与区域划分
## 6. 连续趋近律
## 7. 控制律推导
## 8. 当前优化参数
## 9. 数学符号与 MATLAB 字段映射
## 10. 适用条件
```

列出物理参数、`Qeq`、矩阵元素公式、数值 `A` 与 `B`，并固定状态顺序为
`[小车位置, 小车速度, 摆角, 摆角速度]`。

- [ ] **Step 2: 推导标量虚拟状态和内部动态**

从可控矩阵

```math
\mathcal C=[B\;AB\;A^2B\;A^3B]
```

和变换矩阵 `T` 出发，解释代码如何由目标多项式

```math
(s+p_1)(s+p_2)(s+p_3)
```

构造 `S`，再归一化到 `SB=1`。写清该构造用于保证三维内部动态极点为 `-p1,-p2,-p3`。

- [ ] **Step 3: 写入切换面和四条辅助面**

必须完整给出：

```math
s_1=\sigma+\xi_1\eta,\qquad s_2=\sigma+\xi_2\eta,
```

以及代码对应的四个系数组合：

```math
(\omega_{01},\omega_{02})=(1,\alpha),
\quad(-1,\beta),
\quad(1,-\beta),
\quad(-1,-\alpha).
```

由 `h_i=omega_i1*sigma+omega_i2*eta+m` 展开相平面中 `h0=0` 到 `h3=0` 的四条直线，并与 `run_uas_smc_sim.m` 的绘图公式逐项一致。

- [ ] **Step 4: 推导连续趋近律和控制律**

写出 `epsilon` 的三段定义、`a=-alpha`、`b=beta`、混合项 `phi`、检验量 `N`，并由

```math
\dot\sigma=SAx+SBu
```

推到

```math
u^*=\frac{-SAx+\phi}{SB},\qquad
u=\operatorname{sat}_{[-12,12]}(u^*).
```

说明 `N` 是辅助面导数而不是控制输入，也不是 `phi` 本身。

- [ ] **Step 5: 写入实际参数和代码映射表**

从 MAT 文件填入基线与优化参数、内部极点、`S` 和 `SB`。表格覆盖 `sigma`、`eta`、`s1`、`s2`、`epsilon`、`phi`、`N`、`h`、`region`、`u_unsat` 和 `u`。

- [ ] **Step 6: 验证文档二**

Run:

```powershell
$p='docs/02_UAS-SMC理论推导.md'
$t=Get-Content -Raw -Encoding utf8 $p
@('sigma','eta','s_1','s_2','h_0','h_1','h_2','h_3','epsilon','phi','SB=1','12') | ForEach-Object {
  if (-not $t.Contains($_)) { throw "Missing theory item: $_" }
}
'THEORY_VERIFY: pass'
```

Expected: `THEORY_VERIFY: pass`。

### Task 4: 编写稳定性证明

**Files:**

- Create: `docs/03_UAS-SMC稳定性证明.md`
- Read: `docs/02_UAS-SMC理论推导.md`
- Read: `uas-smc.pdf`
- Read: `uas_smc_controller.m`
- Read: `run_uas_smc_sim.m`

- [ ] **Step 1: 明确假设、命题和结论范围**

使用以下章节：

```markdown
# 当前倒立摆 UAS-SMC 的稳定性证明
## 1. 证明对象与结论
## 2. 假设条件
## 3. 虚拟二阶闭环
## 4. 辅助面单向性与正不变集
## 5. 虚拟状态的渐近收敛
## 6. 内部动态与四维状态收敛
## 7. 连续无抖振性质
## 8. 采样与饱和条件下的结论边界
## 9. 数值证据
## 10. 最终结论
```

结论必须限定为连续、模型准确、全状态可用且控制不饱和的理想闭环。

- [ ] **Step 2: 证明虚拟闭环和辅助面单向性**

代入控制律得到 `dot(sigma)=phi`、`dot(eta)=sigma`，进而逐行得到

```math
\dot h_i=\omega_{i1}\phi+\omega_{i2}\sigma=N_i.
```

引用论文式 (2.57) 对称条件下的 `N_i>=0` 结论，并明确当前程序通过 `minimumN>=-1e-9` 检查数值实现，而该数值检查本身不是证明。

- [ ] **Step 3: 证明虚拟状态与内部状态收敛**

按“进入正不变集”“最大不变零集仅含原点”“`sigma,eta` 渐近收敛”“Hurwitz 内部动态受衰减输入驱动”“原坐标变换可逆”五步展开。不得用 `N>=0` 单独跳步推导全状态稳定。

- [ ] **Step 4: 说明连续性、采样和饱和边界**

分别说明：

- `epsilon` 在边界上的插值如何使 `phi` 连续衔接。
- 数字实现只在采样时刻更新控制量，严格离散证明需要额外误差界。
- 饱和时出现附加误差
  `SB*(sat(u*)-u*)`，因此理想等式 `dot(sigma)=phi` 被破坏。
- 当前初始角和参数摄动测试提供的是局部工作范围内的实用稳定性证据。

- [ ] **Step 5: 验证文档三**

Run:

```powershell
$p='docs/03_UAS-SMC稳定性证明.md'
$t=Get-Content -Raw -Encoding utf8 $p
@('连续','未饱和','N_i','正不变集','内部动态','Hurwitz','采样','饱和','不是证明') | ForEach-Object {
  if (-not $t.Contains($_)) { throw "Missing proof item: $_" }
}
'PROOF_VERIFY: pass'
```

Expected: `PROOF_VERIFY: pass`。

### Task 5: 全局一致性与最终验收

**Files:**

- Verify: `docs/01_仿真使用说明.md`
- Verify: `docs/02_UAS-SMC理论推导.md`
- Verify: `docs/03_UAS-SMC稳定性证明.md`

- [ ] **Step 1: 检查三个文件、占位符和本地链接**

Run:

```powershell
$docs=@('docs/01_仿真使用说明.md','docs/02_UAS-SMC理论推导.md','docs/03_UAS-SMC稳定性证明.md')
foreach($p in $docs){
  if(-not (Test-Path $p)){throw "Missing $p"}
  $t=Get-Content -Raw -Encoding utf8 $p
  if($t -match '待填写|PLACEHOLDER'){throw "Placeholder in $p"}
  if($t.Length -lt 2000){throw "Document too short: $p"}
}
'DOC_FILES_VERIFY: pass'
```

Expected: `DOC_FILES_VERIFY: pass`。

- [ ] **Step 2: 运行 MATLAB 回归测试确保文档工作未影响仿真**

Run:

```powershell
& 'D:\matlab\matlabr2025a\bin\matlab.exe' -batch "r=runtests('test_uas_smc.m'); assert(numel(r)==9 && all([r.Passed])); fprintf('MATLAB_VERIFY: %d/%d passed\n',sum([r.Passed]),numel(r));"
```

Expected: `MATLAB_VERIFY: 9/9 passed`。

- [ ] **Step 3: 人工一致性复核**

逐项确认：

1. 文档一的运行命令能够直接复制。
2. 文档二的四条辅助面方程与相平面绘图代码相同。
3. 文档二的优化参数与 MAT 文件相同。
4. 文档三没有把连续未饱和证明扩展为数字饱和闭环的全局证明。
5. 三份文档互相链接，并能链接到论文、源码、手稿和结果图。

Expected: 五项全部满足后交付。
