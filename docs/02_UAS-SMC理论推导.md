# 倒立摆 UAS-SMC 控制理论推导

> 适用版本：2026-07-17 当前 MATLAB 实现。本文描述“代码实际实现的控制器”，不是对论文符号的机械照抄。

相关文档：[仿真使用说明](01_仿真使用说明.md)｜[稳定性证明](03_UAS-SMC稳定性证明.md)

## 1. 推导目标与资料来源

推导使用以下资料：

1. [`倒立摆状态模型-完整版.pdf`](../倒立摆状态模型-完整版.pdf) 第 3-7 页和第 11 页：倒立摆、电机及状态空间模型；
2. [`uas-smc.pdf`](../uas-smc.pdf) PDF 第 3-5、13-18 页（原文页码 19-21、29-34）：切换面、辅助面、连续无抖振趋近律式 (2.57)；
3. [`uas_smc_controller.m`](../uas_smc_controller.m)：当前控制律的最终定义；
4. [`run_uas_smc_sim.m`](../run_uas_smc_sim.m)：模型、虚拟状态 $S$、积分状态和数字仿真；
5. [`optimize_uas_smc.m`](../optimize_uas_smc.m)：控制参数的优化变量、约束与目标函数；
6. [`推导手稿`](../推导手稿/page1.jpg)：早期的双物理通道推导。

论文中的标量状态 $x_i$ 在本项目中对应虚拟状态 $\sigma$，论文中的 $\int x_i$ 对应

$$
\eta(t)=\int_0^t\sigma(\tau)\,d\tau.
$$

为避免与摆杆质量 $m$ 混淆，本文把论文辅助面常数 $m_i$ 记为 $m_h$；代码字段名为 `m_aux`。

## 2. 倒立摆小角度线性模型

### 2.1 状态与输入

定义

$$
x=\begin{bmatrix}p&\dot p&\theta&\dot\theta\end{bmatrix}^{\mathsf T},
\qquad u=u_a,
$$

其中：

- $p$：小车位置，单位 m；
- $\theta$：摆杆相对直立平衡点的小角度，单位 rad；
- $u_a$：直流电机电枢电压，单位 V。

在直立点附近取 $\sin\theta\approx\theta$、$\cos\theta\approx1$，并采用资料中的电机简化模型，得到

$$
\dot x=Ax+Bu,
$$

$$
A=
\begin{bmatrix}
0&1&0&0\\
0&A_{22}&A_{23}&0\\
0&0&0&1\\
0&A_{42}&A_{43}&0
\end{bmatrix},
\qquad
B=\begin{bmatrix}0&B_{21}&0&B_{41}\end{bmatrix}^{\mathsf T}.
$$

### 2.2 物理参数

| 符号 | 代码字段 | 数值 | 含义 |
|---|---|---:|---|
| $m$ | `m` | `0.0923 kg` | 摆杆质量 |
| $M$ | `M` | `0.1945 kg` | 小车质量 |
| $l$ | `l` | `0.1950 m` | 摆杆质心到转轴距离 |
| $J$ | `J` | `0.0029 kg·m²` | 摆杆绕质心转动惯量 |
| $g$ | `g` | `9.8 m/s²` | 重力加速度 |
| $R$ | `R` | `3.75 Ω` | 电枢电阻 |
| $r$ | `r` | `0.018 m` | 同步带轮半径 |
| $K_m$ | `Km` | `0.183` | 电磁转矩系数 |
| $K_e$ | `Ke` | `0.232` | 反电动势系数 |
| $I$ | `I` | `7.083e-6 kg·m²` | 电机轴及带轮等效转动惯量 |

公共分母为

$$
Q_{\mathrm{eq}}
=mJ+(J+ml^2)\left(M+\frac{I}{r^2}\right)
=0.0016544814366.
$$

矩阵元素为

$$
\begin{aligned}
A_{22}&=-\frac{K_mK_e(J+ml^2)}{Q_{\mathrm{eq}}Rr^2},
&A_{23}&=-\frac{m^2l^2g}{Q_{\mathrm{eq}}},\\
A_{42}&=\frac{mlK_mK_e}{Q_{\mathrm{eq}}Rr^2},
&A_{43}&=\frac{mgl(M+m+I/r^2)}{Q_{\mathrm{eq}}},\\
B_{21}&=\frac{K_m(J+ml^2)}{Q_{\mathrm{eq}}Rr},
&B_{41}&=-\frac{mlK_m}{Q_{\mathrm{eq}}Rr}.
\end{aligned}
$$

代入参数后，代码中的数值模型为

$$
A=\begin{bmatrix}
0&1&0&0\\
0&-135.3751994&-1.918831334&0\\
0&0&0&1\\
0&380.1344331&32.90655397&0
\end{bmatrix},
$$

$$
B=\begin{bmatrix}
0\\10.50324823\\0\\-29.49318877
\end{bmatrix}.
$$

可控矩阵

$$
\mathcal C=\begin{bmatrix}B&AB&A^2B&A^3B\end{bmatrix}
$$

满秩，即 $\operatorname{rank}(\mathcal C)=4$。开环矩阵存在右半平面极点，因此直立点自然不稳定。

## 3. 单输入系统的标量虚拟状态

### 3.1 为什么不直接采用两个物理控制通道

手稿尝试分别为小车误差和摆角误差构造 $H_x$、$H_\theta$，再令 $\dot H_x+\dot H_\theta=N$ 求控制输入。但倒立摆只有一个电压输入 $u_a$，两个物理通道不能被独立指定。直接对两个通道分别求逆会隐含两个独立输入，或只能人为加权合并，无法自动保证未被直接约束的内部动态稳定。

当前实现采用一条标量控制通道：

$$
\sigma=Sx.
$$

控制器只对 $\sigma$ 施加 UAS-SMC，同时在设计 $S$ 时预先保证剩余三维内部动态稳定。

### 3.2 可控标准形构造

令

$$
e_4^{\mathsf T}=\begin{bmatrix}0&0&0&1\end{bmatrix},
\qquad h=e_4^{\mathsf T}\mathcal C^{-1},
$$

构造变换

$$
T=\begin{bmatrix}
h\\hA\\hA^2\\hA^3
\end{bmatrix},
\qquad z=Tx.
$$

由于系统可控，$T$ 可逆。选择三个正数 $p_1<p_2<p_3$，使目标多项式为

$$
(s+p_1)(s+p_2)(s+p_3)
=s^3+c_2s^2+c_1s+c_0.
$$

在坐标 $z$ 中定义

$$
\sigma=c_0z_1+c_1z_2+c_2z_3+z_4
=\begin{bmatrix}c_0&c_1&c_2&1\end{bmatrix}Tx.
$$

所以未归一化向量为

$$
S_0=\begin{bmatrix}c_0&c_1&c_2&1\end{bmatrix}T.
$$

代码进一步归一化：

$$
S=\frac{S_0}{S_0B},
$$

从而

$$
SB=1.
$$

归一化不改变 $\sigma=0$ 对应的内部动态极点。

当前优化极点为

$$
-p_1=-3.781868302,quad
-p_2=-7.887913853,quad
-p_3=-11.91615787.
$$

目标系数为

$$
\begin{bmatrix}c_0&c_1&c_2&1\end{bmatrix}
=\begin{bmatrix}355.4715175&168.8900178&23.58594002&1\end{bmatrix},
$$

最终得到

$$
S=\begin{bmatrix}
-1.229863753&-0.5843272973&-1.237692695&-0.2419994225
\end{bmatrix},
\qquad SB=1.
$$

### 3.3 积分状态

引入

$$
\eta(t)=\int_0^t\sigma(\tau)\,d\tau,
\qquad \dot\eta=\sigma.
$$

因此 UAS-SMC 工作在二维虚拟平面 $(\eta,\sigma)$ 上，而物理系统仍为四维状态。

## 4. 切换面的选择

两条切换面为

$$
s_1=\sigma+\xi_1\eta,
\qquad
s_2=\sigma+\xi_2\eta,
$$

并要求

$$
\xi_1>\xi_2>0.
$$

原因是：每条面单独对应一阶稳定关系 $\dot\eta=-\xi_i\eta$；两条不同斜率的面又能把 $(\eta,\sigma)$ 平面划分成 UAS-SMC 所需的子空间。若 $\xi_1=\xi_2$，两个面重合，分区退化。

当前优化值为

$$
\xi_1=14.39126475,
\qquad
\xi_2=3.53995479.
$$

相平面中两条直线为

$$
s_1=0:\ \sigma=-\xi_1\eta,
\qquad
s_2=0:\ \sigma=-\xi_2\eta.
$$

## 5. 辅助面的选择与区域划分

### 5.1 四个子空间

[`uas_smc_controller.m`](../uas_smc_controller.m) 完全按照 $s_1,s_2$ 的符号选取当前辅助面：

| 区域 | 条件 | $\omega_1$ | $\omega_2$ |
|---:|---|---:|---:|
| 0 | $s_1<0,\ s_2<0$ | $1$ | $\alpha$ |
| 1 | $s_1<0,\ s_2\ge0$ | $-1$ | $\beta$ |
| 2 | $s_1\ge0,\ s_2<0$ | $1$ | $-\beta$ |
| 3 | $s_1\ge0,\ s_2\ge0$ | $-1$ | $-\alpha$ |

### 5.2 四条辅助面

统一写成

$$
h_i=\omega_{i1}\sigma+\omega_{i2}\eta+m_h,
\qquad m_h>0.
$$

当前代码取

$$
\alpha=8,qquad \beta=4.530452048,qquad m_h=0.01.
$$

四个系数对为

$$
\begin{aligned}
(\omega_{01},\omega_{02})&=(1,\alpha),\\
(\omega_{11},\omega_{12})&=(-1,\beta),\\
(\omega_{21},\omega_{22})&=(1,-\beta),\\
(\omega_{31},\omega_{32})&=(-1,-\alpha).
\end{aligned}
$$

这满足论文的中心对称条件：$h_0$ 与 $h_3$ 对称，$h_1$ 与 $h_2$ 对称。展开后：

$$
\begin{aligned}
h_0&=\sigma+\alpha\eta+m_h,
&h_0=0&:\ \sigma=-\alpha\eta-m_h,\\
h_1&=-\sigma+\beta\eta+m_h,
&h_1=0&:\ \sigma=\beta\eta+m_h,\\
h_2&=\sigma-\beta\eta+m_h,
&h_2=0&:\ \sigma=\beta\eta-m_h,\\
h_3&=-\sigma-\alpha\eta+m_h,
&h_3=0&:\ \sigma=-\alpha\eta+m_h.
\end{aligned}
$$

上述六条直线和闭环轨迹绘制在 [`uas_smc_phase_plane.png`](../results/uas_smc_phase_plane.png) 中。

$m_h$ 决定四条辅助面与原点的偏置，使其围成包含原点的凸集；它不直接出现在趋近律 $N$ 中。

## 6. 连续趋近律

### 6.1 连续混合系数

定义 $a=-\alpha$、$b=\beta$。为在切换边界上连续衔接不同表达式，使用

$$
\varepsilon=
\begin{cases}
\dfrac{|s_2|}{|s_1|+|s_2|},
&s_1s_2\le0,\ |s_1|>\mathrm{tol},\\[8pt]
\dfrac{|s_2|}{|s_2|+|\sigma|},
&s_2\sigma\le0,\ |\sigma|>\mathrm{tol},\\[8pt]
1,&\text{其他情形}.
\end{cases}
$$

代码取 `tol=1e-12`，第三段同时处理论文中的 $s_1\sigma\ge0$ 区域和分母退化点。

### 6.2 闭环期望动态

论文式 (2.57) 在本项目中的标量形式为

$$
\phi
=\varepsilon(a\sigma-ks_2)
+(1-\varepsilon)\frac{a+b}{2}\sigma,
$$

其中当前参数为

$$
a=-8,qquad b=4.530452048,qquad k=15.
$$

控制器将 $\phi$ 作为期望的 $\dot\sigma$。

当前辅助面的导数检验量定义为

$$
N=\omega_2\sigma+\omega_1\phi.
$$

因为 $\dot\eta=\sigma$，若实现了 $\dot\sigma=\phi$，则

$$
\dot h
=\omega_1\dot\sigma+\omega_2\dot\eta
=\omega_1\phi+\omega_2\sigma
=N.
$$

所以 $N$ 是当前辅助面的导数，不是控制电压，也不是 $\phi$ 本身。论文证明在相应符号和对称条件下，式 (2.57) 满足 $N\ge0$，且仅在 $(\sigma,\eta)=(0,0)$ 时 $N=0$。

## 7. 控制律推导

由

$$
\sigma=Sx
$$

可得

$$
\dot\sigma=S\dot x=SAx+SBu.
$$

令 $\dot\sigma=\phi$，得到未饱和控制律

$$
u^*=\frac{-SAx+\phi}{SB}.
$$

本项目已经把 $S$ 归一化到 $SB=1$，但代码仍保留一般形式，便于检查输入通道是否奇异。

实际输出加入电压限制：

$$
u=\operatorname{sat}_{[-12,12]}(u^*)
=\min\{12,\max(-12,u^*)\}.
$$

未饱和时，代回可得精确虚拟闭环

$$
\boxed{\dot\sigma=\phi,\qquad\dot\eta=\sigma.}
$$

饱和时该等式会出现附加误差，这也是[稳定性证明](03_UAS-SMC稳定性证明.md)必须单独限定证明范围的原因。

## 8. 当前优化参数

### 8.1 基线和优化结果

| 参数 | 基线 | 当前优化值 |
|---|---:|---:|
| 内部极点 | $[-2,-3,-4]$ | $[-3.781868,-7.887914,-11.916158]$ |
| $\xi_1$ | `2` | `14.39126475` |
| $\xi_2$ | `1` | `3.53995479` |
| $\alpha$ | `1` | `8` |
| $\beta$ | `1` | `4.530452048` |
| $a=-\alpha$ | `-1` | `-8` |
| $b=\beta$ | `1` | `4.530452048` |
| $k$ | `4` | `15` |
| $m_h$ | `0.01` | `0.01` |
| $u_{\max}$ | `12 V` | `12 V` |

### 8.2 优化变量的编码

优化器使用八维变量

$$
q=\begin{bmatrix}
p_1,&p_2-p_1,&p_3-p_2,&\xi_2,&\xi_1-\xi_2,&\alpha,&\beta,&k
\end{bmatrix}.
$$

当前最优向量为

$$
q=\begin{bmatrix}
3.781868302,&4.106045551,&4.028244012,&3.53995479,\\
10.85130996,&8,&4.530452048,&15
\end{bmatrix}.
$$

这种编码天然保证 $p_1<p_2<p_3$ 和 $\xi_1>\xi_2$。优化目标综合调节时间、位置/角度积分误差、峰值位置、电压均方根、电压变化和峰值角度，并对以下违反项施加大罚函数：

- 峰值摆角超过 $10^\circ$；
- 峰值位置超过 $0.20\,\mathrm m$；
- 未饱和控制需求超过 $12\,\mathrm V$；
- 终点未回到允许范围；
- $N<0$。

## 9. 数学符号与 MATLAB 字段映射

| 数学量 | MATLAB 位置 | 说明 |
|---|---|---|
| $x$ | `response.state` | 每行一个采样时刻，四列物理状态 |
| $\eta$ | `response.eta` | $\sigma$ 的积分状态 |
| $\sigma$ | `diagnostic.sigma` | `model.S*x` |
| $s_1,s_2$ | `diagnostic.s1`, `s2` | 两条切换面函数值 |
| $\varepsilon$ | `diagnostic.epsilon` | 连续混合系数 |
| $\phi$ | `diagnostic.phi` | 期望 $\dot\sigma$ |
| $N$ | `diagnostic.N` | 激活辅助面的导数检验量 |
| $\dot N$ | `diagnostic.N_dot` | 仿真后用 `gradient` 求得，不参与控制 |
| $h$ | `diagnostic.h` | 当前区域激活的辅助面函数值 |
| 区域编号 | `diagnostic.region` | `0,1,2,3` |
| $\omega_1,\omega_2$ | `diagnostic.omega1`, `omega2` | 当前辅助面系数 |
| $u^*$ | `diagnostic.u_unsat` | 未饱和电压 |
| $u$ | `response.u` | 限幅后的实际仿真输入 |

## 10. 适用条件

上述推导基于：

1. 直立点附近的小角度线性模型；
2. 四维状态均可用于反馈；
3. $(A,B)$ 可控且 $SB\ne0$；
4. $\xi_1>\xi_2>0$、$\alpha>0$、$\beta>0$、$k>0$；
5. 辅助面满足中心对称条件；
6. 理想等式 $\dot\sigma=\phi$ 只在控制未饱和且模型匹配时严格成立。

当前 MATLAB 仿真使用 `Ts=0.005 s` 的零阶保持数字实现，并加入 $\pm12\,\mathrm V$ 饱和。连续理论如何映射到数字闭环，以及能够严格声称什么，见[稳定性证明](03_UAS-SMC稳定性证明.md)。
