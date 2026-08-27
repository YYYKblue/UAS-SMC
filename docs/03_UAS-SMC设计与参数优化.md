# UAS-SMC 设计与参数优化

本文按当前 MATLAB 实现推导切换面、辅助面、连续趋近律、控制律及参数优化方法。虚拟状态 $\sigma=Sx$ 的构造见[虚拟状态整体设计与推导](02_虚拟状态整体设计与推导.md)。

## 1. 从四阶模型到 UAS-SMC 标量通道

名义模型为

$$
\dot x=Ax+Bu,
$$

定义虚拟状态和积分状态

$$
\sigma=Sx,
\qquad
\eta(t)=\int_0^t\sigma(\tau)\,d\tau.
$$

由此有

$$
\dot\sigma=SAx+SBu,
\qquad
\dot\eta=\sigma.
$$

代码把 $S$ 归一化到 $SB=1$，但保留一般形式。若指定期望虚拟动态 $\dot\sigma=\phi$，则未饱和控制律为

$$
\boxed{
u^*=\frac{-SAx+\phi}{SB}.
}
$$

因此 UAS-SMC 只需在二维平面 $(\eta,\sigma)$ 中设计 $\phi$。在模型匹配且未饱和时，代回可得

$$
\dot\sigma=\phi,
\qquad
\dot\eta=\sigma.
$$

## 2. 两条切换面的选取

按照 [`uas-smc.pdf`](../uas-smc.pdf) 的式 (2.26)，当前标量通道采用

$$
s_1=\sigma+\xi_1\eta,
\qquad
s_2=\sigma+\xi_2\eta,
$$

并要求

$$
\xi_1>\xi_2>0.
$$

这样选取有三点原因：

1. 在 $s_i=0$ 上有 $\dot\eta=\sigma=-\xi_i\eta$，每一条切换面都对应稳定的一阶运动；
2. $\xi_1\ne\xi_2$ 避免两条直线重合，使相平面能够被划分为四个子空间；
3. $\xi_1>\xi_2>0$ 是论文连续无抖振辅助面构造所需的基本几何条件。

相平面中的直线为

$$
s_1=0:\ \sigma=-\xi_1\eta,
\qquad
s_2=0:\ \sigma=-\xi_2\eta.
$$

控制器按 $s_1,s_2$ 的符号定义区域：

| 区域 | 条件 |
|---:|---|
| 0 | $s_1<0$ 且 $s_2<0$ |
| 1 | $s_1<0$ 且 $s_2\ge0$ |
| 2 | $s_1\ge0$ 且 $s_2<0$ |
| 3 | $s_1\ge0$ 且 $s_2\ge0$ |

区域编号完全由控制器的符号判据决定，不能按普通笛卡尔象限编号理解。

## 3. 四条辅助面的结构

当前区域的辅助面统一写成

$$
h_i=\omega_{i1}\sigma+\omega_{i2}\eta+m_h,
\qquad m_h>0.
$$

代码选用中心对称结构

$$
\begin{aligned}
h_0&=\sigma+\alpha\eta+m_h,\\
h_1&=-\lambda\sigma+\lambda\beta\eta+m_h,\\
h_2&=\lambda\sigma-\lambda\beta\eta+m_h,\\
h_3&=-\sigma-\alpha\eta+m_h,
\end{aligned}
$$

其中 $\beta>0$、$\lambda>0$，且 $h_0$ 与 $h_3$、$h_1$ 与 $h_2$ 分别关于原点中心对称。各区域使用的系数为：

| 区域 | $\omega_{i1}$ | $\omega_{i2}$ |
|---:|---:|---:|
| 0 | $1$ | $\alpha$ |
| 1 | $-\lambda$ | $\lambda\beta$ |
| 2 | $\lambda$ | $-\lambda\beta$ |
| 3 | $-1$ | $-\alpha$ |

$m_h$ 使四条直线不经过原点，从而围成包含原点的凸集。它决定集合尺度，但在 $\dot h_i$ 中消失。

## 4. 为什么 $\alpha$ 和 $\lambda$ 不能独立选取

仅满足中心对称只能保证四条辅助面围成一个中心对称四边形，不能保证相邻辅助面的交点落在 $s_1=0$ 或 $s_2=0$ 上。论文式 (2.27) 要求四个顶点交替位于两条切换面上，因此 $\alpha$ 和 $\lambda$ 必须满足额外约束。

### 4.1 在 $s_1=0$ 上约束 $h_1$ 与 $h_3$ 的交点

令 $h_1=h_3=0$，并要求该点满足 $s_1=0$，即 $\sigma=-\xi_1\eta$。代入 $h_3=0$：

$$
(\xi_1-\alpha)\eta+m_h=0,
$$

所以

$$
\eta=-\frac{m_h}{\xi_1-\alpha},
\qquad
\sigma=\frac{m_h\xi_1}{\xi_1-\alpha}.
$$

再代入 $h_1=0$，得到

$$
\lambda(\xi_1+\beta)=\xi_1-\alpha,
$$

即

$$
\lambda=\frac{\xi_1-\alpha}{\xi_1+\beta}.
$$

### 4.2 在 $s_2=0$ 上约束 $h_3$ 与 $h_2$ 的交点

同理，令 $h_3=h_2=0$ 且 $\sigma=-\xi_2\eta$，可得

$$
\lambda=\frac{\alpha-\xi_2}{\xi_2+\beta}.
$$

两式必须同时成立：

$$
\frac{\xi_1-\alpha}{\xi_1+\beta}
=
\frac{\alpha-\xi_2}{\xi_2+\beta}.
$$

解得

$$
\boxed{
\alpha=
\frac{\beta(\xi_1+\xi_2)+2\xi_1\xi_2}
{\xi_1+\xi_2+2\beta},
}
$$

以及

$$
\boxed{
\lambda=
\frac{\xi_1-\xi_2}
{\xi_1+\xi_2+2\beta}.
}
$$

还可以直接验证

$$
\alpha-\xi_2
=\frac{(\xi_1-\xi_2)(\beta+\xi_2)}
{\xi_1+\xi_2+2\beta}>0,
$$

$$
\xi_1-\alpha
=\frac{(\xi_1-\xi_2)(\xi_1+\beta)}
{\xi_1+\xi_2+2\beta}>0,
$$

故

$$
\xi_2<\alpha<\xi_1,
\qquad
\lambda>0.
$$

这就是 [`uas_smc_auxiliary_geometry.m`](../uas_smc_auxiliary_geometry.m) 中几何参数的来源。

## 5. 正不变集及其顶点

四条辅助面围成集合

$$
Q=\left\{(\sigma,\eta)\mid h_i(\sigma,\eta)\ge0,\ i=0,1,2,3\right\}.
$$

按代码的逆时针顶点顺序，有

$$
\begin{aligned}
P_1&=h_1\cap h_3
=\left(
\frac{m_h\xi_1}{\xi_1-\alpha},
-\frac{m_h}{\xi_1-\alpha}
\right),\\
P_2&=h_3\cap h_2
=\left(
-\frac{m_h\xi_2}{\alpha-\xi_2},
\frac{m_h}{\alpha-\xi_2}
\right),\\
P_3&=-P_1,\\
P_4&=-P_2.
\end{aligned}
$$

这里每个坐标对按 $(\sigma,\eta)$ 排列。由构造可知 $P_1,P_3$ 位于 $s_1=0$，$P_2,P_4$ 位于 $s_2=0$。

当前参数下

$$
\begin{aligned}
\xi_1&=17.4331663,&
\xi_2&=5.4331663,\\
\alpha&=9.0510902,&
\beta&=3.6797011,\\
\lambda&=0.3970127,&
m_h&=0.01.
\end{aligned}
$$

四个顶点约为

$$
\begin{bmatrix}
0.0208&-0.0012\\
-0.0150&0.0028\\
-0.0208&0.0012\\
0.0150&-0.0028
\end{bmatrix},
$$

列顺序为 $[\sigma,\eta]$。保存结果中的最大切换面残差约为 $6.94\times10^{-18}$，测试阈值为 $10^{-12}$。

## 6. 连续混合系数与趋近项

令

$$
a=-\alpha,
\qquad
b=\beta.
$$

代码按照论文式 (2.57) 构造连续混合系数

$$
\varepsilon=
\begin{cases}
\dfrac{\lvert s_2\rvert}{\lvert s_1\rvert+\lvert s_2\rvert},
&s_1s_2\le0,\ \lvert s_1\rvert>\mathrm{tol},\\[8pt]
\dfrac{\lvert s_2\rvert}{\lvert s_2\rvert+\lvert\sigma\rvert},
&s_2\sigma\le0,\ \lvert\sigma\rvert>\mathrm{tol},\\[8pt]
1,&\text{其他情形}.
\end{cases}
$$

其中 `tol=1e-12` 用于避免退化分母。期望虚拟动态为

$$
\boxed{
\phi
=\varepsilon(a\sigma-ks_2)
+(1-\varepsilon)\frac{a+b}{2}\sigma,
}
$$

其中 $k>0$。

$\varepsilon$ 在相关边界上按到直线的相对距离插值，使不同区域的表达式连续衔接。控制律没有直接使用 `sign` 函数，因此在连续、未饱和、模型匹配的理想条件下，$u^*$ 是连续函数。

## 7. 趋近律 $N$ 的含义

当前区域的辅助面为

$$
h=\omega_1\sigma+\omega_2\eta+m_h.
$$

在未饱和理想闭环中

$$
\dot h
=\omega_1\dot\sigma+\omega_2\dot\eta
=\omega_1\phi+\omega_2\sigma.
$$

代码定义

$$
\boxed{N=\omega_2\sigma+\omega_1\phi,}
$$

所以 $N=\dot h$。$N$ 不是控制电压，也不是 $\phi$，而是当前激活辅助面的理论导数。

论文在中心对称条件、$\omega_{11}<0$、$\omega_{21}>0$ 和上述 $\phi$ 构造下证明

$$
N(\sigma,\eta)\ge0,
$$

且理想连续模型中只有在 $(\sigma,\eta)=(0,0)$ 时取零。结合论文的分区公共 Lyapunov 函数，可使轨迹趋向原点，并在满足定理条件时使 $Q$ 成为正不变集。

需要区分 $N$ 和 $\dot N$：控制器及理论约束检查的是 $N\ge0$；结果图中的 $\dot N$ 是仿真后对采样数据做数值差分得到的观察量，可以为负。

## 8. 实际电压与饱和

未饱和电压为

$$
u^*=\frac{-SAx+\phi}{SB},
$$

实际仿真输入为

$$
u=\operatorname{sat}_{[-12,12]}(u^*).
$$

若发生饱和，令 $\Delta u=u-u^*$，则

$$
\dot\sigma=\phi+SB\Delta u,
$$

从而 $\dot h=N$ 不再严格成立。因此优化器不仅检查限幅后的电压，还把未饱和需求 $\lvert u^*\rvert\le12\,\mathrm V$ 作为可行性条件。

## 9. 各参数的作用

| 参数 | 设计作用 | 过大或过小的可能影响 |
|---|---|---|
| $p_1,p_2,p_3$ | 决定 $\sigma=0$ 上的三维内部极点 | 太小使内部收敛慢；太大可能增大控制作用和模型敏感性 |
| $\xi_1,\xi_2$ | 决定两条切换面斜率和虚拟平面分区 | 间距太小使几何退化；过大可能使虚拟动态更激进 |
| $\beta$ | 决定侧边辅助面的斜率，并参与派生 $\alpha,\lambda$ | 改变正不变集形状和趋近律斜率 |
| $\alpha$ | 决定 $h_0,h_3$ 斜率 | 由几何约束派生，不独立优化 |
| $\lambda$ | 缩放 $h_1,h_2$ 的状态系数 | 由几何约束派生；不改变对应直线斜率，但会改变直线截距、$Q$ 的形状和 $N$ 的尺度 |
| $k$ | 决定 $a\sigma-ks_2$ 中的趋近强度 | 太小可能收敛慢；太大可能提高电压需求和离散实现压力 |
| $m_h$ | 决定辅助面离原点的偏置及 $Q$ 尺度 | 当前固定为 `0.01`，未纳入优化 |
| $u_{\max}$ | 执行器电压上限 | 当前固定为 $12\,\mathrm V$ |

## 10. 优化变量及约束编码

优化器不直接搜索 $p_1,p_2,p_3,\xi_1,\xi_2,\alpha,\beta,k$。它使用七维正增量向量

$$
q=\begin{bmatrix}
p_1,&p_2-p_1,&p_3-p_2,&\xi_2,&\xi_1-\xi_2,&\beta,&k
\end{bmatrix}.
$$

这样只需限制每一维为正，就能保证

$$
p_1<p_2<p_3,
\qquad
\xi_1>\xi_2>0.
$$

搜索边界为：

| 分量 | 下界 | 上界 |
|---|---:|---:|
| $p_1$ | 0.5 | 5.0 |
| $p_2-p_1$ | 0.25 | 5.0 |
| $p_3-p_2$ | 0.25 | 5.0 |
| $\xi_2$ | 0.2 | 8.0 |
| $\xi_1-\xi_2$ | 0.2 | 12.0 |
| $\beta$ | 0.1 | 8.0 |
| $k$ | 0.2 | 15.0 |

$\alpha$ 和 $\lambda$ 每次都按第 4 节公式解析计算，因此优化过程不会破坏辅助面交点约束。

## 11. 优化工况

默认训练集共 16 个工况：

- 名义参数下的 $+5^\circ,-5^\circ,+8^\circ,-8^\circ$；
- 12 个随机工况，初始摆角在 $[-8^\circ,8^\circ]$ 内均匀抽取；
- 每个随机工况中，$M,m,J,l,R,K_e,K_m,I$ 分别在名义值的 $\pm10\%$ 内独立变化；
- $g$ 保持名义值。

优化结束后再生成 100 个独立随机验证工况，使用不同随机种子。默认种子为 `20260717`，验证种子为 `20260718`，因此结果可重复。

控制器始终使用名义 $A,B,S$ 计算 $u^*$，而被控对象用各摄动工况的实际离散矩阵推进。这对应“名义模型控制器面对参数失配”的数值鲁棒性检查。

## 12. 单工况目标函数

对第 $j$ 个工况定义

$$
\begin{aligned}
J_j={}&0.30\frac{t_{s,j}}{T_{\mathrm{end}}}
+0.15\frac{\mathrm{IAE}_{p,j}}{0.2T_{\mathrm{end}}}
+0.20\frac{\mathrm{IAE}_{\theta,j}}
{\operatorname{rad}(10^\circ)T_{\mathrm{end}}}\\
&+0.10\frac{p_{\max,j}}{0.2}
+0.10\frac{u_{\mathrm{rms},j}}{12}
+0.10\frac{\Delta u_{\mathrm{rms},j}}{12}
+0.05\frac{\theta_{\max,j}}{\operatorname{rad}(10^\circ)}.
\end{aligned}
$$

其中

$$
\mathrm{IAE}_p=\int_0^{T_{\mathrm{end}}}\lvert p(t)\rvert\,dt,
\qquad
\mathrm{IAE}_\theta=\int_0^{T_{\mathrm{end}}}\lvert\theta(t)\rvert\,dt.
$$

目标同时考虑最坏工况和平均工况：

$$
J(q)=0.6\max_jJ_j+0.4\operatorname{mean}_jJ_j+P(q),
$$

其中 $P(q)$ 是约束罚项。

## 13. 可行性条件与罚函数

每个工况必须同时满足：

$$
\begin{aligned}
&\text{全部状态、输入和 }N\text{ 有限},\\
&\theta_{\max}\le10^\circ,\\
&p_{\max}\le0.20\,\mathrm m,\\
&\max_t\lvert u^*(t)\rvert\le12\,\mathrm V,\\
&\lvert p(T_{\mathrm{end}})\rvert<10^{-3}\,\mathrm m,\\
&\lvert\theta(T_{\mathrm{end}})\rvert<0.1^\circ,\\
&\min_tN(t)\ge-10^{-9}.
\end{aligned}
$$

每个不可行工况先增加 $10^6$ 罚值，再按角度、位置、未饱和电压、终值和负 $N$ 的归一化违反量增加 $10^4$ 级罚值。因此收敛图中的群体均值可能比最优可行解高很多个数量级。

## 14. 优化算法

默认首先运行 40 个粒子、60 轮的粒子群算法。第 $j$ 个粒子的速度和位置按

$$
v_j^{(r+1)}
=wv_j^{(r)}
+c_1R_1\left(q_{j,\mathrm{best}}-q_j^{(r)}\right)
+c_2R_2\left(q_{\mathrm{global}}-q_j^{(r)}\right),
$$

$$
q_j^{(r+1)}=q_j^{(r)}+v_j^{(r+1)}
$$

更新，并截断到搜索边界。默认 $w=0.72$、$c_1=c_2=1.49$，速度上限为每一维搜索区间的 $20\%$。

PSO 后执行坐标局部搜索：逐维尝试正负步长，有改进则接受；一轮无改进时步长减半，直到小于区间的 $0.1\%$ 或达到最大轮数。

基线向量

$$
q_{\mathrm{base}}=
\begin{bmatrix}2&1&1&1&1&1&4\end{bmatrix}
$$

被显式放入初始粒子群，并在结束时再次与最优解比较，所以保存的训练目标不会劣于基线目标。

## 15. 当前优化参数

当前保存的最优向量为

$$
q^*\approx
\begin{bmatrix}
4.3453850&2.2695550&4.9080276&5.4331663&12&3.6797011&15
\end{bmatrix}.
$$

解码后：

| 参数 | 当前值 |
|---|---:|
| 内部极点 | $[-4.3453850,-6.6149400,-11.5229676]$ |
| $\xi_1$ | 17.4331663 |
| $\xi_2$ | 5.4331663 |
| $\alpha$ | 9.0510902 |
| $\beta$ | 3.6797011 |
| $\lambda$ | 0.3970127 |
| $a$ | -9.0510902 |
| $b$ | 3.6797011 |
| $k$ | 15 |
| $m_h$ | 0.01 |
| $u_{\max}$ | $12\,\mathrm V$ |

$\xi_1-\xi_2=12$ 和 $k=15$ 正好达到当前搜索上界。这说明“当前边界内的保存解”在这两个方向仍可能希望取更大值，不能据此声称已经找到无约束全局最优参数。若扩大边界，应重新评估电压裕量、采样敏感性和独立验证通过率，而不是只比较训练目标。

## 16. 理论与数字实现的边界

在连续时间、模型匹配、无外扰、全状态反馈且未饱和时，当前设计实现

$$
\dot\sigma=\phi,
\qquad
\dot h=N\ge0,
$$

并可结合论文的分区 Lyapunov 论证和稳定内部动态推出渐近收敛。当前程序实际使用 $200\,\mathrm{Hz}$ 采样、零阶保持、饱和和参数摄动。因此：

- 名义未饱和轨迹与连续理论对应最直接；
- 随机摄动训练与验证是数值鲁棒性证据，不是统一扰动上界证明；
- 发生饱和时，$\dot\sigma=\phi$ 和 $\dot h=N$ 会产生附加误差；
- 小角度模型的结论不能外推到大角度起摆或实机全局稳定性。

## 17. 代码对应关系

| 设计内容 | 实现位置 |
|---|---|
| 切换面、区域、$\varepsilon$、$\phi$、$N$ 和电压 | `uas_smc_controller.m` |
| $\alpha$、$\lambda$、辅助面系数、顶点和残差 | `uas_smc_auxiliary_geometry.m` |
| 名义模型、$S$、采样仿真、指标和绘图 | `run_uas_smc_sim.m` |
| 七维编码、工况、目标函数、PSO 和局部搜索 | `optimize_uas_smc.m` |
| 几何、连续性、对称性和正式结果验收 | `test_uas_smc.m` |
