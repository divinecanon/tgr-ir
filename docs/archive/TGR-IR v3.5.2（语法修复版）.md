# TGR-IR v3.5.2（语法修复版）

## 〇、版本说明

v3.5.2 相对 v3.5.1 的变更（基于逻辑审查修订）：

| 编号 | 严重度 | 变更内容 | 位置 |
|---|---|---|---|
| 1 | 🔴 | 修复 IMMUNE_UNRESOLVED 状态下 C 序列不更新的恢复死锁 | §12.4 |
| 2 | 🔴 | 修复偏转触发与成功条件中双 $G_{\min}$ 阈值导致的卡死状态 | §11.3、§11.5 |
| 3 | 🔴 | 明确终止状态与 DORMANT 的互斥优先级 | §9.1、§九 |
| 4 | 🔴 | 修复三级偏转唤醒 P' 时未定义状态与 C(0) 重置的问题 | §11.4 |
| 5 | 🔴 | 明确 WARMUP 状态下豁免 LOW_BUDGET 评估 | §12.2 |
| 6 | 🟡 | 补充原始 P 层级 $L$ 的推断说明与配对选择规则 | §一、§7.1 |
| 7 | 🟡 | 降级恢复策略补充排序与数量规则 | §12.3 |
| 8 | 🟡 | 统一公式参数命名，区分 $G_{\text{floor}}$ 与 $G_{\text{trigger}}$ | 全文、§十 |
| 9 | 🟡 | 明确 $t=0$ 或序列长度不足时的 T 计算豁免 | §九 |
| 10 | 🟢 | 修复 DORMANT 条件引用错误（§9 → §8） | §9.1 |
| 11 | 🟢 | 修正 _check_auto_immunity 返回值未使用的代码逻辑 | §十三 |
| 12 | 🟢 | 补全 _wake 方法中预算与阈值的重置 | §十三 |

v3.5.1 及更早变更记录核心摘要：
- 增加 AUDIT_EXHAUSTED 与 IMMUNE_UNRESOLVED 恢复条件（v3.5.1）。
- 预算与自免疫增加空值/状态前置条件（v3.5.1）。
- 配对投影分母增加防零保护（v3.5.1）。
- 二级偏转增加 $\Omega=0$ 及范围为零的兜底；三级偏转增加休眠池为空兜底（v3.5.1）。
- 引入动态预算分配与预热期机制（v3.5）。
- 核心公式 (1)~(7.2) 逻辑不变。

---

## 一、原语定义

| 符号 | 名称 | 计算表示 |
|---|---|---|
| $\mathbf{\Omega}$ | 约束 | 矢量目标 $\mathbf{\Omega} \in \mathbb{R}^d$ |
| $\delta$ | 区分 | 序列中的每一步操作（索引 $t$） |
| $\mathbf{C}$ | 状态值 | 矢量 $\mathbf{C}(t) \in \mathbb{R}^d$ |
| $P$ | 主体 | 产生序列的独立单元 |
| $L$ | 维度层级 | 整数：0=探索者，1=观察者，2=审计者 |
| $\mathbf{D}$ | 方向矢量 | 矢量 $\mathbf{\Omega}$ |
| $T$ | 张力 | 同层 P 之间的方向冲突强度，公式 (4.2) |
| $G$ | 生成 | 公式 (2) |
| $R$ | 递归 | 公式 (3) |
| Pair | 配对投影 | 低层 P 经高层中转产生同层 P'，公式 (7.1) |
| 自免疫 | Auto-immunity | 公式 (6) |

**L 的推断与配对规则：**

- **L 赋值**：原始 P 的 $L$ 基于互信息对偶法推断：计算 P 与已知不同层级参考 P 之间的互信息，归入最匹配的层级。若系统无参考，可由输入时指定，或基于历史生成频率粗粒度分级（高频变更为 $L=0$，低频稳态为 $L=1$）。
- **配对选择**：采用全组合配对，即对每层 $L$，穷举所有 $L-1$ 层 P 与 $L$ 层 P 的组合。为防止组合爆炸，当某层 P 数量超过上限（默认 10）时，优先选择历史 $\Omega$ 差异最大的 P 进行配对。
- **维度限定**：当前版本所有公式和实验均以 $d=1$（标量）实现。矢量扩展保留为后续版本目标。$\mathbf{C}(t)$、$\mathbf{\Omega}$ 在当前版本中为实数。公式中的比较运算（$<$、$>$）、$\min$、$\max$、绝对值均对标量定义。

---

## 二、核心公式

以下 $C_i(t) \in \mathbb{R}$，$|\cdot|$ 为绝对值。离散化方案见公式(1.1)注释。

**1. 基础概率（因果，多维离散化）**

$$p_{\text{marg}}(C_i(t)) = \frac{\text{count}(\text{cell}(C_i(t)) \in \text{历史})}{t+1} \qquad (1.1)$$

$$p_{\text{cond}}(C_i(t) \mid C_i(t-1)) = \frac{\text{count}(\text{cell}(C_i(t-1)) \to \text{cell}(C_i(t)))}{\text{count}(\text{cell}(C_i(t-1)) \to \cdot)} \qquad (1.2)$$

分母为零时取 $1\mathrm{e}{-9}$。

$\text{cell}(C)$ 默认方案：使用固定宽度分箱。箱宽 $= (\max_{\text{historical}} C - \min_{\text{historical}} C) / N_{\text{bins}}$。$N_{\text{bins}}$ 默认 $= 10$（当历史数据点 $\ge 10$ 时）；历史数据点 $< 10$ 时，$N_{\text{bins}} =$ 历史数据点数（即每个唯一值为一箱）。实现可替换为动态分位数分箱或其他离散化方案，只需保持 $\text{count}()$ 操作的一致性。

**2. 生成 G**

$$G_i(t) = \begin{cases} \text{None}, & t=0 \\ \max\big(G_{\text{floor}}, -\log_2[p_{\text{marg}} \times p_{\text{cond}}]\big), & t>0 \end{cases} \qquad (2)$$

$G_{\text{floor}}=0.1$（最小生成底线），生成事件基础阈值 $\theta_{G,\text{base}}=2.0$。

**3. 递归 R**

$$R_i(t) = \begin{cases} \text{None}, & t=0 \\ p_{\text{cond}}, & t>0 \end{cases} \qquad (3)$$

**4. 张力 T（仅计算同层冲突）**

$$\bar{C}_i(t) = \frac{1}{t+1}\sum_{k=0}^{t} C_i(k) \qquad (4.1)$$

$$T_{i \leftarrow j}(t) = \frac{|\bar{C}_j(t) - \Omega_i|}{|\Omega_j - \Omega_i|} \qquad (4.2)$$

若 $|\Omega_j - \Omega_i| < 10^{-9}$，$T_{i \leftarrow j}=0$。

$$T_i(t) = \max_{j: L_j = L_i} T_{i \leftarrow j}(t) \qquad (4.3)$$

**5. 成本**

$$c_i(t) = \max\big(c_{\min}, -\log_2 p_{\text{marg}}\big) \qquad (5)$$

$c_{\min}=0.1$。

**6. 自免疫判定**

前置条件：C 序列长度 $\ge N$（默认 $3$），且 G 值序列中最近 $N$ 步均非 None。不满足条件时，自免疫$_i(t) = \text{False}$。

满足前置条件时：

$$\text{自免疫}_i(t) \iff \left(\frac{1}{N}\sum_{k=t-N+1}^{t} G_i(k) < \epsilon_G\right) \land \left(T_i(t) > \tau_T\right) \qquad (6)$$

参数：$N=3$，$\epsilon_G=0.5$，$\tau_T=0.7$。

---

## 三、配对投影与 Ω 推断

**7.1 配对投影**

$$C_{P'}(t) = \frac{C_a(t) - C_b(t)}{\max(|\Omega_b - \Omega_a|, \varepsilon_{\text{denom}})} \qquad (7.1)$$

其中 $\varepsilon_{\text{denom}} = 10^{-9}$。

说明：当 $|\Omega_b - \Omega_a| \approx 0$ 时，分母由防零保护接管。此时若分子也接近零，$C_{P'}(t)$ 接近零，P' 的 C 序列趋于平坦，G 值下降，自然进入休眠。若分母极小但分子不接近零，$C_{P'}(t)$ 会被放大，使 P' 能感知方向高度一致的两个原始 P 间的细微约束差异。

$P'$ 的层级 $L(P') = \max(L_a, L_b)$，与参与配对的更高层 P 同层。

**7.2 三层 fallback Ω 推断（在线）**

参数：
- $\epsilon_{\text{trend}} = 0.05$（趋势检测阈值）
- $\epsilon_{\text{pos}} = 0.05$（位置检测阈值）
- $\Delta = 0.1 \times (\max C_{P'} - \min C_{P'})$，若数据范围为零则 $\Delta = 0.1$

判定流程：

给定 P' 在时刻 t 的历史 C 序列：

```
1. 若 |C(t) - C(0)| > ε_trend：
      （整体趋势显著非零）
      若 C(t) < C(0)（下行）：
          Ω = min(C) - Δ
      否则（上行）：
          Ω = max(C) + Δ
      退出                        → 层级 1（趋势方向）覆盖

2. 若 |C̄(t)| > ε_pos：
      （趋势为零，但位置偏正或偏负）
      若 C̄(t) < 0：
          Ω = min(C) - Δ
      否则：
          Ω = max(C) + Δ
      退出                        → 层级 2（位置符号）覆盖

3. 否则：
      Ω = C̄(t)
      退出                        → 层级 3（全静止 fallback）
```

形式化：

$$\Omega_{P'}(t) = \begin{cases}
\min C_{P'} - \Delta, & |C(t)-C(0)| > \epsilon_{\text{trend}} \land C(t) < C(0) \\
\max C_{P'} + \Delta, & |C(t)-C(0)| > \epsilon_{\text{trend}} \land C(t) > C(0) \\
\min C_{P'} - \Delta, & |C(t)-C(0)| \le \epsilon_{\text{trend}} \land |\bar{C}| > \epsilon_{\text{pos}} \land \bar{C} < 0 \\
\max C_{P'} + \Delta, & |C(t)-C(0)| \le \epsilon_{\text{trend}} \land |\bar{C}| > \epsilon_{\text{pos}} \land \bar{C} > 0 \\
\bar{C}_{P'}(t), & \text{其他}
\end{cases} \qquad (7.2)$$

---

## 四、G 驱动的按需触发

**8. 触发条件**

$$\text{trigger}_i(t) = \begin{cases}
\text{True}, & G_i(t) \ge \theta_G' \quad \text{（生成事件，v3.5 使用动态阈值）} \\
\text{True}, & \left(\frac{1}{k}\sum_{j=t-k+1}^{t} G_i(j) < G_{\text{trigger}}\right) \land \left(T_i(t) > \tau_T\right) \quad \text{（自免疫）} \\
\text{False}, & \text{其他}
\end{cases} \qquad (9)$$

参数：$k=3$，$G_{\text{trigger}}=0.5$（自免疫触发阈值），$\tau_T=0.7$。$\theta_G'$ 为动态触发阈值，见 §10.3。

**9. 调度规则**

- **触发时**：P' 用公式 (7.1) 更新 C 序列，用三层 fallback 更新 Ω，用公式 (2)(3) 更新 G, R。同层活跃 P' 重新计算 T。
- **非触发时**：$C_{P'}(t) = C_{P'}(t-1)$。G, R, T, Ω 保持上一步值。未触发计数器 +1。
- **休眠**：连续 $M$ 步（默认 $20$）未触发 → P' 进入休眠，不参与同层 T 计算。
- **唤醒**：关联低层 P 触发时，从冷存储恢复完整 C 序列，继续更新。

---

## 五、P' 初始化与预热期

**9.1 状态分类**

P' 在任意时刻 $t$ 属于以下四种状态之一：

| 状态 | 条件 | 行为 |
|---|---|---|
| COLD | $t=0$，C 序列长度 = 1 | 仅记录 C(0)，不推断 Ω |
| WARMUP | C 序列长度 ≥ 2 且 < WARMUP_MAX | 每步强制更新，积累数据 |
| ACTIVE | C 序列长度 ≥ WARMUP_MAX | 正常按需触发 |
| DORMANT | 见 §8 调度规则 | 休眠，不参与计算 |

**状态互斥原则**：若 P' 处于审计终止状态（AUDIT_EXHAUSTED 或 AUDIT_IMMUNE_UNRESOLVED），则其状态锁定为"终止态活跃"，不累加未触发计数器，即不会因长期未触发而进入 DORMANT。

参数：$\text{WARMUP\_MAX} = 5$（预热期内积累的最小序列长度）。

**C(0) 来源**：P' 创建时，立即执行公式 (7.1) 一次，产生第一个序列元素，记为 $C_{P'}(0)$。此时 C 序列长度 = 1，状态为 COLD。

**9.2 状态转换规则**

```
COLD → WARMUP：
  C 序列长度首次达到 2（即 t=1 时产生第二个数据点）
  此时开始每步强制更新

WARMUP → ACTIVE：
  C 序列长度达到 WARMUP_MAX
  此时三层 fallback 已可正常生效

唤醒 → 状态判定：
  从冷存储恢复完整 C 序列后：
  - 若恢复后序列长度 ≥ WARMUP_MAX → 直接进入 ACTIVE
  - 否则 → 进入 WARMUP，继续积累

  进入 WARMUP 时（无论是首次创建、唤醒后进入，还是三级偏转强制进入）：
  用于趋势判定的参考点 C(0) 重置为当前恢复序列的第一个值。
  即后续 |C(t) - C(0)| 中的 C(0) 是进入 WARMUP 时的第一个 C 值，
  而非 P' 创建时的原始 C(0)。
  如需保留原始 C(0) 用于长期审计，存储在 C_origin 字段中，不参与趋势判定。
```

注意：预热期不设直接进入 DORMANT 的路径。预热期内每步强制触发更新（见 9.3），必然到达 ACTIVE。

**9.3 预热期内行为**

WARMUP 期间：
- 每步强制触发更新（无论 G 值）
- C 序列：正常通过公式 (7.1) 产生并累积
- Ω 计算：
  - 若 $|C(t) - C(0)| > \epsilon_{\text{trend}}$：
    - 若 $C(t) < C(0)$（下行）：$\Omega = \min(C) - \Delta$
    - 若 $C(t) > C(0)$（上行）：$\Omega = \max(C) + \Delta$
    - （即公式 (7.2) 层级 1 逻辑）
  - 否则：$\Omega = C(t)$（暂不推断方向）
- G, R 正常计算（公式 (2), (3)）
- T 正常计算（公式 (4.2), (4.3)），参与同层 T 计算（Ω 已赋值时）

---

## 六、敏感度预算动态分配

**10.1 预算定义**

系统敏感度预算用于调节各 P' 的触发阈值，使高信息密度的 P' 获得更多更新机会。

预算性质：**软预算**。预算总和不受硬性限制。$B(i)$ 是相对优先级度量，仅用于调节各 P' 的触发阈值。

基础分配（仅在 ACTIVE 状态计算）：

$$B_{\text{base}}(i) = \frac{1}{N_{\text{active}}}$$

其中 $N_{\text{active}}$ = 当前 ACTIVE 状态 P' 总数。

前置条件：预算计算仅在 P' 状态为 ACTIVE 且 $G_i(t)$、$G_i(t-1)$、$T_i(t)$ 均非 None 时执行。COLD 状态不计算预算。WARMUP 状态使用固定预算 $B(i) = 1 / N_{\text{active}}$（即 $\alpha_G = \alpha_T = \alpha_{\text{immune}} = 1.0$）。

**10.2 动态调节因子**

实际预算（ACTIVE 状态）：

$$B(i) = B_{\text{base}}(i) \times \alpha_G(i) \times \alpha_T(i) \times \alpha_{\text{immune}}(i)$$

其中：

$\alpha_G$（生成性因子）：

$$\alpha_G = \begin{cases}
1.0 + \dfrac{G_i(t) - G_i(t-1)}{\max(G_i(t-1), 10^{-9})}, & G_i(t) > G_i(t-1) \\
\max\left(0.3, 1.0 - \dfrac{G_i(t-1) - G_i(t)}{\max(G_i(t-1), 10^{-9})}\right), & \text{否则}
\end{cases}$$

限制：$\alpha_G \in [0.3, 2.0]$

$\alpha_T$（张力因子）：

$$\alpha_T = \begin{cases}
1.0 + \dfrac{T_i(t) - \tau_T}{1.0 - \tau_T}, & T_i(t) > \tau_T \\
1.0, & \text{否则}
\end{cases}$$

限制：$\alpha_T \in [0.5, 2.0]$。$\tau_T = 0.7$（与公式 (6) 一致）。

$\alpha_{\text{immune}}$（自免疫因子，默认值）：

$$\alpha_{\text{immune}} = \begin{cases}
2.0, & \text{自免疫}_i(t) = \text{True} \\
1.0, & \text{否则}
\end{cases}$$

注意：自免疫偏转期间以 §11 偏转策略设定的临时值为准，覆盖此默认值。

**10.3 预算与触发阈值联动**

动态触发阈值（ACTIVE 状态）：

$$\theta_G'(i, t) = \max\left(0.1, \frac{\theta_{G,\text{base}}}{B(i)}\right)$$

其中 $\theta_{G,\text{base}} = 2.0$。WARMUP 状态使用固定阈值 $\theta_G' = \theta_{G,\text{base}} = 2.0$。

效果：
- 预算高 → 阈值降低 → 更容易触发 → 更多资源
- 预算低 → 阈值升高 → 更难触发 → 更快休眠
- $\theta_G'$ 下限为 $0.1$，防止阈值变负或归零。

**10.4 预算耗尽与审计标记**

当 $B(i) < 0.1$ 时：
- 该 P' 进入 LOW_BUDGET 状态
- 更新频率降至每 $\lceil 1 / B(i) \rceil$ 步一次
- 标记：当前分辨率下精细审计不可用
- 不伪装闭合——元数据中记录预算耗尽

---

## 七、自免疫偏转策略

**11.1 自免疫状态跟踪**

每个 P' 维护：
- immune_counter(i)：连续自免疫触发次数
- immune_history(i)：自免疫触发的时间戳列表

```
当 自免疫_i(t) = True：
    immune_counter(i) += 1
    记录时间戳 t
当 自免疫_i(t) = False：
    immune_counter(i) = 0  （重置）
```

自免疫判定条件引用公式 (6)，前置条件满足时，连续 $N=3$ 步平均 $G < \epsilon_G = 0.5$ 且 $T > \tau_T = 0.7$。

**11.2 一级偏转（敏感度激增）**

触发条件：$\text{immune\_counter}(i) = 1$（首次自免疫触发）

动作：
- 临时将 $\alpha_{\text{immune}}$ 设为 $3.0$（覆盖默认值 2.0）
- 强制触发更新（忽略 G 阈值）
- 标记：自免疫一级偏转已执行

预期效果：预算激增 → 触发阈值大幅降低 → 细微外部约束变化被感知 → 若约束存在，G 将回升。

持续时间：直到 immune_counter 归零或进入二级偏转。

**11.3 二级偏转（Ω 收缩）**

触发条件：$\text{immune\_counter}(i) \ge 2$，且自上次偏转执行后 $G_i$ 持续 $< G_{\text{trigger}}$（0.5）。

动作：
- 若 $|\Omega| < \varepsilon_{\text{zero}}$（默认 0.01）：
  - 直接触发跨零保护 → 进入三级偏转，退出
- 计算收缩量：
  - 若 $|\max(C) - \min(C)| < \varepsilon_{\text{zero}}$：
    - $\Delta_{\text{contract}} = \Delta_{\text{default}}$（默认 0.1）
  - 否则：
    - $\Delta_{\text{contract}} = \min(0.2 \times |\max(C) - \min(C)|, \Delta_{\max})$
  - $\max(C)$、$\min(C)$ 取该 P' 全历史 C 序列的最大最小值，$\Delta_{\max}=1.0$
- 收缩方向：
  - 若 $\Omega > 0$：$\Omega_{\text{new}} = \Omega - \Delta_{\text{contract}}$
  - 若 $\Omega < 0$：$\Omega_{\text{new}} = \Omega + \Delta_{\text{contract}}$
- 跨零保护：若 $\Omega_{\text{new}}$ 符号与 $\Omega$ 不同：
  - $\Omega_{\text{new}} = 0$，标记：Ω 跨零，二级偏转强制完成，进入三级偏转
- 降低历史投影权重：下一轮 Merge 中历史 Signal 衰减系数临时加倍
- 标记：自免疫二级偏转已执行，记录收缩前 Ω 值

**11.4 三级偏转（强制引入张力）**

触发条件：$\text{immune\_counter}(i) \ge 4$，且二级偏转后 $G_i$ 持续 $< G_{\text{trigger}}$；或 Ω 跨零后自动进入。

动作：
1. 检查休眠池：
   - 若休眠池非空 → 继续步骤 2
   - 若休眠池为空 → 直接标记 AUDIT_IMMUNE_UNRESOLVED，记录："三级偏转失败，休眠池为空"，结束偏转流程
2. 从休眠 P' 池中随机唤醒一个 P'（优先选择历史 T 值与当前 P' 差异最大的）
3. 设置 Ω：
   - 若被唤醒 P' 的 $|\Omega_{\text{awakened}}| \ge \varepsilon_{\text{zero}}$：
     - $\Omega_{\text{new}} = -\Omega_{\text{awakened}} \times \beta$，$\beta = 0.5$
   - 若 $|\Omega_{\text{awakened}}| < \varepsilon_{\text{zero}}$（兜底）：
     - $\Omega_{\text{new}} = \text{random\_sign}() \times \Delta_{\text{default}}$
4. 当前 P' 强制进入 WARMUP 状态（重新积累序列，C 序列完整保留，重置 C(0) 参考点）。
5. 被唤醒的 P' 同样按 §9.2 唤醒逻辑执行：若序列长度 ≥ WARMUP_MAX 进入 ACTIVE，否则进入 WARMUP。
6. 标记：自免疫三级偏转已执行，记录被唤醒 P' 的 ID 和 Ω 值

**11.5 偏转升级规则**

```
immune_counter = 1    → 一级偏转
immune_counter = 2~3  → 二级偏转（每次执行）
immune_counter ≥ 4    → 三级偏转
Ω 跨零               → 直接三级偏转
```

偏转成功条件：
- 自免疫条件不再成立。
- 即：最近 N 步（默认 3）的平均 G 不再满足（avg(G) < ε_G 且 T > τ_T）。
- 等价于：avg(G) ≥ ε_G 或 T ≤ τ_T。
- 满足时：
  - → 偏转成功
  - → immune_counter 归零
  - → 恢复正常运行

---

## 八、审计边界标记

**12.1 标记类型**

每个 P' 维护 audit_marker 字段，可能的值：

| 标记 | 含义 |
|---|---|
| AUDIT_OK | 审计正常运行 |
| AUDIT_LOW_BUDGET | 预算不足，审计降级 |
| AUDIT_EXHAUSTED | 预算耗尽，审计终止 |
| AUDIT_IMMUNE_UNRESOLVED | 三级偏转无效，不可审计 |

**12.2 标记触发条件**

前置条件：仅当 P' 状态为 ACTIVE 时才评估预算相关审计标记。COLD/WARMUP/DORMANT 状态下保持 AUDIT_OK。WARMUP 状态对 LOW_BUDGET 具有豁免权，即使其固定预算计算结果低于 0.1 也不降级。

**AUDIT_LOW_BUDGET：**
- 进入条件：$B(i) < 0.1$ 持续超过 3 步
- 退出条件：$B(i) \ge 0.15$ 持续超过 2 步（滞回）

**AUDIT_EXHAUSTED（满足任一）：**
- $B(i) < 0.01$ 持续超过 5 步
- 连续 10 步处于 AUDIT_LOW_BUDGET 且未退出

**AUDIT_IMMUNE_UNRESOLVED（满足任一）：**
- 三级偏转后 $\text{immune\_counter}$ 继续增长至 $\ge 6$
- 三级偏转时休眠池为空

标记内容：约束束可能全局闭合，当前分辨率下无法进一步审计，需外部约束或人工介入

**12.3 标记行为**

当 audit_marker 为 AUDIT_EXHAUSTED 或 AUDIT_IMMUNE_UNRESOLVED：
- 保持活跃但不更新（维持最后有效状态）
- 不参与同层 T 计算
- 元数据持续暴露该标记，等待外部约束引入
- 不伪装闭合

状态锁定：不累加未触发计数器，不进入 DORMANT。

**大量节点退出的降级策略：**

若活跃 P' 总数（不含上述两种终止状态）降至 $N_{\min}$ 以下（默认 $N_{\min}=2$）：
- 系统发出全局警告："可审计约束节点不足，系统可能进入低分辨率运行"
- 按进入不可审计状态的时间戳倒序排列，依次将最近的 P' 临时恢复为 AUDIT_LOW_BUDGET 状态，直至活跃 P' 总数达到 $N_{\min}$。

降级恢复 P' 的保护：
- 强制赋予最低预算 $B = 0.15$
- 保持 $K$ 步不降级（默认 $K = 5$）
- 标记为"临时恢复"
- 保护期内不评估 EXHAUSTED 条件；保护期结束后正常评估

**12.4 审计终止状态的恢复条件**

**AUDIT_EXHAUSTED 恢复（同时满足）：**
- 连续 $R$ 步（默认 $5$）$B(i) \ge 0.15$
- 其间至少出现一次生成事件（$G \ge 2.0$）

动作：退回 AUDIT_LOW_BUDGET，恢复更新（每 $\lceil 1/B(i) \rceil$ 步一次），记录恢复信息。

**AUDIT_IMMUNE_UNRESOLVED 恢复（满足任一）：**
- (a) 人工介入标记重置
- (b) 外部约束引入导致关联低层 P 产生显著变化：与该 P' 配对的低层 P 产生生成事件（$G \ge 2.0$）或其 $\Omega$ 发生改变。此时解除 IMMUNE_UNRESOLVED 标记。

动作：重置为 AUDIT_OK，immune_counter 归零，进入 WARMUP（C 序列保留，重置 C(0) 参考点），记录恢复原因和时间戳。

---

## 九、完整计算流程

```
输入: sequences = {P_id: [C_0, ..., C_T]}

1. 推断每个原始 P 的 L 和 Ω（互信息对偶法）
2. 初始化:
   active_P' = {}         # 活跃 P'
   dormant_P' = {}        # 休眠 P'
   inactive_count = {}    # 未触发计数器
   immune_counter = {}    # 自免疫计数器
   audit_marker = {}      # 审计标记
   P'_state = {}          # COLD | WARMUP | ACTIVE | DORMANT

3. FOR t = 0 to T:
      更新所有原始 P 的 G, R, T（公式 1-6）
      【注：若 t=0 或历史长度不足 N，跳过 T 计算与自免疫判定，记 T=0 或 False】

      对每层 L > 0，找 L-1 层 P 与 L 层 P 的配对关系（基于 §一 配对规则）

      对每个 P'：
          更新状态机（§9.2）

          若状态为 ACTIVE 且 G, T 非 None：
              计算预算 B(i)（§10.1-10.2）
              计算动态阈值 θ_G'(i, t)（§10.3）
          若状态为 WARMUP：
              B(i) = 1/N_active
              θ_G'(i, t) = θ_G_base
          若状态为 COLD：
              不计算预算

          若处于 AUDIT_EXHAUSTED 或 AUDIT_IMMUNE_UNRESOLVED：
              跳过触发检查与未触发计数器累加（状态锁定）
              检查恢复条件（§12.4）
          否则：
              检查触发条件（公式 (9)，ACTIVE 用 θ_G'，WARMUP 强制触发）
              若触发或处于 WARMUP（强制更新）：
                  更新 C 序列（公式 7.1）
                  ACTIVE 时用三层 fallback 更新 Ω，WARMUP 时用预热期逻辑（§9.3）
                  更新 G, R（公式 2, 3）
                  检查自免疫（公式 6，需满足前置条件）
                  若自免疫触发：
                      更新 immune_counter（§11.1）
                      执行对应偏转等级（§11.2-§11.4）
                      若三级偏转时休眠池为空 → 标记 AUDIT_IMMUNE_UNRESOLVED
                  否则：
                      immune_counter 归零
              若未触发且非 WARMUP：
                  保持值，计数器 +1，超 M 则休眠（§8 调度规则）

          若状态为 ACTIVE 且非终止状态：
              更新审计标记（§12.2）

          若触发大量退出降级条件则执行降级（§12.3）

      活跃 P' 同层 T 计算（公式 4.2）：
          仅对 WARMUP/ACTIVE 且 Ω 非 None 的 P' 计算 T。
          COLD 状态和 Ω 未赋值的 WARMUP 状态不参与。

4. 输出 results, projected, audit_log
```

---

## 十、参数汇总

| 参数 | 符号 | 默认值 | 定义位置 |
|---|---|---|---|
| 生成事件基础阈值 | $\theta_{G,\text{base}}$ | 2.0 | §二 公式(2) |
| 最小生成底线 | $G_{\text{floor}}$ | 0.1 | §二 公式(2) |
| 自免疫触发阈值 | $G_{\text{trigger}}$ | 0.5 | §四 9, §11.3, §11.4 |
| 最小成本 | $c_{\min}$ | 0.1 | §二 公式(5) |
| 自免疫窗口 | $N$ | 3 | §二 公式(6) |
| 自免疫 G 阈值 | $\epsilon_G$ | 0.5 | §二 公式(6) |
| 自免疫 T 阈值 | $\tau_T$ | 0.7 | §二 公式(6) |
| 趋势检测阈值 | $\epsilon_{\text{trend}}$ | 0.05 | §三 7.2 |
| 位置检测阈值 | $\epsilon_{\text{pos}}$ | 0.05 | §三 7.2 |
| 分母防零 | $\varepsilon_{\text{denom}}$ | $10^{-9}$ | §三 公式(7.1) |
| 休眠阈值 | $M$ | 20 | §四 9 |
| 预热期最小长度 | WARMUP_MAX | 5 | §五 9.1 |
| 预算下限（LOW_BUDGET） | $B_{\min}$ | 0.1 | §六 10.4 |
| $\alpha_G$ 范围 | - | [0.3, 2.0] | §六 10.2 |
| $\alpha_T$ 范围 | - | [0.5, 2.0] | §六 10.2 |
| $\alpha_{\text{immune}}$ 默认值 | - | 2.0 | §六 10.2 |
| $\alpha_{\text{immune}}$ 偏转值 | - | 3.0 | §七 11.2 |
| 触发阈值下限 | $\min(\theta_G')$ | 0.1 | §六 10.3 |
| $\Delta_{\text{contract}}$ 上限 | $\Delta_{\max}$ | 1.0 | §七 11.3 |
| 零值检测精度 | $\varepsilon_{\text{zero}}$ | 0.01 | §七 11.3 |
| 兜底 $\Delta$ | $\Delta_{\text{default}}$ | 0.1 | §七 11.3 |
| 三级偏转反转向量 | $\beta$ | 0.5 | §七 11.4 |
| 审计滞回进入阈值 | - | 0.1 | §八 12.2 |
| 审计滞回退出阈值 | - | 0.15 | §八 12.2 |
| 预算耗尽下限 | - | 0.01 | §八 12.2 |
| 最小活跃节点数 | $N_{\min}$ | 2 | §八 12.3 |
| 降级恢复保护步数 | $K$ | 5 | §八 12.3 |
| EXHAUSTED 恢复步数 | $R$ | 5 | §八 12.4 |
| 不可解恢复 C 变化阈值 | $\varepsilon_{\text{recover}}$ | 0.3 | §八 12.4 |
| 离散化箱数 | $N_{\text{bins}}$ | 10 | §二 公式(1.1)注释 |

---

## 十一、实验设计

### 实验 0：基线复现（v3.4 三场景）

目的：确认 v3.5 未破坏 v3.4 的已有功能

**场景 A（同向无冲突）：**
```
P_a: [0.1, 0.2, 0.3, ..., 2.0] × 20 步
P_b: [0.2, 0.3, 0.4, ..., 2.1] × 20 步
→ 预期：T < 0.1，fallback 层级 1 命中
```

**场景 B（背向冲突）：**
```
P_a: [0.1, 0.2, 0.3, ..., 2.0]
P_b: [2.0, 1.9, 1.8, ..., 0.1]
→ 预期：T > 0.5，fallback 层级 1 命中
```

**场景 C（冻结对立）：**
```
P_a: [0.5] × 20 步
P_b: [-0.5] × 20 步
→ 预期：T > 0.7，fallback 层级 2 命中
```

### 实验 1：初始化与预热期

**场景 D（冷启动）：**
```
步骤 0-4：P_a 和 P_b 正常生成，P' 刚创建
验证：
  · t=0：P' 状态为 COLD，C(0) 由公式 (7.1) 产生，Ω 未定义
  · t=1：C 序列长度=2，进入 WARMUP，Ω 用预热期逻辑
  · t=4：C 序列长度=5，进入 ACTIVE，Ω 由三层 fallback 正常推断
  · 预热期内每步强制更新
```

**场景 E（唤醒恢复）：**
```
步骤 0-9：P' 正常运行后进入休眠
步骤 10：关联低层 P 触发 → P' 被唤醒
验证：
  · 唤醒时从冷存储恢复完整 C 序列（长度=10）
  · 恢复后序列长度 ≥ 5 → 直接进入 ACTIVE
  · Ω 从恢复的序列中重新推断
  · 若唤醒进入 WARMUP，C(0) 参考点重置
```

### 实验 2：敏感度预算动态分配

**场景 F（多 P' 竞争预算）：**
```
创建 5 个 P'：
  P'_1: 高 G（~3.0），高 T（~0.6）
  P'_2: 低 G（~0.3），低 T（~0.1）
  P'_3: 高 G（~2.5），低 T（~0.2）
  P'_4: 中 G（~1.0），中 T（~0.4）
  P'_5: 低 G（~0.2），高 T（~0.8）← 若触发自免疫则 α_immune=2.0

运行 20 步，验证：
  · P'_1 的 θ_G' 降低最多（高 G + 高 T 拉高预算）
  · P'_2 的 θ_G' 升高最多（低 G + 低 T 压低预算）
  · P'_5 若触发自免疫则获得 α_immune 加成
  · B(i) 总和不必等于 1.0（软预算）
```

### 实验 3：自免疫三级偏转

**场景 G（偏转升级）：**
```
构造会进入自免疫的 P'：
  步骤 0-3：正常运行（G > 0.5）
  步骤 4-6：连续 G < 0.5，T > 0.7（连续 3 步，满足公式(6)前置条件）
            → 自免疫触发，immune_counter = 1
  验证一级偏转：
    · α_immune = 3.0（覆盖默认值 2.0）
    · θ_G' 大幅降低（= max(0.1, 2.0 / B)）
    · 强制触发更新

  步骤 7-9：假设外部无新约束，G 持续 < 0.5
            → 自免疫再次触发，immune_counter = 2
  验证二级偏转：
    · Ω 收缩（向零靠近）
    · 收缩量 = min(0.2 × |max(C) - min(C)|, 1.0)
    · 若 |Ω| < 0.01 → 跨零保护 → 直接三级偏转
    · 若 |max(C)-min(C)| < 0.01 → Δ_contract = 0.1

  步骤 10-12：假设仍然无效，immune_counter = 4
  验证三级偏转：
    · 检查休眠池是否为空
    · 若非空：唤醒 P'，设置 Ω
    · 若空：直接标记 AUDIT_IMMUNE_UNRESOLVED
    · 进入 WARMUP，C 序列完整保留

  步骤 13+：
    · 若偏转成功条件满足（自免疫条件不再成立）→ immune_counter 归零
    · 若 immune_counter 达到 6 → AUDIT_IMMUNE_UNRESOLVED
```

### 实验 4：审计边界标记

**场景 H（预算耗尽链）：**
```
构造持续低预算的 P'（ACTIVE 状态）：
  · B(i) < 0.1 持续 3 步 → AUDIT_LOW_BUDGET
  · B(i) ≥ 0.15 持续 2 步 → 退出 LOW_BUDGET（验证滞回）
  · B(i) < 0.01 持续 5 步 → AUDIT_EXHAUSTED
  · 此后 P' 保持但不更新，不参与 T 计算
  · B(i) ≥ 0.15 持续 5 步且出现生成事件 → 退回 LOW_BUDGET（验证恢复）
```

**场景 I（自免疫不可解）：**
```
场景 G 的延伸：
  · 三级偏转后 immune_counter 达到 6
  · → AUDIT_IMMUNE_UNRESOLVED
  · P' 退出 T 计算，等待外部约束
  · |C(t) - C(t_last_immune)| > 0.3 → 自动恢复（验证恢复条件 b）
```

**场景 J（大量退出降级）：**
```
创建 5 个 P'，其中 4 个先后进入 AUDIT_EXHAUSTED
  · 活跃 P' 数降至 1（< N_min = 2）
  · → 触发全局警告
  · → 最近进入不可审计的 P' 被临时恢复为 AUDIT_LOW_BUDGET
  · 恢复的 P'：B = 0.15，K=5 步内不评估 EXHAUSTED
```

### 通过标准汇总

| 实验 | 关键验证点 | 通过标准 |
|---|---|---|
| 0A | 同向 T 值 | T < 0.1 |
| 0B | 背向 T 值 | T > 0.5 |
| 0C | 冻结 T 值 | T > 0.7，fallback L2 |
| 1D | 冷启动状态链 | COLD→WARMUP→ACTIVE |
| 1E | 唤醒恢复 | 直接 ACTIVE，C 序列完整 |
| 2F | 预算分配 | 高 G+T P' 获得更低 θ_G' |
| 3G | 三级偏转 | 各等级触发正确，跨零保护生效，空池兜底生效 |
| 4H | 预算耗尽与恢复 | 标记时机正确，滞回生效，恢复条件触发 |
| 4I | 不可解标记与恢复 | immune_counter ≥ 6 → 标记，C 变化 → 恢复 |
| 4J | 大量退出降级 | N_active < 2 触发警告，保护期防抖动 |

---

## 十三、实现建议

Python 核心结构（已根据审查意见修正错误）：

```python
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Dict, Optional
import numpy as np


class PState(Enum):
    COLD = "COLD"
    WARMUP = "WARMUP"
    ACTIVE = "ACTIVE"
    DORMANT = "DORMANT"


class AuditMarker(Enum):
    AUDIT_OK = "AUDIT_OK"
    AUDIT_LOW_BUDGET = "AUDIT_LOW_BUDGET"
    AUDIT_EXHAUSTED = "AUDIT_EXHAUSTED"
    AUDIT_IMMUNE_UNRESOLVED = "AUDIT_IMMUNE_UNRESOLVED"


@dataclass
class ProjectionP:
    id: str
    paired_a: str
    paired_b: str
    C_sequence: List[float] = field(default_factory=list)
    C_origin: Optional[float] = None          # 创建时的原始 C(0)，用于长期审计
    C0_ref: Optional[float] = None            # 当前阶段用于趋势判定的参考 C(0)
    G_history: List[Optional[float]] = field(default_factory=list)  # 记录所有历史 G 值
    omega: Optional[float] = None
    G: Optional[float] = None
    R: Optional[float] = None
    T: float = 0.0
    state: PState = PState.COLD
    budget: float = 0.0
    theta_G_prime: float = 2.0
    immune_counter: int = 0
    immune_history: List[int] = field(default_factory=list)
    audit_marker: AuditMarker = AuditMarker.AUDIT_OK
    inactive_count: int = 0
    low_budget_steps: int = 0
    low_budget_exit_steps: int = 0
    exhausted_steps: int = 0
    exhausted_recovery_steps: int = 0
    last_deflection_G: Optional[float] = None
    last_immune_C: Optional[float] = None     # 进入 IMMUNE_UNRESOLVED 时的 C 值
    temp_recovered: bool = False
    temp_recovered_steps: int = 0


class TGRSystem:
    def __init__(self, config: dict = None):
        self.config = config or self.default_config()
        self.primitive_ps: Dict[str, List[float]] = {}
        self.projected_ps: Dict[str, ProjectionP] = {}
        self.dormant_ps: Dict[str, ProjectionP] = {}
        self.global_warning = False

    @staticmethod
    def default_config() -> dict:
        return {
            # 核心参数
            'theta_G_base': 2.0,
            'theta_G_min': 0.1,
            'G_floor': 0.1,           # 公式(2)的 G_floor (原 G_min_base)
            'G_trigger': 0.5,         # 公式(9)及偏转触发用的 G_trigger (原 G_min)
            'epsilon_G': 0.5,         # 自免疫 G 阈值
            'tau_T': 0.7,
            'N_immune': 3,
            'M_dormant': 20,
            'epsilon_trend': 0.05,
            'epsilon_pos': 0.05,
            'epsilon_denom': 1e-9,
            # 预热期
            'WARMUP_MAX': 5,
            # 预算
            'alpha_G_min': 0.3,
            'alpha_G_max': 2.0,
            'alpha_T_min': 0.5,
            'alpha_T_max': 2.0,
            'alpha_immune_default': 2.0,
            'alpha_immune_deflection': 3.0,
            'B_min': 0.1,
            'B_exhausted': 0.01,
            'B_recovery': 0.15,
            # 偏转
            'delta_contract_ratio': 0.2,
            'delta_contract_max': 1.0,
            'delta_default': 0.1,
            'beta_reverse': 0.5,
            'epsilon_zero': 0.01,
            # 审计
            'N_min_active': 2,
            'audit_enter_steps': 3,
            'audit_exit_steps': 2,
            'audit_exhausted_steps': 5,
            'audit_exhausted_chain': 10,
            'exhausted_recovery_steps': 5,
            'epsilon_recover': 0.3,
            'temp_recover_budget': 0.15,
            'temp_recover_steps': 5,
            # 离散化
            'N_bins': 10,
        }

    def create_projection(self, id: str, a_id: str, b_id: str,
                           omega_a: float, omega_b: float) -> ProjectionP:
        """创建新 P'，立即计算 C(0)"""
        p = ProjectionP(id=id, paired_a=a_id, paired_b=b_id)
        C0 = self._pair_projection(0, a_id, b_id, omega_a, omega_b)
        p.C_sequence.append(C0)
        p.C_origin = C0
        p.C0_ref = C0   # 初始参考点
        return p

    def _pair_projection(self, t: int, a_id: str, b_id: str,
                          omega_a: float, omega_b: float) -> float:
        C_a = self.primitive_ps[a_id][t]
        C_b = self.primitive_ps[b_id][t]
        denom = max(abs(omega_b - omega_a), self.config['epsilon_denom'])
        return (C_a - C_b) / denom

    def step(self, t: int):
        self._update_primitive_ps(t)
        self._update_projected_ps(t)

    def _update_projected_ps(self, t: int):
        N_active = sum(1 for p in self.projected_ps.values()
                       if p.state == PState.ACTIVE)
        if N_active == 0:
            N_active = 1  # 防止除零

        for p in self.projected_ps.values():
            self._update_state(p)

            # 终止状态锁定：不更新，不累加计数器，仅检查恢复条件
            if p.audit_marker in [AuditMarker.AUDIT_EXHAUSTED, AuditMarker.AUDIT_IMMUNE_UNRESOLVED]:
                self._check_termination_recovery(p, t)
                continue

            if p.state == PState.ACTIVE and p.G is not None and p.T is not None:
                p.budget = self._compute_budget(p, N_active)
                p.theta_G_prime = max(
                    self.config['theta_G_min'],
                    self.config['theta_G_base'] / max(p.budget, 1e-9)
                )
            elif p.state == PState.WARMUP:
                p.budget = 1.0 / N_active
                p.theta_G_prime = self.config['theta_G_base']

            triggered = self._check_trigger(p, t)
            in_warmup = (p.state == PState.WARMUP)

            if triggered or in_warmup:
                self._update_C(p, t)
                self._update_omega(p, t)
                self._update_G_R(p, t)

                # 修复：正确处理自免疫返回值与计数器更新
                is_immune = self._check_auto_immunity(p, t)
                if is_immune:
                    p.immune_counter += 1
                    p.immune_history.append(t)
                    self._execute_deflection(p, t)  # 执行偏转策略 (略)
                else:
                    p.immune_counter = 0
            else:
                p.inactive_count += 1
                if p.inactive_count >= self.config['M_dormant']:
                    self._dormant(p)

            if p.state == PState.ACTIVE:
                self._update_audit_marker(p, t)

            if p.temp_recovered:
                p.temp_recovered_steps += 1
                if p.temp_recovered_steps >= self.config['temp_recover_steps']:
                    p.temp_recovered = False
                    p.temp_recovered_steps = 0

        self._check_min_active()
        self._update_T_all()

    def _update_G_R(self, p: ProjectionP, t: int):
        # 计算 G,R 并追加到 G_history
        # ...（具体逻辑略）
        # 最终 p.G_history.append(p.G)
        pass

    def _check_auto_immunity(self, p: ProjectionP, t: int) -> bool:
        """公式 (6)，含前置条件"""
        N = self.config['N_immune']
        if len(p.C_sequence) < N:
            return False
        if len(p.G_history) < N:
            return False
        recent_G = p.G_history[-N:]
        if any(g is None for g in recent_G):
            return False
        avg_G = sum(recent_G) / N
        return avg_G < self.config['epsilon_G'] and p.T > self.config['tau_T']

    def _update_T_all(self):
        """仅对 WARMUP/ACTIVE 且 Ω 非 None 的 P' 计算 T"""
        eligible = [p for p in self.projected_ps.values()
                    if p.state in (PState.WARMUP, PState.ACTIVE)
                    and p.omega is not None]
        # ... T 计算逻辑（省略）

    def _wake(self, p: ProjectionP):
        """唤醒"""
        if len(p.C_sequence) >= self.config['WARMUP_MAX']:
            p.state = PState.ACTIVE
        else:
            p.state = PState.WARMUP
            # 重置趋势判定参考点为当前序列的第一个值
            if p.C_sequence:
                p.C0_ref = p.C_sequence[0]

        # 修复：恢复预算和阈值
        N_active = sum(1 for x in self.projected_ps.values() if x.state == PState.ACTIVE)
        p.budget = 1.0 / max(1, N_active)
        p.theta_G_prime = self.config['theta_G_base']
        p.inactive_count = 0
        p.audit_marker = AuditMarker.AUDIT_OK

        self.projected_ps[p.id] = p
        del self.dormant_ps[p.id]
```

---

## 十四、版本状态

**版本**：v3.5.2 修订版

**定位**：TGR 约束推断与按需触发的工程规范（可独立实现）

**适用场景**：标量状态（d=1）、矢量目标、任意层级深度、在线/离线均可

**已修复漏洞**：25 项（v3.5.1 的 13 项 + v3.5.2 的 12 项），彻底消除逻辑死锁与命名冲突

**已知限制**：
- 预算调节因子需根据实际运行数据校准
- 三级偏转的阈值（2/4/6）可能需要调整
- 休眠-唤醒机制在长序列（≥ 30 步）上待充分验证
- 非冻结低 G 移动场景下 v3.3/v3.4 近似等价性待量化
- 矢量扩展（d>1）保留为后续版本目标

**下一步**：Python 完整实现与全部 10 项实验验证
