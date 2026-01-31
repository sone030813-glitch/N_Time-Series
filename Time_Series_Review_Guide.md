# Time Series Analysis - Comprehensive Review Guide
# 时间序列分析 - 完整复习指南

---

## Table of Contents | 目录

### Part I: Fundamentals | 基础部分
1. [Weak Stationarity (弱平稳性)](#1-weak-stationarity-弱平稳性)
2. [White Noise (白噪声)](#2-white-noise-白噪声)
3. [Backshift Operator (后移算子)](#3-backshift-operator-后移算子)
4. [Autocovariance and Autocorrelation (自协方差与自相关)](#4-autocovariance-and-autocorrelation-自协方差与自相关)

### Part II: Model Identification | 模型识别
5. [ACF and PACF (自相关函数与偏自相关函数)](#5-acf-and-pacf-自相关函数与偏自相关函数)
6. [Model Identification Strategy (模型识别策略)](#6-model-identification-strategy-模型识别策略)

### Part III: AR, MA, and ARMA Models | AR、MA和ARMA模型
7. [Moving Average (MA) Models (移动平均模型)](#7-moving-average-ma-models-移动平均模型)
8. [Autoregressive (AR) Models (自回归模型)](#8-autoregressive-ar-models-自回归模型)
9. [ARMA Models (ARMA模型)](#9-arma-models-arma模型)

### Part IV: Model Properties | 模型性质
10. [Causality (因果性/平稳性)](#10-causality-因果性平稳性)
11. [Invertibility (可逆性)](#11-invertibility-可逆性)
12. [Relationship Between Causality and Invertibility (因果性与可逆性的关系)](#12-relationship-between-causality-and-invertibility-因果性与可逆性的关系)

### Part V: Non-stationary Models | 非平稳模型
13. [ARIMA Models (ARIMA模型)](#13-arima-models-arima模型)
14. [Differencing (差分)](#14-differencing-差分)
15. [Unit Root and Integration (单位根与单整)](#15-unit-root-and-integration-单位根与单整)

### Part VI: Linear Filters | 线性滤波器
16. [Linear Filters and Stationarity (线性滤波器与平稳性)](#16-linear-filters-and-stationarity-线性滤波器与平稳性)

### Part VII: Practice Problems | 练习题
17. [Worked Examples (详解例题)](#17-worked-examples-详解例题)

---

## Part I: Fundamentals | 基础部分

### 1. Weak Stationarity (弱平稳性)

#### Definition | 定义
A stochastic process $\{X_t\}$ is **weakly stationary** (or covariance stationary) if it satisfies three conditions:

**Condition 1: Constant Mean (均值恒定)**
$$E[X_t] = \mu \quad \text{(constant, independent of } t\text{)}$$

**Condition 2: Finite and Constant Variance (方差有限且恒定)**
$$\text{Var}(X_t) = \sigma^2 < \infty \quad \text{(constant, independent of } t\text{)}$$

**Condition 3: Autocovariance Depends Only on Lag (自协方差只依赖于滞后)**
$$\text{Cov}(X_t, X_{t+h}) = \gamma(h) \quad \text{(depends only on lag } h\text{, not on } t\text{)}$$

#### Key Points | 要点
- All three conditions must be satisfied simultaneously
- Finite variance is crucial: if $\text{Var}(X_t) = \infty$, the process is NOT weakly stationary
- Common mistake: forgetting to check if variance is finite (e.g., Cauchy distribution)

#### Examples | 例子

**Stationary:**
- White noise: $w_t \sim \text{WN}(0, \sigma^2)$
- AR(1) with $|\phi| < 1$: $X_t = \phi X_{t-1} + w_t$
- MA(q) models (always stationary)

**Non-stationary:**
- Random walk: $X_t = X_{t-1} + w_t$ (variance grows with time: $\text{Var}(X_t) = t\sigma^2$)
- Trend: $X_t = \beta t + w_t$ (mean depends on $t$)
- Cauchy i.i.d.: $\text{Var}(X_t) = \infty$

---

### 2. White Noise (白噪声)

#### Definition | 定义
A process $\{w_t\}$ is **white noise** if:
1. $E[w_t] = 0$ for all $t$
2. $\text{Var}(w_t) = \sigma^2$ for all $t$
3. $\text{Cov}(w_t, w_s) = 0$ for all $t \neq s$

Notation: $w_t \sim \text{WN}(0, \sigma^2)$

#### Properties | 性质
- White noise is weakly stationary
- ACF: $\rho(h) = 0$ for all $h \neq 0$
- Building block for time series models

#### Gaussian White Noise | 高斯白噪声
If additionally $w_t \sim N(0, \sigma^2)$ and independent, then it's **Gaussian white noise**

---

### 3. Backshift Operator (后移算子)

#### Definition | 定义
The backshift operator $B$ is defined as:
$$B X_t = X_{t-1}$$
$$B^k X_t = X_{t-k}$$

#### Properties | 性质
- $B^0 = 1$ (identity)
- $B^j B^k = B^{j+k}$
- $(1-B)X_t = X_t - X_{t-1}$ (first difference)
- $(1-B)^d X_t$ represents $d$-th order differencing

#### Polynomial Representation | 多项式表示
- $\phi(B) = 1 - \phi_1 B - \phi_2 B^2 - \cdots - \phi_p B^p$ (AR polynomial)
- $\theta(B) = 1 + \theta_1 B + \theta_2 B^2 + \cdots + \theta_q B^q$ (MA polynomial)

---

### 4. Autocovariance and Autocorrelation (自协方差与自相关)

#### Autocovariance Function (ACF) | 自协方差函数
$$\gamma(h) = \text{Cov}(X_t, X_{t+h}) = E[(X_t - \mu)(X_{t+h} - \mu)]$$

Properties:
- $\gamma(0) = \text{Var}(X_t)$
- $\gamma(h) = \gamma(-h)$ (symmetric)
- $|\gamma(h)| \leq \gamma(0)$

#### Autocorrelation Function (ACF) | 自相关函数
$$\rho(h) = \frac{\gamma(h)}{\gamma(0)} = \frac{\text{Cov}(X_t, X_{t+h})}{\text{Var}(X_t)}$$

Properties:
- $\rho(0) = 1$
- $-1 \leq \rho(h) \leq 1$
- $\rho(h) = \rho(-h)$

#### Partial Autocorrelation Function (PACF) | 偏自相关函数
$$\phi_{hh} = \text{Corr}(X_t, X_{t+h} | X_{t+1}, \ldots, X_{t+h-1})$$

The correlation between $X_t$ and $X_{t+h}$ after removing the linear effect of the intervening variables.

---

## Part II: Model Identification | 模型识别

### 5. ACF and PACF (自相关函数与偏自相关函数)

#### Key Patterns for Model Identification | 识别模型的关键模式

| Model | ACF Pattern | PACF Pattern |
|-------|-------------|--------------|
| **AR(p)** | Tails off (exponential decay or damped oscillation) | Cuts off after lag $p$ |
| **MA(q)** | Cuts off after lag $q$ | Tails off |
| **ARMA(p,q)** | Tails off | Tails off |

#### Terminology | 术语说明

**Cuts off (截尾):** 
- After lag $k$, all values fall within the confidence interval
- Indicates insignificant correlations beyond lag $k$

**Tails off (拖尾):**
- Values gradually decrease (may oscillate)
- Multiple lags remain significant
- Pattern: exponential decay or damped oscillation

#### Confidence Intervals | 置信区间
- Approximate 95% confidence bands: $\pm \frac{1.96}{\sqrt{n}}$
- Values outside these bounds are considered significant

#### Important Notes | 重要注意事项
- A value being negative does NOT mean it's insignificant
- Check if the value exceeds the confidence bounds (regardless of sign)
- First lag is often the most informative

---

### 6. Model Identification Strategy (模型识别策略)

#### Step-by-Step Approach | 逐步方法

**Step 1: Check Stationarity (检查平稳性)**
- Plot the time series
- Check if mean and variance appear constant
- Use statistical tests (e.g., ADF test)
- If non-stationary → consider differencing

**Step 2: Plot ACF and PACF (绘制ACF和PACF图)**
- Look for cutoff patterns
- Identify the dominant pattern

**Step 3: Identify Model Order (识别模型阶数)**
- **If PACF cuts off at lag $p$** → AR(p)
- **If ACF cuts off at lag $q$** → MA(q)
- **If both tail off** → ARMA(p,q)

**Step 4: Fit the Model (拟合模型)**
- Estimate parameters
- Check residuals for white noise

**Step 5: Model Diagnostics (模型诊断)**
- Residual ACF should be insignificant
- Ljung-Box test for residual autocorrelation
- Compare models using AIC/BIC

#### Common Pitfalls | 常见陷阱
- Confusing cutoff vs. tailing off
- Ignoring the magnitude of significant lags
- Not checking if variance is constant over time

---

## Part III: AR, MA, and ARMA Models | AR、MA和ARMA模型

### 7. Moving Average (MA) Models (移动平均模型)

#### MA(1) Model | MA(1)模型
$$X_t = w_t + \theta w_{t-1}$$
where $w_t \sim \text{WN}(0, \sigma^2)$

**Properties:**
- Mean: $E[X_t] = 0$
- Variance: $\text{Var}(X_t) = \sigma^2(1 + \theta^2)$
- ACF: 
  - $\rho(1) = \frac{\theta}{1 + \theta^2}$
  - $\rho(h) = 0$ for $h \geq 2$
- PACF: Tails off exponentially

#### MA(q) Model | MA(q)模型
$$X_t = w_t + \theta_1 w_{t-1} + \theta_2 w_{t-2} + \cdots + \theta_q w_{t-q}$$

Backshift notation:
$$X_t = \theta(B) w_t = (1 + \theta_1 B + \cdots + \theta_q B^q) w_t$$

**Properties:**
- **Always weakly stationary** (no conditions on $\theta$ parameters)
- Mean: $E[X_t] = 0$
- Variance: $\text{Var}(X_t) = \sigma^2(1 + \theta_1^2 + \cdots + \theta_q^2)$
- ACF cuts off after lag $q$
- PACF tails off

**Key Point:** MA models are stationary regardless of parameter values!

---

### 8. Autoregressive (AR) Models (自回归模型)

#### AR(1) Model | AR(1)模型
$$X_t = \phi X_{t-1} + w_t$$

**Properties:**
- **Stationary condition:** $|\phi| < 1$
- Mean: $E[X_t] = 0$
- Variance: $\text{Var}(X_t) = \frac{\sigma^2}{1 - \phi^2}$ (for $|\phi| < 1$)
- ACF: $\rho(h) = \phi^h$ (exponential decay if $\phi > 0$, damped oscillation if $\phi < 0$)
- PACF: $\phi_{11} = \phi$, $\phi_{hh} = 0$ for $h > 1$ (cuts off after lag 1)

#### AR(p) Model | AR(p)模型
$$X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \cdots + \phi_p X_{t-p} + w_t$$

Backshift notation:
$$\phi(B) X_t = w_t$$
where $\phi(B) = 1 - \phi_1 B - \phi_2 B^2 - \cdots - \phi_p B^p$

**Properties:**
- Mean: $E[X_t] = 0$
- ACF: Tails off (exponential decay or damped oscillation)
- PACF: Cuts off after lag $p$

#### AR(2) Model Example | AR(2)模型例子
$$X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + w_t$$

Characteristic equation:
$$1 - \phi_1 z - \phi_2 z^2 = 0$$

**Stationarity condition:** Both roots must lie outside the unit circle (i.e., $|z| > 1$)

---

### 9. ARMA Models (ARMA模型)

#### ARMA(p,q) Model | ARMA(p,q)模型
$$\phi(B) X_t = \theta(B) w_t$$

Expanded form:
$$X_t = \phi_1 X_{t-1} + \cdots + \phi_p X_{t-p} + w_t + \theta_1 w_{t-1} + \cdots + \theta_q w_{t-q}$$

**Properties:**
- Combines AR and MA components
- ACF: Tails off
- PACF: Tails off
- More flexible than pure AR or MA models

#### When to Use ARMA | 何时使用ARMA
- When both ACF and PACF tail off
- Often more parsimonious than pure AR or MA
- Can model a wider range of autocorrelation structures

---

## Part IV: Model Properties | 模型性质

### 10. Causality (因果性/平稳性)

#### Definition | 定义
A process $X_t$ is **causal** (or stationary) if it can be written as a one-sided infinite MA representation:
$$X_t = \sum_{j=0}^{\infty} \psi_j w_{t-j}$$
where $\sum_{j=0}^{\infty} |\psi_j| < \infty$

#### For AR Models | 对于AR模型

**Causality Condition:**
The characteristic equation $\phi(z) = 0$ must have all roots **outside the unit circle** (i.e., $|z| > 1$).

$$\phi(z) = 1 - \phi_1 z - \phi_2 z^2 - \cdots - \phi_p z^p = 0$$

**How to check:**
1. Find the roots of the characteristic equation
2. Compute $|z_i|$ for each root
3. If all $|z_i| > 1$ → Causal ✓
4. If any $|z_i| \leq 1$ → Not causal ✗

#### Example: AR(1)
$$X_t = \phi X_{t-1} + w_t$$

Characteristic equation: $1 - \phi z = 0$ → $z = \frac{1}{\phi}$

Causal if $|z| = \left|\frac{1}{\phi}\right| > 1$ → $|\phi| < 1$ ✓

#### Example: AR(2)
$$X_t = 0.5 X_{t-1} + 0.6 X_{t-2} + w_t$$

Characteristic equation: $1 - 0.5z - 0.6z^2 = 0$
- Roots: $z_1 \approx 0.94$, $z_2 \approx -1.77$
- $|z_1| = 0.94 < 1$ → Not causal ✗

---

### 11. Invertibility (可逆性)

#### Definition | 定义
A process $X_t$ is **invertible** if it can be written as a one-sided infinite AR representation:
$$w_t = \sum_{j=0}^{\infty} \pi_j X_{t-j}$$
where $\sum_{j=0}^{\infty} |\pi_j| < \infty$

In other words, we can recover the white noise $w_t$ from past observations of $X_t$.

#### For MA Models | 对于MA模型

**Invertibility Condition:**
The characteristic equation $\theta(z) = 0$ must have all roots **outside the unit circle** (i.e., $|z| > 1$).

$$\theta(z) = 1 + \theta_1 z + \theta_2 z^2 + \cdots + \theta_q z^q = 0$$

**How to check:**
1. Find the roots of the MA characteristic equation
2. Compute $|z_i|$ for each root
3. If all $|z_i| > 1$ → Invertible ✓
4. If any $|z_i| \leq 1$ → Not invertible ✗

#### Example: MA(1)
$$X_t = w_t + \theta w_{t-1}$$

Characteristic equation: $1 + \theta z = 0$ → $z = -\frac{1}{\theta}$

Invertible if $|z| = \left|\frac{1}{\theta}\right| > 1$ → $|\theta| < 1$ ✓

**If invertible:**
$$w_t = X_t - \theta w_{t-1} = X_t - \theta X_{t-1} + \theta^2 X_{t-2} - \cdots = \sum_{j=0}^{\infty} (-\theta)^j X_{t-j}$$

---

### 12. Relationship Between Causality and Invertibility (因果性与可逆性的关系)

#### Key Distinctions | 关键区别

| Property | What it checks | Applies to | Condition |
|----------|---------------|------------|-----------|
| **Causality** | Stationarity | **AR part** $\phi(B)$ | Roots of $\phi(z)=0$ outside unit circle |
| **Invertibility** | Can express $w_t$ from past $X_t$ | **MA part** $\theta(B)$ | Roots of $\theta(z)=0$ outside unit circle |

#### For Different Models | 对于不同模型

**Pure AR(p) Models:**
- ✓ Need to check: **Causality** (stationarity)
- ✓ Automatically: **Invertible** (trivially, since already in AR form)

**Pure MA(q) Models:**
- ✓ Automatically: **Causal/Stationary** (always stationary)
- ✓ Need to check: **Invertibility**

**ARMA(p,q) Models:**
- ✓ Need to check: **Both causality AND invertibility**
- Check $\phi(z) = 0$ for causality
- Check $\theta(z) = 0$ for invertibility

#### Important Notes | 重要注意事项

1. **Independence:** Causality and invertibility are independent concepts
   - A model can be causal but not invertible
   - A model can be invertible but not causal
   - A model can be both or neither

2. **Stationarity vs. Causality:**
   - For AR models: causality = stationarity
   - "Causal" means the process can be expressed in terms of past shocks only

3. **Why does invertibility matter?**
   - Ensures unique identification of parameters
   - Important for forecasting and model fitting
   - Allows us to compute innovations from observations

---

## Part V: Non-stationary Models | 非平稳模型

### 13. ARIMA Models (ARIMA模型)

#### Definition | 定义
An **ARIMA(p,d,q)** model is defined as:
$$\phi(B)(1-B)^d X_t = \theta(B) w_t$$

Where:
- $p$ = order of autoregression (AR部分的阶数)
- $d$ = degree of differencing (差分次数)
- $q$ = order of moving average (MA部分的阶数)

#### Equivalent Forms | 等价形式

Let $Y_t = (1-B)^d X_t$ (the differenced series), then:
$$\phi(B) Y_t = \theta(B) w_t$$

This is an ARMA(p,q) model on the differenced series.

#### Key Properties | 关键性质
- Original series $X_t$ is **non-stationary** (has unit roots)
- Differenced series $Y_t = \nabla^d X_t$ is **stationary**
- $X_t$ is said to be **integrated of order $d$**, written as $I(d)$

---

### 14. Differencing (差分)

#### First Difference | 一阶差分
$$\nabla X_t = (1-B)X_t = X_t - X_{t-1}$$

#### Second Difference | 二阶差分
$$\nabla^2 X_t = (1-B)^2 X_t = (1-2B+B^2)X_t = X_t - 2X_{t-1} + X_{t-2}$$

#### General d-th Difference | d阶差分
$$\nabla^d X_t = (1-B)^d X_t$$

#### When to Difference | 何时差分
- When the series has a **trend** (non-constant mean)
- When the series has a **unit root**
- When variance increases over time (sometimes log transformation first)

#### How Many Times to Difference | 差分几次
- **Once (d=1):** For linear trend or random walk
- **Twice (d=2):** For quadratic trend (rare in practice)
- **Rule of thumb:** Rarely need $d > 2$

---

### 15. Unit Root and Integration (单位根与单整)

#### Unit Root | 单位根

A process has a **unit root** if the characteristic equation $\phi(z) = 0$ has a root on the unit circle (i.e., $|z| = 1$).

**Example: Random Walk**
$$X_t = X_{t-1} + w_t$$
$$\phi(B) = 1 - B$$
$$\phi(z) = 1 - z = 0 \Rightarrow z = 1$$

Root is exactly 1 (on unit circle) → Has unit root → Non-stationary

#### Integration Order | 单整阶数

A series $X_t$ is **integrated of order $d$**, denoted $X_t \sim I(d)$, if:
- $\nabla^d X_t$ is stationary
- $\nabla^{d-1} X_t$ is non-stationary

**Examples:**
- White noise: $w_t \sim I(0)$ (already stationary)
- Random walk: $X_t = X_{t-1} + w_t \sim I(1)$ (needs one difference)
- $X_t = X_{t-1} + w_t$ where $w_t = w_{t-1} + \varepsilon_t \sim I(2)$ (needs two differences)

#### Testing for Unit Roots | 单位根检验
- **Augmented Dickey-Fuller (ADF) test**
- **Phillips-Perron test**
- **KPSS test**

---

## Part VI: Linear Filters | 线性滤波器

### 16. Linear Filters and Stationarity (线性滤波器与平稳性)

#### Definition | 定义
A **linear filter** applies a weighted sum to a time series:
$$Y_t = \sum_{j=-\infty}^{\infty} \alpha_j X_{t-j}$$

where $\{\alpha_j\}$ is the sequence of filter weights.

#### Stationarity Preservation | 平稳性保持

**Theorem:** If $\{X_t\}$ is stationary and $\sum_{j=-\infty}^{\infty} |\alpha_j| < \infty$ (i.e., $\alpha \in \ell^1(\mathbb{Z})$), then $\{Y_t\}$ is also stationary.

**Proof sketch:**
1. **Mean:** $E[Y_t] = \sum \alpha_j E[X_{t-j}] = \mu_X \sum \alpha_j$ (constant)
2. **Variance:** Convergence ensured by $\sum |\alpha_j| < \infty$
3. **Autocovariance:** $\gamma_Y(h)$ depends only on $h$

#### Special Case | 特殊情况

If the formula is:
$$Y_t = \sum_{i \in \mathbb{Z}} \alpha_i X_i$$
(note: $X_i$ not $X_{t-i}$)

Then $Y_t$ does not depend on $t$ at all! It's a **constant random variable** (trivially stationary).

#### Common Linear Filters | 常见线性滤波器

1. **Moving Average:**
   $$Y_t = \frac{1}{2k+1} \sum_{j=-k}^{k} X_{t-j}$$
   (smoothing filter)

2. **Differencing:**
   $$Y_t = \nabla X_t = X_t - X_{t-1}$$
   ($\alpha_0 = 1, \alpha_1 = -1$, so $\sum |\alpha_j| = 2 < \infty$)

3. **MA(∞) Representation:**
   $$X_t = \sum_{j=0}^{\infty} \psi_j w_{t-j}$$
   (requires $\sum |\psi_j| < \infty$ or $\sum \psi_j^2 < \infty$)

---

## Part VII: Practice Problems | 练习题

### 17. Worked Examples (详解例题)

---

#### **Problem 1: Causality and Invertibility**

**Question:**
Consider the process:
$$X_t = 0.5X_{t-1} + 0.6X_{t-2} + w_t$$

Is this process:
- a) Neither invertible nor causal
- b) Invertible but not causal
- c) Invertible and causal
- d) Not causal but invertible

**Solution:**

**Step 1: Identify the model type**
This is an AR(2) model with no MA part.

**Step 2: Check causality (stationarity)**
Characteristic equation:
$$\phi(z) = 1 - 0.5z - 0.6z^2 = 0$$

Rearrange:
$$0.6z^2 + 0.5z - 1 = 0$$

Using quadratic formula:
$$z = \frac{-0.5 \pm \sqrt{0.25 + 2.4}}{1.2} = \frac{-0.5 \pm \sqrt{2.65}}{1.2} = \frac{-0.5 \pm 1.628}{1.2}$$

Roots:
- $z_1 = \frac{1.128}{1.2} \approx 0.94$
- $z_2 = \frac{-2.128}{1.2} \approx -1.77$

Check:
- $|z_1| = 0.94 < 1$ ✗ (inside unit circle)
- $|z_2| = 1.77 > 1$ ✓ (outside unit circle)

**Conclusion:** NOT causal (one root inside unit circle)

**Step 3: Check invertibility**
Since this is a **pure AR model** (no MA part), it is **automatically invertible** (trivially).

**Answer: b) Invertible but not causal** ✓

**Key lesson:** Pure AR models are always invertible; only need to check causality.

---

#### **Problem 2: ACF and PACF Identification**

**Question:**
Given ACF and PACF plots where:
- ACF: Significant at lags 1 and 5, insignificant elsewhere
- PACF: Significant at lag 1 (large negative value), insignificant at lag 2+

Which model is most appropriate?
- a) AR(1)
- b) MA(6)
- c) MA(2)
- d) ARMA(0,3)

**Solution:**

**Step 1: Analyze PACF**
- PACF has ONE significant lag (lag 1)
- After lag 1, all values within confidence interval
- **PACF cuts off after lag 1** → Suggests **AR(1)**

**Step 2: Analyze ACF**
- ACF significant at lags 1 and 5
- This is **tailing off** (not cutoff)
- Pattern shows oscillation (positive then negative)
- Consistent with AR(1) with negative coefficient

**Step 3: Why not other options?**
- MA(6): Would need PACF to tail off → Doesn't match ✗
- MA(2): Would need ACF to cut off at lag 2 → Doesn't match ✗
- ARMA(0,3) = MA(3): Would need ACF cutoff at lag 3 → Doesn't match ✗

**Answer: a) AR(1)** ✓

**Key lesson:** Prioritize the clearest cutoff pattern. PACF cutoff at lag 1 strongly indicates AR(1).

---

#### **Problem 3: Modeling a Real-World Problem**

**Question:**
The number of people receiving pensions at year $t$ depends on:
- Number of people who retire this year (random white noise)
- Number of people who died (equal to $a$ times the number of retired people)

Write a model. What do you need to check? What can you say about this model?

Note: People can live infinitely long (low mortality rate).

**Solution:**

**Step 1: Set up the model**
Let $P_t$ = number of people receiving pensions at year $t$

$$P_t = P_{t-1} + (\text{new retirees}) - (\text{deaths})$$
$$P_t = P_{t-1} + w_t - a \cdot w_{t-1}$$

where $w_t \sim \text{WN}(0, \sigma^2)$ and $0 < a < 1$ (mortality rate).

**Step 2: Rewrite in standard form**
$$P_t - P_{t-1} = w_t - a \cdot w_{t-1}$$
$$(1-B)P_t = (1-aB)w_t$$

This is an **ARIMA(0,1,1)** or **IMA(1,1)** model.

**Step 3: What to check?**

1. **Unit root:** $\phi(B) = 1-B$ has root $z=1$ (on unit circle)
   - Process is **non-stationary** ✓ (makes sense: pension numbers accumulate)

2. **Invertibility:** MA polynomial $\theta(B) = 1-aB$
   - Root: $z = 1/a$
   - Invertible if $|z| = 1/a > 1$ → **$|a| < 1$** ✓
   - Makes sense: mortality rate should be less than 100%

**Step 4: Interpretation**
- Non-stationary: Pension population grows over time
- MA(1) term: Effect of last year's retirees on this year's deaths
- If $a \approx 0$: Almost no one dies → process like random walk
- First difference $(1-B)P_t$ is stationary

**Answer:** ARIMA(0,1,1) model. Need to check invertibility condition $|a| < 1$. Model is non-stationary but first difference is stationary.

---

#### **Problem 4: Variance of Infinite Sum**

**Question:**
Let $w_t$ be white noise. Define:
$$X_t = \sum_{i \in \mathbb{Z}} \psi_i w_i$$

What property do we need to ensure:
$$\text{Var}(X_t) = \sum_{i \in \mathbb{Z}} \psi_i^2 \text{Var}(w_0) \text{ ?}$$

**Solution:**

**Key requirement:** $\sum_{i \in \mathbb{Z}} \psi_i^2 < \infty$ (square summable, $\ell^2$ condition)

**Why this matters:**

1. **Finite variance:** Ensures $\text{Var}(X_t) < \infty$

2. **Convergence:** For independent random variables:
   $$\text{Var}\left(\sum_{i} \psi_i w_i\right) = \sum_{i} \psi_i^2 \text{Var}(w_i) = \sigma^2 \sum_{i} \psi_i^2$$

3. **Without this condition:** The series may not converge, or variance may be infinite

**Relation to causality:** This is exactly the condition needed for a causal MA(∞) representation!

**Answer:** $\sum \psi_i^2 < \infty$ (square summability / $\ell^2$ condition)

---

#### **Problem 5: Weak Stationarity with Heavy-Tailed Distribution**

**Question:**
Each day in Edinburgh, it rains a height:
$$X_t = \min(0, Y_t)$$
where $Y_t$ are i.i.d. with Cauchy distribution:
$$f(y) = \frac{1}{\pi a(1 + (y/a)^2)}$$

Is the process weakly stationary?

**Solution:**

**Step 1: Check Condition 1 (Constant mean)**
$$E[X_t] = E[\min(0, Y_t)] = \int_{-\infty}^{0} y \cdot f(y) dy$$

For Cauchy distribution, even truncated to $(-\infty, 0]$:
$$\int_{-\infty}^{0} y \cdot \frac{1}{\pi a(1 + (y/a)^2)} dy \text{ diverges}$$

Mean does not exist or is $-\infty$ ✗

**Step 2: Check Condition 2 (Finite variance)**
$$E[X_t^2] = \int_{-\infty}^{0} y^2 \cdot \frac{1}{\pi a(1 + (y/a)^2)} dy$$

As $y \to -\infty$: integrand $\sim y^2/y^2 = 1$ → integral diverges

**$\text{Var}(X_t) = \infty$** ✗✗

**Step 3: Conclusion**
Fails Condition 2 (finite variance), so NOT weakly stationary.

**Answer: False** ✗

**Key lesson:** Heavy-tailed distributions (like Cauchy) often have infinite variance, violating weak stationarity even when the process "looks" stable.

---

#### **Problem 6: Linear Filter with $\ell^1$ Weights**

**Question:**
Let $X_t$ be a stationary process. Let $\alpha \in \ell^1(\mathbb{Z})$ (i.e., $\sum |\alpha_i| < \infty$).

Is it true that:
$$Y_t = \sum_{i \in \mathbb{Z}} \alpha_i X_i$$
is stationary?

**Solution:**

**Key observation:** The formula uses $X_i$ (not $X_{t-i}$), so $Y_t$ does **not depend on $t$**!

$$Y_t = \cdots + \alpha_{-1}X_{-1} + \alpha_0 X_0 + \alpha_1 X_1 + \alpha_2 X_2 + \cdots$$

This is the **same value for all $t$**.

**Analysis:**
- $Y_t = Y$ (a single random variable, independent of $t$)
- Mean: $E[Y_t] = E[Y]$ (constant) ✓
- Variance: $\text{Var}(Y_t) = \text{Var}(Y)$ (constant, finite if $\sum |\alpha_i| < \infty$) ✓
- Autocovariance: $\text{Cov}(Y_t, Y_{t+h}) = \text{Var}(Y)$ for all $h$ (trivially depends only on... well, doesn't depend on $h$, but still stationary in a degenerate sense) ✓

**Answer: True** ✓

This is a **trivially stationary** constant process.

**Note:** If the question were $Y_t = \sum \alpha_i X_{t-i}$, then we'd need to properly verify all three conditions using the linear filter theorem.

---

## Quick Reference Tables | 快速参考表

### Model Comparison | 模型对比

| Model | Form | Always Stationary? | Condition for Stationarity | Condition for Invertibility |
|-------|------|-------------------|---------------------------|----------------------------|
| MA(q) | $X_t = \theta(B)w_t$ | ✓ Yes | None | Roots of $\theta(z)=0$ outside unit circle |
| AR(p) | $\phi(B)X_t = w_t$ | ✗ No | Roots of $\phi(z)=0$ outside unit circle | Always (trivially) |
| ARMA(p,q) | $\phi(B)X_t = \theta(B)w_t$ | ✗ No | Roots of $\phi(z)=0$ outside unit circle | Roots of $\theta(z)=0$ outside unit circle |
| ARIMA(p,d,q) | $\phi(B)(1-B)^d X_t = \theta(B)w_t$ | ✗ No | Has $d$ unit roots | Roots of $\theta(z)=0$ outside unit circle |

### ACF/PACF Pattern Summary | ACF/PACF模式总结

| Model | ACF | PACF | How to Identify |
|-------|-----|------|----------------|
| AR(p) | Tails off (exponential decay or damped oscillation) | Cuts off after lag p | Look for PACF cutoff |
| MA(q) | Cuts off after lag q | Tails off | Look for ACF cutoff |
| ARMA(p,q) | Tails off | Tails off | Both tail off; try different (p,q) combinations |

### Checking Procedures | 检查程序

**For Causality (Stationarity):**
1. Write characteristic equation: $\phi(z) = 1 - \phi_1 z - \cdots - \phi_p z^p = 0$
2. Solve for roots $z_1, z_2, \ldots, z_p$
3. Check: All $|z_i| > 1$? → Causal ✓

**For Invertibility:**
1. Write MA characteristic equation: $\theta(z) = 1 + \theta_1 z + \cdots + \theta_q z^q = 0$
2. Solve for roots $z_1, z_2, \ldots, z_q$
3. Check: All $|z_i| > 1$? → Invertible ✓

**For Weak Stationarity:**
1. Check: $E[X_t] = \mu$ (constant)? ✓
2. Check: $\text{Var}(X_t) = \sigma^2 < \infty$ (finite constant)? ✓
3. Check: $\text{Cov}(X_t, X_{t+h}) = \gamma(h)$ (depends only on $h$)? ✓

---

## Common Mistakes to Avoid | 常见错误

1. **Confusing causality and invertibility**
   - Causality → checks AR part (stationarity)
   - Invertibility → checks MA part
   - They are independent!

2. **Forgetting that pure AR models are always invertible**
   - Don't waste time checking invertibility for AR models

3. **Forgetting that pure MA models are always stationary**
   - MA models don't need causality conditions

4. **Thinking negative values mean "not significant"**
   - Significance is about exceeding confidence bounds (regardless of sign)

5. **Not checking if variance is finite**
   - Heavy-tailed distributions (Cauchy, some Pareto) have infinite variance
   - Infinite variance → NOT weakly stationary

6. **Confusing "cutoff" and "tailing off"**
   - Cutoff: drops to insignificant and stays there
   - Tailing off: gradually decreases but multiple lags remain significant

7. **Applying the wrong representation**
   - $Y_t = \sum \alpha_i X_i$ vs. $Y_t = \sum \alpha_i X_{t-i}$ are very different!

8. **Forgetting unit roots**
   - ARIMA models have unit roots by definition
   - Need to difference before achieving stationarity

---

## Exam Tips | 考试技巧

1. **For identification problems:**
   - Look at PACF first for AR models
   - Look at ACF first for MA models
   - The clearest cutoff is usually the right answer

2. **For causality/invertibility problems:**
   - Identify model type first (AR, MA, or ARMA)
   - Pure AR → check causality only
   - Pure MA → check invertibility only
   - ARMA → check both

3. **For stationarity problems:**
   - Always check all three conditions
   - Don't forget to verify variance is finite
   - Watch out for heavy-tailed distributions

4. **For modeling problems:**
   - Write out the recursive relationship first
   - Convert to standard ARIMA form
   - Identify what needs to be checked

5. **Time management:**
   - ACF/PACF pattern recognition should be quick (30 seconds)
   - Root calculations may take 2-3 minutes
   - Don't spend too long on one problem

---

## Final Checklist Before Exam | 考前最终检查清单

- [ ] Can identify AR, MA, ARMA from ACF/PACF plots
- [ ] Know the three conditions for weak stationarity
- [ ] Can check causality by finding roots of $\phi(z)=0$
- [ ] Can check invertibility by finding roots of $\theta(z)=0$
- [ ] Understand the difference between causality and invertibility
- [ ] Know that pure AR models are always invertible
- [ ] Know that pure MA models are always stationary
- [ ] Can write ARIMA models in standard form
- [ ] Understand when differencing is needed (unit roots)
- [ ] Can recognize when variance is infinite (heavy-tailed distributions)

---

## Glossary | 术语表

| English | 中文 | Definition |
|---------|------|------------|
| Stationarity | 平稳性 | Property where statistical properties don't change over time |
| Causality | 因果性 | Process can be expressed as function of past shocks only |
| Invertibility | 可逆性 | White noise can be recovered from past observations |
| Autocovariance | 自协方差 | Covariance between $X_t$ and $X_{t+h}$ |
| Autocorrelation | 自相关 | Correlation between $X_t$ and $X_{t+h}$ |
| PACF | 偏自相关函数 | Correlation after removing linear effects of intermediate lags |
| White noise | 白噪声 | Uncorrelated random variables with constant variance |
| Unit root | 单位根 | Root of characteristic equation on unit circle |
| Differencing | 差分 | Taking differences to remove trend |
| Integration | 单整 | Number of times differencing is needed to achieve stationarity |
| Backshift operator | 后移算子 | Operator $B$ where $BX_t = X_{t-1}$ |
| Characteristic equation | 特征方程 | Polynomial equation from AR or MA part |
| Cutoff | 截尾 | Pattern where values become insignificant after certain lag |
| Tails off | 拖尾 | Pattern where values gradually decrease |

---

**Good luck with your exam! 考试顺利！** 🎓

