# SwarmEvo

群体智能优化算法库，基于 C11 实现，包含 10 种经典元启发式算法。

## 算法列表

| 缩写 | 算法名称 |
|------|---------|
| GA   | 遗传算法 Genetic Algorithm |
| PSO  | 粒子群算法 Particle Swarm Optimization |
| ACO  | 蚁群算法 Ant Colony Optimization |
| SA   | 模拟退火 Simulated Annealing |
| CSO  | 鸡群算法 Chicken Swarm Optimization |
| CS   | 布谷鸟算法 Cuckoo Search |
| ABC  | 人工蜂群算法 Artificial Bee Colony |
| AIA  | 免疫算法 Artificial Immune Algorithm |
| EDA  | 分布估计算法 Estimation of Distribution Algorithm |
| DE   | 差分进化算法 Differential Evolution |

## 依赖

- GCC（支持 C11）
- [xmake](https://xmake.io/#/getting_started) 构建工具
- algmath 子模块（已包含在 `algmath/` 目录）

## 构建

```bash
# 克隆（含子模块）
git clone --recurse-submodules <repo-url>

# 构建共享库
xmake b alg-project

# 构建并运行测试
xmake b test-alg-project
xmake r test-alg-project
```

## 快速上手

所有算法遵循统一的三步接口：`*_init` → `*_fresh` → `*_free`。

```c
#include "optimization.h"

// 1. 定义目标函数
double my_func(alg_vector *x) {
    double sum = 0;
    for (int i = 0; i < x->size; i++)
        sum += x->vector[i] * x->vector[i];
    return sum;
}

// 2. 初始化优化句柄
optim_handle optim;
double l[6] = {-100, -100, -100, -100, -100, -100};
double r[6] = { 100,  100,  100,  100,  100,  100};
optim_init(&optim, 6, my_func, l, r);

// 3. 运行算法（以 PSO 为例）
pso_handle *pso = pso_init(optim, 50, 0.5, 2.0, 2.0);
pso_fresh(pso, 300);   // 迭代 300 代
optim_print(&optim);   // 打印最优解
pso_free(pso);
optim_free(&optim);
```

### 各算法初始化参数

```c
// GA：种群大小、变异率、交叉率
ga_handle *ga = ga_init(optim, pop_size, mutation_rate, crossover_rate);

// PSO：种群大小、惯性权重 w、学习因子 c1/c2
pso_handle *pso = pso_init(optim, pop_size, w, c1, c2);

// ACO（求解 TSP）：蚂蚁数、城市坐标矩阵、α、β、蒸发率 ρ
aco_handle *aco = aco_init(number, city_coords, alpha, beta, rho);

// SA：初始温度、冷却率
sa_handle *sa = sa_init(optim, temperature, cooling_rate);

// CSO：种群大小、公鸡数 rn、母鸡数 hn、小鸡数 cn、跟随系数 fl
cso_handle *cso = cso_init(optim, pop_size, rn, hn, cn, fl);

// DE：种群大小、缩放因子 F、交叉率 CR
de_handle *de = de_init(optim, pop_size, F, CR);

// ABC：种群大小、最大停滞次数
abc_handle *abc = abc_init(optim, pop_size, max_count);
```

## 项目结构

```
SwarmEvo/
├── algmath/        # 数学基础库（子模块）
├── src/            # 算法实现
│   ├── basic_opti.h/c   # 公共优化基础结构
│   ├── optimization.h   # 统一入口头文件
│   ├── ga/  pso/  aco/  sa/  cso/  cs/  abc/  aia/  eda/  de/
├── test/           # 测试代码
├── docs/           # 算法原理文档
├── matlab/         # MATLAB 参考实现
└── xmake.lua       # 构建配置
```
