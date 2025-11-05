# 🧩 Ranking and Unranking of Discrete Fuzzy Numbers via α-cuts and Interval Orders

---

## 📘 Overview

This repository contains the **reference Python implementation** of the algorithms proposed in the paper:

> **"An Efficient Computational Framework for Discrete Fuzzy Numbers Based on Admissible Orders"**

The goal of this work is to establish a **bijective mapping** between *Discrete Fuzzy Numbers (dfns)* and a contiguous integer index set  
\[
\{0, 1, \dots, N-1\},
\]
where \(N\) is the number of dfns defined over a discrete support \(L_n \to Y_m\).

The mapping is achieved through the **ranking** (`pos`) and **unranking** (`pos⁻¹`) functions based on:
- α-cut representations of dfns,
- a general **interval order** (e.g., *t-inc*, *lexicographic*, or any admissible total order).

This bijection allows dfns to be:
- **Enumerated**, **compared**, and **stored** efficiently,
- **Sampled** uniformly from the discrete fuzzy space,
- **Serialized** and **reconstructed** exactly without information loss.

---

## 🧮 Key Concepts

### 1. Discrete Fuzzy Numbers (dfn)
A dfn is a fuzzy number defined on a discrete finite lattice \(L_n = \{0, \frac{1}{n-1}, \dots, 1\}\).  
Each α-cut of a dfn is represented by a **closed interval** of indices from \(Y_m = \{0, \frac{1}{m-1}, \dots, 1\}\).

### 2. α-cut Representation
Each dfn can be represented as a finite sequence of α-cuts:
\[
x = [ [l_0, r_0], [l_1, r_1], \dots, [l_{n-1}, r_{n-1}] ],
\]
where \(l_i, r_i \in Y_m\) and \(l_i \leq r_i\).

### 3. Canonical Membership Representation (μ-form)

dfns can also be written in the **standard fuzzy-set notation** used in the literature:

\[
x = \{\ dfn(0)/0,\, dfn(1)/1, \dots, dfn(n)/n \}
\]

or equivalently:
\[
x = \{ (y_0 / \alpha_0), (y_1 / \alpha_1), \dots, (y_n / \alpha_n) \}
\]

where each pair \((y_i, \alpha_i)\) expresses the **membership degree** \(\alpha_i \in L_n\) associated to the discrete support value \(y_i \in Y_m\).


Internally, the library automatically converts this μ-representation into its corresponding **α-cut form** to perform ranking and unranking operations.  
Hence, both forms are fully supported:

| Representation | Example | Description |
|----------------|----------|--------------|
| α-cut | `[[0,4],[1,3],[2,2]]` | Internal interval-based form |
| μ-form | `{0/0, 1/1, 2/1}` or JSON: `{"0":0, "1":1, "2":1}` | Standard membership form |

You can input dfns in either form — the system will normalize them internally.


### 4. Admissible Interval Orders
The order relation between α-cuts can be defined in several ways:
- **t-inc** (Total increasing order)
- **lexicographic order**
- **reverse order**
- or any custom admissible total order

The repository allows you to plug in any order that satisfies monotonicity and completeness.

---

## ⚙️ Features

✅ Exact **ranking (pos)** and **unranking (pos⁻¹)** of dfns  
✅ Works with **any admissible interval order** (default: *t-inc*, implemented: *t-inc*, *lex1*, *lex2*)  
✅ Efficient computation via **α-cut decomposition**  
✅ Support for **0-based or 1-based indexing**  
✅ **Bidirectional** consistency guaranteed:  
\[
\text{pos}^{-1}(\text{pos}(x)) = x
\]

---

## 🐍 Installation

This code requires **Python ≥ 3.8** and **NumPy**.

git clone https://github.com/AlejandroMus/dfn-ranking-unranking.git
cd dfn-ranking-unranking
pip install numpy

---

## 🧠 Theoretical Background

The algorithms are based on the following principles:

1. **α-cut decomposition**  
   Each Discrete Fuzzy Number (dfn) is represented as a finite sequence of α-cuts:
   \[
   x = [ [l_0, r_0], [l_1, r_1], \dots, [l_{n-1}, r_{n-1}] ],
   \]
   where \(l_i, r_i \in Y_m\) and \(l_i \leq r_i\).

2. **Interval orders**  
   Each α-cut belongs to a finite set of closed intervals over \(Y_m\).  
   The algorithm allows any *admissible total order* over these intervals (e.g., **t-inc**, **lexicographic**, etc.).  
   These orders must satisfy *monotonicity* and *compatibility* with the fuzzy lattice structure.

3. **Bijection construction**  
   The ranking (`pos`) function defines a unique integer index for every dfn, while  
   the unranking (`pos⁻¹`) function reconstructs the dfn from that index:
   \[
   \text{pos}: DFN(n,m) \rightarrow \{0, \ldots, N-1\}, \quad
   \text{pos}^{-1}: \{0, \ldots, N-1\} \rightarrow DFN(n,m)
   \]
   ensuring that both are perfect inverses:
   \[
   \text{pos}^{-1}(\text{pos}(x)) = x.
   \]

4. **Complexity**  
   The algorithms have theoretical complexity  
   \[
   O(n^2 \log n \, m),
   \]
   with empirical near-linear scaling in \(m\) for fixed \(n\).  
   This makes the method feasible for medium-scale fuzzy lattices.

---

## 🚀 Usage (CLI)

> All commands run from the repo root.  
> Main entry point: `paper_algorithm_full_en_v7.py`  


```bash
python paper_algorithm_full_en_v7.py [ACTION FLAGS] [OPTIONS]
```

---

## ⚙️ Common Options

| Flag | Description |
|------|--------------|
| `--engine {tinc,lex,rev}` | Interval order / engine (default: `tinc`) |
| `--n <int>` | cardinal(L_n) — number of alpha-cuts / levels |
| `--m <int>` | cardinal(Y_m) — discrete grid for endpoints |
| `--index-base {0,1}` | Indexing base for `pos` / `unpos` (default: 0) |
| `--seed <int>` | RNG seed for reproducibility |
| `--quiet` | Less verbose output |

---


### 1️⃣ Rank (pos): dfn → integer index

dfn given as JSON-like α-cuts `[[l0,r0],...,[l_{n-1},r_{n-1}]]`

```bash
python paper_algorithm_full_en_v7.py   --rank   --n 3 --m 5   --dfn "[[0,4],[1,3],[2,2]]"   --engine tinc   --index-base 0
```

**Example Output:**
```
engine=tinc, n=3, m=5, base=0
dfn = [[0, 4], [1, 3], [2, 2]]
pos(dfn) = 17
```

---

### 2️⃣ Unrank (pos⁻¹): integer index → dfn

```bash
python paper_algorithm_full_en_v7.py   --unrank   --n 3 --m 5   --index 17   --engine tinc   --index-base 0
```

**Example Output:**
```
engine=tinc, n=3, m=5, base=0
pos^{-1}(17) = [[0, 4], [1, 3], [2, 2]]
```

---

### 3️⃣ Validate Bijection on a Range

```bash
# Full space (careful if N is large)
python paper_algorithm_full_en_v7.py --validate --n 4 --m 6 --engine tinc

# Sampled validation with fixed seed
python paper_algorithm_full_en_v7.py --validate --n 10 --m 100 --trials 1000 --seed 123
```

Checks that:
```
pos^{-1}(pos(x)) == x
pos(pos^{-1}(i)) == i
```

---

### 4️⃣ Benchmark (micro/milli seconds)

#### Sweep over m
```bash
python paper_algorithm_full_en_v7.py   --bench --engine tinc   --n 10 --m_from 100 --m_to 1000 --m_step 100   --trials 500
```

#### Fixed (n, m) with trials, export CSV
```bash
python paper_algorithm_full_en_v7.py   --bench --engine lex   --n 6 --m 200 --trials 2000   --export bench_n6_m200_lex.csv
```

**Example Output:**
```
engine=tinc, n=10, trials=500
m   N(dfn)   mean_us   p50_us   p95_us
100  ...     42.8      41.9     55.2
200  ...     84.1      82.7     108.3
...
```

If `--export` is set, results are saved as a CSV file.

---

### 5️⃣ Input / Output Helpers

#### Rank a list of dfns from a JSON file
```bash
python paper_algorithm_full_en_v7.py   --rank --n 3 --m 5   --from-file dfn_list.json   --engine tinc
```

#### Unrank a list of indices from a text file
```bash
python paper_algorithm_full_en_v7.py   --unrank --n 3 --m 5   --from-file indices.txt   --engine tinc   --index-base 0   --export unranked.json
```

> `--from-file` automatically detects JSON arrays or newline-separated indices.  
> `--export` writes results as JSON/CSV depending on the action.

---

### 6️⃣ Index Base Conversion (0 ↔ 1)

```bash
python paper_algorithm_full_en_v7.py   --rank --n 3 --m 5   --dfn "[[0,4],[1,3],[2,2]]"   --index-base 1
```

---

### 7️⃣ Reproduce Figures or Tables (optional)

```bash
python paper_algorithm_full_en_v7.py   --repro figure-2   --n 6 --m 60   --engine tinc   --export fig2_data.csv
```

---

### 🆘 Help

```bash
python paper_algorithm_full_en_v7.py --help
```


## 📚 References

If you use this repository, please cite the following work:

> **[Mir, Mus, Riera (In prep.)], "An Efficient Computational Framework for Discrete Fuzzy Numbers Based on Admissible Orders**  
> DOI: *(to be added upon publication)*
> bibtex *(to be added upon publication)*


---

## 🪪 License

This project is released under the **MIT License**:


