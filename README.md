# 🧩 Ranking and Unranking of Discrete Fuzzy Numbers via α-cuts and Interval Orders

---

## 📘 Overview

This repository contains the **reference Python implementation** of the algorithms proposed in the paper:

> **"Ranking and Unranking of Discrete Fuzzy Numbers via α-cuts and Interval Orders"**

The goal of this work is to establish a **bijective mapping** between *Discrete Fuzzy Numbers (DFNs)* and a contiguous integer index set  
\[
\{0, 1, \dots, N-1\},
\]
where \(N\) is the number of DFNs defined over a discrete support \(L_n \to Y_m\).

The mapping is achieved through the **ranking** (`pos`) and **unranking** (`pos⁻¹`) functions based on:
- α-cut representations of DFNs,
- a general **interval order** (e.g., *t-inc*, *lexicographic*, or any admissible total order).

This bijection allows DFNs to be:
- **Enumerated**, **compared**, and **stored** efficiently,
- **Sampled** uniformly from the discrete fuzzy space,
- **Serialized** and **reconstructed** exactly without information loss.

---

## 🧮 Key Concepts

### 1. Discrete Fuzzy Numbers (DFN)
A DFN is a fuzzy number defined on a discrete finite lattice \(L_n = \{0, \frac{1}{n-1}, \dots, 1\}\).  
Each α-cut of a DFN is represented by a **closed interval** of indices from \(Y_m = \{0, \frac{1}{m-1}, \dots, 1\}\).

### 2. α-cut Representation
Each DFN can be represented as a finite sequence of α-cuts:
\[
x = [ [l_0, r_0], [l_1, r_1], \dots, [l_{n-1}, r_{n-1}] ],
\]
where \(l_i, r_i \in Y_m\) and \(l_i \leq r_i\).

### 3. Admissible Interval Orders
The order relation between α-cuts can be defined in several ways:
- **t-inc** (Total increasing order)
- **lexicographic order**
- **reverse order**
- or any custom admissible total order

The repository allows you to plug in any order that satisfies monotonicity and completeness.

---

## ⚙️ Features

✅ Exact **ranking (pos)** and **unranking (pos⁻¹)** of DFNs  
✅ Works with **any admissible interval order** (default: *t-inc*)  
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
   Each Discrete Fuzzy Number (DFN) is represented as a finite sequence of α-cuts:
   \[
   x = [ [l_0, r_0], [l_1, r_1], \dots, [l_{n-1}, r_{n-1}] ],
   \]
   where \(l_i, r_i \in Y_m\) and \(l_i \leq r_i\).

2. **Interval orders**  
   Each α-cut belongs to a finite set of closed intervals over \(Y_m\).  
   The algorithm allows any *admissible total order* over these intervals (e.g., **t-inc**, **lexicographic**, etc.).  
   These orders must satisfy *monotonicity* and *compatibility* with the fuzzy lattice structure.

3. **Bijection construction**  
   The ranking (`pos`) function defines a unique integer index for every DFN, while  
   the unranking (`pos⁻¹`) function reconstructs the DFN from that index:
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

### General form
``
python paper_algorithm_full_en_v7.py [ACTION FLAGS] [OPTIONS]
``

### Common options
--engine {tinc,lex,rev} # Interval order / engine (default: tinc)

--n <int> # |L_n| (number of α-cuts / levels)

--m <int> # |Y_m| (discrete grid for endpoints)

--index-base {0,1} # Indexing base for pos/unpos (default: 0)

--seed <int> # RNG seed for reproducibility

--quiet # Less verbose output

### Examples

1) Rank (pos): DFN → integer index
# DFN given as JSON-like α-cuts [[l0,r0],...,[l_{n-1},r_{n-1}]]
python paper_algorithm_full_en_v7.py \
  --rank \
  --n 3 --m 5 \
  --dfn "[[0,4],[1,3],[2,2]]" \
  --engine tinc \
  --index-base 0

2) Unrank (pos⁻¹): integer index → DFN  
python paper_algorithm_full_en_v7.py \
  --unrank \
  --n 3 --m 5 \
  --index 17 \
  --engine tinc \
  --index-base 0

3) Validate bijection on a range
# Full space (careful if N is large)
python paper_algorithm_full_en_v7.py --validate --n 4 --m 6 --engine tinc

# Sampled validation with fixed seed
python paper_algorithm_full_en_v7.py --validate --n 10 --m 100 --trials 1000 --seed 123

4) Benchmark (micro/milli seconds)
python paper_algorithm_full_en_v7.py \
  --bench --engine tinc \
  --n 10 --m_from 100 --m_to 1000 --m_step 100 \
  --trials 500

# Fixed (n,m) with trials, export CSV
python paper_algorithm_full_en_v7.py \
  --bench --engine lex \
  --n 6 --m 200 --trials 2000 \
  --export bench_n6_m200_lex.csv


## 📚 References

If you use this repository, please cite the following work:

> **[Author(s)], "Ranking and Unranking of Discrete Fuzzy Numbers via α-cuts and Interval Orders," 2024.**  
> *Preprint or conference version, 2024.*  
> DOI: *(to be added upon publication)*  

**Additional related works:**
- [1] R. Fullér and P. Majlender, *On obtaining ranking from fuzzy numbers by maximizing similarity measures*, Fuzzy Sets and Systems, 2001.  
- [2] S. Bustince et al., *Interval-valued fuzzy sets and their applications*, Information Sciences, 2013.  
- [3] K. Atanassov, *Intuitionistic Fuzzy Sets: Theory and Applications*, 1999.  
- [4] [Your Paper/Presentation Title], *FSTA 2024 proceedings*, 2024.

---

## 🪪 License

This project is released under the **MIT License**:


