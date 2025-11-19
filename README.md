# Ranking & Unranking of Discrete Fuzzy Numbers (Discrete Fuzzy Numbers)

Reference Python implementation for **ranking** and **unranking** *Discrete Fuzzy Numbers* (Discrete Fuzzy Numbers) on finite chains. This code is associated to the paper LINK TO BE ADDED ONCE WE HAVE IT

The code includes strict Discrete Fuzzy Number validation, benchmarking, and plotting utilities.

---

## Notation & Conventions

- **Discrete chain**: $L_n=\{0,1,\dots,n\}$ (size $n+1$).
- **Set of membership values**: $Y_m=\{y_1=0<y_2<\cdots<y_m=1\}$.
- **Discrete Fuzzy Number**: will be a sequence $\mu(0),\dots,\mu(n)$ with $\mu(i)\in[0,1]$.
- **$\alpha$-cut** at level $y_j$ ($j\in\{1,\dots,m\}$):
  $$ A^{y_j}=\{\,i\in L_n : \mu(i)\ge y_j\,\}.$$
  For valid Discrete Fuzzy Numbers, each $A^{y_j}$ is a **contiguous interval** in $L_n$ and the family $A^{y_j}$ with $j=1,\ldots,m$ is **nested** (i.e., $A^{y_j}\supseteq A^{y_{j+1}}$).

**Standard display notation**:
$$\{\mu(0)/0,\ \mu(1)/1,\ \dots,\ \mu(n)/n\}.$$

---

## Example: From $\mu$ to $Y_m$-cuts (and intervals)

**Setup.** Let $n=5$ and $m=6$. The chain is $L_n=\{0,1,2,3,4,5\}$ (size $6$)
and the set of membership values is
$$
Y_6=\{y_1=0<0.2<0.4<0.6<0.8<y_6=1\}.
$$
Consider the Discrete Fuzzy Number
$$
A=\big[\,0/0,\ 0.2/1,\ 0.6/2,\ 1.0/3,\ 0.6/4,\ 0.2/5\,\big].
$$

**$Y_m$-cuts and intervals.** For each $j\in\{1,\dots,6\}$,

$$A^{y_j}=\{i\in L_n:\mu(i)\ge y_j\}.$$


| $j$ | $y_j$ | $A^{y_j}$                    | Interval $[a_j,b_j]$ |
|-----:|:-----:|------------------------------|----------------------|
| 6  |  1.0  | $\{3\}$                      | $[3,3]$              |
| 5  |  0.8  | $\{3\}$                      | $[3,3]$              |
| 4  |  0.6  | $\{2,3,4\}$                  | $[2,4]$              |
| 3  |  0.4  | $\{2,3,4\}$                  | $[2,4]$              |
| 2  |  0.2  | $\{1,2,3,4,5\}$              | $[1,5]$              |
| 1  |  0.0  | $\{0,1,2,3,4,5\}$            | $[0,5]$              |

We observe the **nesting** property:

$$A^{y_6}\subseteq A^{y_5}\subseteq A^{y_4}\subseteq A^{y_3}\subseteq A^{y_2}\subseteq A^{y_1}.$$

---

## Requirements

- **Python** 3.9+
- Python packages:
  ```bash
  pip install numpy pandas matplotlib
  ```
- The companion module (for the `module` backend), e.g. `dfn_cuts_rank_unrank.py`,
  must be:
  - in the same directory as `paper_algorithm.py`, or
  - in a directory listed in `PYTHONPATH`.

---


## Command-line interface (CLI)

Show the full help with:

```bash
python paper_algorithm.py --help
```

The main arguments are:

- `--unrank`
  Run **t-inc UNRANK**: from index \(i\) to DFN.

- `--rank`
  Run **t-inc RANK**: from DFN to index \(i\).

- `--n N`
  Chain size parameter: we work on \(L_n = \{0,\dots,n\}\).

- `--m M`
  Number of membership levels (`|Y_m| = m`).

- `--i I`
  Index for unrank, interpreted in the base given by `--index_base`.

- `--mu MU`
  Comma-separated membership vector for rank, e.g.
  `--mu "1,1,1,0.4,0.2,0.2"`.

- `--subindex SUBINDEX`
  Comma-separated integer **subindices** in `{1,...,m}`, of length `n+1`, e.g.
  `--subindex "6,6,6,2,1,1"` for `n=5`, `m=6`.
 

- `--index_base {0,1}`
  Index base used in the CLI:
  - `0` → indices in `{0,...,total-1}` (internal),
  - `1` → indices in `{1,...,total}` (matches the paper; **default**).

- `--verbose`
  If given, show a **detailed step-by-step log** (α-cuts, partial sums, etc.).
  Without `--verbose`, only the **final result** is printed.

- `--bench`
  Run timing benchmarks.

- `--interval_order {lex1|lex2|xy|t-inc|engine}`
  The chosen interval order for benchmarks:
  - `lex1` (lexicographic order),
  - `lex2`,
  - `xy`,
  - `t-inc` (default),
  - `engine` (see below).

- `--m_from`, `--m_to`, `--m_step`
  Range of `m` values for benchmarks.

- `--trials`
  Number of random rank/unrank queries per configuration in the benchmark.

- `--outdir OUTDIR`
  Output directory for logs, CSVs, plots, LaTeX tables (default: `out`).

- `--show_both_indices`
  When printing detailed logs, also display both 0-based and 1-based indices.

- `--plots_from_csv PATH`
  Re-generate plots from an existing CSV file (without re-running the benchmark).

- `--plots_order_label LABEL`
  Label for the "order" in plot titles (e.g. `t-inc/lex1`).

If the script is called without any of `--unrank`, `--rank`, `--bench`, or `--plots_from_csv`,
it will do nothing (you can choose to show the help in that case if desired).


#### Custom interval order with `--interval_order engine`
When the user selects `--interval_order engine`, 
the program expects a file named `engine.py` placed in the same directory as `paper_algorithm.py`.
This file must define the function:


```python
def _engine(I, J) -> bool:
    a, b = I
    c, d = J
    # return True iff [a,b] <= [c,d] in your custom interval order
    return <condition involving a, b, c, d>
```

##### Example: implementing the `lex1` order inside `engine.py`

```python
def _lex1(I, J):
    a, b = I
    c, d = J
    return (a < c) or ((a == c) and (b <= d))
```

## Script Modes

The script supports three general modes:

- **`ranking`**  
  Computes a discrete fuzzy number from a given index.

- **`unranking`**  
  Returns the index of a fuzzy number within its ordered chain.

- **`bench`**  
  Performs a time benchmark by exploring different values for the chain length and membership levels.

---

## Parameters Overview

In the following subsections, we provide a detailed description of all the parameters available in the script.
---

### Help

User can display the help panel (which contains all flags and parameters with their explanation) by running

```bash
python paper_algorithm.py --help
```
will return
```bash
usage: paper_algorithm.py [-h] [--unrank] [--rank] [--n N] [--m M] [--i I] [--mu MU] [--subindex SUBINDEX]
                          [--verbose] [--bench] [--engine {module,tinc}] [--comp {lex1,lex2,tinc}]
                          [--m_from M_FROM] [--m_to M_TO] [--m_step M_STEP] [--trials TRIALS] [--outdir OUTDIR]
                          [--index_base {0,1}] [--show_both_indices] [--plots_from_csv PLOTS_FROM_CSV]
                          [--plots_order_label PLOTS_ORDER_LABEL]

Ranking and unranking of discrete fuzzy numbers on finite chains, with module-based interval orders and benchmarking
utilities.

options:
  -h, --help            show this help message and exit
  --unrank              Run t-inc UNRANK (i -> DFN)
  --rank                Run t-inc RANK (DFN -> i)
  --n N
  --m M
  --i I                 Index for unrank, interpreted in --index_base
  --mu MU               Comma-separated _ for rank
  --subindex SUBINDEX   Comma-separated integer subindices (1...m) for rank
  --verbose             If written, this flag will show detailed step-by-step log
  --bench               Run timing benchmarks
  --interval_order     The chosen interval order: {lex1,lex2, xy, tinc, engine}
  --m_from M_FROM       initial m for the benchmark
  --m_to M_TO           final m for the benchmark
  --m_step M_STEP       m step
  --trials TRIALS       random queries per point
  --outdir OUTDIR       Output directory
  --index_base {0,1}    Index base for demos/logs (1 matches the paper)
  --show_both_indices   Also display the other index base
  --plots_from_csv PLOTS_FROM_CSV
                        Dir to the CSV to regenerate plots without benchmarking
  --plots_order_label PLOTS_ORDER_LABEL
                        'order' name that will appear in the plot title (for instance t-inc/lex1)
```



### Accepted input formats for discrete fuzzy numbers

From the CLI we accept two equivalent ways of specifying the discrete fuzzy number for **ranking**:

#### Real-valued membership vector: `--mu`

A string with real values in `[0,1]`, separated by commas, e.g.:

```bash
--mu "1,1,1,0.4,0.2,0.2"
```

This corresponds to:
$$
  \{1/0,\ 1/1,\ 1/2,\ 0.4/3,\ 0.2/4,\ 0.2/5\}.
$$



#### Discrete membership levels: `--subindex`

A string with **subindices** in `{1,...,m}`, e.g.:

```bash
--subindex "6,6,6,2,1,1"
```

for `n=5`, `m=6`, where subindex 6 represents membership 1, 2 represents 0.2, and 1 represents 0.


## Examples

We demonstrate how to run this code using the examples computed in the paper.
### UNRANK (paper example)

**Goal.** For \(n=5\), \(m=6\), and 1-based index \(i=50\), compute the corresponding discrete fuzzy number in the t-inc order.

Command:

```bash
python paper_algorithm.py --unrank --n 5 --m 6 --i 50 --index_base 1
```
output (without `--verbose`):

```text
DFN (mu-form): { 1.0/0, 1.0/1, 1.0/2, 0.4/3, 0.2/4, 0.2/5 }
```

If you add `--verbose`, the script will also print the detailed log:
$\alpha$-cuts, search intervals, partial sums, and final index in both bases (if `--show_both_indices` is set),
besides writing the log to a text file such as:

- `out/example_tinc_unrank.txt`.

---

### RANK (paper example)

**Goal.** Given the discrete fuzzy number

$$
  A = \{1/0,\ 1/1,\ 1/2,\ 0.2/3,\ 0/4,\ 0/5\},
$$

compute its position \(i\) in the t-inc order on the family of discrete fuzzy numbers for \(n=5\), \(m=6\).

#### 1. Using `--mu`

We specify \(A\) with the membership vector:

```bash
python paper_algorithm.py --rank     --n 5 --m 6     --mu "1,1,1,0.2,0,0"     --index_base 1
```

Typical output:

```text
Index i (1-based) = 55
```

So the position of \(A\) in the total t-inc order (with 1-based indexing, as in the paper) is \(i = 55\).

### 2. Using `--subindex`

For `m = 6`, the discrete fuzzy number

$$
  A = \{1/0,\ 1/1,\ 1/2,\ 0.2/3,\ 0/4,\ 0/5\}
$$

corresponds to the subindices:

- membership 1.0 → subindex 6,
- membership 0.2 → subindex 2,
- membership 0.0 → subindex 1,

so we can also write:

```bash
python paper_algorithm.py --rank     --n 5 --m 6     --subindex "6,6,6,2,1,1"     --index_base 1
```

which produces the **same** index:

```text
Index i (1-based) = 55
```

Again, if you add `--verbose`, a detailed step-by-step log will also be printed and saved to:

- `out/example_tinc_rank.txt`.

---

## Internal validation

Before ranking, the script validates that the input really is a discrete fuzzy number on \(L_n\) with levels in \(Y_m\).


The checks include:

1. **Range**
   - The discrete fuzzy number has exactly `n+1` values.
   - For `--mu`: each $\mu(i)$ lies in `[0,1]`.
   - For `--subindex`: there are `n+1` integers and each subindex lies in `{1,...,m}`.
 
2. **Core**
   - There is at least one point with maximum membership, $\max_i \mu(i) = 1,$ i.e. at least one subindex is equal to `m`.

3. **Membership function behavior**
   - Membership is increasing from the left up to the core,
   - and decreasing from the core to the right.


If any of these conditions fails, a `ValueError` is raised with a readable message, such as:

```text
Invalid DFN (subindex):
- subindex must have length n+1 = 6
- each subindex must be in 1..6
- no-core discrete fuzzy number: there must be at least one point with membership 1 (some subindex = m)
```

---

## Benchmarks and plots

To run timing benchmarks over a range of `m` values, you can use commands like:

```bash
python paper_algorithm.py --bench     --n 10     --m_from 10 --m_to 100 --m_step 10     --trials 100     --interval_order lex1     --outdir out
```

This produces (in `out/`):

- `bench_times_unrank_rank.csv` – average rank/unrank times and standard deviations,
- `fig_times_unrank_rank.png` – runtime vs `m`,
- `fig_loglog_unrank_rank.png` – log–log plot with fitted slopes,
- `table1_results.tex` – LaTeX table with the benchmark results,
- `bench_fit_fromcsv.txt` – slopes/intercepts of the fitted lines.

To regenerate plots from an existing CSV without re-running the benchmarks:

```bash
python paper_algorithm.py --plots_from_csv out/bench_times_unrank_rank.csv  --plots_order_label "t-inc"     --outdir out
```

---

## Troubleshooting

- **`ModuleNotFoundError` for the external module**

  Make sure `dfn_cuts_rank_unrank.py` (or your equivalent) is either:
  - in the same folder as `paper_algorithm.py`, or
  - on your `PYTHONPATH`, e.g.
    ```bash
    export PYTHONPATH="$PYTHONPATH:/path/to/your/module"
    ```

- **Index out of range in UNRANK**

  Check:
  - the value of `--index_base` (0 or 1),
  - and that your `--i` lies in the correct range `{1,...,total}` (for base 1) or
    `{0,...,total-1}` (for base 0), where `total` is the total number of DFNs on `(n,m)`.

- **Invalid DFN**

  Ensure:
  - length = `n+1`,
  - membership values in `[0,1]` (for `--mu`),
  - subindices in `1..m` (for `--subindex`),
  - there is at least one point with membership 1,
  - The membership function is **non-increasing** to the left of the core and **non-decreasing** to the right.
  - $\alpha$-cuts are contiguous and nested.

Use `--verbose` to see a full trace of the α-cuts and partial sums if needed.

---

## Citation

If you use this code for research or teaching, please cite this repository and the associated manuscript. DOI to be included.

---

## License

MIT License
