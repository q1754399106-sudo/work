# Numerical implementation audit

This audit concerns the four `GRAM_LOW_RANK_MG_CACHE_V2` changes. It is not evidence that a candidate outperforms Frozen. Original source remains archived in `candidate/legacy_numerics_v1/`.

Let `L=blkdiag(L_j)` be the Cholesky factors of the block covariance, `V=L\U`, and `K=I+V'*V=T*T'`. With `B=blkdiag(L_j\G_j)` and `Z=B'*V`, the original full-gather normal product is exactly

```
G'*Ceff^(-1)*G*x = B'*B*x - Z*(T'\(T\(Z'*x))).
```

The optimized implementation uses the same block Grams previously used by the preconditioner and exactly this low-rank correction. It retains the original precision-based right-hand side and audits the original normal product after every scale solve. The data residual used to select regularization still uses the original forward/precision functions. A failed solve restarts with the original normal operator. Roundoff can change PCG trajectories, so numerical equivalence must be measured rather than inferred from algebra alone.

The fixed observed pilot, `b=mean((Gdraw-Gref)*pilot)`, the mean correction `D-b`, covariance centering/shrinkage, noise estimates, and all discrepancy targets are unchanged. Caching `MG=M*Gref` uses the full nonstationary reference operator. It does not replace `Gtotal` by an assumed depth/offset factorization. No truth-dependent information was added to inversion.

No scientific estimator bug was found in these four changes. The additional `white_G`/`white_y` context fields consume memory after `Z` has been formed but are not correctness problems. They can be omitted in a later implementation cleanup if memory becomes limiting.

`p1_compare_numerics(outdir)` supplies three checks through the existing MATLAB session:

1. Dense precision/normal-equation checks with low-rank sizes zero and four.
2. Same-observation computational-thread benchmark at the original limit, 1, 2 and 4 threads where allowed. Every timing is written immediately. Settings, RNG and path are restored with `onCleanup` even on failure.
3. Six legacy/optimized candidate/control comparisons sharing one calibration cache and a fixed development observation: DevGradedCompact, Main, SmoothAVO, calibration 310101, observation 410101. Point/bootstrap estimates must match. Selected scales must agree within relative 1e-12; waveform and historical metric differences must be below 1e-5; all original normal residuals must be at most 2e-7.

The computational-thread setting only schedules the same operations. The fastest passing setting is selected using development timing and implementation agreement; no heldout result enters that choice. A single timing is descriptive, so deployment runtime should be estimated conservatively.

Official thread-control semantics: [MathWorks maxNumCompThreads](https://www.mathworks.com/help/matlab/ref/maxnumcompthreads.html). The setter returns the previous thread limit and supports restoring it.

## Actual MATLAB result

The parent runner executed this helper in the existing ResearchMATLAB, MATLAB R2024b Update 6. `preflight/numerical_equivalence_v2/matlab/numerical_equivalence.json` reports PASS. The independent normal-product discrepancies were 3.07e-16 (rank zero) and 5.90e-16 (rank four). The six actual development-case comparisons retained identical selected scales; maximum relative waveform difference was 1.13e-13, maximum historical metric difference was 4.78e-15, and maximum original-precision normal residual was 9.80e-9. Path, RNG and the original computational-thread count were preserved.

| Threads | CombinedUncertainty time (s) | Numerical gate |
|---|---:|---|
| Original 8 | 7.199 | PASS |
| 1 | 8.441 | PASS |
| 2 | 7.928 | PASS |
| 4 | 7.382 | PASS |

Keep the original eight-thread setting. This test found no thread-count acceleration.

| Branch | Legacy (s) | Optimized (s) |
|---|---:|---:|
| NoiseUncertaintyOnly | 12.263 | 7.690 |
| QModelErrorOnly | 13.568 | 7.748 |
| CombinedUncertainty | 13.090 | 7.857 |
| ScalarBudgetDamping | 9.546 | 5.256 |
| NoiseCovarianceOnly | 9.540 | 5.355 |
| NoStructureCombined | 2.157 | 1.707 |

These timings describe one already designated development observation with a shared calibration cache. They do not assess scientific benefit or provide a heldout runtime distribution.
