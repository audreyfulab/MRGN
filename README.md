
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRGN: the Mendelian Randomization Genomic Networks method for causal network inference for genomic trios

<!-- ![](./static/MRGN-logos.jpeg width="800" height="600" style="display: block; margin: 0 auto") -->

<p align="center">
    <img src="static/MRGN-logos.jpeg" width="50%" height="50%" />
</p>
<!--  -->
<!-- badges: start --> <!-- badges: end -->

## Overview

MRGN is a statistical method that leverages the principle of Mendelian randomization 
to infer causal relationship for a genomic trio, which consists of 
a genetic variant (as the instrumental variable) and two molecular phenotypes
(for example, expression of two genes), with or without the presence of 
confounders. Five mutually exclusive causal models may be inferred. 

<p align="center">
    <img src="static/5 causal models.drawio.png" width="60%" height="60%" />
</p>

One of the key challenges in causal inference is accounting for confounding variables $\bf U$, which can significantly impact the results. MRGN addresses this challenge by integrating a regression-based framework for conditional dependence testing that can handle a large number of confounding variables effectively.

$$ T_1 = \alpha_1 +\beta_{11}V+\beta_{12}T_2+{{\bf \Gamma}_1 {\bf U}}+ \epsilon_1; \ \ (1)$$

$$ T_2 = \alpha_2 +\beta_{21}V+\beta_{22}T_1+{{\bf \Gamma}_2} {\bf U}+\epsilon_2  \ \ (2)$$

**Table: Conditional and marginal dependence tests utilized by MRGN**  
*A reference table of the conditional and marginal dependence tests used to determine which model is supported by the data. The first three models each have two possible configurations. Coefficients β₁₁ and β₁₂ are from regression (Eq. 1), and coefficients β₂₁ and β₂₂ from regression (Eq. 2). These coefficients test conditional independence. Correlations $r_{V, T_1}$ and $r_{V, T_2}$ test marginal independence.*

<table>
<thead>
<tr>
<th colspan="2"></th>
<th colspan="4" align="center">Conditional Test</th>
<th colspan="2" align="center">Marginal Test</th>
</tr>
<tr>
<th colspan="2" align="right">Null Hypothesis</th>
<th align="center">T₁ ⫫ V</th>
<th align="center">T₁ ⫫ T₂</th>
<th align="center">T₂ ⫫ V</th>
<th align="center">T₂ ⫫ T₁</th>
<th align="center">V ⫫ T₁</th>
<th align="center">V ⫫ T₂</th>
</tr>
<tr>
<th colspan="2" align="right">Conditioning Set</th>
<td align="center">T₂, <b>U</b></td>
<td align="center">V, <b>U</b></td>
<td align="center">T₁, <b>U</b></td>
<td align="center">V, <b>U</b></td>
<td align="center">--</td>
<td align="center">--</td>
</tr>
<tr>
<th colspan="2" align="right">Parameter</th>
<td align="center">β₁₁</td>
<td align="center">β₁₂</td>
<td align="center">β₂₁</td>
<td align="center">β₂₂</td>
<td align="center">r<sub>V,T₁</sub></td>
<td align="center">r<sub>V,T₂</sub></td>
</tr>
<tr>
<th align="left">Model</th>
<th align="left">Configuration</th>
<th colspan="4"></th>
<th colspan="2"></th>
</tr>
</thead>
<tbody>
<tr>
<td rowspan="2" align="center">M0</td>
<td>V → T₁; T₂ is singleton</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
<td align="center">= 0</td>
<td align="center">= 0</td>
<td align="center">--</td>
<td align="center">--</td>
</tr>
<tr>
<td>V → T₂; T₁ is singleton</td>
<td align="center">= 0</td>
<td align="center">= 0</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
<td align="center">--</td>
<td align="center">--</td>
</tr>
<tr>
<td rowspan="2" align="center">M1</td>
<td>V → T₁ → T₂</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
<td align="center">≠ 0</td>
<td align="center">--</td>
<td align="center">--</td>
</tr>
<tr>
<td>V → T₂ → T₁</td>
<td align="center">= 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">--</td>
<td align="center">--</td>
</tr>
<tr>
<td rowspan="2" align="center">M2</td>
<td>V → T₁ ← T₂</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
</tr>
<tr>
<td>V → T₂ ← T₁</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
<td align="center">≠ 0</td>
</tr>
<tr>
<td align="center">M3</td>
<td>T₁ ← V → T₂</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
<td align="center">≠ 0</td>
<td align="center">= 0</td>
<td align="center">--</td>
<td align="center">--</td>
</tr>
<tr>
<td align="center">M4</td>
<td>T₁ ← V → T₂; T₁ ↔ T₂</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
<td align="center">≠ 0</td>
</tr>
</tbody>
</table>




## Installation

You can install MRGN from GitHub release with

``` r
install.packages("devtools")
devtools::install_github("audreyfulab/MRGN")
```

## Example

``` r
library(MRGN)
#inference on a single trio with an eQTL as the genetic variant
result=infer.trio(M1trio)
print(result)
 
#fast example on 10 eQTL trios from the built in dataset WBtrios
results = sapply(WBtrios[1:10], function(x) infer.trio(x))
print(results)
#return just the inferred model topology
models = sapply(WBtrios[1:10], function(x) infer.trio(x)$Inferred.Model)
print(models)
#fast example on 10 trios with copy number alterations (CNAs) 
#as genetic variants from the built in dataset CNAtrios using permutation
models = sapply(CNAtrios[1:10], function(x) infer.trio(x, is.CNA = TRUE, use.perm = TRUE, nperms = 1000)$Inferred.Model)
print(models)

```

