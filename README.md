
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chemr <img src="man/figures/logo.png" align="right" width="120" />

[![Travis build
status](https://travis-ci.org/paleolimbot/chemr.svg?branch=master)](https://travis-ci.org/paleolimbot/chemr)
[![Coverage
status](https://codecov.io/gh/paleolimbot/chemr/branch/master/graph/badge.svg)](https://codecov.io/github/paleolimbot/chemr?branch=master)

Chemr provides data structures for elements, molecules, and reactions,
to provide a framework for chemical modelling and analysis in R. The
text-based formats for elements, reactions, and molecules broadly
correspond to that used by the
[PHREEQC](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/)
tool, an interface for which is [available in
R](https://cran.r-project.org/package=phreeqc). An extension to the R
library, [tidyphreeqc](https://github.com/paleolimbot/tidyphreeqc) is
under active development.

## Installation

You can install chemr from github with:

``` r
# install.packages("devtools")
devtools::install_github("paleolimbot/chemr")
```

If you can load the package, everything worked!

``` r
library(chemr)
```

## The periodic tibble

The periodic tibble is [Wikipedia’s version of the periodic
table](https://en.wikipedia.org/wiki/List_of_chemical_elements) in
data.frame (tibble) form. You can access these data in a few ways:

``` r
is_element("Hg")
#> [1] TRUE
```

``` r
elmass("Hg")
#>       Hg 
#> 200.5923
```

``` r
elz("Hg")
#> Hg 
#> 80
```

``` r
elsymbol(80)
#> [1] "Hg"
```

You can also access the entire periodic tibble by typing `data(pt)`.

``` r
data(pt)
pt
#> # A tibble: 118 × 7
#>        z symbol name      group period  mass valence  
#>    <int> <chr>  <chr>     <int>  <int> <dbl> <list>   
#>  1     1 H      Hydrogen      1      1  1.01 <int [3]>
#>  2     2 He     Helium       18      1  4.00 <int [1]>
#>  3     3 Li     Lithium       1      2  6.94 <int [2]>
#>  4     4 Be     Beryllium     2      2  9.01 <int [3]>
#>  5     5 B      Boron        13      2 10.8  <int [6]>
#>  6     6 C      Carbon       14      2 12.0  <int [9]>
#>  7     7 N      Nitrogen     15      2 14.0  <int [9]>
#>  8     8 O      Oxygen       16      2 16.0  <int [5]>
#>  9     9 F      Fluorine     17      2 19.0  <int [2]>
#> 10    10 Ne     Neon         18      2 20.2  <int [1]>
#> # … with 108 more rows
```

## Molecules

Molecules are a collection of counted elements (or sub-molecules) with a
charge. While it’s possible to create a molecule by “hand”, it’s much
easier to use the character representation of a molecule, which is
usually what you get when copy/pasting from a source.

``` r
mol("H2O")
#> <mol>
#> [1] H2O
```

And like everything else in R, `mol` objects are vectorized, so you can
serialize an entire column of molecule formulas.

``` r
as_mol(c("H2O", "H+", "Fe(OH)3", "Ca+2"))
#> <mol>
#> [1] H2O     H+      Fe(OH)3 Ca+2
```

You can access the mass, charge, and elemental composition of a molecule
using `mass()`, `charge()`, and `as.data.frame()` or `as.matrix()`

``` r
m <- as_mol(c("H2O", "H+", "Fe(OH)3", "Ca+2"))
mass(m)
#> [1]  18.0150   1.0080 106.8662  40.0784
```

``` r
charge(m)
#> [1] 0 1 0 2
```

``` r
as.data.frame(m)
#>       mol     mass charge H O Fe Ca
#> 1     H2O  18.0150      0 2 1  0  0
#> 2      H+   1.0080      1 1 0  0  0
#> 3 Fe(OH)3 106.8662      0 3 3  1  0
#> 4    Ca+2  40.0784      2 0 0  0  1
```

## Reactions

Reactions are essentially a molecule vector with coefficients (positive
for the left side, negative for the right side). Similar to molecules,
it’s easiest to use the serialized form (conveniently, what is generally
copied/pasted):

``` r
as_reaction("2H2 + O2 = 2H2O")
#> <reaction> 2H2 + O2 = 2H2O   log_k = NA
```

The `is_balanced()` and `balance()` functions will happily balance these
for you, provided you have the correct number of species defined.

``` r
balance("H2 + O2 = H2O")
#> <reaction> 2H2 + O2 = 2H2O   log_k = NA
```

You can access various components of a reaction in the same way as for
molecules:

``` r
r <- as_reaction("2H2 + O2 = 2H2O")
lhs(r)
#> <reaction> 2H2 + O2 =    log_k = NA
```

``` r
rhs(r)
#> <reaction> 2H2O =    log_k = NA
```

``` r
mass(r) # mass balance of the reaction
#> [1] 0
```

``` r
charge(r) # charge balance of the reaction
#> [1] 0
```

``` r
as.data.frame(r)
#>          reaction log_k mol coefficient charge   mass H O
#> 1 2H2 + O2 = 2H2O    NA  H2           2      0  2.016 2 0
#> 2 2H2 + O2 = 2H2O    NA  O2           1      0 31.998 0 2
#> 3 2H2 + O2 = 2H2O    NA H2O          -2      0 18.015 2 1
```

``` r
as.matrix(r)
#>      H  O
#> H2   4  0
#> O2   0  2
#> H2O -4 -2
```

## Molecule and Reaction arithmetic

Various arithmetic operators are available for molecule and reaction
objects, such as `+`, `*` and `==`.

``` r
m <- mol("Fe2O3", "H2O", "NH3", "H+")
m + as_mol("OH-")
#> <mol>
#> [1] Fe2O3OH- H2OOH-   NH3OH-   HOH
```

``` r
m * 2
#> <mol>
#> [1] Fe4O6 H4O2  N2H6  H2+2
```

``` r
m == as_mol(~H2O)
#> [1] FALSE  TRUE FALSE FALSE
```

Reactions have similar arithmetic, with coefficients to various
molecules being added together.

``` r
r1 <- as_reaction("2H2 + O2 = 2H2O", log_k = -46.62)
r1 + as_reaction("H2O = H2O", log_k = 0)
#> <reaction> 2H2 + O2 + H2O = 2H2O + H2O   log_k = -46.62
```

By default the reaction isn’t simplified, but can be using `simplify()`
and `remove_zero_counts()`.

``` r
simplify_reaction(r1 + as_reaction("H2O = H2O", log_k = 0))
#> <reaction> 2H2 + O2 = 2H2O   log_k = -46.62
```

``` r
simplify_reaction(r1 - as_reaction("2H+ + 2OH- = 2H2O", log_k = 14))
#> <reaction> 2H2 + O2 + 0H2O = 2H+ + 2OH-  log_k = -60.62
```

``` r
remove_zero_counts(simplify_reaction(r1 - as_reaction("2H+ + 2OH- = 2H2O", log_k = 14)))
#> <reaction> 2H2 + O2 = 2H+ + 2OH-     log_k = -60.62
```

## The Wish List

There are lots of things missing from this package that should exist,
including the include various parameters to molecules and equations such
as *Δ**H* or aliases (e.g., `CaSO4` as “gypsum”). Additionally, there is
currently no way to indicate hydration in the same way as PHREEQC (e.g.,
`CaSO4:2H2O`). Currently this is possible only as `CaSO4(H2O)2`. Feel
free to [contribute to
development](https://github.com/paleolimbot/chemr) or [submit feature
requests](https://github.com/paleolimbot/chemr/issues) on GitHub.
