#!/usr/bin/env sh

date; echo "expr_emat"
script/expr_emat.R --verbose
date; echo "emat_reduce"
script/emat_reduce.R --verbose
date; echo "emat_norm"
script/emat_norm.R --verbose
date; echo "emat_quality"
script/emat_quality.R --verbose
date; echo "emat_filter"
script/emat_filter.R --verbose
date; echo "emat_cut"
script/emat_cut.R --verbose


# steps below here not used to make data cut
date; echo "load_derived"
script/load_derived.R --verbose
date; echo "emat_derive"
script/emat_derive.R --verbose
date; echo "emat_filter again"
script/emat_filter.R --verbose --pattern.yes='^derive' --best.replicas F
date; echo "done"

