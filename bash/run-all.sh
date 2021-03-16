#!/bin/bash
for seed in {1..100}
do
  Rscript ../simulation/ld_homosc.R $seed
  Rscript ../simulation/ld_heterosc.R $seed
  Rscript ../simulation/hd_homosc.R $seed
  Rscript ../simulation/hd_heterosc.R $seed
done

