# Bioinformatics Workshop
## Informal Bioinformatics Workshop at Smith College: Fall 2019

buoyantWeightAnalysis_base.R = script in base R (no packages loaded) to analyze buoyant weight data

buoyantWeight.csv = buoyant weight data

## Project Description:
Montastraea cavernosa disease susceptibility experiment

Flower Garden Banks National Marine Sanctuary

Sample information:

genet = number to identify unique colonies (genetically distinct individuals)

rep = replicate fragments from the same genet

treat = bacterial treatment (C = control; V = Vibrio challenge)

sam = unique sample ID that includes genet and treatment

bank = bank where that genet was collected (w = West Flower Garden Bank; e = East Flower Garden Bank)

pheno = disease susceptibility phenotype (r = resistant, s = susceptible)

id = description of the replicate fragment used to distinguish fragments of the same genet in photographs

Data types:

weight = buoyant weight in grams

polyp counts = number of polyps counted on the fragment

redness = mean red channel intensity for the micrograph (0-255 scale, no units)

tod = time of death (day post initial infection)

dead = logical descriptor (0 = did not die; 1 = did die)

Gene expression counts matrix:

Columns = samples 

Decoding the sample names...

MC = Montastraea cavernosa; number = genet; letter = replicate (A or B)

Rows = isogroups ("genes")

Cells in the matrix contain the number of times a copy of that isogroup was detected in each sample
