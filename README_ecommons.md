### Data README

Overview of data associated "Synchronized seasonal excretion of multiple coronaviruses in Australian Pteropus spp is associated with co-infections in juvenile and sub-adult bats"


### individual_variant_covariates.csv

Dataset containing metadata and test results from individual bats. Columns include:

- accession_update: key for each sampling instance. This variable is a combination of sampling type, location, and sampling instance at a location
- site: name of roost
- sample_id: individual sample id, unique within accession_update only
- type: type of variant: beta 2d.i, beta 2d.ii, beta 2d.iii, beta 2d.iv, beta 2d.v, or beta 2d.vi
- variant_positive: TRUE or FALSE, for variant testing positive
- bat_age: age of bat, classified as juvenile, sub-adult, or adult
- bat_sex: male or female
- bat_species: bff (for black fying fox) and ghff (for grey-headed flying fox)
- sampling_date: date of sample collection


### combined_out_variant.csv

Dataset containing metadata and test results from pools of bats. Individuals are included as pools of size one. The columns include:

- site: name of roost
- n: pool size
- type: type of variant: beta 2d.i, beta 2d.ii, beta 2d.iii, beta 2d.iv, beta 2d.v, or beta 2d.vi
- variant_positive: TRUE or FALSE, for variant testing positive
- date_end: date of sample collection
- session_type: C for catching and R for underoost fecal collection
