# target_sequencing_analysis
1. (vcf.gz.stats instead of vcf.gz.stat in this rule)we set --latency-wait 13500 instead of 500, since we didnt put output of mutect in rule all and this means rule mutect2 will be considered as success when all output files of mutect2 exits within latency-wait.
2. You cannot directly use if-else in the input directive, but you can use functions or lambda expressions to resolve the input dynamically based on conditions.

3. todo
   1) use a resource config to asign mem, threahs to each rule
   2) pcr status + trim of each reads should be set in unit instead of overall config
