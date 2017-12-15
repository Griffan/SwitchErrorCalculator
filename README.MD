This tool is for haplotype phasing accuracy evaluation. It accepts two sorted(by chr,pos) VCFs as input, without requirements for exactly the same marker set or sample set, as long as your target sample is included in gold_standard.vcf.
perl switch_error_calculator.pl <need_to_evaluate.vcf> <gold_standard.vcf>
