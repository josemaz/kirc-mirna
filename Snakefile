rule all:
	input:
		"data/tables/DEGS/res-all-stage1_ctrl.tsv",
		"data/tables/exp/ctrl-all.tsv"

rule biomart:      
    output:
        "data/tables/biomart.csv"
    shell:
        "Rscript R/biomart.R"

rule clean:
    input:
        "data/tables/biomart.csv"
    output:
        "data/RDS/rnaseq-clean.rds",
        "data/RDS/mirna-clean.rds"
    shell:
        "Rscript R/01-clean-merge.R"

rule norm:
    input:
        "data/RDS/rnaseq-clean.rds",
        "data/RDS/mirna-clean.rds"
    output:
        "data/RDS/rnaseq-norm.rds",
        "data/RDS/mirna-norm.rds"
    shell:
        "Rscript R/02-norm.R"

rule save_exp:
    input:
        "data/RDS/rnaseq-norm.rds",
        "data/RDS/mirna-norm.rds"
    output:
        "data/tables/exp/ctrl-all.tsv"
    shell:
        "Rscript R/03-save-exp.R"

rule degs:
	input:
		"data/RDS/rnaseq-norm.rds",
		"data/RDS/mirna-norm.rds",
		"data/tables/biomart.csv"
	output:
		"data/tables/DEGS/res-all-stage1_ctrl.tsv",
		"data/plots/Volcanos/vol-mir-stage1_ctrl.png",
		"data/plots/Volcanos/vol-rna-stage1_ctrl.png",
		"data/plots/venns/venn-mir-LFC_10-all.png"
	shell:
		"Rscript R/04-deg.R"

