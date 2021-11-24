rule all:
    input:
        "data/RDS/rnaseq-norm.rds",
        "data/RDS/mirna-norm.rds"
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
        "data/tables/exp/ctrl-all.tsv"
    shell:
        "Rscript R/02-norm.R"