rule all:
    input:
        "pipeline/rnaseq-clean.rds",
        "pipeline/mirna-clean.rds"

rule biomart:      
    output:
        "pipeline/biomart.csv"
    shell:
        "Rscript R/biomart.R"

rule clean:
    input:
        "pipeline/biomart.csv"
    output:
        "pipeline/rnaseq-clean.rds",
        "pipeline/mirna-clean.rds"
    shell:
        "Rscript R/01-clean-merge.R"