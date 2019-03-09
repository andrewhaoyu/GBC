rm(list=ls(all=TRUE))
gc()

library(CGEN, lib.loc="/spin1/users/wheelerwi/nilanjan/wga/CGEN")


snp.list      <- list(format="ldat", GLU="", PLINK="", use.GLU=0, use.PLINK=0)
snp.list$dir  <- "/spin1/users/wheelerwi/kai/panscan/data/geno/nov19_merge123/part/" 
snp.list$file <- paste("PanScan123_part", 1:500, ".ldat.gz",  sep="")
snp.list$nsnps.vec <- c(rep(1202, 499), 896)

f <- "/data/wheelerwi/kai/panscan/data/pheno/2015/may29_2015_123/pheno.txt.xls"

pheno.list <- list(file=f, header=1, delimiter="\t", id.var="ID", response.var="case") 
pheno.list$main.vars <- c(
"AGE_LESS_51", "AGE_61_70", "AGE_71_80", "AGE_OVER_80",
"Gender",
"REGION_CNE_noATBC", "REGION_SE",
"Smoke_Stat_Former", "Smoke_Stat_Current", 
"EV1", "EV2", "EV3", "EV4", "EV5",
"STAGE2_EV1", "STAGE2_EV2", "STAGE2_EV3", "STAGE2_EV4", "STAGE2_EV5",
"STAGE3_EV1", "STAGE3_EV2", "STAGE3_EV3", "STAGE3_EV4", "STAGE3_EV5",
"STAGE_1", "STAGE_2")
 
pheno.list$int.vars    <- c("Smoke_Stat_Former", "Smoke_Stat_Current")
pheno.list$strata.vars <- c("Gender",
"REGION_CNE_noATBC", "REGION_SE",
"EV1", "EV2", "EV3", "EV4", "EV5",
"STAGE2_EV1", "STAGE2_EV2", "STAGE2_EV3", "STAGE2_EV4", "STAGE2_EV5",
"STAGE3_EV1", "STAGE3_EV2", "STAGE3_EV3", "STAGE3_EV4", "STAGE3_EV5",
"STAGE_1", "STAGE_2")

pheno.list$keep.vars <- c(pheno.list$id.var, pheno.list$response.var, pheno.list$main.vars, pheno.list$int.vars)

l0 <- list(var="cohort_char", operator="!=", value="ATBC")
pheno.list$subsetData <- list(l0)

out.dir <- "/data/wheelerwi/kai/panscan/results/2015/jul29_scan/out/"

op <- list(n.jobs=1000, out.dir=out.dir)
op$begin.commands.R    <- 'library(CGEN, lib.loc="/data/wheelerwi/nilanjan/wga/CGEN")'
op$begin.commands.qsub <- "module load R"
op$qsub.cmd            <- "bash "

GxE.scan.partition(snp.list, pheno.list, op=op)


