
SI_loc <-'/home/jkb4y/cphgdesk_share/Achilleas/eQTLs_Oct2012/data/si_SNP_NoMHC_20121024_Query.txt'
SI_table <- read.table(SI_loc, T)

metal$SI_annotation <- ifelse(metal$SNP %in% SI_table$conditional_SNP, 'SI','none')
print(metal[metal$SI_annotation=='SI',])

