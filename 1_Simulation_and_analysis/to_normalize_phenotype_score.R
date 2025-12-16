args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please offer hXX list")
}

hXX = args 

simIDList = read.table(paste("./",hXX,"_500.list", sep=""))
ids=scan(paste("./",hXX,"_500.list", sep=""),what="character")


tau_list = c('0.7','0.8','0.9','1.0','1.1','1.2','1.3')
rho_list = c('0.0')
count = 0

for (tau in tau_list){
  for (rho in rho_list){
    for (simID in ids) {
      count = count + 1
      p1=read.table(paste(hXX,"/",simID,"/pheno_neutral_recalculation/gravel.del.only.",hXX,".",simID,".trees.p1.exon.in.smapled.neutral.mutation.rho.",rho,".tau.",tau,".pheno.txt",sep=""))
      p2=read.table(paste(hXX,"/",simID,"/pheno_neutral_recalculation/gravel.del.only.",hXX,".",simID,".trees.p2.exon.in.smapled.neutral.mutation.rho.",rho,".tau.",tau,".pheno.txt",sep=""))
      p3=read.table(paste(hXX,"/",simID,"/pheno_neutral_recalculation/gravel.del.only.",hXX,".",simID,".trees.p3.exon.in.smapled.neutral.mutation.rho.",rho,".tau.",tau,".pheno.txt",sep=""))
      p.sd=sd(c((abs(p1$phenoApos)+abs(p1$phenoBpos)+abs(p1$phenoCpos)+abs(p1$phenoNONEpos)),(abs(p2$phenoApos)+abs(p2$phenoBpos)+abs(p2$phenoCpos)+abs(p2$phenoNONEpos)),(abs(p3$phenoApos)+abs(p3$phenoBpos)+abs(p3$phenoCpos)+abs(p3$phenoNONEpos))))
      p1_normalized = p1
      p2_normalized = p2
      p3_normalized = p3
      cols_to_modify = c("phenoApos", "phenoBpos", "phenoCpos", "phenoNONEpos")
      p1_normalized[cols_to_modify] = p1[cols_to_modify] / p.sd
      p2_normalized[cols_to_modify] = p2[cols_to_modify] / p.sd
      p3_normalized[cols_to_modify] = p3[cols_to_modify] / p.sd
      
      p1_normalized_path = paste(hXX,"/",simID,"/pheno_neutral_recalculation/gravel.del.only.",hXX,".",simID,".trees.p1.exon.in.smapled.neutral.mutation.rho.",rho,".tau.",tau,".pheno.normalized.txt",sep="")
      p2_normalized_path = paste(hXX,"/",simID,"/pheno_neutral_recalculation/gravel.del.only.",hXX,".",simID,".trees.p2.exon.in.smapled.neutral.mutation.rho.",rho,".tau.",tau,".pheno.normalized.txt",sep="")
      p3_normalized_path = paste(hXX,"/",simID,"/pheno_neutral_recalculation/gravel.del.only.",hXX,".",simID,".trees.p3.exon.in.smapled.neutral.mutation.rho.",rho,".tau.",tau,".pheno.normalized.txt",sep="")
      
      write.table(p1_normalized, file = p1_normalized_path, sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
      write.table(p2_normalized, file = p2_normalized_path, sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
      write.table(p3_normalized, file = p3_normalized_path, sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
      print(count)
    }
  }
}

#End

#Normalize phenotype score for deleterious mutation

simIDList = read.table(paste("./",hXX,"_500.list", sep=""))
ids=scan(paste("./",hXX,"_500.list", sep=""),what="character")

tau_list = c('70','80','90','100','110','120','130')
rho_list = c('100')
count = 0

for (tau in tau_list){
  for (rho in rho_list){
    for (simID in ids) {
      count = count + 1
      p1=read.table(paste(hXX,"/",simID,"/pheno/gravel.del.only.",hXX,".",simID,".rho.",rho,".tau.",tau,".p1.pheno.txt",sep=""))
      p2=read.table(paste(hXX,"/",simID,"/pheno/gravel.del.only.",hXX,".",simID,".rho.",rho,".tau.",tau,".p2.pheno.txt",sep=""))
      p3=read.table(paste(hXX,"/",simID,"/pheno/gravel.del.only.",hXX,".",simID,".rho.",rho,".tau.",tau,".p3.pheno.txt",sep=""))
      p.sd=sd(c((abs(p1$phenoApos)+abs(p1$phenoBpos)+abs(p1$phenoCpos)+abs(p1$phenoNONEpos)),(abs(p2$phenoApos)+abs(p2$phenoBpos)+abs(p2$phenoCpos)+abs(p2$phenoNONEpos)),(abs(p3$phenoApos)+abs(p3$phenoBpos)+abs(p3$phenoCpos)+abs(p3$phenoNONEpos))))
      p1_normalized = p1
      p2_normalized = p2
      p3_normalized = p3
      cols_to_modify = c("phenoApos", "phenoBpos", "phenoCpos", "phenoNONEpos")
      p1_normalized[cols_to_modify] = p1[cols_to_modify] / p.sd
      p2_normalized[cols_to_modify] = p2[cols_to_modify] / p.sd
      p3_normalized[cols_to_modify] = p3[cols_to_modify] / p.sd
      
      p1_normalized_path = paste(hXX,"/",simID,"/pheno/gravel.del.only.",hXX,".",simID,".rho.",rho,".tau.",tau,".p1.pheno.normalized.txt",sep="")
      p2_normalized_path = paste(hXX,"/",simID,"/pheno/gravel.del.only.",hXX,".",simID,".rho.",rho,".tau.",tau,".p2.pheno.normalized.txt",sep="")
      p3_normalized_path = paste(hXX,"/",simID,"/pheno/gravel.del.only.",hXX,".",simID,".rho.",rho,".tau.",tau,".p3.pheno.normalized.txt",sep="")
      
      write.table(p1_normalized, file = p1_normalized_path, sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
      write.table(p2_normalized, file = p2_normalized_path, sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
      write.table(p3_normalized, file = p3_normalized_path, sep = " ", row.names = TRUE, col.names = TRUE, quote = FALSE)
      print(count)
    }
  }
}

