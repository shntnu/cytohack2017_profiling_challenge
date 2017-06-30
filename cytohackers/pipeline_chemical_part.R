# convert compounds smiles to structure data format and create cheminer data object
filename='H-BIOA-002.onlycompounds.txt'
ans=fread(filename)
smi=ans$smiles
names(smi)=as.character(ans$broad_sample)
as(smi, "SMIset") 
#cid(smi) <- makeUnique(cid(smi))
saveRDS(smi,paste0(filename,'.molecule.smiles.rds'),compress = TRUE)
sdf <- smiles2sdf(smi)
saveRDS(sdf,paste0(filename,'.molecule.sdf.rds'),compress = TRUE)
write.SDF(sdf, file=paste0(filename,'.molecule.sdf')) 
ap <- sdf2ap(sdf)
saveRDS(ap,paste0(filename,'.molecule.ap.rds'),compress = TRUE)
#list the ones with errors
#which(sapply(as(ap, "list"), length)==1)
fp <- desc2fp(ap, descnames=1024, type="FPset")
saveRDS(fp,paste0(filename,'.molecule.fp.rds'),compress = TRUE)

##############################################################################################################

# normalize the compounds
#ensure using Ambit version>=3.0.1
wget https://sourceforge.net/projects/ambit/files/Ambit2/AMBIT%20applications/ambitcli/ambitcli-3.0.2/ambitcli-3.0.2.jar/download -O ambit.jar --no-check-certificate
java -jar ambit.jar -a standardize -m post -d page=0 -d pagesize=-1       -i H-BIOA-002.onlycompounds.txt.molecule.sdf      -o H-BIOA-002.onlycompounds.txt.molecule.ambit.sdf      -d tautomers=true -d splitfragments=true -d implicith=true      -d smiles=true -d smilescanonical=false -d inchi=true      -d stereo=true -d neutralise=true -d isotopes=true

#data is part of https://github.com/bioinf-jku/project_BBDD

##############################################################################################################
# pairwise similarity matrix

library(data.table)
library(ChemmineR)

#aps = readRDS('H-BIOA-002-3.txt.molecule.ap.rds')
#fps = desc2fp(aps,descnames=1024,type='FPset')

fps = readRDS('H-BIOA-002.onlycompounds.txt.molecule.fp.rds')
L = length(fps)

sims = matrix(nrow=L, ncol=L)
#sims =
for (x in 1:L) {
  for (y in 1:L) {
    sims[x, y] = fpSim(fps[x], fps[y])
    #sims[x] = fpSim(fps[x], fps)
  }
}


colnames(sims) <- rownames(fps@fpma)
rownames(sims) <- rownames(fps@fpma)
saveRDS(sims,'sims.Rds')