eb_crc<-read.table('CRC_longRNAs.txt',header = T)
eb_control<-read.table('Healthy_longRNAs.txt',header = T)
lncRNA_db<-read.csv('lncRNA_db.csv',header = T)
mRNA_db<-read.csv('mRNA_db.csv',header = T)

eb_crc_lnc<-merge(eb_crc,lncRNA_db,by.x = 'Gene.symbol',by.y = 'gene_name',all = F)
eb_crc_m<-merge(eb_crc,mRNA_db,by.x = 'Gene.symbol',by.y = 'gene_name',all = F)
eb_control_lnc<-merge(eb_control,lncRNA_db,by.x = 'Gene.symbol',by.y = 'gene_name',all = F)
eb_control_m<-merge(eb_control,mRNA_db,by.x = 'Gene.symbol',by.y = 'gene_name',all = F)

eb_crc_c_lnc<-merge(eb_crc_lnc,eb_control_lnc,by='Gene.symbol',all=F)
eb_crc_c_m<-merge(eb_crc_m,eb_control_m,by='Gene.symbol',all=F)
rownames(eb_crc_c_lnc)<-eb_crc_c_lnc[,1]
rownames(eb_crc_c_m)<-eb_crc_c_m[,1]
eb_crc_c_lnc<-eb_crc_c_lnc[,-1]
eb_crc_c_m<-eb_crc_c_m[,-1]

eb_crc_c_lnc<-eb_crc_c_lnc[,grep('CRC|Healthy',colnames(eb_crc_c_lnc))]
eb_crc_c_m<-eb_crc_c_m[,grep('CRC|Healthy',colnames(eb_crc_c_m))]

write.csv(eb_crc_c_lnc,'eb_crc_c_lnc.csv')
write.csv(eb_crc_c_m,'eb_crc_c_m.csv')
