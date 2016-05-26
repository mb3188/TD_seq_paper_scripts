resG3_NYCA04 = results(dds1, contrast=c("group","G3","NYCA04"))
resG3_NYCA04_signf<-resG3_NYCA04[which(resG3_NYCA04$padj < 0.1),]
write.table(resG3_NYCA04_signf, "resG3_NYCA04_signf")
write.table(rownames (resG3_NYCA04_signf), "row_resG3_NYCA04_signf")
resG3_NYCA04signif_up<-resG3_NYCA04[which(resG3_NYCA04$padj<0.01 & resG3_NYCA04$log2FoldChange > 0),]
resG3_NYCA04signif_down<-resG3_NYCA04[which(resG3_NYCA04$padj<0.01 & resG3_NYCA04$log2FoldChange < 0),]
write.table(rownames (resG3_NYCA04signif_up), "row_resG3_NYCA04signif_up")
write.table(rownames (resG3_NYCA04signif_down), "row_resG3_NYCA04signif_down")

resG3_NYCB20 = results(dds1, contrast=c("group","G3","NYCB20"))
resG3_NYCB20_signf<-resG3_NYCB20[which(resG3_NYCB20$padj < 0.1),]
write.table(resG3_NYCB20_signf, "resG3_NYCB20_signf")
write.table(rownames (resG3_NYCB20_signf), "row_resG3_NYCB20_signf")
resG3_NYCB20signif_up<-resG3_NYCB20[which(resG3_NYCB20$padj<0.01 & resG3_NYCB20$log2FoldChange > 0),]
resG3_NYCB20signif_down<-resG3_NYCB20[which(resG3_NYCB20$padj<0.01 & resG3_NYCB20$log2FoldChange < 0),]
write.table(rownames (resG3_NYCB20signif_up), "row_resG3_NYCB20signif_up")
write.table(rownames (resG3_NYCB20signif_down), "row_resG3_NYCB20signif_down")

resG3_NYCD15 = results(dds1, contrast=c("group","G3","NYCD15"))
resG3_NYCD15_signf<-resG3_NYCD15[which(resG3_NYCD15$padj < 0.1),]
write.table(resG3_NYCD15_signf, "resG3_NYCD15_signf")
write.table(rownames (resG3_NYCD15_signf), "row_resG3_NYCD15_signf")
resG3_NYCD15signif_up<-resG3_NYCD15[which(resG3_NYCD15$padj<0.01 & resG3_NYCD15$log2FoldChange > 0),]
resG3_NYCD15signif_down<-resG3_NYCD15[which(resG3_NYCD15$padj<0.01 & resG3_NYCD15$log2FoldChange < 0),]
write.table(rownames (resG3_NYCD15signif_up), "row_resG3_NYCD15signif_up")
write.table(rownames (resG3_NYCD15signif_down), "row_resG3_NYCD15signif_down")

resG3_NYCE32 = results(dds1, contrast=c("group","G3","NYCE32"))
resG3_NYCE32_signf<-resG3_NYCE32[which(resG3_NYCE32$padj < 0.1),]
write.table(resG3_NYCE32_signf, "resG3_NYCE32_signf")
write.table(rownames (resG3_NYCE32_signf), "row_resG3_NYCE32_signf")
resG3_NYCE32signif_up<-resG3_NYCE32[which(resG3_NYCE32$padj<0.01 & resG3_NYCE32$log2FoldChange > 0),]
resG3_NYCE32signif_down<-resG3_NYCE32[which(resG3_NYCE32$padj<0.01 & resG3_NYCE32$log2FoldChange < 0),]
write.table(rownames (resG3_NYCE32signif_up), "row_resG3_NYCE32signif_up")
write.table(rownames (resG3_NYCE32signif_down), "row_resG3_NYCE32signif_down")

resG3_GOR = results(dds1, contrast=c("group","G3","GOR/03/69"))
resG3_GORsignf<-resG3_GOR[which(resG3_GOR$padj < 0.1),]
write.table(resG3_GORsignf, "resG3_GORsignf")
write.table(rownames (resG3_GORsignf), "row_resG3_GORsignf")
resG3_GORsignif_up<-resG3_GOR[which(resG3_GOR$padj<0.01 & resG3_GOR$log2FoldChange > 0),]
resG3_GORsignif_down<-resG3_GOR[which(resG3_GOR$padj<0.01 & resG3_GOR$log2FoldChange < 0),]
write.table(rownames (resG3_GORsignif_up), "row_resG3_GORsignif_up")
write.table(rownames (resG3_GORsignif_down), "row_resG3_GORsignif_down")

resG3_NYCF20 = results(dds1, contrast=c("group","G3","NYCF20"))
resG3_NYCF20_signf<-resG3_NYCF20[which(resG3_NYCF20$padj < 0.1),]
write.table(resG3_NYCF20_signf, "resG3_NYCF20_signf")
write.table(rownames (resG3_NYCF20_signf), "row_resG3_NYCF20_signf")
resG3_NYCF20signif_up<-resG3_NYCF20[which(resG3_NYCF20$padj<0.01 & resG3_NYCF20$log2FoldChange > 0),]
resG3_NYCF20signif_down<-resG3_NYCF20[which(resG3_NYCF20$padj<0.01 & resG3_NYCF20$log2FoldChange < 0),]
write.table(rownames (resG3_NYCF20signif_up), "row_resG3_NYCF20signif_up")
write.table(rownames (resG3_NYCF20signif_down), "row_resG3_NYCF20signif_down")


resG3_NYCG31 = results(dds1, contrast=c("group","G3","NYCG31"))
resG3_NYCG31_signf<-resG3_NYCG31[which(resG3_NYCG31$padj < 0.1),]
write.table(resG3_NYCG31_signf, "resG3_NYCG31_signf")
write.table(rownames (resG3_NYCG31_signf), "row_resG3_NYCG31_signf")
resG3_NYCG31signif_up<-resG3_NYCG31[which(resG3_NYCG31$padj<0.01 & resG3_NYCG31$log2FoldChange > 0),]
resG3_NYCG31signif_down<-resG3_NYCG31[which(resG3_NYCG31$padj<0.01 & resG3_NYCG31$log2FoldChange < 0),]
write.table(rownames (resG3_NYCG31signif_up), "row_resG3_NYCG31signif_up")
write.table(rownames (resG3_NYCG31signif_down), "row_resG3_NYCG31signif_down")



resG3_SD2 = results(dds1, contrast=c("group","G3","SD2"))
resG3_SD2_signf<-resG3_SD2[which(resG3_SD2$padj < 0.1),]
write.table(resG3_SD2_signf, "resG3_SD2_signf")
write.table(rownames (resG3_SD2_signf), "row_resG3_SD2_signf")
resG3_SD2signif_up<-resG3_SD2[which(resG3_SD2$padj<0.01 & resG3_SD2$log2FoldChange > 0),]
resG3_SD2signif_down<-resG3_SD2[which(resG3_SD2$padj<0.01 & resG3_SD2$log2FoldChange < 0),]
write.table(rownames (resG3_SD2signif_up), "row_resG3_SD2signif_up")
write.table(rownames (resG3_SD2signif_down), "row_resG3_SD2signif_down")


resG3_NYCC37= results(dds1, contrast=c("group","NYCC37","G3"))
resG3_NYCC37_signf<-resG3_NYCC37[which(resG3_NYCC37$padj < 0.1),]
write.table(resG3_NYCC37_signf, "resG3_NYCC37_signf")
write.table(rownames (resG3_NYCC37_signf), "row_resG3_NYCC37_signf")
resG3_NYCC37signif_up<-resG3_NYCC37[which(resG3_NYCC37$padj<0.01 & resG3_NYCC37$log2FoldChange > 0),]
resG3_NYCC37signif_down<-resG3_NYCC37[which(resG3_NYCC37$padj<0.01 & resG3_NYCC37$log2FoldChange < 0),]
write.table(rownames (resG3_NYCC37signif_up), "row_resG3_NYCC37signif_up")
write.table(rownames (resG3_NYCC37signif_down), "row_resG3_NYCC37signif_down")

resG3_B7268= results(dds1, contrast=c("group","B7268","G3"))
resG3_B7268_signf<-resG3_B7268[which(resG3_B7268$padj < 0.1),]
write.table(resG3_B7268_signf, "resG3_B7268_signf")
write.table(rownames (resG3_B7268_signf), "row_resG3_B7268_signf")
resG3_B7268signif_up<-resG3_B7268[which(resG3_B7268$padj<0.01 & resG3_B7268$log2FoldChange > 0),]
resG3_B7268signif_down<-resG3_B7268[which(resG3_B7268$padj<0.01 & resG3_B7268$log2FoldChange < 0),]
write.table(rownames (resG3_B7268signif_up), "row_resG3_B7268signif_up")
write.table(rownames (resG3_B7268signif_down), "row_resG3_B7268signif_down")

            
resG3_B7268M= results(dds1, contrast=c("group","B7268M","G3"))
resG3_B7268M_signf<-resG3_B7268M[which(resG3_B7268M$padj < 0.1),]
write.table(resG3_B7268M_signf, "resG3_B7268M_signf")
write.table(rownames (resG3_B7268M_signf), "row_resG3_B7268M_signf")
resG3_B7268Msignif_up<-resG3_B7268M[which(resG3_B7268M$padj<0.01 & resG3_B7268M$log2FoldChange > 0),]
resG3_B7268Msignif_down<-resG3_B7268M[which(resG3_B7268M$padj<0.01 & resG3_B7268M$log2FoldChange < 0),]
write.table(rownames (resG3_B7268Msignif_up), "row_resG3_B7268Msignif_up")
write.table(rownames (resG3_B7268Msignif_down), "row_resG3_B7268Msignif_down")             
             
# Count number of occurence of NA in each line, or zero

awk -F\#N/A '{print NF-1}' <fileName>
awk -F\0 '{print NF-1}' <fileName>