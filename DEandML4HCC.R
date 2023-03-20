library(curatedTCGAData)
library(TCGAutils)
library(MultiAssayExperiment)
library(patchwork)

# Get all codes on TCGA
data('diseaseCodes', package = "TCGAutils")
diseaseCodes[,c(1,4)]

# Get sample types on TCGA
data(sampleTypes, package = "TCGAutils")
head(sampleTypes)

liver.rawdata = curatedTCGAData(c("LIHC"), assays =  c("RNASeq2Gene"), 
                                version = '2.0.1',
                                dry.run = FALSE)

rownames(liver.rawdata)
colnames(liver.rawdata)
colData(liver.rawdata)
#getClinicalNames("LIHC")
table(table(sampleMap(liver.rawdata)$primary)) 

# Get how many sample types in each experiment/cancer 
sampleTables(liver.rawdata)
sampleTables(liver.rawdata,vial = TRUE)

# extract assay data

liver.assaydata <- (assay(liver.rawdata[[1]]))

liver.assaydata <- round(liver.assaydata)



##################### DE ##########################


# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked


colnames(liver.assaydata)

sampleInfo <- data.frame("sampleID", "type")

colnames(N) <- c("sampleID", "type")

PT <- data.frame(colnames(liver.assaydata[,which(substr(colnames(liver.assaydata),14,16) == "01A")]),"PT")
N <- data.frame(colnames(liver.assaydata[,which(substr(colnames(liver.assaydata),14,16) == "11A")]),"N")

sampleInfo <- rbind.data.frame (PT,N)
rownames(sampleInfo) <- sampleInfo[,1]

df = subset(liver.assaydata, select = c(colnames(liver.assaydata[,which(substr(colnames(liver.assaydata),14,16) == "01A")]),colnames(liver.assaydata[,which(substr(colnames(liver.assaydata),14,16) == "11A")])) )

all(colnames(df) == rownames(sampleInfo))

dseqData <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(df), colData=sampleInfo, design= ~ type)

dseqData <- DESeq(dseqData)

liver.res <- results(dseqData)

summary(liver.res)
DESeq2::plotMA(liver.res)

############################# ML work #######################

ml.data <- as.data.frame(df[rownames(liver.res[(which((liver.res$log2FoldChange > 3 | liver.res$log2FoldChange < -3) & liver.res$padj < 0.05)) ,]),])

tml <-  as.data.frame(t(ml.data))

rownames(tml) == rownames(sampleInfo)

tml$type <- sampleInfo$type

tml$type == sampleInfo$type

library(caret)
allMLmodels <- unique(caret::modelLookup()$model)


tml$type <- factor(tml$type)
levels(tml$type)

set.seed(123)
trainIndex <- createDataPartition(tml$type, p = 0.70, list = FALSE)
train <- tml[trainIndex, ]
test <- tml[-trainIndex, ]

# Train model using random forest algorithm
rf_model <- train(type ~ ., data = train, method = "rf" )

# Make predictions on test set
predictions <- predict(rf_model, newdata = test )

# Evaluate model performance
cmatrix3 <- confusionMatrix(predictions, test$type)

install.packages("")
BiocManager::install("ggtree")
remotes::install_github('YuLab-SMU/ggtree')

BiocManager::install("GDCRNATools")

library(GDCRNATools)

project <- 'TCGA-CHOL_0406'
rnadir <- paste(project, 'RNAseq', sep='/')

project2 <- 'TARGET-AML'
rnadir2 <- paste(project, 'RNAseq', sep='/')

gdcRNADownload(project.id     = 'TARGET-AML', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir2)


