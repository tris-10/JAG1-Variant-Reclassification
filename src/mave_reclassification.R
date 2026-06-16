# Check for and install dependencies
packages <- c("devtools", "optparse", "dplyr")
options(download.file.method = "wininet")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# Non-CRAN installation
if (!requireNamespace("maveLLR", quietly = TRUE)) {
  #options(download.file.method = "wininet")  # If running on windows
  devtools::install_github("jweile/maveLLR")
}

# Load in libraries
library(optparse)
library(maveLLR)

# Input parameters
option_list <- list(
  make_option(c("-i", "--input"), type="character", 
              help="Input path for TSV file containing positive reference, 
              negative reference, and library MAVE scores", metavar="file"),
  make_option(c("-o", "--output"), type="character", 
              help="Output path for TSV file containing converted LLRp 
              MAVE scores.", metavar="file"),
  make_option(c("-s", "--scoreLabel"), type="character",
              help="Column name for MAVE scores to extract",
              metavar="string"),
  make_option(c("-r", "--refLabel"), type="character",
              help="Column name for reference set categories 
              (benign: -1, library: 0, pathogenic: 1)",
              metavar="string"),
  make_option(c("-b", "--bw"), type="double", default=0.1, 
              help="Bandwidth for LLR KDE [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Stop with an error if missing args
if (is.null(opt$input) || is.null(opt$output) || 
    is.null(opt$scoreLabel) || is.null(opt$refLabel)) {
  stop("--input, --output, --scoreLabel, and --refLabel must all be specified.
       \nUse --help for more info.")
}

inputFile <- opt$input
outputFile <- opt$output
scoreCol <- opt$scoreLabel
refCol <- opt$refLabel

# Load in TSV with MAVE scores
maveTable <- read.table(inputFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter rows: benign (-1), library (0), pathogenic (1) 
mutTable <- subset(maveTable, maveTable[[refCol]] == 1)
wtTable <- subset(maveTable, maveTable[[refCol]] == -1)
libTable <- subset(maveTable, maveTable[[refCol]] == 0)

# Retrieve positive and negative reference variant scores
posScores <- as.numeric(mutTable[[scoreCol]])
negScores <- as.numeric(wtTable[[scoreCol]])
# Retrieve library variant scores
scores <- as.numeric(libTable[[scoreCol]])

# Calculate the LLR using kernel density estimation
llr <- buildLLR.kernel(posScores,negScores,bw=opt$bw,kernel="gaussian")

# Draw a summary plot of the resulting LLR
drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)

# Transform them into their corresponding LLRs
result <- llr$llr(scores)

posResult <- llr$llr(posScores)
negResult <- llr$llr(negScores)

#append associated results to library, positive ref, and negative ref datasets
libTable$llrp_score <- result
mutTable$llrp_score <- posResult
wtTable$llrp_score <- negResult
# Combine into output file
outTable <- rbind(libTable, wtTable, mutTable)

# Categorize variants based off Tavtigian et al., 2018 thresholds
outTable$llrp_functional_consequence <- cut(
  outTable$llrp_score,
  breaks = c(-Inf, -1.27, -0.31, 0.31, 0.63, 1.27, 2.54, Inf),
  labels = c("Norm", "LNorm", "I", "PSu", "PM", "PSt", "PVSt"),
  right = TRUE,      # intervals are (a, b]
  include.lowest = TRUE
)

# Export to TSV
write.table(outTable, file = outputFile, sep = "\t", row.names = FALSE, quote = FALSE)