# =====================================================================
# 21 - sessionInfo template
# Run after the pipeline to record the exact R environment used.
# The original analysis was performed with R 4.2.0.
# =====================================================================

outFile <- "sessionInfo.txt"
si <- sessionInfo()
writeLines(capture.output(si), outFile)
cat("sessionInfo written to", normalizePath(outFile), "\n")
