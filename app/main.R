# Set Working directory via executable
params = commandArgs(trailingOnly=TRUE)

#name pars
upload_folder = params[1]
gmt_folder = params[2]
output_folder = params[3]
by = params[4]


library(EnrichmentAnalysisStepR)

enrichment_workflow_step(upload_folder,
                         gmt_folder = gmt_folder,
                         output_folder = output_folder,
                         by = by)



