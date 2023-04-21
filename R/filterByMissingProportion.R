#' filterByMissingProportion
#' @export filterByMissingProportion
#'
#' @import data.table
filterByMissingProportion <- function(longIntensityDT,
                                  conditionColname,
                                  runIdColname,
                                  idCol = "ID",
                                  imputedColname = "Imputed",
                                  minPropAvail = 0.5){
    # how many replicates in each conditions
    longIntensityDT <- data.table(longIntensityDT)
    filterDT <- unique(longIntensityDT[, .(get(conditionColname), get(runIdColname))])[, .(max_count = .N), by = .(V1)]
    filterDT <- merge(longIntensityDT, filterDT, by.x = conditionColname, by.y = "V1", all.x = T)
    filterDT <- filterDT[get(imputedColname) == 0, .(count_rep = .N,
                                                     max_count = max(max_count, na.rm = T)),
                         by = c(conditionColname, idCol)][, repPC := count_rep/max_count]

    isPresent <- filterDT[repPC >= minPropAvail, unique(get(idCol))]
    finalFiltered <- longIntensityDT[get(idCol) %in% isPresent]
    return(finalFiltered)
}
