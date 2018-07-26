
#' Update clinical data from follow-up table
#' 
#' Survival and vital status are not updated in the clinical data
#' obtained from TCGA, so we need to do it manually.
#' 
#' The data is assumed to be obtained from TCGAbiolinks
#'
#' @param clinical data.frame with clinical information from TCGA project
#' @param follow.up data.frame with follow-up information from TCGA project
#'
#' @return
#' @export
#'
#' @examples
#' library(dplyr)
#' data(clinical, gdc, package = 'skcm.data')
#' clinical.new <- update.survival.from.followup(clinical$all, gdc$follow.up)
#' 
#' problem.barcodes <- c('TCGA-FR-A726', 'TCGA-D3-A8GP', 'TCGA-DA-A1HW', 'TCGA-DA-A1I1', 
#'                       'TCGA-DA-A1I5', 'TCGA-DA-A1I7', 'TCGA-DA-A1IB', 'TCGA-DA-A95W', 
#'                       'TCGA-DA-A95X', 'TCGA-EB-A41A', 'TCGA-EE-A2GJ', 'TCGA-FR-A726', 
#'                       'TCGA-HR-A5NC', 'TCGA-XV-A9VZ')
#'
#' clinical.new %>% filter(bcr_patient_barcode %in% problem.barcodes)
#'
#'
#' #
#' #
#' # Original clinical information from TCGAbiolinks
#'
#' clinical$all %>% filter(bcr_patient_barcode %in% problem.barcodes) %>%
#'  select(bcr_patient_barcode, days_to_last_followup, days_to_death, vital_status,
#'         year_of_form_completion, month_of_form_completion, day_of_form_completion) %>% arrange(bcr_patient_barcode)
#' 
#' #
#' #
#' # Original follow-up information from TCGAbiolinks
#'
#' gdc$follow.up %>% filter(bcr_patient_barcode %in% problem.barcodes) %>%
#'  select(bcr_patient_barcode, days_to_last_followup, days_to_death, vital_status,
#'         year_of_form_completion, month_of_form_completion, day_of_form_completion) %>% arrange(bcr_patient_barcode)
update.survival.from.followup <- function(clinical, follow.up) {
  
  # Replace empty vital status to NA
  levels(follow.up$vital_status) <- levels(follow.up$vital_status) %>% {replace(., . == '', NA)}
  follow.up$vital_status %>% {(all(levels(.) %in% c('Alive', 'Dead')))}
  
  #
  # Build up follow-up information
  f.up.short <- follow.up %>% 
    #
    # Keep individual if it has at least one of the following columns
    filter(!is.na(vital_status) | !is.na(days_to_last_followup) | !is.na(days_to_death)) %>%
    #
    # Group by patient barcode, so that operations are performed per barcode
    group_by(bcr_patient_barcode) %>%
    #
    # Change vital status to 1: Dead 0: Alive
    mutate(vital_status = (vital_status == 'Dead') * 1) %>%
    #
    # Build date for form completion
    mutate(date_form_completion = as.Date(paste0(year_of_form_completion, 
                                                 month_of_form_completion, 
                                                 day_of_form_completion, collapse = ''),
                                          format = '%Y%M%d')) %>%
    #
    # Keep only:
    #  * highest follow-up date
    #  * biggest vital_status (no comming back from dead)
    #  * highest value for days_to_death
    #  * date of form completion (calculated field)
    summarise(days_to_last_followup = max(days_to_last_followup, na.rm = TRUE),
              days_to_death         = max(days_to_death, na.rm = TRUE),
              vital_status          = max(vital_status),
              follow.up.new         = 1,
              date_form_completion  = max(date_form_completion, na.rm = TRUE)) %>%
    #
    # Perform following mutate operations per row
    rowwise %>%
    #
    # Keep only highest value from days to death/follow-up
    mutate(surv_event_time = suppressWarnings(max(days_to_death, days_to_last_followup, na.rm = TRUE))) %>%
    #
    # Replace infinite values by NA (these come from max(NA) = -Inf)
    mutate(surv_event_time = replace(surv_event_time, is.infinite(surv_event_time), NA)) %>%
    #
    # Keep only some columns
    select(bcr_patient_barcode, vital_status, surv_event_time, follow.up.new, date_form_completion)
    
  
  #
  # Prepare clinical data to merge
  #   * vital status 0 or 1 (see above)
  #   * date of form completion (built from day, month and year values)
  clinical.up <- clinical %>% 
    #
    # perform row operations
    rowwise %>%
    #
    # prepare vital status and form completion date
    mutate(vital_status         = (vital_status == 'Dead') * 1,
           date_form_completion = as.Date(paste0(year_of_form_completion, 
                                                 month_of_form_completion, 
                                                 day_of_form_completion, collapse = ''),
                                          format = '%Y%M%d'))
  
  #
  # Merge two tables
  clinical.new <- left_join(clinical.up, 
            f.up.short, 
            by = 'bcr_patient_barcode', 
            suffix = c('__clinical', '__followup')) %>%
    #
    # keep only interesting fields to survival
    select(bcr_patient_barcode, 
           vital_status__clinical,         vital_status__followup, 
           surv_event_time__clinical,      surv_event_time__followup, 
           date_form_completion__clinical, date_form_completion__followup,
           follow.up.new) %>% as.tbl %>%
    # Next mutate operations are performed by row
    rowwise %>%
    #
    # Build:
    #   * vital status with highest value between clinical and follow-up data (again, can't come back from dead)
    #   * same with event time, only biggest one should matter
    #   * calculate days between form_completion
    mutate(vital_status       = suppressWarnings(max(vital_status__clinical, vital_status__followup, na.rm = TRUE)),
           surv_event_time    = suppressWarnings(max(surv_event_time__clinical, surv_event_time__followup, na.rm = TRUE)),
           days_between_forms = as.numeric(date_form_completion__followup - date_form_completion__clinical)) %>%
    #
    # replace infinite values by NA on days_between_form_completion
    mutate(surv_event_time    = replace(surv_event_time, is.infinite(surv_event_time), NA),
           vital_status       = replace(vital_status, is.infinite(vital_status), NA),
           days_between_forms = replace(days_between_forms, is.infinite(days_between_forms), NA)) %>%
    #
    # Select columns
    select(bcr_patient_barcode, 
           vital_status, 
           surv_event_time) %>%
    #
    # sort by patient barcode
    arrange(bcr_patient_barcode)
   
  return(clinical.new)
}