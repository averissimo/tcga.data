
#' Title
#'
#' @param clinical 
#' @param follow.up 
#'
#' @return
#' @export
#'
#' @examples
update.survival.from.followup <- function(clinical, follow.up) {
  
  # Need to add this to main Rmd!
  levels(follow.up$vital_status) <- levels(follow.up$vital_status) %>% {replace(., . == '', NA)}
  
  follow.up$vital_status %>% {(all(levels(.) %in% c('Alive', 'Dead')))}
  
  f.up.short <- follow.up %>% 
    filter(!is.na(vital_status) | !is.na(days_to_last_followup) | !is.na(days_to_death)) %>%
    group_by(bcr_patient_barcode) %>%
    mutate(vital_status = (vital_status == 'Dead') * 1) %>%
    mutate(date_form_completion = as.Date(paste0(year_of_form_completion, 
                                                 month_of_form_completion, 
                                                 day_of_form_completion, collapse = ''),
                                          format = '%Y%M%d')) %>%
    summarise(days_to_last_followup = max(days_to_last_followup, na.rm = TRUE),
              days_to_death         = max(days_to_death, na.rm = TRUE),
              vital_status          = max(vital_status),
              follow.up.new         = 1,
              date_form_completion  = max(date_form_completion, na.rm = TRUE)) %>%
    rowwise %>%
    mutate(surv_event_time = suppressWarnings(max(days_to_death, days_to_last_followup, na.rm = TRUE))) %>%
    mutate(surv_event_time = replace(surv_event_time, is.infinite(surv_event_time), NA)) %>%
    select(bcr_patient_barcode, vital_status, surv_event_time, follow.up.new, date_form_completion)
    
  
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
  
  clinical.new <- left_join(clinical.up, 
            f.up.short, 
            by = 'bcr_patient_barcode', 
            suffix = c('__clinical', '__followup')) %>%
    select(bcr_patient_barcode, 
           vital_status__clinical,
           vital_status__followup, 
           surv_event_time__clinical,
           surv_event_time__followup, 
           date_form_completion__clinical,
           date_form_completion__followup,
           follow.up.new) %>% as.tbl %>%
    rowwise %>%
    mutate(vital_status       = suppressWarnings(max(vital_status__clinical, vital_status__followup, na.rm = TRUE)),
           surv_event_time    = suppressWarnings(max(surv_event_time__clinical, surv_event_time__followup, na.rm = TRUE)),
           days_between_forms = as.numeric(date_form_completion__followup - date_form_completion__clinical)) %>%
    mutate(days_between_forms = replace(days_between_forms, is.infinite(days_between_forms), NA)) %>%
    #
    # Select columns
    select(bcr_patient_barcode, 
           vital_status, 
           surv_event_time) %>%
    mutate(surv_event_time = replace(surv_event_time, is.infinite(surv_event_time), NA),
           vital_status    = replace(vital_status, is.infinite(vital_status), NA)) %>%
    arrange(bcr_patient_barcode)
  
  # clinical.new %>% filter(bcr_patient_barcode %in% c('TCGA-FR-A726', 'TCGA-D3-A8GP', 'TCGA-DA-A1HW', 'TCGA-DA-A1I1', 
  #                                                   'TCGA-DA-A1I5', 'TCGA-DA-A1I7', 'TCGA-DA-A1IB', 'TCGA-DA-A95W', 
  #                                                   'TCGA-DA-A95X', 'TCGA-EB-A41A', 'TCGA-EE-A2GJ', 'TCGA-FR-A726', 
  #                                                   'TCGA-HR-A5NC', 'TCGA-XV-A9VZ'))
   
  return(clinical.new)
}