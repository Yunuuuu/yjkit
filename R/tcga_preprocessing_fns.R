#' remove duplicated samples in TCGA
#'
#' remove duplicated samples in TCGA based on Firehose principle
#'
#' @param barcode a character vector gives barcode of TCGA
#' @details In many instances there is more than one aliquot for a given
#'   combination of individual, platform, and data type. However, only one
#'   aliquot may be ingested into Firehose. Therefore, a set of precedence rules
#'   are applied to select the most scientifically advantageous one among them.
#'
#'   The following precedence rules are applied when the aliquots have differing
#'   analytes. For RNA aliquots, T analytes are dropped in preference to H and R
#'   analytes, since T is the inferior extraction protocol. If H and R are
#'   encountered, H is the chosen analyte. This is somewhat arbitrary and
#'   subject to change, since it is not clear at present whether H or R is the
#'   better protocol. If there are multiple aliquots associated with the chosen
#'   RNA analyte, the aliquot with the later plate number is chosen. For DNA
#'   aliquots, D analytes (native DNA) are preferred over G, W, or X
#'   (whole-genome amplified) analytes, unless the G, W, or X analyte sample has
#'   a higher plate number.
#' @return a tibble with duplicated barcode removed
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @references
#'   \href{http://gdac.broadinstitute.org/runs/stddata__2014_02_15/samples_report/LAML_Replicate_Samples.html}{Replicate Samples - GDAC Firehose}
#' @export
tcga_remove_duplicated_samples <- function(barcode){

  barcode_tibble <- tibble::tibble( barcode = barcode )

  # analyte(H, R, T)
  # plate and portion with highest lexicographical sort value
  barcode_tibble <- dplyr::mutate(
    barcode_tibble,
    bcr_patient_barcode = stringr::str_sub(.data$barcode, 1, 12),
    sample_barcode = stringr::str_sub(.data$barcode, 1, 15),
    sample_vial_barcode = stringr::str_sub(.data$barcode, 1, 16),
    analyte = stringr::str_sub(.data$barcode, 20, 20),
    plate = stringr::str_sub(.data$barcode, 22, 25),
    portion = stringr::str_sub(.data$barcode, 18, 19)
  ) %>%
    dplyr::group_by(.data$sample_barcode) %>%
    dplyr::mutate(
      order = order(.data$analyte, dplyr::desc(.data$plate),
                    dplyr::desc(.data$portion),
                    na.last = TRUE, decreasing = FALSE)
    ) %>%
    dplyr::add_count(
      .data$analyte, .data$plate, .data$portion,
      sort = FALSE, name = "n_smp"
    )

  if (any(barcode_tibble$n_smp > 1)) {
    warning("There remains duplicated samples after applying Firehose TCGA replicated samples preprocessing criteria", call. = FALSE)
  }

  remove_duplicated_tibble <- barcode_tibble %>%
    dplyr::slice_min(order_by = order, n = 1, with_ties = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(c(
      "barcode", "bcr_patient_barcode",
      "sample_barcode", "sample_vial_barcode"
    )))

  message(
    stringr::str_c(
      "A total of",
      nrow(barcode_tibble) -  nrow(remove_duplicated_tibble),
      "samples has been removed",
      sep = " "),
    appendLF = TRUE
  )

  remove_duplicated_tibble
}


# download TCGA Clinical data ---------------------------------------------

#' download TCGA clinical data from indexed data
#' @param project TCGA project ID (see \code{TCGAbiolinks::getGDCprojects()})
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @export
tcga_get_cli_indexed <- function(project) {

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)){
    stop("TCGAbiolinks needed for this function to work. Please install it",
         call. = FALSE)
  }

  cli <- TCGAbiolinks::GDCquery_clinic(
    project, type = "clinical"
  )
  bio <- TCGAbiolinks::GDCquery_clinic(
    project, type = "Biospecimen"
  )
  bio <- dplyr::mutate(
    bio, sample_barcode = stringr::str_sub(.data$submitter_id, 1, 15)
  ) %>%
    dplyr::mutate(
      bcr_patient_barcode = stringr::str_sub(
        .data$submitter_id, 1, 12
      )
    ) %>%
    dplyr::rename(sample_vial_barcode = dplyr::all_of("submitter_id"))

  replicated <- colnames(bio) %in% stringr::str_subset(
    colnames(cli),
    stringr::fixed("bcr_patient_barcode"),
    negate = TRUE
  )

  colnames(bio)[replicated] <- stringr::str_c(
    colnames(bio)[replicated], "_bio", sep = ""
  )

  traits <- dplyr::full_join(
    cli, bio, by = "bcr_patient_barcode"
  ) %>%
    dplyr::select(dplyr::all_of(c("bcr_patient_barcode", "sample_barcode")),
                  dplyr::everything()) %>%
    tibble::as_tibble()

}

#' download TCGA clinical data from xml files ------------------------------
#' @param project TCGA project ID (see \code{TCGAbiolinks::getGDCprojects()})
#' @param path Directory/Folder where the data was downloaded. Default:
#'   \code{here::here("rawdata", "GDCdata")}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @export
tcga_get_cli_xml <- function(project, path = here::here("rawdata", "GDCdata")) {

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks needed for this function to work. Please install it",
         call. = FALSE)
  }

  # Clinical data --------------------------------------------------------

  cli_query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Clinical",
    file.type = "xml"
  )

  TCGAbiolinks::GDCdownload(
    cli_query,
    directory = path,
    files.per.chunk = 8
  )

  clinical <- TCGAbiolinks::GDCprepare_clinic(
    cli_query,
    clinical.info = "patient",
    directory = path
  )

  for (i in c("admin", "radiation", "follow_up", "drug", "new_tumor_event")) {
    message(i, appendLF = TRUE)
    cli_aux <- TCGAbiolinks::GDCprepare_clinic(
      cli_query,
      clinical.info = i,
      directory = path
    )

    if (is.null(cli_aux) || nrow(cli_aux) == 0) next

    # add suffix manually if it already exists
    cli_replicated <- colnames(cli_aux) %in% stringr::str_subset(
      colnames(clinical),
      stringr::fixed("bcr_patient_barcode"),
      negate = TRUE
    )

    colnames(cli_aux)[cli_replicated] <- stringr::str_c(
      colnames(cli_aux)[cli_replicated], "_", i
    )

    # merge clinical data by bcr_patient_barcode
    if (!is.null(cli_aux)) {
      clinical <- dplyr::full_join(
        clinical, cli_aux,
        by = "bcr_patient_barcode"
      )
    }
  }
  clinical <- unique(clinical) %>% tibble::as_tibble()

  # Biospecimen -----------------------------------------------------------

  bio_query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Biospecimen",
    file.type = "xml"
  )

  TCGAbiolinks::GDCdownload(
    bio_query,
    directory = path,
    files.per.chunk = 8
  )

  biospecimen <- TCGAbiolinks::GDCprepare_clinic(
    bio_query,
    clinical.info = "sample",
    directory = path
  )
  biospecimen <- dplyr::select(
    biospecimen,
    dplyr::all_of(c("bcr_patient_barcode", "bcr_sample_barcode", "sample_type"))
  ) %>%
    dplyr::mutate(
      sample_barcode = stringr::str_sub(.data$bcr_sample_barcode, 1, 15)
    ) %>%
    dplyr::rename(sample_vial_barcode = dplyr::all_of("bcr_sample_barcode"))

  biospecimen <- unique(biospecimen) %>% tibble::as_tibble()


  # combine clinical data with biospecimen data ---------------------------

  res <- dplyr::full_join(
    clinical, biospecimen, by = "bcr_patient_barcode"
  ) %>%
    dplyr::select(
      dplyr::all_of(c("bcr_patient_barcode",
                      "sample_barcode",
                      "sample_vial_barcode")),
      dplyr::everything()
    )
}

#' download TCGA clinical data from Biotab files ---------------------------
#' @inheritParams tcga_get_cli_xml
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @export
tcga_get_cli_biotab <- function(project, path = here::here("rawdata", "GDCdata")) {

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks needed for this function to work. Please install it",
         call. = FALSE)
  }

  cli_query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
  )

  TCGAbiolinks::GDCdownload(
    cli_query,
    directory = path,
    files.per.chunk = 8
  )

  clinical.BCRtab.all <- TCGAbiolinks::GDCprepare(
    cli_query,
    directory = path
  )
  clinical <- clinical.BCRtab.all[[ stringr::str_c(
    "clinical_patient_",
    tolower(stringr::str_sub(project, -4, -1)),
    sep = ""
  )]]

  for (
    i in stringr::str_subset(
      names(clinical.BCRtab.all),
      "clinical_patient",
      negate = TRUE
    )
  ){
    message(i, appendLF = TRUE)
    cli_aux <- clinical.BCRtab.all[[i]]
    cli_aux_test <- dplyr::filter(
      cli_aux,
      !(.data$bcr_patient_uuid %in% c("bcr_patient_uuid", "CDE_ID:"))
    )
    if (is.null(cli_aux_test) || nrow(cli_aux_test) == 0) next

    # add suffix manually if it already exists
    cli_replicated <- colnames(cli_aux) %in% stringr::str_subset(
      colnames(clinical),
      "bcr_patient_barcode|bcr_patient_uuid",
      negate = TRUE
    )

    colnames(cli_aux)[cli_replicated] <- stringr::str_c(colnames(cli_aux)[cli_replicated], "_", i, sep = "")

    # merge clinical data by bcr_patient_barcode
    if (!is.null(cli_aux)) {
      clinical <- dplyr::full_join(
        clinical, cli_aux,
        by = c("bcr_patient_uuid", "bcr_patient_barcode")
      )
    }
  }
  clinical <- unique(clinical) %>% tibble::as_tibble()


  # Biospecimen data ------------------------------------------------------

  bio_query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Biospecimen",
    data.type = "Biospecimen Supplement",
    data.format = "BCR Biotab"
  )

  TCGAbiolinks::GDCdownload(
    bio_query,
    directory = path,
    files.per.chunk = 8
  )

  biospecimen.BCRtab.all <- TCGAbiolinks::GDCprepare(
    bio_query, directory = path
  )

  biospecimen <- biospecimen.BCRtab.all[[stringr::str_c(
    "biospecimen_sample_",
    tolower(stringr::str_sub(project, -4, -1)),
    sep = ""
  )]]

  biospecimen <- dplyr::select(
    biospecimen,
    dplyr::all_of(c(
      "bcr_patient_uuid", "bcr_sample_barcode", "sample_type"
    ))
  ) %>%
    dplyr::mutate(
      sample_barcode = stringr::str_sub(.data$bcr_sample_barcode, 1, 15)
    ) %>%
    dplyr::rename(sample_vial_barcode = dplyr::all_of("bcr_sample_barcode"))

  biospecimen <- unique(biospecimen) %>% tibble::as_tibble()


  # Combine clinical data with biospecimen data ---------------------------

  res <- dplyr::full_join(
    clinical, biospecimen, by = "bcr_patient_uuid"
  ) %>%
    dplyr::select(
      dplyr::all_of(c("bcr_patient_uuid",
                      "bcr_patient_barcode",
                      "sample_barcode",
                      "sample_vial_barcode")),
      dplyr::everything()
    )
}



