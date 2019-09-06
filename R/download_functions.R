###############################################################################
#
# Functions for downloading / importing datasets from GEO and ArrayExpress
# - For meta_analysis subproject of drug_markers project
#
###############################################################################

### ======================================================================= ###
#   Helpers: id parsers and validators
### ======================================================================= ###

.build_acc_validator <- function(
                                 regex = c(
                                   "^GPL", "^GSE",
                                   "^2[[:digit:]]{8}$",
                                   "^E-MTAB-[[:digit:]]{4}$"
                                 )) {
  regex <- match.arg(regex)
  function(x) {
    length(x) == 1 && grepl(regex, x)
  }
}

###############################################################################

#' Checks if a vector contains a single GSE accession number (GSE12345).
#'
#' @param        x             A variable for checking to see if it is a GEO
#'   id of an appropriate type.
#' @noRd
#'

.is_gse_acc <- .build_acc_validator("^GSE")

#' Checks if a vector contains a single GPL accession number (GPL12345).
#'
#' @inheritParams   .is_gse_acc
#' @noRd
#'

.is_gpl_acc <- .build_acc_validator("^GPL")

#' Checks that a vector contains a single GDS-UID (200012345)
#'
#' Allows IDs to be specified as strings, factors or numbers. Returns false if
#' more than one uid is passed in, or if the uid is not of GDS format
#' (200012345, ie, 9 characters, starting with 2).
#'
#' @inheritParams   .is_gse_acc
#' @noRd
#'

.is_gds_uid <- .build_acc_validator("^2[[:digit:]]{8}$")

#' Checks that a vector contains a single ArrayExpress UID
#'
#' Only considers IDs of the form "E-MTAB-2893"
#'
#' @param        x             A variable that should be checked to see if it
#'   is an Array-Express accession number
#' @noRd
#'

.is_aryx_acc <- .build_acc_validator("^E-MTAB-[[:digit:]]{4}$")

###############################################################################

#' Constructs the filename for a GPL dataset
#'
#' @param        gpl_acc       String. A single GPL accession number, as used
#'   at Gene-Expression Omnibus. Of the form GPL12345.
#' @param        annot_gpl     Logical. Should the function return the URL for
#'   the most recent (NCBI-constructed) annotations for the platform, or the
#'   user/vendor-submitted 'soft' file.
#' @param        for_url       Logical. Is the required filename for use when
#'   searching the NCBI ftp site, or when storing the data locally?
#'
#'

.make_gpl_filename <- function(
                               gpl_acc,
                               annot_gpl = TRUE,
                               for_url = FALSE) {
  # *.soft files are <ID>_family.soft.gz at NCBI, but stored as <ID>.soft
  #   locally.
  # *.annot.gz files are <ID>.annot.gz at NCBI and when stored locally.

  stopifnot(.is_gpl_acc(gpl_acc))
  stopifnot(length(annot_gpl) == 1 && is.logical(annot_gpl))
  stopifnot(length(for_url) == 1 && is.logical(for_url))

  if (annot_gpl) {
    paste0(gpl_acc, ".annot.gz")
  } else if (for_url) {
    paste0(gpl_acc, "_family.soft.gz")
  } else {
    paste0(gpl_acc, ".soft")
  }
}

###############################################################################

#' Converts a GPL accession number into the FTP URL for that annotation
#' dataset.
#'
#' @inheritParams   .make_gpl_filename
#'
#' @param        ftp_prefix    Prefix for the URL.
#'
#' @importFrom   stringr       str_replace
#'

.make_gpl_url <- function(
                          gpl_acc,
                          annot_gpl = TRUE,
                          ftp_prefix = paste0(
                            "ftp://ftp.ncbi.nlm.nih.gov/",
                            "geo/platforms"
                          )) {
  # Validity testing performed in .make_gpl_filename
  gpl_filename <- .make_gpl_filename(gpl_acc, annot_gpl, for_url = TRUE)
  gpl_nnn <- str_replace(gpl_acc, "[[:digit:]]{1,3}$", "nnn")

  infix <- if (annot_gpl) {
    "annot"
  } else {
    "soft"
  }

  paste(ftp_prefix, gpl_nnn, gpl_acc, infix, gpl_filename, sep = "/")
}

###############################################################################

#' Extracts GDS, GSE, and GPL ids/accessions from an NCBI-eSummary document
#'
#' Takes an XML string from NCBI::eUtils::eSummary over the Gene-Expression
#' Omnibus dataset. Extracts the GPL-accession (GPL987), GSE accession
#' (GSE12345) and the GDS-uid (200012345) for each GSE-dataset that is present
#' in the XML string. If there are no DocSum entries in the XML string, the
#' function exits. If any of the datasets have more than one entry for their
#' uid, accession or GPL value, the function exits.
#'
#' @param        xml_string    An XML string, as returned by httr::GET or
#'   RCurl::getURL. For the function to work, the XML should have arisen from
#'   eUtils search of GEO using eSummary and should contain a DocSum node for
#'   each dataset in the search.
#'
#' @importFrom   dplyr         bind_rows
#' @importFrom   magrittr      %>%
#' @importFrom   tibble        as_tibble
#' @importFrom   XML           getNodeSet   xmlParse   xmlValue
#' @noRd
#'

.parse_gds_to_gpl <- function(
                              xml_string) {
  # Checks validity for all XML DocSums - is a single GPL id mentioned?
  # Parses the GPL uid out of the DocSum for each GPL id
  # Converts to GPL accessions (format: GPL98765)
  # Returns GDS uid to GPL-accession mapping in a tibble

  doc <- XML::xmlParse(xml_string)
  nodes <- XML::getNodeSet(doc, "//DocSum")

  if (length(nodes) == 0) {
    stop("No dataset nodes in xml_string")
  }

  # Extract the Id, and values for Accession and GPL from each node
  predicates <- list(
    uid = "./Id",
    gse = "./Item[@Name = 'Accession']",
    gpl = "./Item[@Name = 'GPL']"
  )

  # Prefix the GPL value with "GPL"
  modifiers <- list(
    uid = identity,
    gse = identity,
    gpl = function(x) paste0("GPL", x)
  )

  # Make a data-frame that contains the GDS-uid, the GSE accession and the
  # GPL accession for each dataset; then bind the dataframes all together
  dfs <- lapply(
    nodes,
    function(n) {
      .extract_gse_info_from_xml <- function(
                                                   pred,
                                                   mod) {
        vals <- vapply(XML::getNodeSet(n, pred), XML::xmlValue, character(1))
        if (length(vals) != 1) {
          stop(paste0("Multiple values for predicate: ", pred))
        }
        mod(vals)
      }

      Map(
        .extract_gse_info_from_xml, predicates, modifiers
      ) %>%
        tibble::as_tibble()
    }
  )

  dplyr::bind_rows(dfs)
}

###############################################################################

#' Obtains the GPL accession(s) for a set of input GDS uids from NCBI-eUtils
#'
#' Uses NCBI-eUtils. Takes GEO dataset UIDs (of the format 200012345) and
#' creates a mapping from these UIDs to the (default) GEO platform accession
#' number used by GEO to annotate the GEO dataset. GPL accessions will be of
#' the form "GPL1234" and GDS uids will be returned in the same for as input
#' (ie, 200012345). Does not make an effort to ensure input order is the same
#' as output order.
#'
#' @param        ids           ABC
#' @param        eutils_prefix   DEF
#'
#' @importFrom   RCurl         getURL
#' @importFrom   tibble        tibble
#'
#' @export
#'
dl_gds_to_gpl <- function(
                          ids,
                          eutils_prefix = paste0(
                            "https://eutils.ncbi.nlm.nih.gov/",
                            "entrez/eutils/esummary.fcgi?db=gds"
                          )) {
  # gds uids are of the form "200012345", "200031365"
  # ie, all have 9 characters and start with 200...
  # gds uid "200012345" corresponds to GSE accession "GSE12345"

  if (length(ids) == 0) {
    return(
      tibble::tibble(
        uid = character(0),
        gse = character(0),
        gpl = character(0)
      )
    )
  }

  # Die if any invalid GDS uids are passed in
  stopifnot(
    all(vapply(ids, .is_gds_uid, logical(1)))
  )

  # Obtain the .xml data from esummary for all input GDS uids.
  # - Note that the xml-url does not pass RCurl::url.exists since it is a
  # query URL; so I'm not sure how best to ensure the URL is legal / safe
  # before calling getURL.
  gds_query <- paste0(eutils_prefix, "&id=", paste(ids, collapse = ","))
  xml_string <- RCurl::getURL(gds_query)

  # .parse_gds_to_gpl is used to check validity of the XML document, and to
  # reformat it into an appropriate format for the user
  .parse_gds_to_gpl(xml_string)
}

### ======================================================================= ###
#   Downloaders: Functions
### ======================================================================= ###

#' Download GPL annotations from GEO
#'
#' TODO: docs.
#' Note that invalid GPL id kills the function.
#' Note that invalid combination of annot.gpl and gpl.acc kill the function.
#' Wrapper around getGEO for GPLs that checks that an GPL dataset is available
#'   before attempting to download it.
#'
#' @inheritParams   .make_gpl_url
#'
#' @param        dest_dir      Directory in which the GPL annotation dataset
#'   should be saved. Default is \code{"./data/ext/GEOquery"} and this
#'   directory should exist before calling this function.
#'
#' @importFrom   GEOquery      getGEO
#' @importFrom   RCurl         url.exists
#'
dl_gpl_annotations <- function(
                               gpl_acc,
                               dest_dir = file.path(
                                 "data", "ext", "GEOquery"
                               ),
                               annot_gpl = TRUE) {
  # TODO: ? Return a list containing the dir and filename for the downloaded
  # / existing annotations.
  #
  stopifnot(.is_gpl_acc(gpl_acc))
  stopifnot(dir.exists(dest_dir))
  stopifnot(
    length(annot_gpl) == 1 &&
      is.logical(annot_gpl)
  )

  # The filepath into which getGEO would store the GPL dataset
  gpl_fname <- .make_gpl_filename(gpl_acc, annot_gpl, for_url = FALSE)
  gpl_filepath <- file.path(dest_dir, gpl_fname)

  if (file.exists(gpl_filepath)) {
    message(sprintf("GPL filename %s already exists", gpl_filepath))
  } else {
    gpl_url <- .make_gpl_url(gpl_acc, annot_gpl)

    if (!RCurl::url.exists(gpl_url)) {
      stop(
        sprintf(
          "GPL URL (%s) does not exist for GPL id: %s", gpl_url, gpl_acc
        )
      )
    }

    gpl <- GEOquery::getGEO(
      GEO = gpl_acc,
      destdir = dest_dir,
      AnnotGPL = annot_gpl
    )
  }

  stopifnot(file.exists(gpl_filepath))

  list(
    gpl_dir = dest_dir,
    gpl_files = gpl_fname,
    annot_gpl = annot_gpl
  )
}

###############################################################################

#' Download processed data from GEO for a GSE accession number
#'
#' Given a GSE accession number (GSE12345), this function downloads the
#' processed dataset from the Gene Expression Omnibus as a matrix. Processed
#' data are stored as \code{dest_dir/<gse_acc>_series_matrix.txt.gz} and any
#' annotation data is stored in \code{dest_dir/}.
#'
#' @param        acc           A single GSE accession number of the form
#'   "GSE12345".
#' @param        dest_dir      An existing directory. The function dies if this
#'   directory is missing. Data related to gse-id are stored either in this
#'   directory (the processed data matrix and the annotation data) or in a
#'   subdirectory <gse.id> of this directory (for the raw data).
#'
#' @return       A list giving the directories and filenames for all downloaded
#'   files. The list contains some entries that are only relevant to
#'   array-express downloads, rather than GEO downloads.
#'
#' @export
#'
dl_geo_processed <- function(
                             acc,
                             dest_dir = file.path("data", "ext", "GEOquery")) {
  stopifnot(dir.exists(dest_dir))
  stopifnot(.is_gse_acc(acc))

  gse_processed_fname <- paste0(acc, "_series_matrix.txt.gz")
  gse_processed_path <- file.path(dest_dir, gse_processed_fname)

  # Check to see if the processed dataset has already been downloaded for
  # this GSE accession, and then downloads it if not.
  if (file.exists(gse_processed_path)) {
    message(sprintf(
      "GSE dataset %s already exists", gse_processed_path
    ))
  } else {
    gq_processed <- GEOquery::getGEO(
      GEO = acc,
      destdir = dest_dir,
      getGPL = FALSE
    )
  }

  # Make list of files based on the downloaded stuff
  list(
    processed_dir = dest_dir,
    processed_files = gse_processed_fname
  )
}

###############################################################################

#' Download raw data from GEO for a GSE accession number
#'
#' Given a GSE identifier, this function downloads the raw data for that
#' dataset from the Gene Expression Omnibus. Raw data is stored in
#' <dest_dir>/<gse_acc>/*_RAW.tar.
#'
#' The function checks whether the raw data has already been
#' downloaded before pulling it from GEO again.
#'
#' @inheritParams   dl_geo_processed
#'
#' @return       A list giving the directories and filenames for all downloaded
#'   files. The list contains some names that are relevant to array-express
#'   downloads, rather than GEO downloads.
#'
#' @importFrom   GEOquery      getGEO   getGEOSuppFiles
#'
#' @export
#'
dl_geo_raw <- function(
                       acc,
                       dest_dir = file.path("data", "ext", "GEOquery")) {
  # TODO: Higher level calling function that decides how to download
  # datasets, based on config options.
  # TODO: Check if a 'samples' file is downloaded

  # ? should annot.gpl be an argument
  # ? should the function determine the default GPL id for a GSE, and return
  # that in the list
  #
  # - We store the matrix files and platform annotations in <dest_dir>/
  # - We store the GEO-dataset-specific .cel files etc in <dest_dir>/<id>/
  # - This is so that the same plaform annotation data is not downloaded
  # multiple times

  stopifnot(dir.exists(dest_dir))
  stopifnot(.is_gse_acc(acc))

  # Setting up GEO-specific subdir of extdata.path is done by getGEOSuppFiles
  gse_destdir <- file.path(dest_dir, acc)
  gse_raw_path <- file.path(gse_destdir, paste0(acc, "_RAW.tar"))

  # Note that processed-data must also be downloaded, so that phenotypic data
  # can be parsed out of it
  gse_processed_fname <- paste0(acc, "_series_matrix.txt.gz")
  gse_processed_path <- file.path(dest_dir, gse_processed_fname)

  # Checks to see if the .tar of the raw data has already been downloaded for
  # this dataset, and then downloads it if not.
  if (dir.exists(gse_destdir) && file.exists(gse_raw_path)) {
    message(sprintf("GSE dataset %s already exists", gse_raw_path))
  } else {
    # Download a tar of the CEL files
    gq_files <- GEOquery::getGEOSuppFiles(
      GEO = acc,
      makeDirectory = TRUE,
      baseDir = dest_dir
    )

    gp_files <- dl_geo_processed(
      acc = acc,
      dest_dir = dest_dir
    )
  }

  stopifnot(file.exists(gse_raw_path))
  stopifnot(file.exists(gse_processed_path))

  # Make list of files based on the downloaded stuff
  list(
    raw_dir = gse_destdir,
    raw_archive = dir(gse_destdir, pattern = "RAW"),
    processed_dir = dest_dir,
    processed_files = gse_processed_fname
  )
}

###############################################################################

#' Download raw data from ArrayExpress
#'
#' @param        acc           An Array-Express accession number, eg,
#'   "E-MTAB-2893"
#' @param        dest_dir      The parent directory in which all ArrayExpress
#'   datasets are downloaded to. Note that the datasets for \code{acc} all get
#'   stored into \code{`<dest_dir>/<acc>`}
#'
#' @importFrom   ArrayExpress   getAE
#' @importFrom   magrittr      %>%
#' @importFrom   purrr         discard
#'
#' @include      list_manipulation.R
#' @export
#'
dl_aryx_raw <- function(
                        acc,
                        dest_dir = tempdir()) {
  stopifnot(dir.exists(dest_dir))
  stopifnot(.is_aryx_acc(acc))

  dataset_path <- file.path(dest_dir, acc)
  if (!dir.exists(dataset_path)) {
    dir.create(dataset_path)
  }

  ae <- ArrayExpress::getAE(
    accession = acc,
    path = dataset_path,
    type = "raw"
  ) %>%
    purrr::discard(is.null)

  list(
    raw_dir = dataset_path,
    raw_files = ae[["rawFiles"]],
    raw_archive = ae[["rawArchive"]],
    processed_dir = ae[["path"]],
    processed_files = get_or_else(
      ae, "processedFiles", as.character(NA)
    ),
    processed_archive = get_or_else(
      ae, "processedArchive", as.character(NA)
    ),
    sdrf = ae[["sdrf"]],
    idf = ae[["idf"]],
    adf = ae[["adf"]]
  )
}

###############################################################################

# Delete downloaded datasets

###############################################################################
