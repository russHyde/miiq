###############################################################################
#
###############################################################################

#' import_gpl
#'
#' Reads in the microarray annotation data stored in a local copy of a GPL
#' dataset (as downloaded from NCBI::GEO).
#'
#' @param        gpl_file      Filename for a GPL dataset (either a *.soft or a
#'   annot.gz file).
#' @param        gpl_dir       Dirname corresponding to the GPL dataset in
#'   \code{gpl_file}.
#'
#' @importFrom   GEOquery      parseGEO
#'
import_gpl <- function(
                       gpl_file,
                       gpl_dir) {
  stopifnot(dir.exists(gpl_dir))

  gpl_filepath <- file.path(gpl_dir, gpl_file)
  stopifnot(file.exists(gpl_filepath))

  GEOquery::parseGEO(gpl_filepath)
}

###############################################################################

#' add_gpl_to_eset
#'
#' Replaces the featureData slot of an ExpressionSet with the corresponding
#' rows of a "GEOquery::GPL" object.
#'
#' @param        eset          An \code{ExpressionSet}
#' @param        gpl_data      A \code{GPL} from \code{GEOquery} containing
#'   annotation data for the \code{ExpressionSet} \code{eset}
#'
#' @importClassesFrom   Biobase   ExpressionSet   AnnotatedDataFrame
#'
#' @importFrom   Biobase       featureData<-
#' @importFrom   GEOquery      Columns   Table
#' @importFrom   methods       is
#'
add_gpl_to_eset <- function(
                            eset,
                            gpl_data) {
  stopifnot(methods::is(eset, "ExpressionSet"))
  stopifnot(methods::is(gpl_data, "GPL"))

  gpl_table <- GEOquery::Table(gpl_data)
  stopifnot(all(rownames(eset) %in% gpl_table[, "ID"]))

  annot_ord <- match(
    rownames(eset),
    gpl_table[, "ID"]
  )

  fd <- new(
    "AnnotatedDataFrame",
    data = gpl_table[annot_ord, ],
    varMetadata = GEOquery::Columns(gpl_data)
  )
  rownames(fd) <- rownames(eset)

  Biobase::featureData(eset) <- fd
  eset
}

###############################################################################

#' import_geo_processed
#'
#' @param        mic           A MicroarrayImportConfig object
#'
#' @importFrom   GEOquery      getGEO
#'
#' @export
#'
import_geo_processed <- function(
                                 mic) {
  stopifnot(is(mic, "MicroarrayImportConfig"))

  path <- mic@processed_dir
  gse_matrix <- file.path(path, mic@processed_files)

  stopifnot(
    dir.exists(path) &&
      file.exists(gse_matrix)
  )

  if (is.na(mic@annot_gpl)) {
    stop("annot_gpl should be TRUE/FALSE in geo_processed_import")
  }

  GEOquery::getGEO(
    filename = gse_matrix,
    destdir = path,
    getGPL = TRUE,
    GSEMatrix = TRUE,
    AnnotGPL = mic@annot_gpl
  )
}

###############################################################################

# TODO: rewrite to use info stored in MicroarrayImportConfig
# TODO: should path/archive/proc_path/gse_matrix be checked for validity as
# files/dirs? - or should this be done in setValidity for mic?

#' import_geo_affy_raw
#'
#' @inheritParams   import_geo_processed
#'
#' @param        tar           Shell command for tar/untar-ing files. Default
#' is 'tar' since I can't use absolute paths to /bin/tar.
#'
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @importFrom   Biobase       pData
#' @importFrom   methods       is
#' @importFrom   oligo         read.celfiles
#' @importFrom   oligoClasses   list.celfiles
#' @importFrom   stringr       str_extract
#' @importFrom   utils         untar
#'
#' @export
#'
import_geo_affy_raw <- function(
                                mic,
                                tar = "tar") {
  stopifnot(methods::is(mic, "MicroarrayImportConfig"))

  path <- mic@raw_dir

  # ? Can archive be more than one entry?
  archive <- file.path(path, mic@raw_archive)

  stopifnot(
    dir.exists(path) &&
      file.exists(archive)
  )

  if (is.na(mic@annot_gpl)) {
    stop("annot_gpl should be TRUE/FALSE in geo_affy_raw_import")
  }

  # untar the archive
  utils::untar(archive, exdir = path, tar = tar)

  # import all the cel files
  cels <- oligoClasses::list.celfiles(
    path, listGzipped = TRUE
  )

  gsms <- stringr::str_extract(cels, "GSM\\d+")

  eset_raw <- oligo::read.celfiles(
    filenames = file.path(path, cels),
    sampleNames = gsms
  )

  # TODO: import the sample annotations from the .soft file
  # TODO: drop the GPL/AnnotGPL stuff - annotate the genes at a later stage
  eset_processed <- import_geo_processed(mic)

  # Add phenotype data from the processed to the raw dataset
  sample_ord <- match(gsms, colnames(eset_processed))
  Biobase::pData(eset_raw) <- Biobase::pData(eset_processed)[sample_ord, ]

  eset_raw
}

###############################################################################

#' import_aryx_affy_raw
#'
#' TODO
#'
#' @param        mic           A MicroarrayImportConfig
#'
#' @export
#'
import_aryx_affy_raw <- function(
                                 mic) {
  ArrayExpress::ae2bioc(
    mageFiles = list(
      rawFiles = mic@raw_files,
      sdrf = mic@sdrf,
      idf = mic@idf,
      adf = mic@adf,
      path = mic@raw_dir
    )
  )
}

###############################################################################
