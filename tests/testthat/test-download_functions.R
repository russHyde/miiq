###############################################################################

context("Tests for downloading of GEO and ArrayExpress Datasets")

###############################################################################

test_that(".is_gpl_acc", {
  expect_true(
    .is_gpl_acc("GPL12345"),
    info = "Genuine GPL id"
  )
  expect_false(
    .is_gpl_acc(c("GPL12345", "GPL987676")),
    info = ".is_gpl_acc only works for singleton GPL ids"
  )
  expect_false(
    .is_gpl_acc("NOT A GPL"),
    info = "Invalid GPL id"
  )
})

###############################################################################

test_that(".make_gpl_filename", {
  expect_equal(
    object = .make_gpl_filename(
      "GPL12345",
      annot_gpl = TRUE,
      for_url = TRUE
    ),
    expected = "GPL12345.annot.gz",
    info = "Filename at NCBI for GPL12345 annot_gpl"
  )
  expect_equal(
    object = .make_gpl_filename(
      "GPL12345",
      annot_gpl = FALSE,
      for_url = TRUE
    ),
    expected = "GPL12345_family.soft.gz",
    info = "Filename at NCBI for GPL12345 standard GPL"
  )
  expect_equal(
    object = .make_gpl_filename(
      "GPL12345",
      annot_gpl = TRUE,
      for_url = FALSE
    ),
    expected = "GPL12345.annot.gz",
    info = "Local filename for GPL12345 annot_gpl"
  )
  expect_equal(
    object = .make_gpl_filename(
      "GPL12345",
      annot_gpl = FALSE,
      for_url = FALSE
    ),
    expected = "GPL12345.soft",
    info = "Local filename for GPL12345 standard GPL"
  )
})

###############################################################################

test_that(".make_gpl_url", {
  ftp_platform_prefix <- "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/"

  expect <- paste0(
    ftp_platform_prefix,
    "GPL12nnn/GPL12345/annot/GPL12345.annot.gz"
  )
  expect_equal(
    object = .make_gpl_url(
      "GPL12345",
      annot_gpl = TRUE
    ),
    expected = expect,
    info = "annot_gpl url for a GPL"
  )

  expect <- paste0(
    ftp_platform_prefix,
    "GPL12nnn/GPL12345/soft/GPL12345_family.soft.gz"
  )
  expect_equal(
    object = .make_gpl_url(
      "GPL12345",
      annot_gpl = FALSE
    ),
    expected = expect,
    info = ".soft url for a GPL"
  )

  expect <- paste0(
    ftp_platform_prefix,
    "GPLnnn/GPL12/soft/GPL12_family.soft.gz"
  )
  expect_equal(
    object = .make_gpl_url(
      "GPL12",
      annot_gpl = FALSE
    ),
    expected = expect,
    info = ".soft url for a GPL with less than 3 digits"
  )
})

###############################################################################

test_that(".is_gds_uid", {
  expect_false(
    object = .is_gds_uid("NOT A GDS"),
    info = "Doesn't start with '2'"
  )
  expect_false(
    object = .is_gds_uid("200"),
    info = "Has less than 9 characters"
  )
  expect_false(
    object = .is_gds_uid("2001234560"),
    info = "Has more than 9 characters"
  )
  expect_false(
    object = .is_gds_uid("20012345x"),
    info = "Contains a non-digit"
  )
  expect_true(
    object = .is_gds_uid(200123456),
    info = "GDS uid as a number"
  )
  expect_true(
    object = .is_gds_uid(factor("200123456")),
    info = "GDS uid as a number"
  )
  # Note that .is_gds_uid returns false if multiple entries are passed in
  expect_false(
    object = .is_gds_uid(c("200012345", "200098765")),
    info = "Only currently allows 1 UID to be passed at once"
  )
})

###############################################################################

#' @importFrom   tibble        tibble
#'
test_that(".parse_gds_to_gpl", {

  # No entries
  xml <- '<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN"
"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">
<eSummaryResult>
<ERROR>UID=200999999: cannot get document summary</ERROR>

</eSummaryResult>'
  expect_error(
    object = .parse_gds_to_gpl(xml),
    info = "Valid GDS xml, but there's no entries inside it"
  )

  # Single entry
  xml <- '<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN"
"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">
<eSummaryResult>
<DocSum>
    <Id>200031365</Id>
    <Item Name="Accession" Type="String">GSE31365</Item>
    <Item Name="GPL" Type="String">6244</Item>
    <Item Name="GSE" Type="String">31365</Item>
</DocSum>
</eSummaryResult>'
  expect <- tibble::tibble(
    uid = "200031365",
    gse = "GSE31365",
    gpl = "GPL6244"
  )
  expect_equal(
    object = .parse_gds_to_gpl(xml),
    expected = expect,
    info = "Valid GDS xml, and valid map from GSE to GPL for a single
 dataset"
  )

  # Two entries
  xml <- '<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN"
"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">
<eSummaryResult>
<DocSum>
    <Id>200031365</Id>
    <Item Name="Accession" Type="String">GSE31365</Item>
    <Item Name="GPL" Type="String">6244</Item>
    <Item Name="GSE" Type="String">31365</Item>
</DocSum>
<DocSum>
    <Id>200011578</Id>
    <Item Name="Accession" Type="String">GSE11578</Item>
    <Item Name="GPL" Type="String">571</Item>
    <Item Name="GSE" Type="String">11578</Item>
</DocSum>
</eSummaryResult>'
  expect <- tibble::tibble(
    uid = c("200031365", "200011578"),
    gse = c("GSE31365", "GSE11578"),
    gpl = c("GPL6244", "GPL571")
  )
  expect_equal(
    object = .parse_gds_to_gpl(xml),
    expected = expect,
    info = "Valid GDS xml, and valid map from GSE to GPL for a two
datasets"
  )
})

###############################################################################
test_that("dl_gds_to_gpl", {
  expect_equal(
    object = dl_gds_to_gpl(character(0)),
    expected = tibble::tibble(
      uid = character(0), gse = character(0), gpl = character(0)
    ),
    info = "No GEO ids in the input"
  )

  expect_error(
    object = dl_gds_to_gpl("NOT A GDS"),
    info = "Invalid GDS uid passed to .dl_gds_to_gpl"
  )

  expect_equal(
    object = dl_gds_to_gpl("200011578"),
    expected = tibble::tibble(
      uid = "200011578",
      gse = "GSE11578",
      gpl = "GPL571"
    ),
    info = "Single GDS uid, mapped to a GPL id"
  )

  expect_equal(
    object = dl_gds_to_gpl(c("200011578", "200031365")),
    expected = tibble::tibble(
      uid = c("200011578", "200031365"),
      gse = c("GSE11578", "GSE31365"),
      gpl = c("GPL571", "GPL6244")
    ),
    info = "Single GDS uid, mapped to a GPL id"
  )
})

###############################################################################
