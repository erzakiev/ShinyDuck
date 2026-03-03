make_state <- function() {
  rv <- reactiveValues(
    app_state = "idle",   # idle | loading | ready | error
    multiplegroups = 0,   # 0 or 1, set to 1 if more than 2 groups in design
    # project
    referenceGenomeChoice = 2,
    projFolderFull = NULL,

    # inputs / design
    files = NULL,
    colData = NULL,
    multiple_groups = 0,

    # stage flags
    fastqc_finished = 0,
    bbduk_finished  = 0,
    salmon_finished = 0,

    # core results
    tpms_tbl = NULL,
    deseq = NULL,
    deg_tbl = NULL,
    txi = NULL,
    txi_tpms = NULL,
    res_txi_deseq = NULL,
    res_DEGs_txi_deseq = NULL,
    txi_deseq_deseq = NULL

  )
  rv
}
