make_state <- function() {
  rv <- reactiveValues(
    app_state = "idle",   # idle | loading | ready | error

    # project
    referenceGenomeChoice = 2,
    projFolderFull = NULL,

    # inputs / design
    files = NULL,
    colData = NULL,
    multiple_groups = 0, # 0 or 1, set to 1 if more than 2 groups in design

    # stage flags
    fastqc_finished = 0,
    bbduk_finished  = 0,
    salmon_finished = 0,

    # core results
    txi = NULL,
    txi_tpms = NULL,
    txi_deseq = NULL,
    res_txi_deseq = NULL,
    res_DEGs_txi_deseq = NULL,
    txi_deseq_deseq = NULL,
    GO_result = NULL,
    vst_data = NULL
  )
  rv
}
