#' @export
require(clusterProfiler)
mp_go <- function(
    gl,
    ont = "All",
    pvalueCutoff = 0.1,
    TERM2GENE = read.table("~/wkdir/reference/Mpolymorpha/TERM2GENE.txt", header = TRUE),
    TERM2NAME = read.table("~/wkdir/reference/Mpolymorpha/TERM2NAME.txt", header = TRUE)) {
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)

    if (ont == "BP") {
        sub_terms <- names(goterms[goterms == "BP"])
    } else if (ont == "MF") {
        sub_terms <- names(goterms[goterms == "MF"])
    } else if (ont == "CC") {
        sub_terms <- names(goterms[goterms == "CC"])
    } else {
        sub_terms <- names(goterms)
    }

    TERM2GENE <- TERM2GENE[TERM2GENE$term %in% sub_terms, ]
    TERM2NAME <- TERM2NAME[TERM2NAME$term %in% sub_terms, ]

    res <- clusterProfiler::enricher(
        gene = gl,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        pvalueCutoff = 0.1,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2
    )

    res <- as.data.frame(res)
    if (nrow(res) != 0) {
        res$ont <- ont
        return(res)
    } else {
        return(NULL)
    }
}
