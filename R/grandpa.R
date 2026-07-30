#' @importFrom igraph vcount ecount is_igraph vertex_attr_names vertex_attr V
#' @importFrom igraph as_edgelist make_empty_graph set_vertex_attr is_directed
#' @importFrom igraph degree add_edges are_adjacent as_undirected cluster_louvain
#' @importFrom igraph membership degree_distribution
#' @importFrom dplyr count mutate arrange desc select group_by summarise across
#' @importFrom dplyr contains everything rowwise filter %>%
#' @importFrom tidyr uncount expand_grid
#' @importFrom stats approxfun
#' @importFrom igraph vcount ecount is_igraph vertex_attr_names vertex_attr V "V<-"
#' @importFrom utils globalVariables
NULL

# Register non-standard-evaluation column names referenced inside
# dplyr/tidyr pipelines so R CMD check does not flag them as undefined
# global variables.
utils::globalVariables(c("p", "n", "addToBag", "ix", "V1", "V2", "areConnected"))

# ============================================================================
# PATCH NOTES (vs original GRANDPA from CarlyBobak/grandpa)
# ----------------------------------------------------------------------------
# Bug fixes:
#   * Line 49 (orig): `CMSSub` -> `G`. Was a leftover dev reference.
#   * Line 197 (orig): explicit column-1 access in 1-label branch.
#   * Line 254/255 (orig) in sampleEdges and 156/157 in sampleVertices:
#     row-count after uncount() can be off by +/- 1 from `mt`/`nt` due to
#     rounding. Now padded/trimmed to exactly mt/nt before index assignment.
#   * Line 315 (orig): parenthesis bug
#       `length(C2V$map[[lab2]]==1)` -> `length(C2V$map[[lab2]]) == 1`.
#   * Line 371 (orig): undefined `resolution`; now an explicit argument
#     with a default of 1.
#   * Line 417 (orig): `nrow(df)==F` -> `nrow(df) == 0`.
#   * createMap/sampleVertices: force `ix` to integer to prevent character
#     coercion through dplyr pipelines that would break degree() lookups.
#   * addAugment: coerce factor labels to character immediately to avoid
#     downstream cbind/paste/count surprises.
#
# Deprecation fixes:
#   * `list.vertex.attributes()` -> `vertex_attr_names()` (lines 98, 158).
#   * `count(across())` -> `count(across(everything()))` (lines 112, 120, 124).
#   * `across(!c(ix))` -> `across(-ix)` (line 224).
#   * `as.undirected()` -> `as_undirected()` (line 369).
#   * `are.connected()` -> `are_adjacent()` (lines 315, 323, 415).
#
# Style / safety:
#   * `T`/`F` -> `TRUE`/`FALSE` throughout function defaults and predicates.
#   * `1:n` -> `seq_len(n)` to avoid the n=0 edge case.
#   * `summarise()` gains `.groups = "drop"` to silence dplyr's grouping
#     warning.
#
# Functional changes: none. The output of grandpa() with identical seeds
# and inputs should be identical to the original (modulo cases where the
# original would have errored or warned about deprecations).
# ============================================================================


#' Conducts a GRANDPA procedure
#'
#' Uses the GRANDPA framework to generate a network with matching attribute
#' probabilities.
#'
#' @param G an igraph object containing the original graph. Attributes that
#'   will be used in the modeling procedure must contain the word 'Label'
#'   (case-insensitive) in their attribute name.
#' @param nt the number of desired nodes in the generated graph.
#' @param mt the number of desired edges in the generated graph.
#' @param preventSelf logical; should self connections be prevented?
#' @param preventDups logical; should duplicate connections be prevented?
#' @param augment logical; should degree augmentation occur?
#' @param augmentCommunity logical; should communities be augmented?
#' @param nbins the number of degree bins used in the degree augmentation.
#' @param degType string; the type of degree to be used in the degree
#'   augmentation. See ?degree for details.
#' @param seed a random seed set for reproducibility.
#'
#' @return an igraph object with attribute and augmentation labels.
#' @export
grandpa <- function(G, nt = vcount(G), mt = ecount(G),
                    preventSelf = TRUE, preventDups = TRUE,
                    augment = FALSE, augmentCommunity = FALSE,
                    nbins = 3, degType = "all", seed = 0) {

  if (!is_igraph(G)) {
    stop("Graph is not an igraph object")
  }

  # Add augmentation label
  if (augment) {
    G <- addAugment(G, nbins, degType)
  }

  if (augmentCommunity) {
    G <- addCommunityAugment(G)
  }

  # Check that a vertex attribute named "label" exists.
  if (!any(grepl("label", vertex_attr_names(G), ignore.case = TRUE))) {
    stop("At least one vertex attribute must contain the word 'label'")
  }

  # Calculate observed probabilities of attributes for vertices and edges.
  obsProb <- calcObsProbs(G)
  observedLabels     <- obsProb[[1]]
  observedEdgeLabels <- obsProb[[2]]

  # Sample vertices and edges.
  sampled_v <- sampleVertices(observedLabels, nt, seed)
  Gt        <- initializeTargetG(sampled_v, nt, d = is_directed(G))
  C2V       <- createMap(sampled_v)
  sampled_e <- sampleEdges(observedEdgeLabels, mt, seed)

  # Build the graph.
  Gt <- simGraph(Gt, sampled_e, C2V, seed, preventSelf, preventDups)

  return(Gt)
}


#' Calculates the observed edge and vertex probabilities
#'
#' This function is typically called by the grandpa function rather than
#' directly.
#'
#' @param G An igraph object containing vertex attributes whose names
#'   include the substring 'Label' (case-insensitive).
#'
#' @return a list of tidy tables containing the vertex probabilities and
#'   edge probabilities respectively.
#' @export
calcObsProbs <- function(G) {

  # Identify label-bearing vertex attributes.
  attr_names <- vertex_attr_names(G)
  Labels <- attr_names[grepl("label", attr_names, ignore.case = TRUE)]
  if (length(Labels) == 0) {
    stop("No vertex attributes containing 'label' were found.")
  }

  # Collect all label values into a single data frame.
  observedLabels <- vertex_attr(G, Labels[1])
  if (length(Labels) > 1) {
    for (i in 2:length(Labels)) {
      observedLabels <- cbind(observedLabels, vertex_attr(G, Labels[i]))
    }
  }
  observedLabels <- data.frame(observedLabels, stringsAsFactors = FALSE)
  colnames(observedLabels) <- Labels
  rownames(observedLabels) <- V(G)$name

  observedLabels_summary <- observedLabels %>%
    count(across(everything())) %>%
    mutate(p = n / sum(n)) %>%
    arrange(desc(p))

  # Edge probabilities.
  Es <- as_edgelist(G, names = FALSE)

  if (length(Labels) == 1) {
    observedEdgeLabels <- data.frame(observedLabels[Es[, 1], ],
                                     observedLabels[Es[, 2], ],
                                     stringsAsFactors = FALSE) %>%
      count(across(everything())) %>%
      mutate(p = n / sum(n)) %>%
      arrange(desc(p))
  } else {
    composite <- do.call(paste, c(observedLabels, sep = ""))
    observedEdgeLabels <- data.frame(composite[Es[, 1]],
                                     composite[Es[, 2]],
                                     stringsAsFactors = FALSE) %>%
      count(across(everything())) %>%
      mutate(p = n / sum(n)) %>%
      arrange(desc(p))
  }

  return(list(observedLabels_summary, observedEdgeLabels))
}


#' Creates the vertices that will be used in the generated graph
#'
#' @param observedLabels a data frame containing observed label
#'   combinations and their probabilities from the original graph.
#' @param nt the number of desired nodes in the generated graph.
#' @param seed the random seed used for reproducibility.
#'
#' @return a data frame containing node attributes.
#' @export
sampleVertices <- function(observedLabels, nt, seed) {
  set.seed(seed)
  observedLabels$addToBag <- round(observedLabels$p * nt, 0)
  sampled_v <- observedLabels %>%
    select(contains("label", ignore.case = TRUE), addToBag) %>%
    uncount(addToBag)

  # uncount() can produce a row count slightly off from nt due to rounding
  # of `addToBag`. Pad or trim so the result has exactly nt rows.
  if (nrow(sampled_v) > nt) {
    sampled_v <- sampled_v[seq_len(nt), , drop = FALSE]
  } else if (nrow(sampled_v) < nt) {
    pad_n  <- nt - nrow(sampled_v)
    pad_ix <- sample(seq_len(nrow(sampled_v)), pad_n, replace = TRUE)
    sampled_v <- rbind(sampled_v, sampled_v[pad_ix, , drop = FALSE])
  }
  sampled_v <- as.data.frame(sampled_v[sample(seq_len(nt)), , drop = FALSE])
  colnames(sampled_v) <- colnames(observedLabels)[
    grepl("label", colnames(observedLabels), ignore.case = TRUE)
  ]
  sampled_v$ix <- as.integer(seq_len(nt))
  return(sampled_v)
}


#' Initializes an empty graph with the generated vertices
#'
#' @param sampled_v data frame of sampled vertices and their labels.
#' @param nt number of desired nodes.
#' @param d logical; is the graph directed?
#'
#' @return an empty igraph object with the desired vertices.
#' @export
initializeTargetG <- function(sampled_v, nt, d) {
  Gt <- make_empty_graph(nt, directed = d)

  label_cols <- colnames(sampled_v)[-ncol(sampled_v)]  # drop `ix`
  for (cn in label_cols) {
    Gt <- set_vertex_attr(Gt, cn, seq_len(nt),
                          value = as.character(sampled_v[, cn]))
  }

  if (length(label_cols) > 1) {
    V(Gt)$SummaryLabel <- do.call(
      paste,
      c(sampled_v[, label_cols, drop = FALSE], sep = "")
    )
  } else {
    Gt <- set_vertex_attr(Gt, "SummaryLabel", seq_len(nt),
                          value = as.character(sampled_v[, label_cols[1]]))
  }
  return(Gt)
}


#' Creates the map of label combinations to vertex IDs
#'
#' @param sampled_v data frame of sampled vertices with an `ix` column.
#' @return a list mapping label combinations to integer vertex IDs.
#' @export
createMap <- function(sampled_v) {
  # Ensure `ix` is integer.
  sampled_v$ix <- as.integer(sampled_v$ix)

  C2V <- suppressMessages(
    sampled_v %>%
      group_by(across(-ix)) %>%
      summarise(map = list(ix), .groups = "drop")
  )
  names(C2V$map) <- do.call(paste, c(C2V[, -ncol(C2V)], sep = ""))
  return(C2V)
}


#' Creates the edges that will be used in the generated graph
#'
#' @param observedEdgeLabels data frame of edge label combinations and
#'   their probabilities.
#' @param mt number of desired edges.
#' @param seed random seed.
#'
#' @return a data frame of sampled edges.
#' @export
sampleEdges <- function(observedEdgeLabels, mt, seed) {
  set.seed(seed)
  observedEdgeLabels$addToBag <- round(observedEdgeLabels$p * mt, 0)
  sampled_e <- observedEdgeLabels %>%
    select(1:2, addToBag) %>%
    uncount(addToBag)

  # Force row count to exactly mt.
  if (nrow(sampled_e) > mt) {
    sampled_e <- sampled_e[seq_len(mt), , drop = FALSE]
  } else if (nrow(sampled_e) < mt) {
    pad_n  <- mt - nrow(sampled_e)
    pad_ix <- sample(seq_len(nrow(sampled_e)), pad_n, replace = TRUE)
    sampled_e <- rbind(sampled_e, sampled_e[pad_ix, , drop = FALSE])
  }
  sampled_e <- sampled_e[sample(seq_len(mt)), , drop = FALSE]
  sampled_e$ix <- seq_len(mt)
  colnames(sampled_e) <- c(paste("Label", 1:2, sep = ""), "ix")
  return(sampled_e)
}


#' Procedure to add edges to the generated graph
#'
#' @param Gt the empty graph.
#' @param sampled_e data frame of sampled edges with label columns.
#' @param C2V mapping of label combinations to vertex IDs.
#' @param seed random seed.
#' @param preventSelf logical; prevent self connections?
#' @param preventDups logical; prevent duplicate connections?
#'
#' @return a generated igraph object.
#' @export
simGraph <- function(Gt, sampled_e, C2V, seed, preventSelf, preventDups) {
  set.seed(seed)
  first <- TRUE

  for (i in seq_len(nrow(sampled_e))) {
    lab1 <- as.character(sampled_e[i, 1])
    lab2 <- as.character(sampled_e[i, 2])

    if (length(C2V$map[[lab1]]) > 1) {
      v1 <- sample(C2V$map[[lab1]], 1,
                   prob = 1 / (degree(Gt, v = C2V$map[[lab1]]) + 1e-4))
    } else {
      v1 <- C2V$map[[lab1]]
    }

    if (length(C2V$map[[lab2]]) > 1) {
      v2 <- sample(C2V$map[[lab2]], 1,
                   prob = 1 / (degree(Gt, v = C2V$map[[lab2]]) + 1e-4))
    } else {
      v2 <- C2V$map[[lab2]]
    }

    if (preventSelf && v1 == v2) {
      v_ix <- preventSelfConnections(v1, v2, lab1, lab2, C2V)
      v1 <- v_ix[1]
      v2 <- v_ix[2]
    }

    # Avoid drawing the same edge twice when both labels resolve to a
    # single vertex each.
    if (length(C2V$map[[lab1]]) == 1 && length(C2V$map[[lab2]]) == 1 &&
        are_adjacent(Gt, v1, v2) && preventDups) {
      if (first) {
        warning("Could not prevent duplicate connections. ",
                "Consider simplifying to unique edges")
        first <- FALSE
      }
      next
    }

    if (preventDups && are_adjacent(Gt, v1, v2)) {
      v_ix <- preventDuplicateConnections(v1, v2, lab1, lab2, C2V, Gt)
      v1 <- v_ix[1]
      v2 <- v_ix[2]
    }

    Gt <- add_edges(Gt, c(v1, v2))
  }
  return(Gt)
}


#' Add degree augmentation label
#'
#' @param G an igraph object.
#' @param nbins the number of degree bins.
#' @param degType string; type of degree (see ?degree).
#'
#' @return an igraph object with `augLabel` vertex attribute.
#' @export
addAugment <- function(G, nbins, degType) {
  deg <- degree(G, mode = degType)
  # Coerce factor to character immediately to avoid downstream
  # cbind/paste/count surprises with factor levels.
  V(G)$augLabel <- as.character(ggplot2::cut_interval(deg, nbins))
  return(G)
}


#' Add community augmentation label
#'
#' This function is typically called by the `grandpa` function rather than
#' directly. By default it implements Louvain clustering to detect
#' community membership. For best results we recommend fitting a community
#' algorithm outside of GRANDPA and attaching the result as the
#' `CommunityLabel` vertex attribute prior to calling `grandpa()`.
#'
#' @param G an igraph object.
#' @param resolution Louvain resolution parameter (default 1).
#' @param seed random seed for reproducibility (default NULL = no seed).
#'
#' @return an igraph object with `CommunityLabel` vertex attribute.
#' @export
addCommunityAugment <- function(G, resolution = 1, seed = NULL) {
  if (is_directed(G)) {
    warning("Communities are detected on an undirected graph")
    tempG <- as_undirected(G, mode = "collapse")
  } else {
    tempG <- G
  }
  if (!is.null(seed)) set.seed(seed)
  Gclust <- cluster_louvain(tempG, resolution = resolution)
  V(G)$CommunityLabel <- as.character(membership(Gclust))
  return(G)
}


#' Prevents self connections during graph generation
#'
#' @param v1 a vertex id.
#' @param v2 a second vertex id.
#' @param lab1 the label corresponding to vertex id 1.
#' @param lab2 the label corresponding to vertex id 2.
#' @param C2V mapping of label categories to vertex IDs.
#'
#' @return updated vertex IDs.
#' @export
preventSelfConnections <- function(v1, v2, lab1, lab2, C2V) {
  while (v1 == v2 && length(C2V$map[[lab1]]) > 1) {
    v1 <- sample(C2V$map[[lab1]], 1)
  }
  while (v1 == v2 && length(C2V$map[[lab2]]) > 1) {
    v2 <- sample(C2V$map[[lab2]], 1)
  }
  if (v1 == v2) {
    warning("Could not prevent self connection")
  }
  return(c(v1, v2))
}


#' Prevents duplicate connections during graph generation
#'
#' @param v1 a vertex id.
#' @param v2 a second vertex id.
#' @param lab1 the label corresponding to vertex id 1.
#' @param lab2 the label corresponding to vertex id 2.
#' @param C2V mapping of label categories to vertex IDs.
#' @param Gt the graph being generated.
#'
#' @return updated vertex IDs.
#' @export
preventDuplicateConnections <- function(v1, v2, lab1, lab2, C2V, Gt) {
  df <- expand_grid(C2V$map[[lab1]], C2V$map[[lab2]])
  colnames(df) <- c("V1", "V2")
  df <- df %>%
    rowwise() %>%
    mutate(areConnected = are_adjacent(Gt, V1, V2)) %>%
    filter(!areConnected)

  if (nrow(df) == 0) {
    warning("Could not prevent duplicate connections. ",
            "Consider simplifying to unique edges")
    return(c(v1, v2))
  } else {
    ix <- sample(seq_len(nrow(df)), 1)
    return(as.numeric(df[ix, 1:2]))
  }
}


#' Calculates the NRMSE between the degree distributions of two graphs
#'
#' @param G the original graph.
#' @param Gt the generated graph.
#'
#' @return numeric; the NRMSE between the degree distributions.
#' @export
degree_error <- function(G, Gt) {
  sourceCCDF <- degree_distribution(G,  cumulative = TRUE)
  targetCCDF <- degree_distribution(Gt, cumulative = TRUE)

  source_f <- approxfun(sourceCCDF)
  target_f <- approxfun(targetCCDF)

  overlap_x <- min(length(sourceCCDF), length(targetCCDF))
  i <- seq(1, overlap_x, length.out = 500)

  h1 <- source_f(i) - target_f(i)
  error <- sqrt(sum(h1^2) / length(h1)) / mean(target_f(i))
  return(error)
}
