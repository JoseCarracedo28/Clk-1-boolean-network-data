###############################################################################
# End-to-end reproducible verification for Reviewer 3 (comments C1, C3, C4)
# clk-1 mutant Boolean network (clk1 = 0), Configurations 1 and 2.
#
# BoolNet pipeline mirroring clk1_async_verification.py. Produces:
#   1. Synchronous attractors (both configs)            BoolNet
#   2. Strategy A: reduced AMPK-ROS-HIF-1 submodule, sync vs async (both configs)
#   3. Strategy B: asynchronous updates from biological initial conditions
#   4. mpbn most-permissive attractors (both configs)   via reticulate -> mpbn
# Figure S7 is produced by the companion Python script.
#
# Verified results:
#   Config 1: sync 16/10/11-state (basin 78/14/8%); async = single complex
#             attractor >12.5M states (cf. Reviewer 1); mpbn = 1 attractor, 24
#             free nodes -> identical partition across the three schemes.
#             Strategy A submodule: no oscillation under either scheme.
#   Config 2: sync = 10 cyclic attractors; mpbn = 1 attractor, only the 4-node
#             UPRmt loop free -> the 10 cycles are synchronous artifacts.
#             Strategy A submodule: oscillates under sync, converges under async
#             (textbook synchronous artifact).
#
# Requires: install.packages("BoolNet"); reticulate + Python "mpbn" (pip install mpbn)
###############################################################################

library(BoolNet)

## ---------------------------------------------------------------------------
## Build a config-specific network file (clk1 = 0 for the mutant).
## Configs differ ONLY in ros, hif1, ampk (the AMPK-ROS-HIF-1 loop).
## ---------------------------------------------------------------------------
build_network <- function(cfg) {
  common <- c(
    "targets, factors",
    "clk1, 0",
    "ETC, clk1",
    "pink1, !ETC",
    "skn1, ros | !ETC",
    "taf4, !clk1",
    "met, !hif1",
    "hlh11, !atfs1 & !mtor",
    "unc51, ampk & !mtor",
    "hlh30, !mtor",
    "creb, taf4 & !crtc1",
    "UP, !clk1 & (ETC2 | !hsp60)",
    "ATG, unc51 & hlh30",
    "lipl4, hlh30",
    "MTG, pink1 & skn1 & unc51",
    "fzo1, creb",
    "clpp1, UP",
    "LDs, ATG",
    "ETC2, fzo1 & MTG",
    "atfs1, clpp1",
    "hsp60, atfs1",
    "atgl1, !hlh11",
    "lip, lipl4 & LDs & atgl1",
    "betaox, lip",
    "ATP, (ETC2 & betaox) | ETC",
    "crtc1, !ampk",
    "mtor, !ampk & !unc51")
  loop <- if (cfg == 1)
    c("ros, ETC2 | (!ampk & !clk1)", "hif1, ros", "ampk, !ATP | (ros & !hif1)")
  else
    c("ros, !ETC & (!ampk | hif1)", "hif1, ros & !ampk", "ampk, ros & !hif1")
  path <- sprintf("clk1_mutant_cfg%d.txt", cfg)
  writeLines(c(common, loop), path)
  loadNetwork(path)
}

## Decode an attractor's states and report oscillating / fixed nodes.
osc_partition <- function(attr, idx) {
  genes <- attr$stateInfo$genes
  nwords <- ceiling(length(genes) / 32)
  st <- matrix(attr$attractors[[idx]]$involvedStates, nrow = nwords)
  mat <- t(apply(st, 2, function(col) {
    bits <- unlist(lapply(col, function(x) as.integer(intToBits(x))))
    bits[seq_along(genes)]
  }))
  colnames(mat) <- genes
  list(osc  = genes[apply(mat, 2, function(c) length(unique(c)) > 1)],
       fix1 = genes[apply(mat, 2, function(c) all(c == 1))],
       fix0 = genes[apply(mat, 2, function(c) all(c == 0))])
}

## ---------------------------------------------------------------------------
## Strategy A: reduced AMPK-ROS-HIF-1 submodule (ETC2, betaox as clamped inputs)
## ---------------------------------------------------------------------------
strategy_A <- function(cfg) {
  cat(sprintf("\n[Strategy A] AMPK-ROS-HIF-1 submodule (cfg %d):\n", cfg))
  ## NOTE: in cfg 2, ros = !ETC & (!ampk | hif1); with ETC = 0 this is (!ampk | hif1).
  loop <- if (cfg == 1)
    c("ampk, !ATP | (ros & !hif1)", "ros, ETC2 | !ampk",
      "hif1, ros", "ATP, ETC2 & betaox")
  else
    c("ampk, ros & !hif1", "ros, !ampk | hif1",
      "hif1, ros & !ampk", "ATP, ETC2 & betaox")
  writeLines(c("targets, factors", "ETC2, ETC2", "betaox, betaox", loop),
             sprintf("ampk_ros_hif1_cfg%d.txt", cfg))
  mod <- loadNetwork(sprintf("ampk_ros_hif1_cfg%d.txt", cfg))
  s <- getAttractors(mod, type = "synchronous",  method = "exhaustive")
  a <- getAttractors(mod, type = "asynchronous", method = "random", startStates = 5000)
  cat("   synchronous:\n");  print(s)
  cat("   asynchronous:\n"); print(a)
  ## Expected: cfg1 -> only fixed points (both schemes); cfg2 -> synchronous cycle
  ## that disappears under asynchronous updating (spurious synchronous oscillation).
}

## ---------------------------------------------------------------------------
## mpbn most-permissive attractors via reticulate (order-independent ground truth)
## ---------------------------------------------------------------------------
mpbn_attractors <- function(cfg) {
  library(reticulate)
  mpbn <- import("mpbn")
  rules <- list(
    clk1="0", ETC="clk1", pink1="!ETC", skn1="ros | !ETC", taf4="!clk1",
    met="!hif1", hlh11="!atfs1 & !mtor", unc51="ampk & !mtor", hlh30="!mtor",
    creb="taf4 & !crtc1", UP="!clk1 & (ETC2 | !hsp60)", ATG="unc51 & hlh30",
    lipl4="hlh30", MTG="pink1 & skn1 & unc51", fzo1="creb", clpp1="UP", LDs="ATG",
    ETC2="fzo1 & MTG", atfs1="clpp1", hsp60="atfs1", atgl1="!hlh11",
    lip="lipl4 & LDs & atgl1", betaox="lip", ATP="(ETC2 & betaox) | ETC",
    crtc1="!ampk", mtor="!ampk & !unc51")
  if (cfg == 1)
    rules <- c(rules, list(ros="ETC2 | (!ampk & !clk1)", hif1="ros",
                           ampk="!ATP | (ros & !hif1)"))
  else
    rules <- c(rules, list(ros="!ETC & (!ampk | hif1)", hif1="ros & !ampk",
                           ampk="ros & !hif1"))
  bn <- mpbn$MPBooleanNetwork(rules)
  attrs <- iterate(bn$attractors())
  cat(sprintf("\n[mpbn] cfg %d: %d most-permissive attractor(s)\n", cfg, length(attrs)))
  for (a in attrs) {
    free <- names(a)[sapply(a, function(v) identical(v, "*"))]
    cat("   free (oscillating):", length(free), "->", paste(sort(free), collapse=", "), "\n")
  }
}

## ===========================================================================
for (cfg in 1:2) {
  cat("\n", strrep("=", 70), "\nCONFIGURATION ", cfg, "\n", strrep("=", 70), "\n", sep="")
  net <- build_network(cfg)

  ## 1. Synchronous attractors (random sampling reproduces basin proportions;
  ##    exhaustive is feasible but heavy because clk1 is a constant -> 2^28 states)
  attr_sync <- getAttractors(net, type = "synchronous",
                             method = "random", startStates = 200000)
  cat(sprintf("[Synchronous] %d attractor(s) found\n", length(attr_sync$attractors)))
  nwords <- ceiling(length(net$genes) / 32)
  sizes <- sapply(attr_sync$attractors,
                  function(a) length(a$involvedStates) / nwords)
  for (i in order(-sizes)) {
    p <- osc_partition(attr_sync, i)
    cat(sprintf("   |A|=%d  #oscillating=%d\n", sizes[i], length(p$osc)))
  }

  ## 2. Strategy A
  strategy_A(cfg)

  ## 3. Strategy B: asynchronous updates from biological seeds (dominant attractor)
  dom <- which.max(sizes)
  seed_mat <- matrix(attr_sync$attractors[[dom]]$involvedStates, nrow = nwords)
  seed_list <- lapply(seq_len(ncol(seed_mat)), function(j) {
    bits <- unlist(lapply(seed_mat[, j], function(x) as.integer(intToBits(x))))
    bits[seq_along(net$genes)]
  })
  cat("\n[Strategy B] asynchronous from biological seeds (this can be slow; the\n",
      "  asynchronous attractor is a single very large complex attractor,\n",
      "  >12.5M states for cfg 1, cf. Reviewer 1):\n", sep="")
  attr_async <- getAttractors(net, type = "asynchronous",
                              method = "chosen", startStates = seed_list)
  print(attr_async)

  ## 4. mpbn ground truth
  mpbn_attractors(cfg)
}

cat("\nFor Figure S7 and a fast (R-free) reproduction, run clk1_async_verification.py\n")
