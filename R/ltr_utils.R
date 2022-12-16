add_coordinates_of_closest_neighbor <- function(gff) {
  gff <- gff[order(seqnames(gff), as.vector(start(gff)))]
  # split to chromosomes:
  gff_parts <- split(gff, as.vector(seqnames(gff)))
  upstreams <- c(sapply(gff_parts, function(x) c(1, head(end(x), -1))))
  downstreams <- c(sapply(gff_parts, function(x) c(start(x)[-1], seqlengths(x)[runValue(seqnames(x[1]))])))
  gff_updated <- unlist(gff_parts)
  gff_updated$upstream_domain <- unlist(upstreams)
  gff_updated$downstream_domain <- unlist(downstreams)
  names(gff_updated) <- NULL
  return(gff_updated)
}

get_domain_clusters_alt <- function(gff, dist_models, threshold=0.99){
  # gff <- sort(gff, by = ~ seqnames * start)
  ## it must be already sorted by seqnames and start
  ## it must take into account strand
  strand1 <- head(strand(gff) == "+", -1)
  strand2 <- strand(gff)[-1] == "+"
  domain_pairs <- data.frame(
    D1 = paste0(head(gff$Name,-1),"_S"),
    D2 = paste0(gff$Name[-1],"_S"),
    C1 =  head(gff$Final_Classification,-1),
    C2 =  gff$Final_Classification[-1],
    S1 =  head(strand(gff),-1),
    S2 =  strand(gff)[-1],
    start1 =  ifelse(strand1, head(start(gff),-1), head(end(gff),-1)),
    start2 = ifelse(strand2, start(gff)[-1], end(gff)[-1]),
    chrpos= paste0(seqnames(gff)[-1],"_",start(gff)[-1])
  )
  domain_pairs$D_A <- with(domain_pairs, ifelse(D1 < D2, D1, D2))
  domain_pairs$D_B <- with(domain_pairs, ifelse(D1 > D2, D1, D2))
  domain_pairs$quantile <- 1
  same_element_cluster <- domain_pairs$S1 == domain_pairs$S2 & domain_pairs$C1 == domain_pairs$C2
  domain_pairs$same_element_cluster <- same_element_cluster

  q_plus_function <- apply(domain_pairs[same_element_cluster,],
                                         1,
                                         function(x) dist_models$plus[[x["C1"]]][[x["D_A"]]][[x["D_B"]]])

  D <- abs(domain_pairs[same_element_cluster,]$start1 - domain_pairs[same_element_cluster,]$start2)

  q_plus_value <- mapply(function(x, D)if(is.null(x)){0}else{x(D)}, x= q_plus_function, D = D)

  domain_pairs$quantile[same_element_cluster] <- q_plus_value

  domain_pairs$split_position <- !same_element_cluster | domain_pairs$quantile > threshold
  # combine clusters based on distances and original cluster
  clusters <- paste(cumsum(c(TRUE, domain_pairs$split_position)), get_domain_clusters(gff))
  return(clusters)
}


get_domain_clusters <- function(gff) {
  # get consecutive domains from same linage
  # must be sorted!
  plus_order_split <- c(0, as.numeric(cumsum(head(gff$domain_order, -1) >= gff$domain_order[-1] & strand(gff)[-1] == '+')))
  minus_order_split <- rev(as.numeric(cumsum(rev(c(gff$domain_order[-1] >= head(gff$domain_order,-1) & strand(gff)[-1] == '-', FALSE)))))
  # split based on classification - must be same and consecutive
  x <- rle(gff$Final_Classification)
  # split on strand change
  s <- rep(seq_along(runLength(strand(gff))), runLength(strand(gff)))
  domain_cluster <- paste0(rep(seq_along(x$lengths), x$lengths), "_", seqnames(gff),
                           "_", plus_order_split, "_", minus_order_split, "_", s)
  clusters <- factor(domain_cluster, levels = unique(domain_cluster))
  return(clusters)
}

# create partial TE from clusters
get_partial_te_from_cluster_of_domains <- function (gpart){
  ID <- paste("TE_partial_", gpart$Cluster[1], sep="")
  te_partial <- GRanges(type="transposable_element_partial",
                        strand=strand(gpart)[1],
                        ID = ID,
                        source = "dante_ltr",
                        seqnames=seqnames(gpart)[1],
                        Final_Classification=gpart$Final_Classification[1],
                        Name=gpart$Final_Classification[1],
                        IRanges(min(start(gpart)), max(end(gpart))))
  gpart$Parent <- ID
  return(c(te_partial,gpart))
}

# create partial TE from clusters
get_partial_te_from_cluster_of_domains <- function(gpart, ID) {
  te_partial <- makeGRangesFromDataFrame(
    data.frame(type = "transposable_element_partial",
               strand = gpart$strand[1],
               ID = ID,
               source = "dante_ltr",
               seqnames = gpart$seqnames[1],
               Final_Classification =
                 gpart$Final_Classification[1],
               Name = gpart$Final_Classification[1],
               IRanges(min(gpart$start), max(gpart$end))),
    keep.extra.columns = TRUE)


  gpart$Parent <- ID
  gpart_gr <- makeGRangesFromDataFrame(gpart, keep.extra.columns = TRUE)
  return(c(te_partial, gpart_gr))
}

# function to count for each element number of occurences:
count_occurences_for_each_element <- function(x) {
  counts_unique <- table(x)
  counts <- counts_unique[x]
  counts
}


clean_domain_clusters <- function(gcl, lineage_domain_span, min_domains) {
  ## remove clusters wich does not have enough domains or domains
  ## are on different strand
  N_domains <- sapply(gcl, nrow)
  N_unique_domains <- sapply(gcl, function(x)length(unique(x$Name)))
  S <- sapply(gcl, function(x)paste(sort(unique(x$strand)), collapse = " "))
  S_OK <- S %in% c("+", "-")
  max_span <- lineage_domain_span[sapply(gcl, function(x)  x$Final_Classification[1])]
  # set to zero if lineage is not covered in lineage domain span
  max_span[is.na(max_span)] <- 0
  span <- sapply(gcl, function(x)max(x$end) - min(x$start))
  cond1 <- S_OK &
    N_unique_domains == N_domains &
    N_domains >= min_domains &
    span <= max_span
  return(gcl[cond1])
}

check_ranges <- function(gx, s, offset = OFFSET) {
  # check is range is not out of sequence length
  START <- sapply(gx, function(x)min(x$start)) - offset
  END <- sapply(gx, function(x)max(x$end)) + offset
  MAX <- seqlengths(s)[sapply(gx, function(x)as.character(x$seqnames[1]))]
  good_ranges <- (START > 0) & (END <= MAX)
  return(good_ranges)
}

get_ranges <- function(gx, offset = OFFSET) {
  S <- sapply(gx, function(x)min(x$start))
  E <- sapply(gx, function(x)max(x$end))
  gr <- GRanges(seqnames = sapply(gx, function(x)x$seqnames[1]), IRanges(start = S - offset, end = E + offset))
}

get_ranges_left <- function(gx, offset = OFFSET, offset2 = 300) {
  ## offset2 - how many nt cen LTR extend to closes protein domain
  ## this is necassary as some detected proteins domains does not have correct bopundaries
  ## if LTR retrotransposons insters to other protein domain.
  S <- sapply(gx, function(x)min(x$start))
  max_offset <- S - sapply(gx, function(x)min(x$upstream_domain)) + 10
  offset_adjusted <- ifelse(max_offset < offset, max_offset, offset)
  gr <- GRanges(seqnames = sapply(gx, function(x)x$seqnames[1]), IRanges(start = S - offset_adjusted, end = S + offset2))
  return(gr)
}

get_ranges_right <- function(gx, offset = OFFSET, offset2 = 300) {
  E <- sapply(gx, function(x)max(x$end))
  max_offset <- sapply(gx, function(x)max(x$downstream_domain)) - E + 10
  offset_adjusted <- ifelse(max_offset < offset, max_offset, offset)
  gr <- GRanges(seqnames = sapply(gx, function(x)x$seqnames[1]), IRanges(start = E - offset2, end = E + offset_adjusted))
  return(gr)
}

firstTG <- function(ss) {
  x <- matchPattern("TG", ss)
  if (length(x) == 0) {
    return(0)
  }else {
    return(min(start(x)))
  }
}

lastCA <- function(ss) {
  x <- matchPattern("CA", ss)
  if (length(x) == 0) {
    return(0)
  }else {
    return(max(start(x)))
  }
}

domain_distance <- function(d_query, d_reference){
  if (d_query == d_reference){
    return (0)
  }
  d_query_p <-  strsplit(d_query," ")[[1]]
  d_reference_p <-  strsplit(d_reference," ")[[1]]
  d <- length(d_reference_p) - sum(d_query_p == d_reference_p[d_reference_p %in% d_query_p])
  return(d)
}


trim2TGAC <- function(bl) {
  for (i in 1:nrow(bl)) {
    cons <- consensusString(c(bl$qseq[i], bl$sseq[i]), ambiguityMap="?")
    tg_P <- firstTG(cons)
    ca_P <- lastCA(cons)
    e_dist <- bl$length[i] - ca_P
    max_dist <- 50 # was 25
    # count mismatches in the trimming region (?)
    S_mismatch <- sum(strsplit(substring(cons,1, tg_P),"")[[1]] == "?")
    E_mismatch <- sum(strsplit(substring(cons, ca_P, nchar(cons)),"")[[1]] == "?")
    no_match <- any(tg_P == 0, ca_P == 0)
    if (!no_match &
      tg_P - S_mismatch < max_dist &
      e_dist - E_mismatch < max_dist) {
      ## trim alignment
      bl[i,] <- trim_blast_table(bl[i,], tg_P, e_dist - 1)
    }

  }
  return(bl)
}

trim_blast_table <- function(b, T1, T2) {
  b$qstart <- b$qstart + T1
  b$sstart <- b$sstart + T1
  b$qend <- b$qend - T2
  b$send <- b$send - T2
  b$sseq <- substring(b$sseq, T1, b$length - T2)
  b$qseq <- substring(b$qseq, T1, b$length - T2)
  b$length <- nchar(b$sseq)
  return(b)
}

blast_all2all <- function(seqs, db=NULL, ncpus=1, word_size=20, perc_identity=90, max_target_seq = 5000, task = "megablast", additional_params= ""){
  if (ncpus == 1){
    query <- list(seqs)
  }else{
    query <-split(seqs, round(seq(1,ncpus,length.out = length(seqs))))
    if (length(query) < ncpus){
      ncpus <- length(query)
    }
  }
  if(is.null(db)){
    # search against itself
    db <- seqs
  }
  qf <-tempfile(fileext=paste0("_",1:ncpus,".fasta"))
  outf <-tempfile(fileext=paste0("_",1:ncpus,".csv"))
  dbf <- tempfile()
  script <-  tempfile()
  writeXStringSet(db, dbf)
  mapply(query, qf, FUN = writeXStringSet)
  cols <- "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
  cmd_db <-  paste("makeblastdb -dbtype nucl -in ", dbf)
  cmd_blast <-  paste("blastn -task ", task, " -word_size", word_size,
                    "-outfmt \"6 ", cols, "\" ",
                    "-perc_identity", perc_identity, " -min_raw_gapped_score 500",
                    "-max_target_seqs", max_target_seq, additional_params,
                    "-query", qf, "-db", dbf, "-out", outf,
                    "&"
  )

  # TODO - inspect only forward strand??
  system(cmd_db)
  cmd_all <- paste(paste(cmd_blast,collapse="\n"),"\nwait")
  cat(cmd_all, file = script)
  system(paste("sh ", script))

  bl_list <- lapply(outf, read.table, stringsAsFactors = FALSE, col.names = unlist(strsplit(cols, " ")), sep="\t", comment.char = "")
  bl_table <- do.call(rbind, bl_list)
  unlink(qf)
  #unlink(outf)
  unlink(dbf)
  unlink(script)
  return(bl_table)
}

identify_conflicts <- function (bl){
  QL <- gsub(".+[|]", "", bl$qaccver)
  SL <- gsub(".+[|]", "", bl$saccver)
  id_with_conflict <- unique(c(bl$qaccver[QL != SL],bl$saccver[QL != SL]))
  id_ok <- setdiff(unique(c(bl$qaccver,bl$saccver)), id_with_conflict)
  single_hit <- sapply(
    sapply(
      split(bl$qaccver, bl$saccver), unique
    ), length) == 1
  id_with_no_hits <- names(single_hit)[single_hit] # except hit to itself
  return(list(
    ok = id_ok,
    conflict = id_with_conflict,
    no_hit = id_with_no_hits)
  )
}


analyze_TE <- function(seqs, ncpus = 10, word_size = 20, perc_identity = 90, ...){
  blt <- blast_all2all(seqs, ncpus = ncpus, word_size = word_size, perc_identity = perc_identity, ...)
  te_conflict_info <- identify_conflicts(blt)
  blt_te_ok <- blast_table_subset(blt, te_conflict_info$ok)
  te_ok_lineages <- split(blt_te_ok,
                                   gsub(
                                     ".+[|]",
                                     "",
                                     blt_te_ok$qaccver))
  gr_representative <- GRangesList(mclapply(te_ok_lineages,
                             FUN = get_representative_ranges,
                             mc.cores = ncpus
  ))
  seqs_representative <- getSeq(seqs, Reduce(c, gr_representative))
  names(seqs_representative) <- seqnames(Reduce(c, gr_representative))

  return(
    list(
      seqs_representative = seqs_representative,
      te_conflict_info = te_conflict_info,
      gr_representative = gr_representative,
      blast = blt
    )
  )
}

query_coverage <- function(blt){
  blt <- blt[blt$qaccver != blt$saccver,]
  Q_lengths <-  blt$qlen[!duplicated(blt$qaccver)]
  names(Q_lengths) <- blt$qaccver[!duplicated(blt$qaccver)]
  gr <- GRanges(seqnames = blt$qaccver, seqlengths = Q_lengths,
               IRanges(start = blt$qstart, end = blt$qend))
  return(coverage(gr))
}

multiplicity_of_te <- function(blt){
  # exclude self to self hits and calculate coverage + mean_multiplicity of TE
  # assuption is that TE which are 'identical' to other TE from the same lineage are
  # likely correct
  blt_no_self <- blt[blt$qaccver != blt$saccver,]
  cvr <- query_coverage(blt_no_self)
  L <- sapply(cvr, function(x) sum(width(x)))
  C1 <- sapply(cvr, function(x) sum(as.numeric(runValue(x) >= 1) * runLength(x)))
  multiplicity <- sapply(cvr, function(x) sum(as.numeric(runValue(x)) * runLength(x)))/L
  data.frame(L = L, C1 = C1,  multiplicity = multiplicity )
}

verify_based_on_multiplicity <- function(TE_info, min_coverage=0.99, min_multiplicity=3){
  blt <- TE_info$blast[TE_info$blast$qaccver %in% TE_info$te_conflict_info$ok,]
  mp <- multiplicity_of_te(blt)
  id_ok_mp_verified <- rownames(mp)[mp$C1/mp$L > min_coverage & mp$multiplicity >= min_multiplicity]
  return(list(multiplicity = mp,
              id_ok_mp_verified = id_ok_mp_verified))

}

compare_TE_datasets <- function(Q, S, word_size = 20, min_coverage = 0.95, ncpus=10, ...){
  blt <- blast_all2all(seqs = Q, db = S, ncpus = ncpus, word_size = word_size, ...)
  QL <- gsub(".+[|]", "", blt$qaccver)
  SL <- gsub(".+[|]", "", blt$saccver)
  id_with_conflict <- unique(c(blt$qaccver[QL != SL]))
  id_ok <- setdiff(unique(blt$qaccver), id_with_conflict)
  # check coverage hits
  blt_ok <- blt[blt$qaccver %in% id_ok,]
  Q_lengths <-  blt_ok$qlen[!duplicated(blt_ok$qaccver)]
  names(Q_lengths) <- blt_ok$qaccver[!duplicated(blt_ok$qaccver)]
  gr <- GRanges(seqnames = blt_ok$qaccver, seqlengths = Q_lengths,
               IRanges(start = blt_ok$qstart, end = blt_ok$qend))
  cvr <- coverage(gr)
  L <- sapply(cvr, function(x) sum(width(x)))
  C1 <- sapply(cvr, function(x) sum(as.numeric(runValue(x) >= 1) * runLength(x)))
  Max_uncovered <- sapply(cvr, function(x){
    if(any(runValue(x)==0)){
      return(max(runLength(x)[runValue(x) == 0]))
    }else{
      return(0)
    }
  })

  # verified based on hit to reference - S
  C1_prop <- C1/L
  pass <-  (C1_prop >= min_coverage & Max_uncovered < 500)
  if (any(pass)){
    id_ok_verified <-  names(C1_prop)
  }else {
    id_ok_verified <- NULL
  }
  return(list(id_with_conflict = id_with_conflict,
              id_ok = id_ok,
              id_ok_verified = id_ok_verified
  ))
}



blast_table_subset <- function(bl,id){
  return(bl[bl$qaccver %in% id & bl$saccver %in% id,, drop = FALSE])
}

get_representative_ranges <-  function(bl, min_length = 200, min_identity = 98){
  bl <- bl[bl$pident>=min_identity, , drop=FALSE]
  bl <- bl[bl$pident>=min_identity & bl$length >= min_length, , drop=FALSE]
  score <- sort(unlist(by(bl$bitscore, bl$qaccver, sum, simplify = FALSE)),
               decreasing = TRUE)
  L <-  bl$qlen[!duplicated(bl$qaccver)]
  names(L) <- bl$qaccver[!duplicated(bl$qaccver)]
  gr <- GRanges(seqnames = bl$qaccver,
               IRanges(start = bl$qstart, end = bl$qend),
               seqlengths = L,
               subject = bl$saccver,
               sstart = ifelse(bl$send < bl$sstart, bl$send, bl$sstart),
               send = ifelse(bl$send > bl$sstart, bl$send, bl$sstart))
  SN <-  levels(seqnames(gr))
  inc <- rep(TRUE, length(gr))
  MSK <- GRangesList()
  for (i in names(score)){
    inc[gr$subject %in% i] <- FALSE
    gr_part <- gr[seqnames(gr) %in% i & inc]
    MSK[[i]] <- GRanges(seqnames = factor(gr_part$subject, levels = SN),
                       IRanges(start = gr_part$sstart, end = gr_part$send),
                       seqlengths = L
    )
  }
  gout <- unlist(MSK)

  full_gr <- GRanges(seqnames = factor(SN, levels = SN),
                     IRanges(start = rep(1,length(L)),
                            end = L)
  )
  unmasked_gr <- GenomicRanges::setdiff(full_gr, gout)
  return(unmasked_gr[width(unmasked_gr) >= min_length])
}

expected_diversity <- function(seqs, niter=100, km = 6){
  L <- nchar(seqs)
  R <- matrix(ncol = niter, nrow = length(seqs))
  for (i in 1:niter){
    seqs_rnd <- DNAStringSet(sapply(L, function(n) paste(sample(c("A", "C", "T", "G"), n, replace=TRUE), collapse="")))
    R[,i] <- seq_diversity(seqs_rnd, km = km)$richness
  }
  R

}

seq_diversity <- function (seqs, km=6){
  K <- oligonucleotideFrequency(seqs, width=km)>0
  P <- t(K)/rowSums(K)
  # shannon index
  SI <- apply(P, 2, function(x) {x1 <- x[x>0]; -sum(x1*log(x1))})
  # richness
  R <- rowSums(K)
  list(richness=R, shannon_index=SI)
}

mask_tandem_repeats <- function(fasta_file){
  # use tidehunte to detect tandem repeats
  # this function require tidehunter and bedtools to be installed
  bed_mask <- tempfile()
  fasta_masked <- tempfile()
  system(paste("TideHunter -f 2" , fasta_file, " |cut -f 1,4,5 >", bed_mask))
  system(paste("bedtools maskfasta -fi", fasta_file, "-bed", bed_mask, "-fo", fasta_masked, "-soft"))
  bed <- read.table(bed_mask, header = FALSE, col.names = c("seqname", "start", "end"), sep = "\t")
  unlist(bed_mask)
  return(list(fasta=fasta_masked, bed=bed))
}

overlap_size_perc <- function(interval1, interval2){
  if (interval1[1] > interval2[2] | interval1[2] < interval2[1]){
    return(0)
  }else{
    interval_size <- interval1[2] - interval1[1] + 1
    ovl_size <- (min(interval1[2], interval2[2]) - max(interval1[1], interval2[1]))
    return(ovl_size/interval_size)
  }
}

blast <- function(s1, s2, expected_aln_lenght=NULL, min_identity=70){
  tmp1 <- tempfile()
  tmp2 <- tempfile()
  tmp_out <- tempfile()
  writeXStringSet(DNAStringSet(s1), tmp1)
  writeXStringSet(DNAStringSet(s2), tmp2)
  # masking - test
  masks1 <- mask_tandem_repeats(tmp1)
  masks2 <- mask_tandem_repeats(tmp2)



  # alternative blast:
  cmd <- paste("blastn -task blastn -lcase_masking -word_size 7 -dust no -gapextend 1 -gapopen 2 -reward 1 -penalty -1",
               " -query ", masks1$fasta, ' -subject ', masks2$fasta, ' -strand plus ', '-perc_identity ', min_identity,
               '-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"',
               '  -out', tmp_out)

  system(cmd)
  out_raw <- read.table(tmp_out, as.is = TRUE, sep = "\t",
                        col.names = strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq", split = ' ')[[1]])
  if (nrow(out_raw) == 0) {
    return(out_raw)
  }
  # calculate form bed what percent of the blast hits is masked
  # if more than 90% is masked remove the hit
  out_raw$masked_perc_query <- 0
  out_raw$masked_perc_subject <- 0
  for (bl in 1:nrow(out_raw)){
    int_query <- c(out_raw$qstart[bl], out_raw$qend[bl])
    int_subject <- c(out_raw$sstart[bl], out_raw$send[bl])
    # check query
    q_ovl <- 0
    for (m in seq_along(masks1$bed$seqname)){
      int_mask <- c(masks1$bed$start[m], masks1$bed$end[m])
      q_ovl <-  q_ovl + overlap_size_perc(int_query, int_mask)
    }
    out_raw$masked_perc_query[bl] <- q_ovl
    # check subject
    s_ovl <- 0
    for (m in seq_along(masks2$bed$seqname)){
      int_mask <- c(masks2$bed$start[m], masks2$bed$end[m])
      s_ovl <- s_ovl + overlap_size_perc(int_subject, int_mask)
    }
    out_raw$masked_perc_subject[bl] <- s_ovl
  }

  out_raw <- out_raw[out_raw$masked_perc_query < 0.9 & out_raw$masked_perc_subject < 0.9, , drop = FALSE]

  if (nrow(out_raw) == 0) {
    return(out_raw)
  }
  out <- trim2TGAC(out_raw)
  # remove alingment shorted that
  # out <- out[out$length > 100,]
  # alngment must be at least 80% of expected LTRmin95
  out <- out[out$length > (expected_aln_lenght *0.8),]
  if (nrow(out) == 0) {
    return(out)
  }
  ## filter for TGCA
  TG_L <- substring(out$qseq, 1, 2) == "TG"
  TG_R <- substring(out$sseq, 1, 2) == "TG"
  CA_L <- substring(out$qseq, out$length - 1, out$length) == "CA"
  CA_R <- substring(out$sseq, out$length - 1, out$length) == "CA"
  cond <- TG_L & TG_R & CA_L & CA_R
  out <- out[cond, , drop = FALSE]
  unlink(tmp1)
  unlink(tmp2)
  unlink(tmp_out)
  return(out)
}

get_correct_coordinates <- function(b) {
  do.call(rbind, strsplit(b$qaccver, split = "_"))
}

evaluate_ltr <- function(GR, GRL, GRR, blast_line, Lseq, Rseq) {
  LTR_L <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRL) + blast_line$qstart - 2,
                                               end = start(GRL) + blast_line$qend - 1))
  LTR_R <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRR) + blast_line$sstart - 2,
                                               end = start(GRR) + blast_line$send - 1))

  TSD_L <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRL) + blast_line$qstart - 3:10,
                                               end = start(GRL) + blast_line$qstart - 3))
  TSD_R <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames(GR),
                                               start = start(GRR) + blast_line$send,
                                               end = start(GRR) + blast_line$send + 0:7))

  TSD_L_seq <- DNAStringSet(substring(Lseq, blast_line$qstart - 2:9, blast_line$qstart - 2))
  TSD_R_seq <- DNAStringSet(substring(Rseq, blast_line$send + 1, blast_line$send + 1:8))

  matching_tsd <- TSD_R_seq == TSD_L_seq
  matching_tsd[1:3] <- FALSE # remove short tsd
  p <- which(matching_tsd)
  if (length(p) > 0) {
    TSD_Length <- max(p)
    TSD_sequence <- TSD_L_seq[TSD_Length]
    TSD_position <- append(TSD_R[TSD_Length], TSD_L[TSD_Length])
  }else{
    # check if there is a TSD, allow for 1 mismatch, min lenght must be 5
    N_mismatch <- sapply(mapply("!=", strsplit(as.character(TSD_L_seq), ""),
                               strsplit(as.character(TSD_R_seq), "")), sum)
    p <- which(N_mismatch <= 1 & nchar(TSD_L_seq) >= 5)
    if (length(p) > 0) {
      TSD_Length <- max(p)
      TSD_sequence <- paste0(TSD_R_seq[TSD_Length],"/", TSD_L_seq[TSD_Length])
      TSD_position <- append(TSD_R[TSD_Length], TSD_L[TSD_Length])
    }else {
    TSD_Length <- 0
    TSD_sequence <- ""
    TSD_position <- NA
    }
  }

  TE_Length <- end(LTR_R) - start(LTR_L)
  LTR_Identity <- blast_line$pident
  out <- list(TSD_position = TSD_position, TSD_sequence = TSD_sequence, TSD_Length = TSD_Length,
              LTR_R_position = LTR_R, LTR_L_position = LTR_L, TE_Length = TE_Length, LTR_Identity = LTR_Identity)
  return(out)
}

get_best_ltr <- function(x) {
  tsd_ok <- sapply(x, function(k)k$TSD_Length > 3)
  te_length_ok <- sapply(x, function(k)k$TE_Length < 30000)
  ltr_length_ok <- sapply(x, function(k)width(k$LTR_R_position) >= 100 & width(k$LTR_L_position) >= 100)
  print("--------------------------------------------------------------------")
  if (sum(tsd_ok & te_length_ok & ltr_length_ok) >= 1) {
    # return the first one (best bitscore)
    print("number of OK ltr after filtering")
    print(length(x[tsd_ok & te_length_ok]))
    return(x[tsd_ok & te_length_ok][1])
  }
  if (any(te_length_ok & ltr_length_ok)) {
    print("number of OK ltr after filtering")
    print(length(x[te_length_ok & ltr_length_ok]))
    return(x[te_length_ok & ltr_length_ok][1])
  }else {
    return(NULL)
  }
}

get_te_gff3 <- function(g, ID) {
  D <- makeGRangesFromDataFrame(g$domain, keep.extra.columns = TRUE)
  sn <- seqnames(D)[1]
  S <- strand(D)[1]
  TE <- GRanges(seqnames = sn,
                IRanges(start = start(g$ltr_info[[1]]$LTR_L_position),
                        end = end(g$ltr_info[[1]]$LTR_R_position)), strand = S)
  TE$type <- "transposable_element"

  TE$Ndomains <- nrow(g$domain)
  TE$ID <- ID

  if (as.character(S) == "+") {
    LTR_5 <- g$ltr_info[[1]]$LTR_L
    LTR_3 <- g$ltr_info[[1]]$LTR_R
  }else {
    LTR_3 <- g$ltr_info[[1]]$LTR_L
    LTR_5 <- g$ltr_info[[1]]$LTR_R
  }
  LTR_5$LTR <- '5LTR'
  LTR_3$LTR <- '3LTR'
  LTR_5$type <- "long_terminal_repeat"
  LTR_3$type <- "long_terminal_repeat"
  strand(LTR_3) <- S
  strand(LTR_5) <- S
  LTR_3$Parent <- ID
  LTR_5$Parent <- ID
  LTR_3$Final_Classification <- D$Final_Classification[1]
  LTR_5$Final_Classification <- D$Final_Classification[1]
  LTR_5$LTR_Identity <- g$ltr_info[[1]]$LTR_Identity
  LTR_3$LTR_Identity <- g$ltr_info[[1]]$LTR_Identity

  TE$LTR_Identity <- g$ltr_info[[1]]$LTR_Identity
  TE$LTR5_length <- width(LTR_5)
  TE$LTR3_length <- width(LTR_3)

  if (is.na(g$ltr_info[[1]]$TSD_position)[1]) {
    # no TSD found
    TSD <- NULL
    TE$TSD <- 'not_found'
  }else {
    TSD <- g$ltr_info[[1]]$TSD_position
    TSD$type <- "target_site_duplication"
    TSD$Parent <- ID
    TE$TSD <- as.character(g$ltr_info[[1]]$TSD_sequence)
  }

  TE$Name <- TE$Final_Classification <- D$Final_Classification[1]

  D$Parent <- ID
  out <- c(TE, LTR_3, LTR_5, D, TSD)
  return(out)
}

get_TE <- function(Lseq, Rseq, domains_gff, GR, GRL, GRR,LTR_length) {
  xx <- blast(Lseq, Rseq, LTR_length)
  # sort by ltr position - first are the ltr closest to the domains
  ord <-  order(xx$qend - xx$sstart, decreasing = TRUE)
  xx <- xx[ord, , drop = FALSE]
  if (nrow(xx) == 0) {
    return(NULL)
  }else {
    ltr_tmp <- list()
    for (j in 1:nrow(xx)) {
      ltr_tmp[[j]] <- evaluate_ltr(GR, GRL, GRR, xx[j, , drop = FALSE], Lseq, Rseq)
    }
    print("number of ltr_tmp found")
    print(length(ltr_tmp))
    if (length(ltr_tmp) == 2) {
      print(ltr_tmp)
      print("~~~~~~~~~~~~~~~~~~~~~~~~~~")
    }
    ltr <- get_best_ltr(ltr_tmp)
    print("====================================================================")
    if (length(ltr) == 0) {
      return(NULL)
      ## add good ltr
    }else {
      return(list(domain = domains_gff, ltr_info = ltr, blast_out = xx))
    }
  }
}

add_pbs <- function(te, s, trna_db) {
  ltr5 <- te[which(te$LTR == "5LTR")]
  STRAND <- as.character(strand(te)[1])
  if (STRAND == "+") {
    pbs_gr <- GRanges(seqnames(ltr5), IRanges(start = end(ltr5) + 1, end = end(ltr5) + 31))
    pbs_s <- reverseComplement(getSeq(s, pbs_gr))
  }else {
    pbs_gr <- GRanges(seqnames(ltr5), IRanges(end = start(ltr5) - 1, start = start(ltr5) - 30))
    pbs_s <- getSeq(s, pbs_gr)
  }

  names(pbs_s) <- "pbs_region"
  # find trna match
  tmp <- tempfile()
  tmp_out <- tempfile()
  writeXStringSet(DNAStringSet(pbs_s), tmp)
  # alternative blast:
  cmd <- paste("blastn -task blastn -word_size 7 -dust no -perc_identity 100",
               " -query ", tmp, ' -db ', trna_db, ' -strand plus ',
               '-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"',
               '  -out', tmp_out)

  system(cmd)
  pbs_exact_gr <- NULL
  out_raw <- read.table(tmp_out, as.is = TRUE, sep = "\t",
                        col.names = strsplit(
                          "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
                          split = ' ')[[1]])
  if (nrow(out_raw) > 0) {
    cca <- grepl("CCA$", out_raw$qseq)
    to_end <- out_raw$send == 23 # align to end of sequence
    max_dist <- out_raw$qend > 25 # max 5 bp from ltr
    min_length <- out_raw$length >= 10
    out_pass <- out_raw[cca & to_end & max_dist & min_length,]
    if (nrow(out_pass) > 0) {
      trna_id <- out_pass$saccver[1]
      if (STRAND == "+") {
        S <- end(ltr5) + 32 - out_pass$qend[1]
        E <- end(ltr5) + 32 - out_pass$qstart[1]
      }else {
        S <- start(ltr5) - 31 + out_pass$qstart[1]
        E <- start(ltr5) - 31 + out_pass$qend[1]
      }
      pbs_exact_gr <- GRanges(seqnames(ltr5), IRanges(start = S, end = E))
      pbs_exact_gr$trna_id <- trna_id
      pbs_exact_gr$Length <- out_pass$Length
      strand(pbs_exact_gr) <- STRAND
      pbs_exact_gr$type <- 'primer_binding_site'
      pbs_exact_gr$Parent <- te[1]$ID
      te$trna_id <- c(trna_id, rep(NA, length(te) - 1))

    }
  }
  te <- append(te, pbs_exact_gr)
  unlink(tmp)
  unlink(tmp_out)
  return(te)
}


get_te_rank <- function (gr){
  DL_id <- gr$ID[gr$type == "transposable_element" &
                  gr$TSD == "not_found" &
                  is.na(gr$trna_id)]
  DLT_id <- gr$ID[gr$type == "transposable_element" &
                          gr$TSD != "not_found" &
                          is.na(gr$trna_id)]
  DLTP_id <- gr$ID[gr$type == "transposable_element" &
                              gr$TSD != "not_found" &
                              !is.na(gr$trna_id)]
  DLP_id <- gr$ID[gr$type == "transposable_element" &
                          gr$TSD == "not_found" &
                          !is.na(gr$trna_id)]
  Rank <- character(length(gr))
  ID <- unlist(ifelse(is.na(gr$ID), gr$Parent, gr$ID))

  Rank[ID %in% DL_id] <- "DL"
  Rank[ID %in% DLT_id] <- "DLT"
  Rank[ID %in% DLP_id] <- "DLP"
  Rank[ID %in% DLTP_id] <- "DLTP"
  return(Rank)
}

dante_filtering <- function(dante_gff, min_similarity=0.4,
                            min_identity=0.2, Relative_Length=0.6,
                            min_relat_interuptions=8) {
  include <- as.numeric(dante_gff$Similarity) >= min_similarity &
    as.numeric(dante_gff$Identity) >= min_identity &
    as.numeric(dante_gff$Relat_Length) >= Relative_Length &
    as.numeric(dante_gff$Relat_Interruptions) <= min_relat_interuptions
  # alternative threshold - if similarity is high, allow loew relative length
  min_rel_length_sim <- 0.35
  include2 <- as.numeric(dante_gff$Similarity) * as.numeric(dante_gff$Relat_Length) > min_rel_length_sim &
    as.numeric(dante_gff$Identity) >= min_identity &
    as.numeric(dante_gff$Relat_Interruptions) <= min_relat_interuptions

  include[is.na(include)] <- FALSE
  include2[is.na(include2)] <- FALSE
  return(dante_gff[include | include2,])
}

get_te_statistics <- function(gr, RT){
  Ranks <- c("D", "DL", "DLT", "DLP", "DLTP")
  all_class <- names(sort(table(RT$Final_Classification), decreasing = TRUE))
  RT_domain <- as.integer(table(factor(RT$Final_Classification, levels = all_class)))
  rank_table <- list()
  for (i in 1:length(Ranks)) {
    gr_part <- gr[gr$type == "transposable_element" & gr$Rank == Ranks[i]]
    rank_table[[Ranks[i]]] <- as.integer(table(factor(gr_part$Final_Classification, levels = all_class)))

  }
  out <- cbind(do.call(cbind, rank_table), RT_domain=RT_domain)
  out <- rbind(out, Total = colSums(out))
  rownames(out) <- c(all_class, "Total")
  return(out)
}


get_te_statistics_old <- function(gr, RT) {
  DOMAINS <- gr[gr$type == "transposable_element" & gr$Rank == "D"]
  DOMAINS_LTR <- gr[gr$type == "transposable_element" & gr$Rank == "DL"]
  DOMAINS_LTR_TSD <- gr[gr$type == "transposable_element" & gr$Rank == "DLT"]
  DOMAINS_LTR_TSD_PBS <- gr[gr$type == "transposable_element" &  gr$Rank == "DLTP"]
  DOMAINS_LTR_PBS <- gr[gr$type == "transposable_element"  & gr$Rank == "DLP"]

  all_class <- names(sort(table(RT$Final_Classification), decreasing = TRUE))
  RT_domain <- as.integer(table(factor(RT$Final_Classification, levels = all_class)))
  D <- as.integer(table(factor(DOMAINS$Final_Classification, levels = all_class)))
  DL <- as.integer(table(factor(DOMAINS_LTR$Final_Classification, levels = all_class)))
  DLT <- as.integer(table(factor(DOMAINS_LTR_TSD$Final_Classification, levels = all_class)))
  DLTP <- as.integer(table(factor(DOMAINS_LTR_TSD_PBS$Final_Classification, levels = all_class)))
  DLP <- as.integer(table(factor(DOMAINS_LTR_PBS$Final_Classification, levels = all_class)))
  out <- data.frame(RT_domain = RT_domain,
                    DOMAINS = D,
                    DOMAINS_LTR = DL,
                    DOMAINS_LTR_TSD = DLT,
                    DOMAINS_LTR_PBS = DLP,
                    DOMAINS_LTR_TSD_PBS = DLTP,
                    row.names = all_class
  )
  total <- colSums(out)
  out <- rbind(out, Total = total)
  return(out)
}

getSeqNamed <- function(s, gr, name = NULL) {
  spart <- getSeq(s, gr)
  if (is.null(name)){
    id1 <- paste0(seqnames(gr), '_', start(gr), "_", end(gr))
  }else{
    id1 <- mcols(gr)[,name]
  }
  id2 <- gr$Final_Classification
  names(spart) <- paste0(id1, "#", id2)
  spart
}

get_TE_id <- function (gr){
  gr_te <- gr[gr$type == "transposable_element"]
  ID <- gr_te$ID
  A <- paste0(seqnames(gr_te), '_', start(gr_te), "_", end(gr_te))
  B <- gr_te$Final_Classification
  names(ID) <- paste0(A, "#", B)

}


get_te_sequences <- function (gr, s) {
  Ranks <- c("D", "DL", "DLT", "DLP", "DLTP")
  s_te <- list()
    for (i in 1:length(Ranks)) {
        gr_te <- gr[gr$type == "transposable_element" & gr$Rank == Ranks[i]]
        if (length(gr_te) > 0) {
        s_te[[Ranks[i]]] <- getSeqNamed(s, gr_te)
        }
    }
  return(s_te)
}


cd_hit_est <- function(seqs, min_identity = 0.9, word_size = 10, ncpu = 2){
  # runs cd-hi-est and return table with cluster membership, and size and if reads was repesentative
  # input sequences must be in the same orientation!!!
  sfile <- tempfile()
  fasta_out <- tempfile()
  clstr <- paste0(fasta_out,".clstr")
  # cdhit is triming names!!
  ori_names <-  names(seqs)
  names(seqs) <- seq_along(seqs)
  writeXStringSet(seqs, sfile)
  cmd <- paste("cd-hit-est",
               "-i", sfile,
               "-o", fasta_out,
               "-c", min_identity,
               "-n", word_size,
               "-T", ncpu,
               "-r", 0)
  system(cmd)
  cls_raw <-  grep("^>", readLines(clstr), invert = TRUE, value = TRUE)
  unlink(fasta_out)
  unlink(clstr)
  index <-  gsub("\t.+","",cls_raw)
  id <-  as.numeric(gsub("[.].+","",
                       gsub(".+>", "", cls_raw))
  )
  is_representative <-  id %in% id[grep("[*]$",cls_raw)]
  membership <-  cumsum(index=="0")
  cluster_size <-  tabulate(membership)[membership]
  # reorder
  ord <- order(id)
    cls_info <- data.frame(
      seq_id = ori_names,
      membership = membership[ord],
      cluster_size = cluster_size[ord],
      is_representative = is_representative[ord]
    )
  return(cls_info)
}


dante2dist <-  function(dc){
  # dc - cluster of domains, granges
  dpos <- get_dom_pos(dc)
  D <-  c(dist(dpos))
  NN <-  combn(names(dpos), 2)
  # order feature in pair alphabeticaly
  N1 <-  ifelse(NN[1,]<NN[2,], NN[1,], NN[2, ])
  N2 <-  ifelse(NN[1,]>NN[2,], NN[1,], NN[2, ])
  dfd <-  data.frame(F1 = N1, F2 = N2, distance = D)
}



dante_to_quantiles <- function(dc, model_function, lineage=NULL){
  dfd <- dante2dist(dc)
  if (is.null(lineage)){
    lineage <- dc$Final_Classification[1]
  }
  # check if lineage has defined ecdf
  if (lineage %in% names(model_function$plus)){
    fn <- model_function$plus[[lineage]]
    fnm <- model_function$minus[[lineage]]
    dstat1 <- mapply(mapply(dfd$F1, dfd$F2, FUN=function(a, b)fn[[a]][[b]]), dfd$distance, FUN=function(f, x)f(x))
    dstat2 <- mapply(mapply(dfd$F1, dfd$F2, FUN=function(a, b)fnm[[a]][[b]]), dfd$distance, FUN=function(f, x)f(-x))
    dout <- cbind(dfd, dstat1, dstat2, dstat12 = ifelse(dstat1>dstat2, dstat2, dstat1))
    return(dout)
  }else{
    return(NULL)
  }
}

get_dom_pos <- function(g){
  if (length(g) == 0){
    return(NULL)
  }
  gdf <-  data.frame(rexdb = as.character(seqnames(g)),
                   domain = g$Name,
                   S = start(g),
                   E = end(g),
                   stringsAsFactors = FALSE)
  gSmat <- xtabs(S  ~ rexdb + domain, data = gdf)
  colnames(gSmat) <-  paste0(colnames(gSmat),"_S")
  gEmat <- xtabs(E  ~ rexdb + domain, data = gdf)
  colnames(gEmat) <-  paste0(colnames(gEmat),"_E")
  SEmat <- cbind(gSmat,gEmat)
  # reorder
  dom_position <- colMeans(SEmat)
  SEmat_sorted <-  SEmat[,order(dom_position)]
  return(SEmat_sorted)
}
