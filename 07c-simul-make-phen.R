## sum PGS for different chr and add transf. + nois

library(data.table)
library(ggplot2)
library(hash)

set.seed(283619)


n_causal_snps_str  <- c("100", "1k", "10k")
architecture       <- c("rg", "cl")
architecture_hash  <- hash(architecture, c("regular", "clustered"))

transf_type        <- c("nt", "sg")  # no transformation / sigmoid
transf_timing      <- c("bn", "an")  # before/after noise
transf_timing_hash <- hash(transf_timing, c("before add. noise", "after add. noise"))
heritability       <- c(0.3, 0.6)


for (ncs in n_causal_snps_str) {
    for (arch in architecture) {

        ## load chr1 score
        pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                               ncs, "_", arch, "-chr1.sscore"),
                        data.table = FALSE)

        ## rescale components by allele count
        pgs_df[, 5] <- pgs_df[, 5] * pgs_df$ALLELE_CT

        ## add remaining components
        for (chr in 2:22) {
            pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                                       ncs, "_", arch, "-chr", chr, ".sscore"),
                               data.table = FALSE)
            pgs_df_chr[, 5] <- pgs_df_chr[, 5] * pgs_df_chr$ALLELE_CT
            pgs_df[, 5] <- pgs_df[, 5] + pgs_df_chr[, 5]
        }
        pgs_df <- pgs_df[, c(1, 2, 5)]
        

        ## add transformation + noise
        for (h2 in heritability) {
            for (tr in transf_type) {
                if (tr == "sg") {
                    for (tr_t in transf_timing) {

                        phen_df <- pgs_df

                        if (tr_t == "bn") {

                            ## transform
                            phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3]))

                            ## add environmental noise
                            colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                            phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                                rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                        } else {

                            ## add environmental noise
                            colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                            phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                                rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                            ## transform
                            phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3]))

                        }

                        ## export
                        phen_df_out <- phen_df[, 2:3]
                        colnames(phen_df_out)[1] <- "#IID"
                        phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                        fwrite(phen_df_out,
                               file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2, ".tab"),
                               sep = '\t', na = 'NA', quote = FALSE)
                    }        
                } else {
                    
                    phen_df <- pgs_df

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, "_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                    ## export
                    phen_df_out <- phen_df[, 2:3]
                    colnames(phen_df_out)[1] <- "#IID"
                    phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                    fwrite(phen_df_out,
                           file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, "_", h2, ".tab"),
                           sep = '\t', na = 'NA', quote = FALSE)
                }
            }
        }
    }
}




## scaled sigmoids
set.seed(948293)

for (ncs in n_causal_snps_str) {
    for (arch in architecture) {

        ## load chr1 score
        pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                               ncs, "_", arch, "-chr1.sscore"),
                        data.table = FALSE)

        ## rescale components by allele count
        pgs_df[, 5] <- pgs_df[, 5] * pgs_df$ALLELE_CT

        ## add remaining components
        for (chr in 2:22) {
            pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                                       ncs, "_", arch, "-chr", chr, ".sscore"),
                                data.table = FALSE)
            pgs_df_chr[, 5] <- pgs_df_chr[, 5] * pgs_df_chr$ALLELE_CT
            pgs_df[, 5] <- pgs_df[, 5] + pgs_df_chr[, 5]
        }
        pgs_df <- pgs_df[, c(1, 2, 5)]
        

        ## add transformation + noise
        for (h2 in heritability) {
            tr <- "sgsc"

            for (tr_t in transf_timing) {

                phen_df <- pgs_df

                if (tr_t == "bn") {

                    score <- phen_df[, 3]  # store original score for plotting

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                    ## add environmental noise
                    transform <- as.numeric(scale(phen_df[, 3]))  # store transformation before noise but after rescaling for plotting
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                    ## plot
                    p_df <- cbind(phen_df, score, transform)
                    colours <- c("Initial" = "black", "Before noise" = "#D35F27", "Final" = "#0073AD")

                    p <- ggplot(p_df) +
                        geom_jitter(aes(x = score, y = 0, color = "Initial"),
                                    size = 1, shape = ".") +
                        geom_point(aes(x = score, y = .data[[colnames(phen_df)[3]]], color = "Final"),
                                   size = 1, shape = ".", alpha = 0.5) +
                        geom_line(aes(x = score, y = transform, color = "Before noise")) +
                        ggtitle(colnames(phen_df)[3]) +
                        labs(x = "Simulated phenotype before any scaling/transformation",
                             y = "Simulated phenotype after scaling/transformation",
                             color = "Legend") +
                        scale_colour_manual(values = colours) +
                        guides(colour = guide_legend(override.aes = list(shape = 19))) +
                        theme_bw() +
                        theme(plot.title = element_text(size = 9),
                              axis.title = element_text(size = 8),
                              axis.text  = element_text(size = 7),
                              legend.position = c(0.5, 0.08),
                              legend.direction = "horizontal",
                              legend.title = element_blank(),
                              legend.text = element_text(size = 7),
                              legend.background = element_rect(fill = "white", color = "black", size = 0.1),
                              legend.margin = margin(c(1, 10, 1, 1)))

                    ggsave(p, file = paste0("../data/phenotypes/simulations-prep/figs/", colnames(phen_df)[3],
                                            "-transform.png"),
                           type = 'cairo-png', width = 1200/120, height = 800/120, units = 'in', dpi = 120)

                } else {

                    score <- phen_df[, 3]  # store original score for plotting

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))
                    noise <- phen_df[, 3]

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                    ## plot
                    p_df <- cbind(phen_df, score, noise)
                    colours <- c("Initial" = "black", "After noise" = "#0073AD", "Final" = "#D35F27")

                    p <- ggplot(p_df) +
                        geom_jitter(aes(x = score, y = 0, color = "Initial"),
                                    size = 1, shape = ".") +
                        geom_point(aes(x = score, y = noise, color = "After noise"),
                                   size = 1, shape = ".", alpha = 0.5) +
                        geom_point(aes(x = score, y = .data[[colnames(phen_df)[3]]], color = "Final"),
                                   size = 1, shape = ".", alpha = 0.5) +
                        ggtitle(colnames(phen_df)[3]) +
                        labs(x = "Simulated phenotype before any scaling/transformation",
                             y = "Simulated phenotype after scaling/transformation",
                             color = "Legend") +
                        scale_colour_manual(values = colours) +
                        guides(colour = guide_legend(override.aes = list(shape = 19))) +
                        theme_bw() +
                        theme(plot.title = element_text(size = 9),
                              axis.title = element_text(size = 8),
                              axis.text  = element_text(size = 7),
                              legend.position = c(0.5, 0.08),
                              legend.direction = "horizontal",
                              legend.title = element_blank(),
                              legend.text = element_text(size = 7),
                              legend.background = element_rect(fill = "white", color = "black", size = 0.1),
                              legend.margin = margin(c(1, 10, 1, 1)))

                    ggsave(p, file = paste0("../data/phenotypes/simulations-prep/figs/", colnames(phen_df)[3],
                                            "-transform.png"),
                           type = 'cairo-png', width = 1200/120, height = 800/120, units = 'in', dpi = 120)

                }

                ## export
                phen_df_out <- phen_df[, 2:3]
                colnames(phen_df_out)[1] <- "#IID"
                phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                fwrite(phen_df_out,
                       file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2, ".tab"),
                       sep = '\t', na = 'NA', quote = FALSE)
            }        
        }
    }
}




## make versions with scaled betas
set.seed(583920)

for (ncs in n_causal_snps_str) {
    for (arch in architecture) {

        ## load chr1 score
        pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                               ncs, "_", arch, "-chr1.sscore"),
                        data.table = FALSE)

        ## rescale components by allele count
        pgs_df[, 6] <- pgs_df[, 6] * pgs_df$ALLELE_CT

        ## add remaining components
        for (chr in 2:22) {
            pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                                       ncs, "_", arch, "-chr", chr, ".sscore"),
                               data.table = FALSE)
            pgs_df_chr[, 6] <- pgs_df_chr[, 6] * pgs_df_chr$ALLELE_CT
            pgs_df[, 6] <- pgs_df[, 6] + pgs_df_chr[, 6]
        }
        pgs_df <- pgs_df[, c(1, 2, 6)]
        

        ## add transformation + noise
        for (h2 in heritability) {

            ## no transformation
            tr <- "nt"
                   
            phen_df <- pgs_df

            ## add environmental noise
            colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, "_a5_", h2)
            phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

            ## export
            phen_df_out <- phen_df[, 2:3]
            colnames(phen_df_out)[1] <- "#IID"
            phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

            fwrite(phen_df_out,
                   file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, "_a5_", h2, ".tab"),
                   sep = '\t', na = 'NA', quote = FALSE)


            ## scaled sigmoids
            tr <- "sgsc"

            for (tr_t in transf_timing) {

                phen_df <- pgs_df

                if (tr_t == "bn") {

                    score <- phen_df[, 3]  # store original score for plotting

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                    ## add environmental noise
                    transform <- as.numeric(scale(phen_df[, 3]))  # store transformation before noise but after rescaling for plotting
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_a5_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                    ## plot
                    p_df <- cbind(phen_df, score, transform)
                    colours <- c("Initial" = "black", "Before noise" = "#D35F27", "Final" = "#0073AD")

                    p <- ggplot(p_df) +
                        geom_jitter(aes(x = score, y = 0, color = "Initial"),
                                    size = 1, shape = ".") +
                        geom_point(aes(x = score, y = .data[[colnames(phen_df)[3]]], color = "Final"),
                                   size = 1, shape = ".", alpha = 0.5) +
                        geom_line(aes(x = score, y = transform, color = "Before noise")) +
                        ggtitle(colnames(phen_df)[3]) +
                        labs(x = "Simulated phenotype before any scaling/transformation",
                             y = "Simulated phenotype after scaling/transformation",
                             color = "Legend") +
                        scale_colour_manual(values = colours) +
                        guides(colour = guide_legend(override.aes = list(shape = 19))) +
                        theme_bw() +
                        theme(plot.title = element_text(size = 9),
                              axis.title = element_text(size = 8),
                              axis.text  = element_text(size = 7),
                              legend.position = c(0.5, 0.08),
                              legend.direction = "horizontal",
                              legend.title = element_blank(),
                              legend.text = element_text(size = 7),
                              legend.background = element_rect(fill = "white", color = "black", size = 0.1),
                              legend.margin = margin(c(1, 10, 1, 1)))

                    ggsave(p, file = paste0("../data/phenotypes/simulations-prep/figs/", colnames(phen_df)[3],
                                            "-transform.png"),
                           type = 'cairo-png', width = 1200/120, height = 800/120, units = 'in', dpi = 120)

                } else {

                    score <- phen_df[, 3]  # store original score for plotting

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_a5_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))
                    noise <- phen_df[, 3]

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                    ## plot
                    p_df <- cbind(phen_df, score, noise)
                    colours <- c("Initial" = "black", "After noise" = "#0073AD", "Final" = "#D35F27")

                    p <- ggplot(p_df) +
                        geom_jitter(aes(x = score, y = 0, color = "Initial"),
                                    size = 1, shape = ".") +
                        geom_point(aes(x = score, y = noise, color = "After noise"),
                                   size = 1, shape = ".", alpha = 0.5) +
                        geom_point(aes(x = score, y = .data[[colnames(phen_df)[3]]], color = "Final"),
                                   size = 1, shape = ".", alpha = 0.5) +
                        ggtitle(colnames(phen_df)[3]) +
                        labs(x = "Simulated phenotype before any scaling/transformation",
                             y = "Simulated phenotype after scaling/transformation",
                             color = "Legend") +
                        scale_colour_manual(values = colours) +
                        guides(colour = guide_legend(override.aes = list(shape = 19))) +
                        theme_bw() +
                        theme(plot.title = element_text(size = 9),
                              axis.title = element_text(size = 8),
                              axis.text  = element_text(size = 7),
                              legend.position = c(0.5, 0.08),
                              legend.direction = "horizontal",
                              legend.title = element_blank(),
                              legend.text = element_text(size = 7),
                              legend.background = element_rect(fill = "white", color = "black", size = 0.1),
                              legend.margin = margin(c(1, 10, 1, 1)))

                    ggsave(p, file = paste0("../data/phenotypes/simulations-prep/figs/", colnames(phen_df)[3],
                                            "-transform.png"),
                           type = 'cairo-png', width = 1200/120, height = 800/120, units = 'in', dpi = 120)

                }

                ## export
                phen_df_out <- phen_df[, 2:3]
                colnames(phen_df_out)[1] <- "#IID"
                phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                fwrite(phen_df_out,
                       file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, tr_t, "_a5_", h2, ".tab"),
                       sep = '\t', na = 'NA', quote = FALSE)
            }        
        }
    }
}
