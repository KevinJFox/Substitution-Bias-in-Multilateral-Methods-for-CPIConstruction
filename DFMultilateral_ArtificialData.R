#########################################################################
###### R Code for Multilateral Methods with Scanner Data ################
###### Reference: Diewert, W.E. and K.J. Fox (2020), 
###### "Subsititution Bias in Multilatreal Methods for CPI Construction Using Scanner Data"
#########################################################################
###### Contact: Kevin Fox (K.Fox@unsw.edu.au)
###### 30 June 2020

########################################################
#### Multilateral Methods Code in R by Shipei Zeng #####
########################################################


###### Outputs from this script: Results using artificial data
###### 
###### Table 1 - Table 4: df_1, df_2, df_3, df_4
###### Table A1 -Table A6: df_a1, df_a2, df_a3, df_a4, df_a5, df_a6
###### Figure 1 - Figure 3: see plots
######
###### Latex code for making tables is available by write files for
###### tbl_xx_latex (e.g. write(tbl_1_latex, "tbl_1_latex.txt"))


#########################################################################
########################### Loading Packages ############################
#########################################################################
library(xtable)
library(rstudioapi) # extract path for setting up the working directory
#########################################################################



#########################################################################
########################### Function Set-up #############################
#########################################################################

### --------------------- General Form (begin) -----------------------###

### a warp-up function for bilateral price indexes
# df: an ordered dataframe structure by period
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
# type: "f" for "Fisher" (default), "t" for "Tornqvist",
#     "l" for "Laspeyres", "p" for "Paasche"
#     vector is allowed, e.g., c("f", "l", "p")
# seq: "ch" for "chained indexes" (default),
#     "fb" for "fixed base indexes
#     vector is allowed, e.g., c("ch", "fb")

pb_index <- function(df, p, qty, id, tm, type="f", seq="ch") {
        # data preparation
        df <- df[order(as.factor(df[,tm]), df[,id]),]
        tm_list <- unique(df[,tm])
        val <- df[,p]*df[,qty]
        # value shares
        df[,"s"] <- unlist(by(val, as.factor(df[,tm]), function(x) {x/sum(x)}))
        # run basic_index function
        pb_out <- basic_index(df, p, qty, id, tm, tm_list)
        # user options
        use_opt <-  unlist(lapply(type, function(x) paste(x, seq, sep="_"))) 
        return(pb_out[,c(tm, use_opt)])
}

### a warp-up function for multilateral price indexes
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
# type: "geks" for "GEKS" (default), "ccdi" for "CCDI",
#     "wtpd" for "weighted time product dummy", 
#     "gk" for "Geary-Khamis"
# len: window length for linked indexes, single window
#     as len=NULL (default)
# link: the linking observation position in the overlapping, 
#     effective when len is not NULL, no linking as link
#     =NULL (default)
#     mean splice as link="mean",
#     options include "mean" and numbers

pm_index <- function(df, p, qty, id, tm, type="geks", len=NULL, link=NULL) {
        # data preparation
        df <- df[order(as.factor(df[,tm]), df[,id]),]
        tm_list <- unique(df[,tm])
        val <- df[,p]*df[,qty]
        # value shares
        df[,"s"] <- unlist(by(val, as.factor(df[,tm]), function(x) {x/sum(x)}))
        # user options
        use_list <- c("geks", "ccdi", "wtpd", "gk")
        use_id <- which(use_list %in% type)
        # construct a window
        if (is.null(len)) {
                len <- length(tm_list)
        }
        wd_tick <- length(tm_list)-len+1
        df_wd <- df[df[,tm] %in% tm_list[wd_tick:(wd_tick+len-1)],]
        tm_list_wd <- unique(df_wd[,tm])
        pindex_out <- lapply(use_id, function(fun_id) {
                p_fun <- switch(fun_id, geks_index, ccdi_index, 
                                wtpd_index, gk_index)
                p_fun_out <- p_fun(df_wd, p, qty, id, tm, tm_list_wd)
        })
        pindex_out <- data.frame(tm_list_wd, do.call(cbind, pindex_out))
        colnames(pindex_out) <- c(tm, use_list[use_id])
        # iterative update
        while (nrow(pindex_out) < length(tm_list)) {
                wd_tick <- wd_tick-1
                df_wd <- df[df[,tm] %in% tm_list[wd_tick:(wd_tick+len-1)],]
                tm_list_wd <- unique(df_wd[,tm])
                pindex_temp <- lapply(use_id, function(fun_id) {
                        p_fun <- switch(fun_id, geks_index, ccdi_index, 
                                        wtpd_index, gk_index)
                        p_fun_out <- p_fun(df_wd, p, qty, id, tm, tm_list_wd)
                })
                pindex_temp <- data.frame(tm_list_wd, do.call(cbind, pindex_temp))
                colnames(pindex_temp) <- c(tm, use_list[use_id])
                # append the linked indexes in a reversed way
                if (is.numeric(link)) {
                        # for link as numbers
                        link_obs_out <- pindex_out[link, -1]
                        link_obs_temp <- pindex_temp[link+1, -1]
                        link_rate <- link_obs_temp/link_obs_out
                } else if (link=="mean") {
                        # for link as "mean"
                        link_obs_out <- pindex_out[1:(len-1), -1]
                        link_obs_temp <- pindex_temp[2:len, -1]
                        link_rate <- as.data.frame(link_obs_temp/link_obs_out)
                        link_rate <- apply(link_rate, 2, function(x) {exp(mean(log(x)))})
                }
                link_fill <- nrow(pindex_out)+1-len
                link_obs_fill <- as.data.frame(tail(pindex_out, link_fill)[, -1])
                link_obs_fill <- rbind(link_rate, link_obs_fill)
                link_obs_fill <- as.data.frame(apply(link_obs_fill, 2, function(x) {x*x[1]}))
                link_obs_fill <- data.frame(tail(pindex_out, link_fill)[, 1], link_obs_fill[-1,])
                colnames(link_obs_fill) <- c(tm, use_list[use_id])
                pindex_temp <- rbind(pindex_temp, link_obs_fill)
                pindex_out <- pindex_temp
        }
        rownames(pindex_out) <- NULL
        pindex_out <- pindex_out[, c(tm, type)]
        return(pindex_out)
}



### ---------------------- General Form (end) ------------------------###

### -------------------- Specific Form (begin) -----------------------###

### fixed base and chained price indexes
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods

basic_index <- function(df, p, qty, id, tm, tm_list) {
        # fixed base indexes
        p_0 <- df[df[,tm]==tm_list[1], p]
        q_0 <- df[df[,tm]==tm_list[1], qty]
        s_0 <- df[df[,tm]==tm_list[1], "s"]
        pindex_fb <- lapply(tm_list, function(x) {
                p_1 <- df[df[,tm]==x, p]
                q_1 <- df[df[,tm]==x, qty]
                s_1 <- df[df[,tm]==x, "s"]
                las_fb <- crossprod(p_1, q_0)/crossprod(p_0, q_0)
                paa_fb <- crossprod(p_1, q_1)/crossprod(p_0, q_1)
                fis_fb <- sqrt(las_fb*paa_fb)
                tor_fb <- prod((p_1/p_0)^(0.5*(s_0+s_1)))
                output_fb <- c(fis_fb, tor_fb, las_fb, paa_fb)
                return(output_fb)
        })
        pindex_fb <- do.call(rbind, pindex_fb)
        # chained indexes
        pindex_ch_g <- lapply(tm_list[-1], function(y) {
                p_1 <- df[df[,tm]==y, p]
                q_1 <- df[df[,tm]==y, qty]
                s_1 <- df[df[,tm]==y, "s"]
                y_tick <- which(tm_list==y)
                p_0 <- df[df[,tm]==tm_list[y_tick-1], p]
                q_0 <- df[df[,tm]==tm_list[y_tick-1], qty]
                s_0 <- df[df[,tm]==tm_list[y_tick-1], "s"]
                las_ch <- crossprod(p_1, q_0)/crossprod(p_0, q_0)
                paa_ch <- crossprod(p_1, q_1)/crossprod(p_0, q_1)
                fis_ch <- sqrt(las_ch*paa_ch)
                tor_fb <- prod((p_1/p_0)^(0.5*(s_0+s_1)))
                output_ch <- c(fis_ch, tor_fb, las_ch, paa_ch)
                return(output_ch)
        })
        pindex_ch_g <- do.call(rbind, pindex_ch_g)
        pindex_ch_g <- rbind(rep(1,4), pindex_ch_g)
        pindex_ch <- apply(pindex_ch_g, 2, cumprod)
        # put fixed base indexes and chained indexed together
        pindex_out <- data.frame(tm_list, pindex_fb, pindex_ch)
        colnames(pindex_out) <- c(tm,"f_fb","t_fb","l_fb","p_fb",
                              "f_ch","t_ch","l_ch","p_ch")
        return(pindex_out)
}

### GEKS and CCDI price indexes
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
geks_index <- function(df, p, qty, id, tm, tm_list) {
        # GEKS indexes
        p_loop <- lapply(tm_list, function(y) {
                # fixed base indexes
                p_0 <- df[df[,tm]==y, p]
                q_0 <- df[df[,tm]==y, qty]
                s_0 <- df[df[,tm]==y, "s"]
                pindex_fb <- lapply(tm_list, function(x) {
                        p_1 <- df[df[,tm]==x, p]
                        q_1 <- df[df[,tm]==x, qty]
                        s_1 <- df[df[,tm]==x, "s"]
                        las_fb <- crossprod(p_1, q_0)/crossprod(p_0, q_0)
                        paa_fb <- crossprod(p_1, q_1)/crossprod(p_0, q_1)
                        fis_fb <- sqrt(las_fb*paa_fb)
                        return(fis_fb)
                })
                pindex_fb <- do.call(rbind, pindex_fb)
        })
        p_loop <- Reduce("*", p_loop)
        p_loop <- p_loop^(1/length(tm_list))
        # normalised
        p_loop <- p_loop/p_loop[1]
        return(p_loop)
}

ccdi_index <- function(df, p, qty, id, tm, tm_list) {
        # GEKS indexes
        p_loop <- lapply(tm_list, function(y) {
                # fixed base indexes
                p_0 <- df[df[,tm]==y, p]
                s_0 <- df[df[,tm]==y, "s"]
                pindex_fb <- lapply(tm_list, function(x) {
                        p_1 <- df[df[,tm]==x, p]
                        s_1 <- df[df[,tm]==x, "s"]
                        tor_fb <- prod((p_1/p_0)^(0.5*(s_0+s_1)))
                        return(tor_fb)
                })
                pindex_fb <- do.call(rbind, pindex_fb)
        })
        p_loop <- Reduce("*", p_loop)
        p_loop <- p_loop^(1/length(tm_list))
        # normalised
        p_loop <- p_loop/p_loop[1]
        return(p_loop)
}

### weighted time product dummy method
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
wtpd_index <- function(df, p, qty, id, tm, tm_list) {
        # value share matrix
        s_mat <- matrix(df[,"s"], byrow=TRUE, nrow=length(tm_list))
        w_tnj_list <- lapply(as.data.frame(t(s_mat)), function(x) {
                s_nj <- x %*% t(x)
                diag(s_nj) <- 0
                return(s_nj)
        })
        w_nj_mat <- Reduce("+", w_tnj_list)
        f_nj_mat <- t(apply(w_nj_mat, 1, function(x) {x/sum(x)}))
        w_nk_sum <- apply(w_nj_mat, 1, function(x) {sum(x)})
        f_tnj_list <- lapply(w_tnj_list, function(x) {
                apply(x, 2, function(y) {y/w_nk_sum})
        })
        # I_N, F and f used to solve equations
        # I_N: diag_mat
        # F: f_nj_mat
        # f: f_n_vec
        n_id <- ncol(s_mat)
        diag_mat <- diag(nrow=n_id)
        y_mat <- log(matrix(df[, p], byrow=TRUE, nrow=length(tm_list)))
        y_tnj_list <- lapply(as.data.frame(t(y_mat)), function(x) {
                y_nj <- matrix(rep(x, n_id), nrow=n_id)
                y_jn <- matrix(rep(x, n_id), nrow=n_id, byrow=TRUE)
                temp <- y_nj-y_jn
                return(temp)
        })
        fy_tnj_list <- Map("*", f_tnj_list, y_tnj_list)
        f_n_vec <-  apply(Reduce("+", fy_tnj_list), 1, function(x) {sum(x)})
        # check if the inverse of I_N-F exists
        i_n_f <- diag_mat-f_nj_mat
        inv_check <- try(solve(i_n_f, f_n_vec), TRUE)
        if (inherits(inv_check, "try-error")) {
                beta <- c(solve(i_n_f[-n_id, -n_id], f_n_vec[-n_id]),0)
        } else {
                beta <- solve(i_n_f, f_n_vec)
        }
        # weighted time product dummy price level
        b <- exp(beta)
        p_temp <- lapply(tm_list, function(x) {
                p_t <- df[df[,tm]==x, p]
                s_t <- df[df[,tm]==x, "s"]
                a_t <- prod((p_t/b)^(s_t))
        })
        pindex_wtpd <- do.call(rbind, p_temp)
        # normalised
        pindex_wtpd <- pindex_wtpd/pindex_wtpd[1]
        return(pindex_wtpd)
}

### Geary-Khamis multilateral indexes
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
gk_index <- function(df, p, qty, id, tm, tm_list) {
        # matrices required to solve the equation system
        q_mat <- matrix(df[,qty], byrow=TRUE, nrow=length(tm_list))
        q_hat_inv <-  solve(diag(apply(q_mat, 2, sum)))
        s_mat <- matrix(df[,"s"], byrow=TRUE, nrow=length(tm_list))
        q_list <- as.list(as.data.frame(t(q_mat)))
        s_list <- as.list(as.data.frame(t(s_mat)))
        sq_list <- Map(function(x,y) {x %*% t(y)}, s_list, q_list)
        sq_sum <- Reduce("+", sq_list)
        # I_N, C and 0_N used to solve equations
        # I_N: diag_mat
        # C: c_mat
        # 0_N: vec_0
        n_id <- ncol(s_mat)
        diag_mat <- diag(nrow=n_id)
        c_mat <- q_hat_inv %*% sq_sum
        i_n_c <- diag_mat-c_mat
        vec_part <- c_mat[-nrow(c_mat), ncol(c_mat)]
        b_vec <- c(solve(i_n_c[-n_id, -n_id], vec_part), 1)
        # GK indexes
        p_temp <- lapply(tm_list, function(x) {
                p_1 <- df[df[,tm]==x, p]
                q_1 <- df[df[,tm]==x, qty]
                gk_p <- crossprod(p_1, q_1)/crossprod(b_vec, q_1)
                return(gk_p)
        })
        pindex_gk <- do.call(rbind, p_temp)
        # normalised
        pindex_gk <- pindex_gk/pindex_gk[1]
        return(pindex_gk)
}

### CES indexes
# df: a dataframe structure
# p: variable name for prices
# val: variable name for expenditures
# id: variable name product id
# tm: variable name for periods
# sigma: elasticity
# alpha: a vector of parameters following the 
#     same order as ascending categories of "id"
# len: window length for linked indexes, single window
#     as len=NULL (default)
# link: the linking observation position in the overlapping, 
#     effective when len is not NULL, 
#     mean splice as link="mean" (default),
#     options include "mean" and numbers

ces_index <- function(df, p, val, id, tm, sigma, alpha, len=NULL, link="mean") {
        # data preparation
        df <- df[order(as.factor(df[,tm]), df[,id]),]
        tm_list <- unique(df[,tm])
        alpha <- rep(alpha, times=length(tm_list))
        s_temp <- alpha*df[,p]^(1-sigma)
        s <- unlist(by(s_temp, as.factor(df[,tm]), function(x) {x/sum(x)}))
        df[, "q"] <- df[, val]*s/df[, p]
        # construct a window
        if (is.null(len)) {
                len <- length(tm_list)
        }
        wd_tick <- length(tm_list)-len+1
        df_wd <- df[df[,tm] %in% tm_list[wd_tick:(wd_tick+len-1)],]
        s_temp_wd <- s_temp[df[,tm] %in% tm_list[wd_tick:(wd_tick+len-1)]]
        tm_list_wd <- unique(df_wd[,tm])
        # index output
        pindex_ces <- as.vector(unlist(by(s_temp_wd, as.factor(df_wd[,tm]), function(x) {
                (sum(x))^(1/(1-sigma))})))
        # normalised
        pindex_ces <- pindex_ces/pindex_ces[1]
        pindex_out <- data.frame(tm_list_wd, pindex_ces)
        colnames(pindex_out) <- c(tm, "ces")
        # iterative update
        while (nrow(pindex_out) < length(tm_list)) {
                wd_tick <- wd_tick-1
                df_wd <- df[df[,tm] %in% tm_list[wd_tick:(wd_tick+len-1)],]
                s_temp_wd <- s_temp[df[,tm] %in% tm_list[wd_tick:(wd_tick+len-1)]]
                tm_list_wd <- unique(df_wd[,tm])
                pindex_temp <- as.vector(unlist(by(s_temp_wd, as.factor(df_wd[,tm]), function(x) {
                        (sum(x))^(1/(1-sigma))})))
                pindex_temp <- pindex_temp/pindex_temp[1]
                pindex_temp <- data.frame(tm_list_wd, pindex_temp)
                colnames(pindex_temp) <- c(tm, "ces")
                # append the linked indexes in a reversed way
                if (is.numeric(link)) {
                        # for link as numbers
                        link_obs_out <- pindex_out[link, -1]
                        link_obs_temp <- pindex_temp[link+1, -1]
                        link_rate <- link_obs_temp/link_obs_out
                } else if (link=="mean") {
                        # for link as "mean"
                        link_obs_out <- pindex_out[1:(len-1), -1]
                        link_obs_temp <- pindex_temp[2:len, -1]
                        link_rate <- as.data.frame(link_obs_temp/link_obs_out)
                        link_rate <- apply(link_rate, 2, function(x) {exp(mean(log(x)))})
                }
                link_fill <- nrow(pindex_out)+1-len
                link_obs_fill <- as.data.frame(tail(pindex_out, link_fill)[, -1])
                link_obs_fill <- rbind(link_rate, link_obs_fill)
                link_obs_fill <- as.data.frame(apply(link_obs_fill, 2, function(x) {x*x[1]}))
                link_obs_fill <- data.frame(tail(pindex_out, link_fill)[, 1], link_obs_fill[-1,])
                colnames(link_obs_fill) <- c(tm, "ces")
                pindex_temp <- rbind(pindex_temp, link_obs_fill)
                pindex_out <- pindex_temp
        }
        # data output
        rownames(pindex_out) <- NULL
        df_out <- df[order(as.numeric(rownames(df))), c(tm, id, p, "q")]
        out_list <- list(pindex_out, df_out)
        names(out_list) <- c("CES Indexes", "Data")
        return(out_list)
}

### relative price dissimilarity
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
# type: "f" for "Fisher" (default), "t" for "Tornqvist",
# sim: "al" for "weighted asymptotically linear index of  
#     relative price dissimilarity" (default), "lq" for 
#     "weighted log quadratic index of relatve dissimilarity"

p_cmp <- function(df, p, qty, id, tm, type="f", sim="al") {
        # data preparation
        df <- df[order(as.factor(df[,tm]), df[,id]),]
        tm_list <- unique(df[,tm])
        val <- df[,p]*df[,qty]
        # value shares
        df[,"s"] <- unlist(by(val, as.factor(df[,tm]), function(x) {x/sum(x)}))
        # a matrix for dissimilarity measures
        rel_ij <- lapply(tm_list, function(perd_i) {
                rel_i <- lapply(tm_list, function(perd_j) {
                        df_temp <- df[df[,tm]==perd_i | df[,tm]==perd_j,]
                        if (perd_i==perd_j) {p_temp <- 1
                        } else if (perd_i>perd_j){
                                p_temp <- 1/pb_index(df_temp, p, qty, id, tm, type, "fb")[-1,-1]
                        } else {
                                p_temp <- pb_index(df_temp, p, qty, id, tm, type, "fb")[-1,-1]
                        }
                        s_0 <- df_temp[df_temp[,tm]==perd_i,"s"]
                        s_1 <- df_temp[df_temp[,tm]==perd_j,"s"]
                        p_0 <- df_temp[df_temp[,tm]==perd_i, p]
                        p_1 <- df_temp[df_temp[,tm]==perd_j, p]
                        if (sim=="lq") {
                                rel_temp <- sum(((log(p_1/(p_temp*p_0)))^2)*(s_1+s_0)*0.5)
                        } else if (sim=="al") {
                                rel_temp <- sum((p_1/(p_temp*p_0)+(p_temp*p_0)/p_1-2)*(s_1+s_0)*0.5)
                        }
                        return(rel_temp)
                })
                rel_i <- unlist(rel_i)
                return(rel_i)
        })
        rel_ji <- t(do.call(rbind, rel_ij))
        p_rel_out <- data.frame(tm_list, rel_ji)
        colnames(p_rel_out) <- c(tm, tm_list)
        return(p_rel_out)
}

### relative price linking
# df: a dataframe structure
# p: variable name for prices
# qty: variable name for quantities
# id: variable name product id
# tm: variable name for periods
# type: "f" for "Fisher" (default), "t" for "Tornqvist",
# sim: "al" for "weighted asymptotically linear index of  
#     relative price dissimilarity" (default), "lq" for 
#     "weighted log quadratic index of relatve dissimilarity"

sim_index <- function(df, p, qty, id, tm, type="f", sim="al") {
        sim_mat <- p_cmp(df, p, qty, id, tm, type, sim)
        # data preparation
        df <- df[order(as.factor(df[,tm]), df[,id]),]
        tm_list <- unique(df[,tm])
        # recursive linking
        p_list <- list()
        p_list[[1]] <- 1
        for (i in 2:length(tm_list)) {
                link_id <- which.min(sim_mat[1:(i-1),i+1])
                df_temp <- df[df[,tm]==tm_list[i] | df[,tm]==tm_list[link_id],]
                p_temp <- pb_index(df_temp, p, qty, id, tm, type, "fb")[-1,-1]
                p_link <- p_temp*p_list[[link_id]]
                p_list[[i]] <- p_link
        }
        pindex_out <- data.frame(tm_list, unlist(p_list))
        colnames(pindex_out) <- c(tm, sim)
        return(pindex_out)
}
  
  
  

### --------------------- Specific Form (end) ------------------------###

#########################################################################



#########################################################################
################################# Tables ################################
#########################################################################

### set up working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Table 1:
### "Price and Expenditure Data for the Artificial Data Set"
# originial data: df_1
t <- 1:12
p_1 <- c(2,1.75,1.6,1.5,1.45,1.4,1.35,1.3,1.25,1.2,1.15,1.1)
p_2 <- c(1,0.5,1.05,1.10,1.12,1.15,1.18,0.6,1.2,1.25,1.28,1.3)
p_3 <- c(1,0.95,0.9,0.85,0.4,0.8,0.75,0.72,0.7,0.4,0.7,0.65)
p_4 <- c(0.5,0.55,0.6,0.65,0.7,0.75,0.7,0.65,0.7,0.75,0.75,0.8)
e <- c(10,13,11,12,15,13,14,17,15,18,16,17)
df_1 <- data.frame(t, p_1, p_2, p_3, p_4, e)
# data into a tidy form
p_list <- list(p_1, p_2, p_3, p_4)
df_temp <- lapply(1:4, function(x) {cbind(t, p_list[[x]], e)})
df_1_tidy <- do.call(rbind, df_temp)
p_type <- rep(1:4, each=12)
df_1_tidy <- data.frame(df_1_tidy, p_type)
colnames(df_1_tidy) <- c("t", "p", "e", "id")
df_1_tidy$id <- as.character(df_1_tidy$id)

### Table 2:
### "Differences at Period 12, D(t, sigma), between the Single Window
### CCDI Price Levels and the Linked CCDI Price Levels as Functions of
### the Linking Period t and the Elasticity of Substitution"
# price difference output: df_2
# parameter: sigma_vec, alpha_vec, link_list (linking position), len_val 
#     (window length)
sigma_vec <- c(0, 0.5, 1.001, 2, 4, 10)
alpha_vec <- c(0.2,0.2,0.2,0.4)
link_list <- as.list(1:10)
link_list[[11]] <- "mean"
len_val <- 11
p_diff_list <- lapply(sigma_vec, function(sigma_val) {
        df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                                   sigma=sigma_val, alpha=alpha_vec)[[2]]
        # sales adjusted
        q2_adjust <- (df_1_tidy_ces[,"t"] %in% c(3, 9)) &
                df_1_tidy_ces[,"id"]=="2"
        df_1_tidy_ces[q2_adjust, "q"] <- df_1_tidy_ces[q2_adjust, "q"]/2
        q3_adjust <- (df_1_tidy_ces[,"t"] %in% c(6, 11)) &
                df_1_tidy_ces[,"id"]=="3"
        df_1_tidy_ces[q3_adjust, "q"] <- df_1_tidy_ces[q3_adjust, "q"]/2
        # single window CCDI
        p_sgl <- pm_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t",
                            type="ccdi")
        p_sgl_tail <- tail(p_sgl, 1)[,-1]
        # rolling windows CCDI
        p_wd_tail <- lapply(link_list, function(link_val) {
                p_wd <- pm_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t",type="ccdi",
                         len=len_val, link=link_val)
                p_temp <- tail(p_wd, 1)[,-1]
                return(p_temp)
        })
        p_wd_tail <- unlist(p_wd_tail)
        p_diff <- p_wd_tail-p_sgl_tail
        return(p_diff)
})
df_2 <- data.frame(c(2:11,"mean"), do.call(cbind, p_diff_list))
colnames(df_2) <- c("t", paste("alpha=",sigma_vec, sep=""))

### Table 3:
### "Biases at Period 12, B(t, sigma), as Functions of the linking period
### and the Elasticity of Substitution"
# price difference output: df_3
# parameter: sigma_vec, alpha_vec, link_list (linking position), len_val 
#     (window length)
sigma_vec <- c(0, 0.5, 1.001, 2, 4, 10)
alpha_vec <- c(0.2,0.2,0.2,0.4)
link_list <- as.list(1:10)
link_list[[11]] <- "mean"
len_val <- 11
p_diff_list <- lapply(sigma_vec, function(sigma_val) {
        df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                                   sigma=sigma_val, alpha=alpha_vec)[[2]]
        # sales adjusted
        q2_adjust <- (df_1_tidy_ces[,"t"] %in% c(3, 9)) &
                df_1_tidy_ces[,"id"]=="2"
        df_1_tidy_ces[q2_adjust, "q"] <- df_1_tidy_ces[q2_adjust, "q"]/2
        q3_adjust <- (df_1_tidy_ces[,"t"] %in% c(6, 11)) &
                df_1_tidy_ces[,"id"]=="3"
        df_1_tidy_ces[q3_adjust, "q"] <- df_1_tidy_ces[q3_adjust, "q"]/2
        # single window CES
        p_sgl <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                           sigma=sigma_val, alpha=alpha_vec)[[1]]
        p_sgl_tail <- tail(p_sgl, 1)[,-1]
        # rolling windows CCDI
        p_wd_tail <- lapply(link_list, function(link_val) {
                p_wd <- pm_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t",type="ccdi",
                                 len=len_val, link=link_val)
                p_temp <- tail(p_wd, 1)[,-1]
                return(p_temp)
        })
        p_wd_tail <- unlist(p_wd_tail)
        p_diff <- p_wd_tail-p_sgl_tail
        return(p_diff)
})
df_3 <- data.frame(c(2:11,"mean"), do.call(cbind, p_diff_list))
colnames(df_3) <- c("t", paste("alpha=",sigma_vec, sep=""))

### Table 4:
### "Mean Absolute Differences in Percentage Points between CES and 
### Ten Approximating Indexes as Functions of the Elasticity of
### Substitution"
# price difference output: df_4
# parameters: sigma_vec, alpha_vec
sigma_vec <- c(0, 0.5, 1.001, 2, 4, 10)
alpha_vec <- c(0.2,0.2,0.2,0.4)
p_diff_mat <- lapply(sigma_vec, function(sigma_val) {
        # CES indexes
        p_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                           sigma=sigma_val, alpha=alpha_vec)[[1]]
        df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                                   sigma=sigma_val, alpha=alpha_vec)[[2]]
        # sales adjusted
        q2_adjust <- (df_1_tidy_ces[,"t"] %in% c(3, 9)) &
                df_1_tidy_ces[,"id"]=="2"
        df_1_tidy_ces[q2_adjust, "q"] <- df_1_tidy_ces[q2_adjust, "q"]/2
        q3_adjust <- (df_1_tidy_ces[,"t"] %in% c(6, 11)) &
                df_1_tidy_ces[,"id"]=="3"
        df_1_tidy_ces[q3_adjust, "q"] <- df_1_tidy_ces[q3_adjust, "q"]/2
        # chained Fisher and Tornqvist indexes
        p_ft_ch <- pb_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t", 
                            type=c("f","t"), seq=c("ch","fb"))[, c(1, 2, 4, 3, 5)]
        # GESK, CCDI, weighted time product dummy, and GK indexes
        p_multi <- pm_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t",
                            type=c("wtpd", "gk", "geks", "ccdi"))
        # AL indexes
        p_al <- sim_index(df_1_tidy_ces, "p", "q", "id", "t", type="f", sim="al")
        # LQ indexes
        p_lq <- sim_index(df_1_tidy_ces, "p", "q", "id", "t", type="t", sim="lq")
        # exclude periods 3, 6, 9, 11
        p_mat <- data.frame(p_ces, p_ft_ch, p_multi, p_al, p_lq)[c(-3,-6,-9,-11),c(-3,-8,-13,-15)]
        p_diff <- abs(p_mat[, -1] -p_mat[,"ces"])
        p_mean_diff <- apply(p_diff, 2, mean)[-1]
        # convert to percentage points
        p_mean_diff <- p_mean_diff*100
        return(p_mean_diff)
})
df_4 <- data.frame(sigma_vec, do.call(rbind, p_diff_mat))

### Table A1:
### "Price and Quantity Data for Two Products for Four Periods"
# orignial data: df_a1
mat_a1 <- matrix(c(1,2,3,4,1,0.5,1,1,1,1,1,1,10,5000,1,10,100,100,100,100),
                 nrow = 4)
colnames(mat_a1) <- c("t","p_1","p_2","q_1","q_2")
df_a1 <- as.data.frame(mat_a1)
# data into a tidy form
mat_a1_1 <- mat_a1[,c("t", "p_1", "q_1")]
mat_a1_1 <- cbind(mat_a1_1, rep(1, 4))
mat_a1_2 <- mat_a1[,c("t", "p_2", "q_2")]
mat_a1_2 <- cbind(mat_a1_2, rep(2, 4))
mat_a1_tidy <- rbind(mat_a1_1, mat_a1_2)
colnames(mat_a1_tidy) <- c("t", "p", "q", "id")
df_a1_tidy <- as.data.frame(mat_a1_tidy)
df_a1_tidy$id <- as.character(df_a1_tidy$id)

### Table A2:
### "Fixed Base and Chained Fisher, Tornqvist, Laspeyres and Paasche Indexes"
# price index output: df_a2
df_a2 <- pb_index(df_a1_tidy, p="p", qty="q", id="id", tm="t",
                  type=c("f", "t", "l", "p"), seq=c("fb", "ch"))[,c(1,2,6,8,3,5,7,9)]

### Table A3:
### " Asymptotic Linear Measures of Price Dissimilarity for sigma = 0"
# relative price dissimilarity output: df_a3
# Fisher base for asymptotic linear measures
# parameters: sigma_val, alpha_vec 
sigma_val <- 0
alpha_vec <- c(0.2,0.2,0.2,0.4)
df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                           sigma=sigma_val, alpha=alpha_vec)[[2]]
# sales adjusted
q2_adjust <- (df_1_tidy_ces[,"t"] %in% c(3, 9)) &
        df_1_tidy_ces[,"id"]=="2"
df_1_tidy_ces[q2_adjust, "q"] <- df_1_tidy_ces[q2_adjust, "q"]/2
q3_adjust <- (df_1_tidy_ces[,"t"] %in% c(6, 11)) &
        df_1_tidy_ces[,"id"]=="3"
df_1_tidy_ces[q3_adjust, "q"] <- df_1_tidy_ces[q3_adjust, "q"]/2
# measures output
df_a3 <- p_cmp(df_1_tidy_ces, "p", "q", "id", "t", type="f", sim="al")

### Table A4:
### "Log Quadratic Measures of Price Dissimilarity for sigma = 0"
# relative price dissimilarity output: df_a4
# Tornqvist base for log quadratic measures
# parameters: follow table a3 
# measures output
df_a4 <- p_cmp(df_1_tidy_ces, "p", "q", "id", "t", type="t", sim="lq")

### Table A5:
### "Alternative Price Levels for Different Methods and Elasticities of 
### Substitution"
# price index output: df_a5
# parameters: sigma_vec, alpha_vec
# sigma=1.001 is used for "sigma=1" case
sigma_vec <- c(0, 0.5, 1.001, 2, 4, 10, 20)
alpha_vec <- c(0.2,0.2,0.2,0.4)
# loop over sigma_val
df_a5 <- lapply(sigma_vec, function(sigma_val) {
        # CES price indexes
        p_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                           sigma=sigma_val, alpha=alpha_vec)[[1]]
        df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                                   sigma=sigma_val, alpha=alpha_vec)[[2]]
        # chained Fisher and Tornqvist indexes
        p_ft_ch <- pb_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t", 
                            type=c("f","t"), seq="ch")
        # GESK, CCDI, weighted time product dummy, and GK indexes
        p_multi <- pm_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t",
                            type=c("wtpd", "gk", "geks", "ccdi"))
        df_temp <- data.frame(p_ces, p_ft_ch, p_multi)[, c(-3, -6)]
        return(df_temp)
})
names(df_a5) <- paste("sigma = ", sigma_vec, sep="")

### Table A6:
### "Price Levels for the Sales Adjusted Data"
# price index output: df_a6
# parameters: sigma_vec, alpha_vec
# sigma=1.001 is used for "sigma=1" case
sigma_vec <- c(0, 0.5, 1.001, 2, 4, 10)
alpha_vec <- c(0.2,0.2,0.2,0.4)
# loop over sigma_val
df_a6 <- lapply(sigma_vec, function(sigma_val) {
        # CES price indexes
        p_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                           sigma=sigma_val, alpha=alpha_vec)[[1]]
        df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                                   sigma=sigma_val, alpha=alpha_vec)[[2]]
        # sales adjusted
        q2_adjust <- (df_1_tidy_ces[,"t"] %in% c(3, 9)) &
                df_1_tidy_ces[,"id"]=="2"
        df_1_tidy_ces[q2_adjust, "q"] <- df_1_tidy_ces[q2_adjust, "q"]/2
        q3_adjust <- (df_1_tidy_ces[,"t"] %in% c(6, 11)) &
                df_1_tidy_ces[,"id"]=="3"
        df_1_tidy_ces[q3_adjust, "q"] <- df_1_tidy_ces[q3_adjust, "q"]/2
        # chained Fisher and Tornqvist indexes
        p_ft_ch <- pb_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t", 
                            type=c("f","t"), seq=c("ch","fb"))[, c(1, 2, 4, 3, 5)]
        # GESK, CCDI, weighted time product dummy, and GK indexes
        p_multi <- pm_index(df_1_tidy_ces, p="p",qty="q",id="id",tm="t",
                            type=c("wtpd", "gk", "geks", "ccdi"))
        # AL indexes
        p_al <- sim_index(df_1_tidy_ces, "p", "q", "id", "t", type="f", sim="al")
        # LQ indexes
        p_lq <- sim_index(df_1_tidy_ces, "p", "q", "id", "t", type="t", sim="lq")
        # combine columns
        df_temp <- data.frame(p_ces, p_ft_ch, p_multi, p_al, p_lq)[, c(-3,-8,-13,-15)]
        return(df_temp)
})
names(df_a6) <- paste("sigma = ", sigma_vec, sep="")



#########################################################################



#########################################################################
################################ Figures ################################
#########################################################################

### Figure 1:
### "Alternative Price Levels for Different Methods and Elasticities of 
### Substitution"
# data source: df_a5
pdf("Figure_1.pdf")
op<-par(oma=c(4,0,0,0), # Room for the legend
        mfrow = c(2, 3))
sigma_name <- c(0, 0.5, 1, 2, 4, 10)
# loop over df_a5
for (i in 1:6) {
        plot(df_a5[[i]][,c("t","ccdi")], xlab="Period", 
             ylab="Index", main = bquote(sigma == .(sigma_name[i])),
             ylim=c(0.75, 1.6), type="l", col="gold",lwd=2, cex.axis=1.2, cex.main=2, cex.lab=1.5)
        lines(df_a5[[i]][,c("t","f_ch")], col=2, lwd=2)
        lines(df_a5[[i]][,c("t","t_ch")], col=3, lwd=2)
        lines(df_a5[[i]][,c("t","wtpd")], col=4, lwd=2)
        lines(df_a5[[i]][,c("t","gk")], col=6, lwd=2)
        lines(df_a5[[i]][,c("t","geks")], col="violetred", lwd=2)
        lines(df_a5[[i]][,c("t","ces")], col=1, lwd=2)
}
# add a title
# mtext("Figure 1: Alternative Price Levels for Different Methods and 
#      Elasticities of Substitution", side = 3, outer =
#              T, cex = 1)
legend_texts = c(
        as.expression(bquote(CES))
        , as.expression(bquote(FCH))
        , as.expression(bquote(TCH))
        , as.expression(bquote(WTPD))
        , as.expression(bquote(GK))
        , as.expression(bquote(GEKS))
        , as.expression(bquote(CCDI))
)
# add a legend
par(op) # Leave the last plot
op <- par(usr=c(0,1,0,1), # Reset the coordinates
          xpd=NA)         # Allow plotting outside the plot region
legend(-.1,-0.12 # Find suitable coordinates by trial and error
       , legend = legend_texts
       , col = c(1,2,3,4,6,"violetred","gold")
       , lty = 1, lwd=3, seg.len=0.9, x.intersp=0.1, cex=1.2, bty="n", horiz=T)
dev.off()

### Figure 2:
### "Sales Adjusted Data, Commodity 2"
# data source: df_1_tidy
# sigma = 2 as an illustration
sigma_val <- 2
df_1_tidy_ces <- ces_index(df_1_tidy, "p", "e", "id", "t", 
                           sigma=sigma_val, alpha=alpha_vec)[[2]]
# sales adjusted
q2_adjust <- (df_1_tidy_ces[,"t"] %in% c(3, 9)) &
        df_1_tidy_ces[,"id"]=="2"
df_1_tidy_ces[q2_adjust, "q"] <- df_1_tidy_ces[q2_adjust, "q"]/2
q3_adjust <- (df_1_tidy_ces[,"t"] %in% c(6, 11)) &
        df_1_tidy_ces[,"id"]=="3"
df_1_tidy_ces[q3_adjust, "q"] <- df_1_tidy_ces[q3_adjust, "q"]/2
# order and select commodity 2
df_1_tidy_ces <- df_1_tidy_ces[order(as.factor(df_1_tidy_ces[,"t"])),]
df_1_tidy_ces[, "val"] <- df_1_tidy_ces[, "p"]*df_1_tidy_ces[, "q"]
df_sales_2 <- df_1_tidy_ces[df_1_tidy_ces[,"id"]=="2", ]
# grid layout
layout(matrix(c(1,2,3),nrow=1), widths=c(1,1,1), heights=c(1,1,1), TRUE)
par(mar = c(4,4,2,1))
# loop over columns
col_names <- c("Price", "Quantity", "Expenditure")
ylab_names <- c("Index", "Index", "Value")
ylim_list <- list(c(0.5, 1.4), c(.58, 7.2), c(0.6, 4))
pdf("Figure_2.pdf")
layout(matrix(c(1,2,3),nrow=1), widths=c(1,1,1), heights=c(1,1,1), TRUE)
par(mar = c(4,4,2,1))
for (i in 1:3) {
        plot(df_sales_2[, c(1, i+2)], xlab="Period", 
             ylab=ylab_names[i], main = col_names[i],
             ylim=ylim_list[[i]], type="l", lwd=2, cex.main=1.5, cex.lab=1.2)
}
dev.off()

### Figure 3:
### "Alternative Price Levels for Sales Adjusted Data"
# data source: df_a6
pdf("Figure_3.pdf")
op<-par(oma=c(4,0,0,0), # Room for the legend
        mfrow = c(2, 3))
sigma_name <- c(0, 0.5, 1, 2, 4, 10)
# loop over df_a6
for (i in 1:6) {
  plot(df_a5[[i]][,c("t","ccdi")], xlab="Period", 
       ylab="Index", main = bquote(sigma == .(sigma_name[i])),
       ylim=c(0.75, 1.6), type="l", col="gold",lwd=2, cex.axis=1.2, cex.main=2, cex.lab=1.5)
  lines(df_a6[[i]][,c("t","f_ch")], col=2, lwd=2)
  lines(df_a6[[i]][,c("t","f_fb")], col="cadetblue", lwd=2)
  lines(df_a6[[i]][,c("t","wtpd")], col=4, lwd=2)
  lines(df_a6[[i]][,c("t","gk")], col=6, lwd=2)
  lines(df_a6[[i]][,c("t","geks")], col="violetred", lwd=2)
  lines(df_a6[[i]][,c("t","ces")], col=1, lwd=2)
}
# add a title
# mtext("Figure 3: Alternative Price Levels for Sales Adjusted Data", 
#      side = 3, outer = T, cex = 1)
legend_texts = c(
  as.expression(bquote(CES))
  , as.expression(bquote(FCH))
  , as.expression(bquote(FFB))
  , as.expression(bquote(WTPD))
  , as.expression(bquote(GK))
  , as.expression(bquote(GEKS))
  , as.expression(bquote(CCDI))
)

# add a legend
par(op) # Leave the last plot
op <- par(usr=c(0,1,0,1), # Reset the coordinates
          xpd=NA)         # Allow plotting outside the plot region
legend(-.1,-0.12 # Find suitable coordinates by trial and error
       , legend = legend_texts
       , col = c(1,2,"cadetblue",4,6,"violetred","gold")
       , lty = 1, lwd=3, seg.len=0.9, x.intersp=0.1, cex=1.2, bty="n", horiz=T)
dev.off()

#########################################################################


#########################################################################
############################ Export to Latex ############################
#########################################################################

### Table 1:
### "Price and Expenditure Data for the Artificial Data"
# general framework
tbl_1 <- xtable(df_1, caption = "Price and Expenditure Data for the 
       Artificial Data Set", floating=TRUE, latex.environments="center")
# column alignment
align(tbl_1) <- rep("r", 7)
# column names
names(tbl_1) <- c("$t$", "$p_{t1}$", "$p_{t2}$", "$p_{t3}$", 
                  "$p_{t4}$", "$e^{t}$")
# number of digits
digits(tbl_1) <- matrix(c(0, 0, 2, 2, 2, 2, 0), byrow = TRUE, ncol=7,
                        nrow=12)
# bold values (not compatible with number of digits)
# tbl_1[2,3] <- paste0("{\\bfseries", tbl_1[2,3], "}")
# suppress row names, column sanitization, caption position
addtorow <- list()
addtorow$pos <- list(nrow(df_1))
addtorow$command <- c("\\hline \n \\multicolumn{6}{c} {Note: Sales prices are in 
                      bold type.} \n")
tbl_1_latex <- capture.output(print(tbl_1, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top",
      add.to.row = addtorow,
      hline.after = c(-1, 0)))
write(tbl_1_latex, "tbl_1_latex.txt")

### Table 2:
### "Differences at Period 12, D(t, sigma), between the Single Window
### CCDI Price Levels and the Linked CCDI Price Levels as Functions of
### the Linking Period t and the Elasticity of Substitution"
# adjust characters
df_2[,1] <- c(as.character(2:11),"Mean")
# general framework
tbl_2 <- xtable(df_2, caption = "Differences at Period 12, $D(t, \\sigma)$, 
                between the Single Window CCDI Price Levels and the
                Linked CCDI Price Levels as Functions of the Linking 
                Period t and the Elasticity of Substitution", 
                floating=TRUE, latex.environments="center")
# column alignment
align(tbl_2) <- rep("r", 8)
# column names
names(tbl_2) <- c("$t$", "$D(t,0)$", "$D(t,0.5)$", "$D(t,1)$", 
                  "$D(t,2)$", "$D(t,4)$", "$D(t,10)$")
# number of digits
digits(tbl_2) <- matrix(c(0, 0, 5, 5, 5, 5, 5, 5), byrow = TRUE, ncol=8,
                        nrow=11)
# suppress row names, column sanitization, caption position
tbl_2_latex <- capture.output(print(tbl_2, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top",
      hline.after = c(-1, 0, 10, 11)))
write(tbl_2_latex, "tbl_2_latex.txt")

### Table 3:
### "Biases at Period 12, B(t, sigma), as Functions of the linking period
### and the Elasticity of Substitution"
# adjust characters
df_3[,1] <- c(as.character(2:11),"Mean")
# general framework
tbl_3 <- xtable(df_3, caption = "Biases at Period 12, $B(t, \\sigma)$, as
                Functions of the Linking Period t and the Elasticity of
                Substitution", 
                floating=TRUE, latex.environments="center")
# column alignment
align(tbl_3) <- rep("r", 8)
# column names
names(tbl_3) <- c("$t$", "$B(t,0)$", "$B(t,0.5)$", "$B(t,1)$", 
                  "$B(t,2)$", "$B(t,4)$", "$B(t,10)$")
# number of digits
digits(tbl_3) <- matrix(c(0, 0, 5, 5, 5, 5, 5, 5), byrow = TRUE, ncol=8,
                        nrow=11)
# suppress row names, column sanitization, caption position
tbl_3_latex <- capture.output(print(tbl_3, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top",
      hline.after = c(-1, 0, 10, 11)))
write(tbl_3_latex, "tbl_3_latex.txt")

### Table 4:
### "Mean Absolute Differences in Percentage Points between CES and 
### Ten Approximating Indexes as Functions of the Elasticity of
### Substitution"
# adjust characters
df_4[,1] <- as.character(c(0, 0.5, 1, 2, 4, 10))
# general framework
tbl_4 <- xtable(df_4, caption = "Mean Absolute Differences in 
                Percentage Points between CES and Ten Approximating
                Indexes as Functions of the Elasticity of Substitution", 
                floating=TRUE, latex.environments="center")
# column alignment
align(tbl_4) <- rep("r", 12)
# column names
names(tbl_4) <- c("$\\sigma$", "$B_{FCH}$", "$B_{TCH}$", "$B_{FFB}$", 
                  "$B_{TFB}$", "$B_{WTPD}$", "$B_{GK}$", "$B_{GEKS}$",
                  "$B_{CCDI}$", "$B_{AL}$", "$B_{LQ}$")
# number of digits
digits(tbl_4) <- matrix(c(0, 0, rep(2, 10)), byrow = TRUE, ncol=12,
                        nrow=6)
# suppress row names, column sanitization, caption position
tbl_4_latex <- capture.output(print(tbl_4, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top"))
write(tbl_4_latex, "tbl_4_latex.txt")

### Table A1:
### "Price and Quantity Data for Two Products for Four Periods"
# general framework
tbl_a1 <- xtable(df_a1, caption = "Price and Quantity Data for
                 Two Products for Four Periods", 
                floating=TRUE, latex.environments="center")
# column alignment
align(tbl_a1) <- rep("r", 6)
# column names
names(tbl_a1) <- c("$t$", "$p_{1}^{t}$", "$p_{2}^{t}$", "$q_{1}^{t}$", 
                  "$q_{2}^{t}$")
# number of digits
digits(tbl_a1) <- matrix(c(0, 0, 1,1,0,0), byrow = TRUE, ncol=6,
                        nrow=4)
# suppress row names, column sanitization, caption position
tbl_a1_latex <- capture.output(print(tbl_a1, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top"))
write(tbl_a1_latex, "tbl_a1_latex.txt")

### Table A2:
### "Fixed Base and Chained Fisher, Tornqvist, Laspeyres and Paasche Indexes"
# general framework
tbl_a2 <- xtable(df_a2, caption = "Fixed Base and Chained Fisher, T\\\"{o}rnqvist,
                 Laspeyres and Paasche Indexes", 
                 floating=TRUE, latex.environments="center")
# column alignment
align(tbl_a2) <- c("r","r", rep("l", 7))
# column names
names(tbl_a2) <- c("$t$", "$P_{F(FB)}$", "$P_{L(FB)}$", "$P_{P(FB)}$", 
                   "$P_{F(CH)}$", "$P_{T(CH)}$", "$P_{L(CH)}$", "$P_{P(CH)}$")
# number of digits
digits(tbl_a2) <- matrix(c(0, 0, rep(3,7)), byrow = TRUE, ncol=9,
                         nrow=4)
# suppress row names, column sanitization, caption position
tbl_a2_latex <- capture.output(print(tbl_a2, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top"))
write(tbl_a2_latex, "tbl_a2_latex.txt")

### Table A3:
### " Asymptotic Linear Measures of Price Dissimilarity for sigma = 0"
# general framework
tbl_a3 <- xtable(df_a3, caption = "Asymptotic Linear Measures of Price
                 Dissimilarity for $\\sigma$ = 0", 
                 floating=TRUE, latex.environments="center")
# column alignment
align(tbl_a3) <- rep("r", 14)
# column names
names(tbl_a3) <- c("$t$", paste("$\\Delta_{AL}^{", 1:12, "t}$", sep=""))
# number of digits
digits(tbl_a3) <- matrix(c(0, 0, rep(3,12)), byrow = TRUE, ncol=14,
                         nrow=12)
# suppress row names, column sanitization, caption position
tbl_a3_latex <- capture.output(print(tbl_a3, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top"))
write(tbl_a3_latex, "tbl_a3_latex.txt")

### Table A4:
### "Log Quadratic Measures of Price Dissimilarity for sigma = 0"
# general framework
tbl_a4 <- xtable(df_a4, caption = "Log Quadratic Measures of Price
                 Dissimilarity for $\\sigma$ = 0", 
                 floating=TRUE, latex.environments="center")
# column alignment
align(tbl_a4) <- rep("r", 14)
# column names
names(tbl_a4) <- c("$t$", paste("$\\Delta_{LQ}^{", 1:12, "t}$", sep=""))
# number of digits
digits(tbl_a4) <- matrix(c(0, 0, rep(3,12)), byrow = TRUE, ncol=14,
                         nrow=12)
# suppress row names, column sanitization, caption position
tbl_a4_latex <- capture.output(print(tbl_a4, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top"))
write(tbl_a4_latex, "tbl_a4_latex.txt")

### Table A5:
### "Alternative Price Levels for Different Methods and Elasticities of 
### Substitution"
# general framework
df_a5_list <- do.call(rbind, df_a5)
tbl_a5 <- xtable(df_a5_list, caption = "Alternative Price Levels for 
                 Different Methods and Elasticities of Substitution",
                 latex.environments="center")
# column alignment
align(tbl_a5) <- rep("r", 9)
# column names
names(tbl_a5) <- c("$t$", "$\\pi_{CES}^{t}$", "$\\pi_{FCH}^{t}$",
                   "$\\pi_{TCH}^{t}$", "$\\pi_{WTPD}^{t}$",
                   "$\\pi_{GK}^{t}$", "$\\pi_{GEKS}^{t}$",
                   "$\\pi_{CCDI}^{t}$")
# number of digits
digits(tbl_a5) <- c(0,0,rep(4,7))
# suppress row names, column sanitization, caption position
addtorow <- list()
addtorow$pos <- list(0, 12, 24, 36, 48, 60, 72)
addtorow$command <- paste("\\multicolumn{8}{c}{\\bfseries Alternative Price Levels when $\\sigma$ = ", 
                           c(0.5, 1, 2, 4, 10, 20),
                           "}\\\\ \n", sep="")
addtorow$command <- c("\\hline \n \\multicolumn{8}{c}{\\bfseries Alternative Price Levels when $\\sigma$ = 0}\\\\ \n",
                      addtorow$command)
tbl_a5_latex <- capture.output(print(tbl_a5, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top",
      tabular.environment = "longtable",
      floating = FALSE,
      hline.after=c(-1),
      add.to.row = addtorow))
write(tbl_a5_latex, "tbl_a5_latex.txt")

### Table A6:
### "Price Levels for the Sales Adjusted Data"
# general framework
df_a6_list <- do.call(rbind, df_a6)
tbl_a6 <- xtable(df_a6_list, caption = "Price Levels for the Sales
                 Adjusted Data",
                 latex.environments="center")
# column alignment
align(tbl_a6) <- rep("r", 13)
# column names
names(tbl_a6) <- c("$t$", "$\\pi_{CES}^{t}$", "$\\pi_{FCH}^{t}$",
                   "$\\pi_{TCH}^{t}$", "$\\pi_{FFB}^{t}$",
                   "$\\pi_{TFB}^{t}$", "$\\pi_{WTPD}^{t}$",
                   "$\\pi_{GK}^{t}$", "$\\pi_{GEKS}^{t}$",
                   "$\\pi_{CCDI}^{t}$", "$\\pi_{AL}^{t}$",
                   "$\\pi_{LQ}^{t}$")
# number of digits
digits(tbl_a6) <- c(0,0,rep(4,11))
# suppress row names, column sanitization, caption position
addtorow <- list()
addtorow$pos <- list(0, 12, 24, 36, 48, 60)
addtorow$command <- paste("\\multicolumn{12}{c}{\\bfseries Alternative Price Levels when $\\sigma$ = ", 
                          c(0.5, 1, 2, 4, 10),
                          "}\\\\ \n", sep="")
addtorow$command <- c("\\hline \n \\multicolumn{12}{c}{\\bfseries Alternative Price Levels when $\\sigma$ = 0}\\\\ \n",
                      addtorow$command)
tbl_a6_latex <- capture.output(print(tbl_a6, include.rownames = FALSE,
      sanitize.text.function=function(x){x},
      caption.placement = "top",
      tabular.environment = "longtable",
      floating = FALSE,
      hline.after=c(-1),
      add.to.row = addtorow))
write(tbl_a6_latex, "tbl_a6_latex.txt")

#########################################################################
