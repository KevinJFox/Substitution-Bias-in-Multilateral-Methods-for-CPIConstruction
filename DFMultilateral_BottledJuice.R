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


###### Exercises on scanner data: Bottled Juice
###### From Dominick's supermarket Database
###### https://www.chicagobooth.edu/research/kilts/datasets/dominicks


###### Outputs from this script: Results using Bottled Juice data
###### 
###### Table A7: df_a7
###### Latex code for making table: tbl_a7_latex.txt
###### Data for plotting using in Multilateral_Fig4.R:df_a7_bottled_juice.xlsx
######
###### Figures: 
###### Figure_bottled_juice: Plots FCH TCH WTPD GK GEKS CCDI
###### Figure_bottled_juice_sim: Plots FCH WTPD GK AL LQ CCDI


#########################################################################
########################### Loading Packages ############################
#########################################################################
library(rstudioapi) # extract path for setting up the working directory
library(dplyr) # create unit values by month
library(openxlsx)
library(xtable) # for creating tables


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
################# Running Functions on Real Data ########################
#########################################################################

### set up working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the data
load("bottled_juice.RData")

# Some preliminary cleaning up
# get rid of any products with zero prices, these are probably errors:
bottled_juice <- bottled_juice[!(bottled_juice$PRICE == 0),]
# in case there's some NA week observations that we can't use:
bottled_juice <- bottled_juice[!is.na(bottled_juice$WEEK),]
# in case there are any variables flagged as dodgy:
bottled_juice <- bottled_juice[bottled_juice$OK == 1,]
# sort the products by week
bottled_juice <- bottled_juice[order(bottled_juice$WEEK),]
# create a date variable. Week one finished 20/09/1989 according to the Dominicks manual
bottled_juice$date <- as.Date("13-09-1989","%d-%m-%Y") + 7*bottled_juice$WEEK
# create a month variable (using 4 weeks as a month)
bottled_juice$month <- (bottled_juice$WEEK-1) %/% 4 + 1

# To compute index numbers we need the month/quarter/year indices to be
# continuous with no gaps, so let's check if that's true
# If so, the following will return a value of "TRUE". 
all(min(bottled_juice$month):max(bottled_juice$month) %in% unique(bottled_juice$month))

# Recalling the definitions of the MOVE and QTY variables, we now create our Quantity variable
bottled_juice$Quantity <- bottled_juice$MOVE/bottled_juice$QTY

# Compute the unit values by month
bottled_juice_monthly <- bottled_juice %>% mutate(value=PRICE*Quantity) %>% group_by(month, UPC) %>%
        mutate(value_sum=sum(value), q_sum=sum(Quantity), unit_value=value_sum/q_sum) %>%
        filter(!duplicated(UPC)) %>% ungroup() %>% select(UPC, month, q_sum, unit_value) %>%
        as.data.frame()
# make a matched sample
count_UPC <- table(bottled_juice_monthly$UPC)
match_UPC <- count_UPC[count_UPC==100]
bottled_juice_monthly_match <- bottled_juice_monthly[as.character(bottled_juice_monthly$UPC) %in% names(match_UPC),]

### -------------------- Applications to Table --------------------- ###

### Table A7:
### "Alternative Price Levels for Different Methods, Bottled Juice"
# price index output: df_a7
p_ft_ch <- pb_index(bottled_juice_monthly_match, p="unit_value", qty="q_sum", id="UPC", tm="month", 
                    type=c("f","t"), seq=c("ch","fb"))[, c(1, 2, 4, 3, 5)]
# GEKS, CCDI, weighted time product dummy, and GK indexes
p_multi <- pm_index(bottled_juice_monthly_match, p="unit_value", qty="q_sum", id="UPC", tm="month",
                    type=c("wtpd", "gk", "geks", "ccdi"))
# AL indexes
p_al <- sim_index(bottled_juice_monthly_match, p="unit_value", qty="q_sum", id="UPC", tm="month", type="f", sim="al")
# LQ indexes
p_lq <- sim_index(bottled_juice_monthly_match, p="unit_value", qty="q_sum", id="UPC", tm="month", type="t", sim="lq")
# combine columns
df_a7 <- data.frame(p_ft_ch, p_multi, p_al, p_lq)[, c(-6,-11, -13)]


###--------------------- LaTex Tables ------------------------------- ###
### Table A7:
### "Alternative Price Levels for Different Methods, Bottled Juice"
# general framework
tbl_a7 <- xtable(df_a7, caption = "Alternative Price Levels for 
                 Different Methods, Bottled Juice",
                 latex.environments="center")
# column alignment
align(tbl_a7) <- rep("r", 12)
# column names
names(tbl_a7) <- c("$t$", "$\\pi_{FCH}^{t}$", "$\\pi_{TCH}^{t}$", 
                   "$\\pi_{FFB}^{t}$", "$\\pi_{TFB}^{t}$", 
                   "$\\pi_{WTPD}^{t}$", "$\\pi_{GK}^{t}$", "$\\pi_{GEKS}^{t}$",
                   "$\\pi_{CCDI}^{t}$", "$\\pi_{AL}^{t}$", "$\\pi_{LQ}^{t}$")
# number of digits
digits(tbl_a7) <- c(0,0,rep(4,10))
tbl_a7_latex <- capture.output(print(tbl_a7, include.rownames = FALSE,
                                     sanitize.text.function=function(x){x},
                                     caption.placement = "top",
                                     tabular.environment = "longtable",
                                     floating = FALSE,
                                     hline.after=c(-1)))
write(tbl_a7_latex, "tbl_a7_latex.txt")


### -------------------- Applications to Plots ---------------------- ###

# data source: df_a7
#%%%%%%%%%%
pdf("Figure_bottled_juice.pdf")
plot(df_a7[,c("month","f_ch")], ylim=c(0.8, 1.32),
     ylab="Index", xlab="Month", type="l", col=2, lwd=2, cex.main=1.5)
lines(df_a7[,c("month","t_ch")], col=3, lwd=2)
lines(df_a7[,c("month","wtpd")], col=4, lwd=2)
lines(df_a7[,c("month","gk")], col=6, lwd=2)
lines(df_a7[,c("month","geks")], col="violetred", lwd=2)
lines(df_a7[,c("month","ccdi")], col="gold", lwd=2)

legend_texts = c(as.expression(bquote(FCH))
                 , as.expression(bquote(TCH))
                 , as.expression(bquote(WTPD))
                 , as.expression(bquote(GK))
                 , as.expression(bquote(GEKS))
                 , as.expression(bquote(CCDI))
)

legend("bottom", legend = legend_texts
       , col = c(2:4,6, "violetred", "gold")
       , lty = 1, lwd=3, seg.len=0.9, cex=0.9, bty="n", horiz=T)
dev.off()
##########

### Similar to above, but plotting similarity indexes instead of GEKS and t_ch
pdf("Figure_bottled_juice_sim.pdf")
plot(df_a7[,c("month","f_ch")], col=2, ylim=c(0.8, 1.32),
     ylab="Index", xlab="Month", type="l", lwd=2, cex.main=1.42)
lines(df_a7[,c("month","wtpd")], col=4, lwd=2)
lines(df_a7[,c("month","gk")], col=6, lwd=2)
lines(df_a7[,c("month","al")], col="violetred", lwd=2)
lines(df_a7[,c("month","lq")], col=3, lwd=2)
lines(df_a7[,c("month","ccdi")], col="gold", lwd=2)


legend_texts = c(as.expression(bquote(FCH))
                 , as.expression(bquote(WTPD))
                 , as.expression(bquote(GK))
                 , as.expression(bquote(AL))
                 , as.expression(bquote(LQ))
                 , as.expression(bquote(CCDI))
)

legend("bottom", legend = legend_texts
       , col = c(2,4,6, "violetred", 3, "gold")
       , lty = 1, lwd=3, seg.len=0.9, cex=0.9, bty="n", horiz=T)
dev.off()
#########################################################################

### WRITE OUT DATA FOR USE IN A COMBINED PLOT: DFMultilateral_Fig4.R
write.xlsx(df_a7, "df_a7_bottled_juice.xlsx")

