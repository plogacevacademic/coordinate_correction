


# trial = 03
# p <-
#   samples %>% subset(trial_id==trial) %>% ggplot(aes(gxR, -gyR)) + 
#   geom_point(color = "red", alpha = 0.1) +
#   geom_point(aes(gxL, -gyL), color = "blue", alpha = 0.1)
# p





## detect fixations ##
######################


library(saccades)

detect_fixations_binocular <- function(samples)
{
    fixations <- list()
    fixations$L <- samples %>% 
                      dplyr::select(time, trial=trial_id, x=gxL, y=gyL) %>%
                      saccades::detect.fixations()
    fixations$R <- samples %>%
                      dplyr::select(time, trial=trial_id, x=gxR, y=gyR) %>%
                      saccades::detect.fixations()
  
    iL <- with(fixations$L, Intervals(cbind(start, end)))
    iR <- with(fixations$R, Intervals(cbind(start, end)))
    
    overlap_LR <- interval_overlap(iL, iR) %>% sapply(function(x) { length(x) > 0 })
    overlap_RL <- interval_overlap(iR, iL) %>% sapply(function(x) { length(x) > 0 })
    
    list(L = fixations$L[overlap_LR,], R = fixations$R[overlap_RL,])
}


# to-do: Find out why some fixations have NAs for x and/or y coordinates
fixations <- detect_fixations_binocular(samples)
fixations$L %<>% subset(dur > 30) 
fixations$R %<>% subset(dur > 30) 



# The algorithm needs to handle the following problems:
# - Fixation sequence is shifted upward or downward (do some people maybe really read under or over the line - and if so, by how much can they? - and of so - does the probability of that happening depend on processing ease?)
#   -> But if the data need to corrected, by would I trust the x-coordinates, except maybe the relative order of coordinates?
# - Fixation sequence is shifted and distorted (ideally, linearly)
# => Both of those don't need any estimation of the SD, but rather a maximal range (the range of foveal vision) [though, ]

library(rethinking)

y_min <- min(aois[[cur_trial_id]]$y1)
y_max <- min(aois[[cur_trial_id]]$y2)
y_mid <- (y_min + y_max)/2
x_min <- min(aois[[cur_trial_id]]$x1)
x_max <- max(aois[[cur_trial_id]]$x2)
x_mid <- (x_min + x_max)/2

transform_par <- function(p) {
  names(p) <- c("y_delta", "y_drift", "y_buffer", "p_text")
  #print(p)
  
  p %<>% as.list()
  p$y_buffer %<>% exp()
  p$p_text <- plogis(p$p_text)/2 + 1/2
  p
}

# map_to_word_pos <- function(x, y, aoi, y_buffer)
# {
#   stopifnot(length(x) == 1 && length(y) == 1)
#   
#   is_in_aoi <- with(aoi, x > x1 & 
#                          x < x2 & 
#                          y > (y1 - y_buffer) &
#                          y < (y2 + y_buffer) )
#   
#   n_matching_positions = sum(is_in_aoi)
#   if (n_matching_positions == 0) {
#     return (NA)
#     
#   } else if (n_matching_positions == 1) {
#     return (aoi$word_no[is_in_aoi])
#     
#   } else {
#     stop("Overlapping AOIs. Fixation located in multiple AOIs.")
#   }
# }

p_fixation_on_pos <- function(x, plateau_width = 40) {
  normalizing_factor <- (function(x) 2*x - 1)( plogis(plateau_width/2) ) 
  (plogis(x+plateau_width/2) - plogis(x-plateau_width/2)) / normalizing_factor
}
# plot(p_fixation_on_line, xlim = c(-40, 40))


map_to_word_pos_prob <- function(x, y, aoi, y_width)
{
  stopifnot(length(x) == 1 && length(y) == 1)
  
  # to-do: modify the AOIs, in order to map fixation that happen to fall 'between' words onto something   
  
  is_in_aoi <- with(aoi, x > x1 & 
                         x < x2 )
  
  n_matching_positions = sum(is_in_aoi)
  if (n_matching_positions == 0) {
      return( c(NA, NA) )
  }
  
  if ( n_matching_positions > 1 ) {
    stop("Overlapping AOIs. Fixation located in multiple AOIs.")
  }
  
  idx_fixated_x_coord <- which(is_in_aoi)
  p = p_fixation_on_pos(y-aoi$y[idx_fixated_x_coord], y_width)

  return (c(aoi$word_no[is_in_aoi], p))
}

map_to_word_positions <- function(p, df, aoi)
{
    aoi$y <- with(aoi, (y1+y2)/2) 
    
    with(df, sapply(1:nrow(df), function(idx) { 
        map_to_word_pos(x = x[idx], 
                        y = y[idx] - p$y_delta - (x[idx]-x_mid) * p$y_drift, 
                        aoi, 
                        p$y_buffer) 
    }))
}


fn <- function(p, df, aoi) {
  p %<>% transform_par()
  
  df$word_no <- map_to_word_positions(p, df, aoi = aois[[cur_trial_id]])
  total <- df %>% group_by(word_no) %>% dplyr::summarise(total = sum(dur))
  total <- aoi %>% left_join(total, by = "word_no")
  
  logLik <- ifelse(is.na(total$total), 1-p$p_text, p$p_text) %>% log() %>% sum()
  #print(logLik)
  
  prior <- dnorm(p$y_delta, mean = 0, sd = 40, log = T) +
             dnorm(p$y_drift, mean = 0, sd = 40, log = T) +
             dnorm(p$y_buffer, mean = 0, sd = 40, log = T) +
             dbeta( (p$p_text-1/2)*2, 5, 1)
             #dnorm(p$p_text, mean = qlogis(.9), sd = 1, log = T)

  res <- sum(logLik) + sum(prior)
  
  print(res)
  res
}



start_par <- c(y_delta = 0, y_drift = 0,
               y_buffer = log(90),
               p_text = qlogis(.8))
p <- start_par %>% transform_par()



cur_trial_id = 13
df <- fixations$L %>% subset(trial == cur_trial_id) %>% subset(!is.na(x) | !is.na(y))
fn(start_par, df = df, aoi = aois[[cur_trial_id]])


# library(GA)
# res <- GA::ga(type = "real-valued", fn, 
#                df = df, aoi = aois[[cur_trial_id]],
#                lower = rep(-5, 4), 
#                upper = rep(5, 4)
#               #,
#                #pmutation = .2
#                )
# apply(res@solution, MARGIN = 1, function(x) transform_par(x) %>% t)
#
# p <- res@solution[1,] %>% transform_par
# df$word_no <- map_to_word_positions(p, df, aoi = aois[[cur_trial_id]])


res <- optim(start_par, fn, df = df, aoi = aois[[cur_trial_id]], method = "SANN",
             control = list(fnscale = -1, maxit = 500))
res


p <- res$par %>% transform_par()
df$word_no <- map_to_word_positions(p, df, aoi = aois[[cur_trial_id]])





df$word_no[is.na(df$word_no)] <- "NA"
df$fixation_id <- 1:nrow(df)

df %>% ggplot(aes(x, -y)) +
        geom_path(data = df, color="blue") +
        geom_text(aes(label = sprintf("%s [%d]", word_no, fixation_id) ), 
                  size = 4, color = "black", alpha = 1)
p + scale_x_continuous(limits = c(0, 1500)) +
    scale_y_continuous(limits = c(250,750))


# cbind(xL$y,
# fn(start_par, df = xL)
# )




P(seq of fixation on words and outside|seq of recorded fixations, AOIs) 

P(seq of fixation on words and outside|seq of recorded fixations, AOIs) 




error_pix = 3*18
cur_trial_id = 13
{
xL <- fixations$L %>% subset(trial == cur_trial_id) %>% subset(!is.na(x) | !is.na(y))
xR <- fixations$R %>% subset(trial == cur_trial_id) %>% subset(!is.na(x) | !is.na(y))


p <- xL %>% ggplot(aes(x, y)) +
        geom_point(data = xL, aes(x, y, size = dur), color = "blue", alpha = 0.5) + geom_path(data = xL, color="blue") 
        #+ # alpha = -start+min(start)
        #geom_point(data = xR, aes(x, y, size = dur), color = "red", alpha = 0.5) + geom_path(data = xR, color="red")
p <- p + with(aois[[cur_trial_id+1]], geom_rect(xmin=min(x1)-error_pix, xmax=max(x2)+error_pix, 
                                                ymin=min(y1)-error_pix, ymax=max(y2)+error_pix, alpha = 0.01))
p + scale_x_continuous(limits = c(0, 1500)) +
          scale_y_continuous(limits = c(250,750))
}


