dir_data <- '/Users/user/prj/dialect-loss'
setwd(dir_data)
dir_results <- '/Users/user/prj/dialect-loss/results_bootstrap/'
bootstrap_iterations <- 1000
do_bootstrap <- F
library('ggplot2')
library('lme4')
library(data.table)
yearmin <- 1800
yearmax <- 2000
speakers_table <- read.csv('speakers_needed.csv')
speakers <- speakers_table$speaker

variables <- c('var1_adj', 'var2_ae', 'var3_sh', 'var4_n', 
               'var5_sja_ek', 'var6_sja_v', 
               'var7_ot', 'var8_tu', 'var10_one', 'var11_ei', 'var12_kto')

read_table <- function(var) {
  filename <- paste(var,'csv', sep='.')
  table <- read.csv(filename)
  table <- table[table$speaker %in% speakers, ]
  table <- table[table$total > 0, ]
  table$prop <- table$cons/table$total
  table$age <- 2017 - table$year
  table$year20 <- table$year - 1920
  table
} 

unaggregate <- function(table) {
  cons <- as.vector(table$cons)
  inn <- as.vector(table$inn)
  data <- data.frame(speaker=rep(table$speaker, 2), year = rep(table$year, 2), 
                     gender = rep(table$gender, 2), cons = rep(c(1, 0), each=nrow(table)), 
                     value = c(cons, inn))
  data <- data[rep(seq_len(nrow(data)), data$value), 1:4]
  data$age <- 2017 - data$year
  data$year20 <- data$year - 1920
  data
}

reg <- function(data) {
  glmer(cons ~ year20 + (1|speaker), data = data, family='binomial')
}

bag <- function(table) {
  indicies = sample(1:nrow(table), replace = T, size=nrow(table))
  table[indicies,]
}

confband <- function(model) {
  vcov.m <- vcov(model)
  
  confband <- data.frame(year=seq(yearmin, yearmax, length.out = 10000))
  confband$year20 <- confband$year - 1920
  
  mm <- model.matrix(~ year20, confband)
  vars <- mm %*% vcov.m %*% t(mm)
  sds <- sqrt(diag(vars))
  z.val <- qnorm(1 - (1 - 0.95)/2)
  
  confband$pred_ <- predict(model, confband, type="link", re.form = NA)
  confband$lower_ <- confband$pred_ - z.val * sds
  confband$upper_ <- confband$pred_ + z.val * sds
  confband$pred <- plogis(confband$pred_)
  confband$lower <- plogis(confband$lower_)
  confband$upper <- plogis(confband$upper_)
  confband
}

midpoint_confint <- function(band) {
  list(left=band[min(which(band$lower_<0)), "year"],
       right=band[min(which(band$upper_<0)),"year"],
       midpoint=band[min(which(band$pred_<0)),"year"])
}

midpoints <- list()

for (var in variables) {
  
  table <- read_table(var)
  data <- unaggregate((table))
  
  fit_re <- reg(data)
  summary(fit_re)
  band <- confband(fit_re)
  if(do_bootstrap) {
    bootstrap_predicts <- vector("list", bootstrap_iterations)
    
    for (i in 1:bootstrap_iterations) {
      bs_data <- unaggregate(bag(table))
      bootstrap_predicts[[i]] <- confband(reg(bs_data), bs_data)$pred
      print(i)
    }
    bs_confint <- t(sapply(transpose(bootstrap_predicts), 
                           function(x){quantile(x, probs=c(0.05/2, 1-0.05/2), 
                                                na.rm=TRUE,
                                                type=8)}))
    colnames(bs_confint) <- c("bs_lower", "bs_upper")
    band <- cbind(band, bs_confint)
  }

  if (F) {
    fig <- ggplot() +
      geom_point(data=table, aes(x=year, y=prop, size=total), shape=21, fill='grey', color='black', alpha=0.6) +
      geom_line(size=1, alpha=0.8, aes(year, pred), data=band) +
      geom_ribbon(aes(ymin=lower, ymax=upper, x=year), alpha=0.3, data=band, 
                  color='steelblue2', fill='steelblue2') +
      theme_bw() +
      theme(plot.title = element_text(lineheight=1, size=12, family="serif",
                                      margin=margin(0,0,10,0)), 
            axis.text=element_text(size=10, family="serif"), 
            axis.title=element_text(size=12, family="serif"),
            axis.title.y=element_text(margin=margin(0,10,0,0)), 
            axis.title.x=element_text(margin=margin(10,0,0,0)), 
            legend.text=element_text(size=10, family="serif"),
            legend.title=element_text(size=12, family="serif"),
            plot.margin=unit(c(0.4,0.15,0.15,0.15), 'cm')) +
      scale_x_continuous('year of birth', breaks=seq(yearmin, yearmax, 10)) +
      scale_y_continuous('probability', breaks=seq(0,1,0.1), limits=c(0, 1)) +
      scale_size_continuous(range=c(1,8)) +
      guides(size = guide_legend(title = 'number of\nobservations'))
      if(do_bootstrap) {
        fig <- fig + geom_ribbon(aes(ymin=bs_lower, ymax=bs_upper, x=year), alpha=0.3, data=band, 
                    color='orange1', fill='orange1')
      }
    }
  midpoints[[var]] <- midpoint_confint(band)
  
#  ggsave(paste(dir_results, paste(var, "png", sep="."), sep='/'), fig, height=4, width = 7)
}

fig
