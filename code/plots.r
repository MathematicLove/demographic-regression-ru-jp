SAVE_FILES <- TRUE
OUT_ROOT <- "plots_demography"   

open_dev <- function(filename=NULL, width=1920, height=1080) {  
  if (!is.null(filename) && SAVE_FILES) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    png(filename, width=width, height=height, res=120)  
    assign(".dev_is_png", TRUE, envir=.GlobalEnv)
  } else {
    sys <- tolower(Sys.info()[["sysname"]])
    if (grepl("windows", sys)) windows(width=16, height=9)      
    else if (grepl("darwin", sys)) quartz(width=16, height=9)   
    else x11(width=16, height=9)                                
    assign(".dev_is_png", FALSE, envir=.GlobalEnv)
  }
}
close_dev <- function() {
  if (exists(".dev_is_png", envir=.GlobalEnv) && get(".dev_is_png", envir=.GlobalEnv)) dev.off()
}

safe_name <- function(x) {
  x <- gsub(" ", "_", x); x <- gsub("ё","е",x); gsub("[^A-Za-zА-Яа-я0-9_]+","",x)
}
pretty_title <- function(country, var) {
  paste(country, "—", switch(var,
    "Общая_численность"="Общая численность",
    "Родившиеся"="Рождаемость",
    "Умершие"="Смертность",
    "Приток_иностранцев"="Приток иностранцев",
    var))
}

col_linear <- "#1f77b4"; col_quad <- "#ff7f0e"; col_exp <- "#2ca02c"; col_mean <- "#9467bd"; col_fact <- "#000000"

ru <- data.frame(
  Год = 2014:2024,
  Общая_численность = c(144025334,146743989,147580009,147182316,147797071,147840696,147959284,147455745,146980061,146447424,146028325),
  Родившиеся = c(1942683,1940579,1888729,1690307,1604344,1481074,1436514,1398253,1306200,1264938,1222408),
  Умершие    = c(1912347,1908541,1891015,1826125,1828910,1798307,2138586,2441594,1905800,1760200,1818635),
  Гражданство_получено = c(157800,209799,265319,257800,269400,497800,656300,735385,691045,378500,209000),
  Страна="Россия", stringsAsFactors=FALSE
)
jp <- data.frame(
  Год = 2014:2024,
  Общая_численность = c(127237150,127094745,126932772,126706210,126443180,126166948,126146099,125502000,124947000,124500000,123753040),
  Родившиеся = c(1003609,1005721,977242,946146,918400,865239,840832,831000,770747,727000,720988),
  Умершие    = c(1273025,1290510,1308158,1340567,1362470,1381093,1372648,1475000,1568961,1590503,1618684),
  Иностр_резиденты = c(3300000,2232189,2382822,2561848,2731093,2933137,2887116,2823565,3075213,3410992,3760000),
  Страна="Япония", stringsAsFactors=FALSE
)

ru$Естественный_прирост <- ru$Родившиеся - ru$Умершие
ru$Приток_иностранцев  <- ru$Гражданство_получено
ru$Интегральный_баланс  <- ru$Естественный_прирост + ru$Приток_иностранцев

jp$Приток_иностранцев  <- c(NA, diff(jp$Иностр_резиденты))
jp$Естественный_прирост <- jp$Родившиеся - jp$Умершие
jp$Интегральный_баланс  <- jp$Естественный_прирост + jp$Приток_иностранцев

safe_log <- function(y) ifelse(y<=0 | is.na(y), NA, log(y))

fit_models <- function(x_year, y) {
  df <- data.frame(t=x_year, y=y)
  m_lin  <- lm(y ~ t, data=df)
  m_quad <- lm(y ~ t + I(t^2), data=df)
  ly <- safe_log(y)
  if (all(is.na(ly))) { m_exp <- NULL; aic_exp <- NA; bic_exp <- NA
  } else { m_exp <- lm(ly ~ t, data=data.frame(t=x_year, ly=ly)); aic_exp <- AIC(m_exp); bic_exp <- BIC(m_exp) }
  m_mean <- lm(y ~ 1, data=df)
  aics <- c(AIC(m_lin), AIC(m_quad), aic_exp, AIC(m_mean))
  bics <- c(BIC(m_lin), BIC(m_quad), bic_exp, BIC(m_mean))
  names(aics) <- names(bics) <- c("Linear","Quadratic","Exponential","Mean")
  list(models=list(Linear=m_lin, Quadratic=m_quad, Exponential=m_exp, Mean=m_mean), AIC=aics, BIC=bics)
}
predict_model <- function(model, t_new, type=c("linear","quadratic","exponential","mean")){
  type <- match.arg(type)
  if (type=="exponential") { if (is.null(model)) return(rep(NA_real_, length(t_new))); return(exp(predict(model, newdata=data.frame(t=t_new)))) }
  if (type=="mean") return(rep(coef(model)[1], length(t_new)))
  predict(model, newdata=data.frame(t=t_new))
}
choose_winner <- function(aics,bics){
  a <- names(which.min(aics)); b <- names(which.min(bics)); v <- table(c(a,b)); list(AIC=a, BIC=b, Final=names(which.max(v)))
}

years_future <- 2025:2030
fit_and_forecast_series <- function(df, yname, drop_na=FALSE){
  d <- df[,c("Год",yname)]; names(d) <- c("Год","Y"); if (drop_na) d <- d[!is.na(d$Y), ]
  x <- d$Год; y <- d$Y
  fm <- fit_models(x,y); win <- choose_winner(fm$AIC,fm$BIC)
  fc <- data.frame(
    Год=years_future,
    Linear      = predict_model(fm$models$Linear,      years_future, "linear"),
    Quadratic   = predict_model(fm$models$Quadratic,   years_future, "quadratic"),
    Exponential = predict_model(fm$models$Exponential, years_future, "exponential"),
    Mean        = predict_model(fm$models$Mean,        years_future, "mean"),
    stringsAsFactors=FALSE)
  list(hist=list(t=x, y=y), models=fm$models, AIC=fm$AIC, BIC=fm$BIC, winners=win, forecast=fc)
}

targets <- c("Общая_численность","Родившиеся","Умершие","Приток_иностранцев")
RU <- JP <- list()
for (nm in targets) {
  RU[[nm]] <- fit_and_forecast_series(ru, nm, drop_na=FALSE)
  JP[[nm]] <- fit_and_forecast_series(jp, nm, drop_na=(nm=="Приток_иностранцев"))
}

build_choice_table <- function(res_list, country){
  out <- data.frame(Страна=character(), Показатель=character(), AIC_лучший=character(), BIC_лучший=character(), Итоговая_модель=character(), stringsAsFactors=FALSE)
  for (nm in names(res_list)) {
    w <- res_list[[nm]]$winners
    out <- rbind(out, data.frame(Страна=country, Показатель=nm, AIC_лучший=w$AIC, BIC_лучший=w$BIC, Итоговая_модель=w$Final))
  }
  out
}
choice_tbl <- rbind(build_choice_table(RU,"Россия"), build_choice_table(JP,"Япония"))
choice_tbl <- choice_tbl[order(choice_tbl$Страна, choice_tbl$Показатель), ]
cat("\n ВЫБОР МОДЕЛЕЙ (AIC/BIC) \n")
print(choice_tbl, row.names = FALSE)

get_final_forecast <- function(res_list){
  acc <- list()
  for (nm in names(res_list)) {
    fin <- res_list[[nm]]$winners$Final
    fc  <- res_list[[nm]]$forecast
    acc[[nm]] <- data.frame(Год=fc$Год, Показатель=nm, Значение=fc[[fin]])
  }
  do.call(rbind, acc)
}
to_wide <- function(df){
  u <- sort(unique(df$Год)); out <- data.frame(Год=u, stringsAsFactors=FALSE)
  for (p in unique(df$Показатель)) {
    sub <- df[df$Показатель==p, c("Год","Значение")]; names(sub) <- c("Год",p)
    out <- merge(out, sub, by="Год", all.x=TRUE)
  }
  out
}

ru_fc_all <- get_final_forecast(RU); jp_fc_all <- get_final_forecast(JP)
ru_wide <- to_wide(ru_fc_all); jp_wide <- to_wide(jp_fc_all)
ru_wide$Естественный_прирост <- ru_wide$Родившиеся - ru_wide$Умершие
ru_wide$Интегральный_баланс <- ru_wide$Естественный_прирост + ru_wide$Приток_иностранцев
jp_wide$Естественный_прирост <- jp_wide$Родившиеся - jp_wide$Умершие
jp_wide$Интегральный_баланс <- jp_wide$Естественный_прирост + jp_wide$Приток_иностранцев

plot_all_models_single <- function(res, title_text, ylab="чел.", filename=NULL) {
  xh <- res$hist$t; yh <- res$hist$y
  xr <- range(c(xh, years_future))
  yr <- range(c(yh, res$forecast$Linear, res$forecast$Quadratic, res$forecast$Exponential, res$forecast$Mean), na.rm=TRUE)
  open_dev(filename)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar=c(5,6,4,2), cex.axis=1.3, cex.lab=1.4, cex.main=1.5)  
  plot(xh, yh, type="b", pch=16, col=col_fact, lwd=3, xlim=xr, ylim=yr, 
       xlab="Год", ylab=ylab, main=paste0(title_text, " — все модели"), cex=1.2)
  grid(col="gray80", lty=2)  
  xx <- seq(min(xh), max(xh), by=1)
  lines(xx, predict(res$models$Linear,    newdata=data.frame(t=xx)), col=col_linear, lwd=4)
  lines(xx, predict(res$models$Quadratic, newdata=data.frame(t=xx)), col=col_quad,   lwd=4)
  if (!is.null(res$models$Exponential)) lines(xx, exp(predict(res$models$Exponential, newdata=data.frame(t=xx))), col=col_exp, lwd=4)
  lines(xx, rep(coef(res$models$Mean)[1], length(xx)), col=col_mean, lwd=4)
  points(res$forecast$Год, res$forecast$Linear, col=col_linear, pch=17, cex=1.5)
  points(res$forecast$Год, res$forecast$Quadratic, col=col_quad, pch=17, cex=1.5)
  if (!all(is.na(res$forecast$Exponential))) points(res$forecast$Год, res$forecast$Exponential, col=col_exp, pch=17, cex=1.5)
  points(res$forecast$Год, res$forecast$Mean, col=col_mean, pch=17, cex=1.5)
  legend("topleft", legend=c("Факт","Linear","Quadratic","Exponential","Mean"),
         col=c(col_fact,col_linear,col_quad,col_exp,col_mean), 
         lty=c(1,1,1,1,1), pch=c(16,NA,NA,NA,NA), lwd=c(3,4,4,4,4), 
         pt.cex=1.2, cex=1.3, bty="n", bg="white")
  close_dev()
}

plot_best_model_single <- function(res, title_text, ylab="чел.", filename=NULL) {
  xh <- res$hist$t; yh <- res$hist$y
  xr <- range(c(xh, years_future))
  yr <- range(c(yh, res$forecast$Linear, res$forecast$Quadratic, res$forecast$Exponential, res$forecast$Mean), na.rm=TRUE)
  fin <- res$winners$Final
  col_fin <- switch(fin, Linear=col_linear, Quadratic=col_quad, Exponential=col_exp, Mean=col_mean)
  open_dev(filename)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar=c(5,6,4,2), cex.axis=1.3, cex.lab=1.4, cex.main=1.5)  
  plot(xh, yh, type="b", pch=16, col=col_fact, lwd=3, xlim=xr, ylim=yr, 
       xlab="Год", ylab=ylab, main=paste0(title_text, " — лучшая по AIC/BIC: ", fin), cex=1.2)
  grid(col="gray80", lty=2)  
  xx <- seq(min(xh), max(xh), by=1)
  if (fin=="Linear") { 
    lines(xx, predict(res$models$Linear, newdata=data.frame(t=xx)), col=col_fin, lwd=6)
    lines(res$forecast$Год, res$forecast$Linear, col=col_fin, lwd=6, type="b", pch=17, cex=1.5)
  } else if (fin=="Quadratic") { 
    lines(xx, predict(res$models$Quadratic, newdata=data.frame(t=xx)), col=col_fin, lwd=6)
    lines(res$forecast$Год, res$forecast$Quadratic, col=col_fin, lwd=6, type="b", pch=17, cex=1.5)
  } else if (fin=="Exponential" && !is.null(res$models$Exponential)) { 
    lines(xx, exp(predict(res$models$Exponential, newdata=data.frame(t=xx))), col=col_fin, lwd=6)
    lines(res$forecast$Год, res$forecast$Exponential, col=col_fin, lwd=6, type="b", pch=17, cex=1.5)
  } else { 
    lines(xx, rep(coef(res$models$Mean)[1], length(xx)), col=col_fin, lwd=6)
    lines(res$forecast$Год, res$forecast$Mean, col=col_fin, lwd=6, type="b", pch=17, cex=1.5)
  }
  legend("topleft", legend=c("Факт", paste("Прогноз (", fin, ")")), 
         col=c(col_fact, col_fin), lty=1, lwd=c(3,6), pch=c(16,17), 
         pt.cex=1.2, cex=1.3, bty="n", bg="white")
  close_dev()
}

plot_balance_single <- function(hist_df, fore_df, title_text, filename=NULL) {
  rng <- range(c(hist_df$Интегральный_баланс, fore_df$Интегральный_баланс), na.rm=TRUE)
  open_dev(filename)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar=c(5,6,4,2), cex.axis=1.3, cex.lab=1.4, cex.main=1.5)  
  plot(hist_df$Год, hist_df$Интегральный_баланс, type="b", pch=16, col="#2E4057", lwd=4, 
       ylim=rng, xlab="Год", ylab="чел.", main=title_text, cex=1.3)
  grid(col="gray80", lty=2)  
  lines(fore_df$Год, fore_df$Интегральный_баланс, type="b", pch=17, col="#C44536", lwd=5, lty=1)
  abline(h=0, lty=3, col="gray40", lwd=2)
  legend("topleft", legend=c("Факт","Прогноз (по лучшим моделям)"), 
         col=c("#2E4057","#C44536"), pch=c(16,17), lty=c(1,1), lwd=c(4,5), 
         pt.cex=1.3, cex=1.3, bty="n", bg="white")
  close_dev()
}

var_labels <- c("Общая_численность","Родившиеся","Умершие","Приток_иностранцев")
for (v in var_labels) {
  fn1 <- file.path(OUT_ROOT, "all_models", paste0("ALLMODELS_Russia_", safe_name(v), ".png"))
  plot_all_models_single(RU[[v]], pretty_title("Россия", v), "чел.", if (SAVE_FILES) fn1 else NULL)

  fn2 <- file.path(OUT_ROOT, "all_models", paste0("ALLMODELS_Japan_", safe_name(v), ".png"))
  plot_all_models_single(JP[[v]], pretty_title("Япония", v), "чел.", if (SAVE_FILES) fn2 else NULL)

  fn3 <- file.path(OUT_ROOT, "best_models", paste0("BEST_Russia_", safe_name(v), ".png"))
  plot_best_model_single(RU[[v]], pretty_title("Россия", v), "чел.", if (SAVE_FILES) fn3 else NULL)

  fn4 <- file.path(OUT_ROOT, "best_models", paste0("BEST_Japan_", safe_name(v), ".png"))
  plot_best_model_single(JP[[v]], pretty_title("Япония", v), "чел.", if (SAVE_FILES) fn4 else NULL)
}

ru_hist_bal <- ru[, c("Год","Интегральный_баланс")]
jp_hist_bal <- jp[!is.na(jp$Интегральный_баланс), c("Год","Интегральный_баланс")]
ru_fore_bal <- ru_wide[, c("Год","Интегральный_баланс")]
jp_fore_bal <- jp_wide[, c("Год","Интегральный_баланс")]

fn5 <- file.path(OUT_ROOT, "balance", "BALANCE_Russia.png")
plot_balance_single(ru_hist_bal, ru_fore_bal, "Россия — Интегральный баланс (факт и прогноз)", if (SAVE_FILES) fn5 else NULL)

fn6 <- file.path(OUT_ROOT, "balance", "BALANCE_Japan.png")
plot_balance_single(jp_hist_bal, jp_fore_bal, "Япония — Интегральный баланс (факт и прогноз)", if (SAVE_FILES) fn6 else NULL)

if (SAVE_FILES) {
  cat("\nPNG сохранены в: ", normalizePath(OUT_ROOT), "\n", sep="")
} else {
  cat("\nГрафики выведены на экран. Чтобы СНИМАТЬ PNG-файлы, поставь SAVE_FILES <- TRUE вверху скрипта.\n")
}

if (!exists("SAVE_FILES")) SAVE_FILES <- TRUE
if (!exists("OUT_ROOT")) OUT_ROOT <- "plots_demography"
dir.create(file.path(OUT_ROOT, "hypotheses"), showWarnings = FALSE, recursive = TRUE)

future_bal <- rbind(
  data.frame(Страна="Россия", Год=ru_wide$Год,  Интегральный_баланс=ru_wide$Интегральный_баланс),
  data.frame(Страна="Япония", Год=jp_wide$Год, Интегральный_баланс=jp_wide$Интегральный_баланс)
)
future_bal_262030 <- future_bal[ future_bal$Год %in% 2026:2030, ]
future_bal_262030 <- future_bal_262030[ order(future_bal_262030$Год, future_bal_262030$Страна), ]

classify_sign <- function(x) ifelse(x > 0, "Улучшение", "Ухудшение")
sign_tbl <- future_bal_262030
sign_tbl$Статус <- classify_sign(sign_tbl$Интегральный_баланс)

cat("\n СТАТУСЫ (2026–2030)\n")
print(sign_tbl, row.names = FALSE)

years <- sort(unique(sign_tbl$Год))
combo <- character(length(years))
for (i in seq_along(years)) {
  y  <- years[i]
  stR <- sign_tbl$Статус[sign_tbl$Страна=="Россия" & sign_tbl$Год==y]
  stJ <- sign_tbl$Статус[sign_tbl$Страна=="Япония" & sign_tbl$Год==y]
  if (stR=="Улучшение" && stJ=="Улучшение") combo[i] <- "C1_оба_улучш"
  else if (stR=="Ухудшение" && stJ=="Ухудшение") combo[i] <- "C2_оба_ухудш"
  else if (stR=="Улучшение" && stJ=="Ухудшение") combo[i] <- "C3_РФ_улучш_JP_ухудш"
  else combo[i] <- "C4_JP_улучш_РФ_ухудш"
}

cats <- c("C1_оба_улучш","C2_оба_ухудш","C3_РФ_улучш_JP_ухудш","C4_JP_улучш_РФ_ухудш")
O <- sapply(cats, function(cn) sum(combo==cn))
obs_counts <- data.frame(Категория=cats, Наблюдено=O, stringsAsFactors = FALSE)

cat("\n частоты за 5 лет\n")
print(obs_counts, row.names = FALSE)

## 4) Сложные гипотезы (ожидания)
## H1: улучшение в обеих странах (каждый год)          -> (5,0,0,0)
## H2: ухудшение в обеих странах (каждый год)          -> (0,5,0,0)
## H3: РФ улучш, Япония ухудш (каждый год)             -> (0,0,5,0)
## H4: Япония улучш, РФ ухудш (каждый год)             -> (0,0,0,5)
E_H1 <- c(5,0,0,0); E_H2 <- c(0,5,0,0); E_H3 <- c(0,0,5,0); E_H4 <- c(0,0,0,5)

smooth_expect <- function(E) { E1 <- E + 0.5; E1 * 5 / sum(E1) }

H <- list(H1=E_H1, H2=E_H2, H3=E_H3, H4=E_H4)
chi_res <- data.frame(Гипотеза=character(), Chi2=numeric(), df=integer(), p_value=numeric(), stringsAsFactors = FALSE)

for (hn in names(H)) {
  E <- smooth_expect(H[[hn]])
  chi2 <- sum( (O - E)^2 / E )
  pval <- pchisq(chi2, df=3, lower.tail = FALSE)
  chi_res <- rbind(chi_res, data.frame(Гипотеза=hn, Chi2=chi2, df=3, p_value=pval))
}
chi_res <- chi_res[order(-chi_res$p_value), ]

cat("\n ПО ГИПОТЕЗАМ (от большего p-value) \n")
print(chi_res, row.names = FALSE)

best <- chi_res[1, , drop=FALSE]
best$Решение <- ifelse(best$p_value >= 0.05, "Не отвергается на уровне 5%", "Отвергается (но ближе к данным по min χ^2)")
print(best, row.names = FALSE)

plot_hypo_bars <- function() {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar=c(8,7,5,2), cex.axis=1.4, cex.lab=1.5, cex.main=1.6)  
  
  colors <- c("#4CAF50", "#F44336", "#2196F3", "#9C27B0")
  colors_alpha <- adjustcolor(colors, alpha.f = 0.8)
  
  bp <- barplot(height=O, 
                names.arg=c("C1: оба улучш","C2: оба ухудш","C3: РФ+ / JP-","C4: JP+ / РФ-"),
                col=colors_alpha, 
                las=2, 
                ylab="Частота за 2026–2030",
                main="Комбинации исходов демографического баланса",
                border=colors,
                lwd=2,
                ylim=c(0, max(O)*1.2))
  
  text(x=bp, y=O, labels=O, pos=3, cex=1.8, col=colors, font=2)
  
  legend_text <- paste0(chi_res$Гипотеза, ": p = ", 
                       formatC(chi_res$p_value, format="f", digits=4))
  
  legend("topright", bty="n", cex=1.4,
         legend=legend_text,
         fill=colors_alpha,
         border=colors,
         title="Гипотезы (p-value)")
  
  grid(NA, NULL, col="gray80", lty=3)
}

plot_hypo_bars()

if (SAVE_FILES) {
  fn <- file.path(OUT_ROOT, "hypotheses", "hypotheses_chi2_bars.png")
  png(fn, width=1920, height=1080, res=120); plot_hypo_bars(); dev.off() 
  cat("\nГрафик гипотез сохранён в: ", normalizePath(fn), "\n", sep="")
}