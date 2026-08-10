library(readxl)
library(stringi)
library(lavaan)
library(boot)
library(psych)
library(writexl)

rm(list = ls())
options(width = 340, scipen = 999, stringsAsFactors = FALSE)
set.seed(9186)

## ---- Paths ---------------------------------------------------------------
## Run this script with the repository root as the working directory.
## kok     : derived tables that ship with the repository (data/derived/).
## ham_kok : participant-level source files (analizlik.xlsx, "ham veriler.csv").
##           These are NOT distributed with the repository. Point SVT_RAW_DIR
##           at the folder that holds them, or edit the default below.
kok     <- "data/derived/"
ham_kok <- Sys.getenv("SVT_RAW_DIR", unset = "raw-data/")
if (!dir.exists(kok)) dir.create(kok, recursive = TRUE)

cizgi <- function(b) cat("\n", strrep("=", 130), "\n", b, "\n", strrep("=", 130), "\n", sep = "")

tr_from <- c("\u0131","\u0130","\u015F","\u015E","\u011F","\u011E","\u00E7","\u00C7",
             "\u00F6","\u00D6","\u00FC","\u00DC","\u00E2","\u00C2","\u00EE","\u00CE","\u00FB","\u00DB")
tr_to   <- c("i","I","s","S","g","G","c","C","o","O","u","U","a","A","i","I","u","U")

asciile <- function(x) {
  x <- stri_replace_all_regex(x, "[\\uFEFF\\u200B\\u00A0]", " ")
  x <- stri_replace_all_fixed(x, tr_from, tr_to, vectorize_all = FALSE)
  x <- stri_trans_general(x, "Latin-ASCII")
  stri_replace_all_regex(x, "[^\\x20-\\x7E\\t]", "")
}
stopifnot(asciile("\u015Eehir \u0130stanbul y\u0131l\u0131") == "Sehir Istanbul yili")

cizgi("BOLUM 1: VERI YUKLEME VE MADDE DUZEYI YENIDEN URETIM")

d <- as.data.frame(read_xlsx(paste0(ham_kok, "analizlik.xlsx")))
cat("analizlik.xlsx :", nrow(d), "x", ncol(d), "  eksik:", sum(is.na(d)), "\n")

sat <- asciile(readLines(paste0(ham_kok, "ham veriler.csv"), warn = FALSE, encoding = "UTF-8"))
gec <- file.path(tempdir(), "hv.csv")
writeLines(sat, gec, useBytes = TRUE)
n_alan <- max(count.fields(gec, sep = ",", quote = "\""))

h <- read.csv(gec, sep = ",", header = FALSE, skip = 1, check.names = FALSE,
              na.strings = c(""," ","NA","N/A","NULL","#NULL!"),
              col.names = paste0("V", seq_len(n_alan)))
pilot <- rowSums(!is.na(h[, 86:ncol(h), drop = FALSE])) > 0
cat("pilot satir dislandi :", sum(pilot), "\n")
h <- h[!pilot, 1:85, drop = FALSE]
kr <- vapply(h, is.character, logical(1))
h[kr] <- lapply(h[kr], function(x) stri_trim_both(stri_replace_all_regex(x, "\\s+", " ")))

kisa <- c("ZamanDamgasi","YasHam","Cinsiyet","Egitim","MedeniDurum","Gelir","Meslek",
          paste0("ERS",1:4), paste0("KAB",1:8), paste0("UYG",1:8), paste0("ULS",1:3),
          paste0("GUV",1:17), paste0("MEM",1:26), paste0("SOY",1:12))
stopifnot(length(kisa) == ncol(h))
names(h) <- kisa

h5a <- c("Kesinlikle Katilmiyorum"=1,"Katilmiyorum"=2,"Kararsizim"=3,"Katiliyorum"=4,"Kesinlikle Katiliyorum"=5)
h5b <- c("Hic Katilmiyorum"=1,"Katilmiyorum"=2,"Kararsizim"=3,"Katiliyorum"=4,"Tamamen Katiliyorum"=5)
h5c <- c("Hic Katilmiyorum"=1,"Katilmiyorum"=2,"Ne Katiliyorum Ne Katilmiyorum"=3,"Katiliyorum"=4,"Tamamen Katiliyorum"=5)
h4  <- c("Cok zor"=1,"Oldukca zor"=2,"Oldukca kolay"=3,"Cok kolay"=4)
er_m <- c(paste0("ERS",1:4), paste0("KAB",1:8), paste0("UYG",1:8), paste0("ULS",1:3))
gv_m <- paste0("GUV",1:17); mm_m <- paste0("MEM",1:26); sy_m <- paste0("SOY",1:12)
for (a in er_m) h[[a]] <- unname(h5a[h[[a]]])
for (a in gv_m) h[[a]] <- unname(h5b[h[[a]]])
for (a in mm_m) h[[a]] <- unname(h5c[h[[a]]])
for (a in sy_m) h[[a]] <- unname(h4[h[[a]]])
madde_hepsi <- c(er_m, gv_m, mm_m, sy_m)
cat("cevrim NA :", sum(is.na(h[, madde_hepsi])), "\n")

yas <- suppressWarnings(as.numeric(h$YasHam))
yas[h$YasHam == "15.12.2002"] <- 22
yas[h$YasHam == "1985"] <- 40
h <- h[!is.na(yas) & yas >= 18, ]
stopifnot(nrow(h) == nrow(d))

M <- as.matrix(h[, madde_hepsi])
h$Longstring <- apply(M, 1, function(r) max(rle(as.vector(r))$lengths))
h$IRV <- apply(M, 1, sd)
h$FarkliSecenek <- apply(M, 1, function(r) length(unique(r)))
h$DuzBlok <- rowSums(vapply(list(er_m, gv_m, mm_m, sy_m), function(a)
  apply(h[, a, drop = FALSE], 1, function(r) sd(r) == 0), logical(nrow(h))))
h$Kabullenici <- rowMeans(h[, c("GUV4","GUV15")]) - rowMeans(h[, setdiff(gv_m, c("GUV4","GUV15"))])
h$GUV4 <- 6 - h$GUV4; h$GUV15 <- 6 - h$GUV15

stil <- c("Longstring","IRV","Kabullenici","FarkliSecenek","DuzBlok")
for (v in c(stil, madde_hepsi)) d[[v]] <- h[[v]]

cat("\nYENIDEN URETIM KONTROLU (hepsi 1 olmali)\n")
kon <- c(ERISIM = cor(rowMeans(d[, er_m]), d$ERISIM),
         MEMNUNIYET = cor(rowMeans(d[, mm_m]), d$MEMNUNIYET),
         GUVEN = cor(rowSums(d[, gv_m]), d$GUVEN),
         SOY = cor(rowMeans(d[, sy_m]), d$SOY_ort),
         MTT = cor(rowMeans(d[, paste0("MEM",1:9)]), d$MTT),
         ULS = cor(rowMeans(d[, paste0("ULS",1:3)]), d$ULS))
print(round(kon, 6))
stopifnot(all(kon > 0.999999))

d_tam <- d

cizgi("BOLUM 2: KORELASYON YAPISI VE ORTAK YONTEM VARYANSI")

skor <- c("ERS","KAB","UYG","ULS","ERISIM","GUV_Saglayici","GUV_Odeyici","GUVEN",
          "MTT","FCV","RNS","ILT","MEMNUNIYET","SOY_index","Yas")
cat("Spearman korelasyon matrisi:\n")
print(round(cor(d[, skor], method = "spearman"), 3))

cat("\nORTAK YONTEM VARYANSI GOSTERGELERI (tum orneklem)\n")
cmb_gost <- function(i) {
  eE <- eigen(cor(d[i, er_m]))$values
  eM <- eigen(cor(d[i, mm_m]))$values
  eG <- eigen(cor(d[i, gv_m]))$values
  data.frame(n = sum(i),
             alfa_E = round(psych::alpha(d[i, er_m], warnings=FALSE)$total$raw_alpha, 3),
             alfa_M = round(psych::alpha(d[i, mm_m], warnings=FALSE)$total$raw_alpha, 3),
             alfa_G = round(psych::alpha(d[i, gv_m], warnings=FALSE)$total$raw_alpha, 3),
             ev1ev2_E = round(eE[1]/eE[2], 1), ev1ev2_M = round(eM[1]/eM[2], 1),
             ilkfak_E = round(100*eE[1]/length(eE), 1), ilkfak_M = round(100*eM[1]/length(eM), 1),
             harman_ilkfak = round(100*eigen(cor(d[i, c(er_m, mm_m, gv_m, sy_m)]))$values[1]/78, 1),
             longstring = round(mean(d$Longstring[i]), 1), duzblok = round(mean(d$DuzBlok[i]), 2))
}
print(cmb_gost(rep(TRUE, nrow(d))), row.names = FALSE, right = FALSE)

cat("\nHarman tek faktor testi (78 madde, tek faktorlu EFA):\n")
hf <- eigen(cor(d[, madde_hepsi]))$values
cat("ilk ozdeger :", round(hf[1], 2), "  ikinci :", round(hf[2], 2),
    "  oran :", round(hf[1]/hf[2], 1), "  ilk faktor varyansi : %",
    round(100*hf[1]/78, 1), "\n", sep = "")

cizgi("BOLUM 3: SVT MAKINESI")

gost <- list(ERISIM = c("ERS","KAB","UYG","ULS"), MEMNUNIYET = c("MTT","FCV","RNS","ILT"))
aile <- c(ERS="E", KAB="E", UYG="E", ULS="E", ERISIM="E",
          MTT="M", FCV="M", RNS="M", ILT="M", MEMNUNIYET="M",
          GUV_Saglayici="G", GUV_Odeyici="G", GUVEN="G",
          SOY_index="S", Yas="Y", setNames(rep("T", length(stil)), stil))
tum_mod <- c("Cinsiyet","YasGrup","EgitimGrup","MedeniDurum","Gelir","MeslekGrup")
mod_sec <- function(z) if (z == "Yas") setdiff(tum_mod, "YasGrup") else tum_mod
esik <- c(cfi = 0.010, tli = 0.010, rmsea = 0.015, srmr = 0.030)

fit4 <- function(mdl, dat, g, esitlik = NULL)
  tryCatch(cfa(mdl, data = dat, group = g, group.equal = esitlik,
               std.lv = TRUE, estimator = "MLR", missing = "fiml"),
           error = function(e) NULL, warning = function(w) NULL)

fm4 <- function(f) if (is.null(f)) rep(NA, 4) else
  as.numeric(fitMeasures(f, c("cfi.scaled","tli.scaled","rmsea.scaled","srmr")))

mi_ham <- function(x, g) {
  mdl <- paste0("F =~ ", paste(gost[[x]], collapse = " + "))
  dd <- d; dd[[g]] <- factor(dd[[g]])
  f0 <- fit4(mdl, dd, g)
  f1 <- fit4(mdl, dd, g, "loadings")
  f2 <- fit4(mdl, dd, g, c("loadings","intercepts"))
  f3 <- fit4(mdl, dd, g, c("loadings","intercepts","residuals"))
  if (any(vapply(list(f0,f1,f2,f3), is.null, logical(1)))) return(NULL)
  m <- rbind(fm4(f0), fm4(f1), fm4(f2), fm4(f3))
  colnames(m) <- c("cfi","tli","rmsea","srmr")
  dl <- rbind(
    metric = c(m[1,1]-m[2,1], m[1,2]-m[2,2], m[2,3]-m[1,3], m[2,4]-m[1,4]),
    scalar = c(m[2,1]-m[3,1], m[2,2]-m[3,2], m[3,3]-m[2,3], m[3,4]-m[2,4]),
    strict = c(m[3,1]-m[4,1], m[3,2]-m[4,2], m[4,3]-m[3,3], m[4,4]-m[3,4]))
  colnames(dl) <- names(esik)
  if (any(is.na(dl))) return(NULL)
  list(delta = dl, fitleri = m,
       invaryant = all(dl <= matrix(esik, nrow = 3, ncol = 4, byrow = TRUE)))
}

adaptif_agirlik <- function(x, modl) {
  D <- do.call(rbind, lapply(modl, function(g) {
    r <- mi_cache[[paste(x, g, sep = "|")]]
    if (is.null(r)) return(NULL)
    apply(r$delta, 2, function(v) max(pmax(0, v))) }))
  inv <- vapply(modl, function(g) {
    r <- mi_cache[[paste(x, g, sep = "|")]]
    if (is.null(r)) NA else r$invaryant }, logical(1))
  if (is.null(D) || nrow(D) < 3) return(c(0.40, 0.10, 0.30, 0.20))
  N <- pmin(sweep(D, 2, esik, "/"), 1)
  VS <- 1/(1 + apply(N, 2, function(v) sd(v)/(mean(v)+0.001)))
  R <- suppressWarnings(cor(N)); R[is.na(R)] <- 0
  RP <- 1/(1 + apply(abs(R), 1, function(v) mean(v[v < 1])))
  DP <- vapply(seq_len(ncol(N)), function(j)
    if (sum(!inv, na.rm=TRUE) > 0 && sum(inv, na.rm=TRUE) > 0)
      max(0, mean(N[!inv, j]) - mean(N[inv, j])) else sd(N[, j]), numeric(1))
  if (all(DP == 0)) DP <- apply(N, 2, sd)
  w <- DP * RP * VS
  if (sum(w) == 0) return(c(0.40, 0.10, 0.30, 0.20))
  as.numeric(w/sum(w))
}

mi_skor <- function(x, g, w) {
  r <- mi_cache[[paste(x, g, sep = "|")]]
  if (is.null(r)) return(NULL)
  nrm <- pmin(apply(r$delta, 2, function(v) max(pmax(0, v))) / esik, 1)
  list(skor = sum(w * nrm), invaryant = r$invaryant)
}

stab_metrik <- function(x, y, z, g, min_n = 25) {
  ind <- gost[[x]]
  m0 <- paste0("F =~ ", paste(ind, collapse=" + "), "\n", y, " ~ c*F")
  m1 <- paste0("F =~ ", paste(ind, collapse=" + "), "\n", z, " ~ a*F\n", y, " ~ cprime*F + b*", z)
  m2 <- paste0("F =~ ", paste(ind, collapse=" + "), "\n", z, " ~ xz*F")
  dd <- d[!is.na(d[[g]]), ]
  b0 <- b1 <- b2 <- nn <- numeric(0); ad <- character(0)
  for (l in levels(factor(dd[[g]]))) {
    sub <- dd[dd[[g]] == l, ]
    if (nrow(sub) < min_n) next
    f0 <- tryCatch(sem(m0, sub, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    f1 <- tryCatch(sem(m1, sub, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    f2 <- tryCatch(sem(m2, sub, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(f0) || is.null(f1) || is.null(f2)) next
    q0 <- standardizedSolution(f0); q1 <- standardizedSolution(f1); q2 <- standardizedSolution(f2)
    v0 <- q0$est.std[q0$label == "c"][1]
    v1 <- q1$est.std[q1$label == "cprime"][1]
    v2 <- q2$est.std[q2$label == "xz"][1]
    if (any(is.na(c(v0,v1,v2)))) next
    b0 <- c(b0,v0); b1 <- c(b1,v1); b2 <- c(b2,v2); nn <- c(nn,nrow(sub)); ad <- c(ad,l)
  }
  if (length(b0) < 2) return(NULL)
  wm <- function(v) weighted.mean(v, nn)
  ws <- function(v) sqrt(sum(nn*(v-wm(v))^2)/sum(nn))
  eps <- 1e-4; e <- 1e-8
  mu0 <- wm(b0); mu1 <- wm(b1); mu2 <- wm(b2)
  sd0 <- ws(b0); sd1 <- ws(b1); sd2 <- ws(b2)
  cv0 <- 100*sd0/max(abs(mu0), eps)
  cv1 <- 100*sd1/max(abs(mu1), eps)
  cv2 <- 100*sd2/max(abs(mu2), eps)
  dlog <- log(cv0+e) - log(cv1+e)
  dmu <- log(abs(mu1)+e) - log(abs(mu0)+e)
  dsg <- log(sd0+e) - log(sd1+e)
  list(delta_log = dlog, cv0 = cv0, cv1 = cv1, cv_xz = cv2,
       pct = 100*(1-exp(-dlog)),
       ocr = mean(as.integer(abs(b1-mu1) < abs(b0-mu0))),
       os = max(0,dmu)/(max(0,dmu)+max(0,dsg)+e),
       dmu = dmu, dsg = dsg, K = length(b0), n = sum(nn),
       gruplar = data.frame(grup = ad, n = nn, b0 = b0, b1 = b1))
}

svt_tam <- function(x, y, z, R = 1000, ayrinti = FALSE) {
  modl <- mod_sec(z)
  w <- adaptif_agirlik(x, modl)
  res <- list()
  for (g in modl) {
    mi <- mi_skor(x, g, w); st <- stab_metrik(x, y, z, g)
    if (is.null(mi) || is.null(st)) next
    res[[g]] <- c(list(moderator = g, mi_skor = mi$skor, invaryant = mi$invaryant), st)
  }
  if (length(res) < 3) return(NULL)
  dv <- vapply(res, function(r) r$delta_log, numeric(1))
  wv <- vapply(res, function(r) r$n, numeric(1))
  ov <- vapply(res, function(r) r$ocr, numeric(1))
  sv <- vapply(res, function(r) r$os, numeric(1))
  mv <- vapply(res, function(r) r$mi_skor, numeric(1))
  iv <- vapply(res, function(r) r$invaryant, logical(1))
  cv0v <- vapply(res, function(r) r$cv0, numeric(1))
  cv1v <- vapply(res, function(r) r$cv1, numeric(1))
  bo <- boot(dv, function(dd, idx) sum(dd[idx] * wv[idx]/sum(wv[idx])), R = R)
  md <- bo$t0; se <- sd(bo$t)
  pb <- pnorm(md/se, lower.tail = FALSE)
  ci <- tryCatch(boot.ci(bo, type = "perc")$percent[4:5], error = function(e) c(NA,NA))
  wvar <- sum(wv*(dv-md)^2)/sum(wv)
  dc <- if (is.finite(wvar) && wvar > 0) md/sqrt(wvar) else NA
  kr2 <- (dv > 0) | (ov >= 0.5)
  pbin <- binom.test(sum(kr2), length(kr2), 0.5, "greater")$p.value
  mos <- weighted.mean(sv, wv)
  k1 <- (pb < 0.05 && md > 0); k2 <- (pbin < 0.05)
  out <- data.frame(
    X = x, Y = y, Z = z, K = length(res),
    delta_l = round(md, 3), se = round(se, 3), d = round(dc, 3),
    CI_alt = round(ci[1], 3), CI_ust = round(ci[2], 3),
    p_boot = round(pb, 4), p_binom = round(pbin, 4),
    poz_dl = sum(dv > 0), ocr50 = sum(ov >= 0.5), kriter2 = sum(kr2),
    CV0 = round(weighted.mean(cv0v, wv), 1), CV1 = round(weighted.mean(cv1v, wv), 1),
    OS = round(mos, 3), OCR = round(mean(ov), 3),
    mi_ihlal = sum(!iv), r_MI_dl = round(suppressWarnings(cor(mv, dv)), 3),
    boot_gecti = k1, binom_gecti = k2,
    karar = ifelse(k1 && k2, "STABILIZATOR", "Degil"),
    mekanizma = if (!(k1 && k2)) "Yok" else if (mos < 0.3) "Type A" else if (mos > 0.7) "Type B" else "Type AB",
    stringsAsFactors = FALSE)
  attr(out, "dv") <- dv; attr(out, "ov") <- ov
  attr(out, "cv0v") <- cv0v; attr(out, "cv1v") <- cv1v
  attr(out, "modad") <- names(res)
  if (ayrinti) attr(out, "detay") <- res
  out
}

grid_calistir <- function(gr, etiket) {
  cat("\n>>>", etiket, ":", nrow(gr), "uclu\n")
  t0 <- Sys.time(); tut <- list(); mod_dl <- list()
  for (i in seq_len(nrow(gr))) {
    r <- tryCatch(svt_tam(gr$X[i], gr$Y[i], gr$Z[i]), error = function(e) NULL)
    if (!is.null(r)) {
      tut[[length(tut)+1]] <- r
      mod_dl[[length(mod_dl)+1]] <- data.frame(
        etiket = etiket, X = r$X, Y = r$Y, Z = r$Z,
        moderator = attr(r, "modad"), dl = attr(r, "dv"), ocr = attr(r, "ov"),
        cv0 = attr(r, "cv0v"), cv1 = attr(r, "cv1v"), stringsAsFactors = FALSE)
    }
    if (i %% 20 == 0) cat("   ", i, "/", nrow(gr), " ",
                          round(difftime(Sys.time(), t0, units="min"),1), "dk\n")
  }
  S <- do.call(rbind, tut)
  S$etiket <- etiket
  S$r_XZ <- round(mapply(function(a,b) cor(d[[a]], d[[b]], method="spearman"), S$X, S$Z), 3)
  S$r_ZY <- round(mapply(function(a,b) cor(d[[a]], d[[b]], method="spearman"), S$Z, S$Y), 3)
  S$maxZ <- round(pmax(abs(S$r_XZ), abs(S$r_ZY)), 3)
  S$etk_p <- round(mapply(function(x,y,z)
    summary(lm(scale(d[[y]]) ~ scale(d[[x]]) * scale(d[[z]])))$coefficients[4,4], S$X, S$Y, S$Z), 4)
  S$C2 <- ifelse(S$maxZ < .30 & S$etk_p > .05, "gecti", "kaldi")
  S$taban_ok <- S$CV0 >= 10
  cat("   tamamlandi:", round(difftime(Sys.time(), t0, units="min"),1), "dk\n")
  list(S = S, mod = do.call(rbind, mod_dl))
}

mi_kur <- function() {
  mc <- list()
  for (x in names(gost)) for (g in tum_mod) {
    k <- paste(x, g, sep = "|")
    mc[[k]] <- mi_ham(x, g)
    cat(sprintf("  %-12s %-12s -> %s\n", x, g,
                ifelse(is.null(mc[[k]]), "BASARISIZ", ifelse(mc[[k]]$invaryant, "invaryant", "IHLAL"))))
  }
  mc
}

cizgi("BOLUM 4: TAM ORNEKLEM MI DURUMU")
mi_cache <- mi_kur()

cizgi("BOLUM 5: GRID 1 - OLCEK SKORLARI Z OLARAK")
sonuc_ad <- c("MTT","FCV","RNS","ILT","MEMNUNIYET","ERS","KAB","UYG","ULS","ERISIM",
              "GUV_Saglayici","GUV_Odeyici","GUVEN","SOY_index")
stab_ad <- c("ERS","KAB","UYG","ULS","ERISIM","MTT","FCV","RNS","ILT","MEMNUNIYET",
             "GUV_Saglayici","GUV_Odeyici","GUVEN","SOY_index","Yas")
g1 <- expand.grid(X = names(gost), Y = sonuc_ad, Z = stab_ad, stringsAsFactors = FALSE)
g1 <- g1[aile[g1$X] != aile[g1$Y] & aile[g1$X] != aile[g1$Z] & aile[g1$Y] != aile[g1$Z], ]
R1 <- grid_calistir(g1, "olcek_Z")
print(R1$S[order(-R1$S$delta_l), c("X","Y","Z","K","delta_l","p_boot","p_binom","poz_dl",
                                   "CV0","CV1","OS","OCR","mi_ihlal","karar","maxZ","C2","taban_ok")],
      row.names = FALSE, right = FALSE)

cizgi("BOLUM 6: GRID 2 - YANIT STILI GOSTERGELERI Z OLARAK")
cat("Stil gostergelerinin odak degiskenlerle korelasyonu (C2):\n")
print(round(cor(d[, c(stil,"ERISIM","MEMNUNIYET","GUVEN","SOY_index","Yas")],
                method="spearman")[stil, c("ERISIM","MEMNUNIYET","GUVEN","SOY_index","Yas")], 3))
cat("\nStil gostergelerinin moderatorlerle eta karesi (C1):\n")
print(do.call(rbind, lapply(stil, function(z)
  data.frame(gosterge = z, t(round(vapply(tum_mod, function(g) {
    a <- summary(aov(d[[z]] ~ factor(d[[g]])))[[1]]; a[1,2]/sum(a[,2]) }, numeric(1)), 3))))),
  row.names = FALSE, right = FALSE)

g2 <- expand.grid(X = names(gost), Y = sonuc_ad, Z = stil, stringsAsFactors = FALSE)
g2 <- g2[aile[g2$X] != aile[g2$Y], ]
R2 <- grid_calistir(g2, "stil_Z")
print(R2$S[order(-R2$S$delta_l), c("X","Y","Z","K","delta_l","p_boot","p_binom","poz_dl",
                                   "CV0","CV1","OS","OCR","mi_ihlal","karar","maxZ","C2","taban_ok")],
      row.names = FALSE, right = FALSE)

cizgi("BOLUM 7: ALT ORNEKLEM CMB TANISI")
alt <- list(
  "Tum orneklem"     = rep(TRUE, nrow(d_tam)),
  "Saglik personeli" = d_tam$MeslekGrup %in% c("Hekim-hemsire-ebe","Saglik teknisyeni","Tibbi sekreter"),
  "Saglik klinik"    = d_tam$MeslekGrup %in% c("Hekim-hemsire-ebe","Saglik teknisyeni"),
  "Saglik disi"      = !(d_tam$MeslekGrup %in% c("Hekim-hemsire-ebe","Saglik teknisyeni","Tibbi sekreter")),
  "Dikkatli"         = d_tam$Longstring < 30 & d_tam$DuzBlok <= 1 & d_tam$FarkliSecenek >= 3,
  "Dikkatsiz"        = !(d_tam$Longstring < 30 & d_tam$DuzBlok <= 1 & d_tam$FarkliSecenek >= 3))
print(do.call(rbind, lapply(names(alt), function(a)
  cbind(altkume = a, cmb_gost(alt[[a]])))), row.names = FALSE, right = FALSE)

cat("\nODAK YOL CV'LERI (SVT hammaddesi)\n")
cv_alt <- function(i, x, y, min_n = 20)
  vapply(tum_mod, function(g) {
    b <- vapply(unique(d_tam[[g]][i]), function(l) {
      j <- i & d_tam[[g]] == l
      if (sum(j) < min_n) return(NA_real_)
      coef(lm(scale(d_tam[[y]][j]) ~ scale(d_tam[[x]][j])))[2] }, numeric(1))
    b <- b[!is.na(b)]
    if (length(b) < 2) return(NA_real_)
    100*sd(b)/max(abs(mean(b)), 1e-4) }, numeric(1))
ciftler <- list(c("ERISIM","MEMNUNIYET"), c("ERISIM","MTT"), c("ERISIM","FCV"),
                c("ERISIM","SOY_index"), c("MEMNUNIYET","SOY_index"))
print(do.call(rbind, lapply(names(alt), function(a)
  do.call(rbind, lapply(ciftler, function(p) {
    cv <- cv_alt(alt[[a]], p[1], p[2])
    data.frame(altkume = a, n = sum(alt[[a]]), yol = paste(p[1], "->", p[2]),
               kullanilabilir = sum(!is.na(cv)), ort_CV = round(mean(cv, na.rm=TRUE), 1),
               maks_CV = round(max(cv, na.rm=TRUE), 1)) })))),
  row.names = FALSE, right = FALSE)

cizgi("BOLUM 8: DIKKATLI ALT ORNEKLEMINDE SVT")
d <- d_tam[alt[["Dikkatli"]], ]
cat("N =", nrow(d), "\n")
for (g in tum_mod) {
  tb <- table(d[[g]])
  cat(sprintf("  %-12s duzey=%d min=%2d  %s\n", g, length(tb), min(tb),
              paste(paste0(names(tb), ":", as.integer(tb)), collapse="  ")))
}
cat("\nMI yeniden hesaplaniyor:\n")
mi_cache <- mi_kur()
g3 <- expand.grid(X = names(gost), Y = sonuc_ad, Z = stab_ad, stringsAsFactors = FALSE)
g3 <- g3[aile[g3$X] != aile[g3$Y] & aile[g3$X] != aile[g3$Z] & aile[g3$Y] != aile[g3$Z], ]
R3 <- grid_calistir(g3, "dikkatli")
print(R3$S[order(-R3$S$delta_l), c("X","Y","Z","K","delta_l","p_boot","p_binom","poz_dl",
                                   "CV0","CV1","OS","OCR","mi_ihlal","karar","maxZ","C2","taban_ok")],
      row.names = FALSE, right = FALSE)

d <- d_tam
mi_cache <- mi_kur()
cat("\nd tam ornekleme geri yuklendi, N =", nrow(d), "\n")

cizgi("BOLUM 9: ONERME 3 AMPIRIK TESTI")

TUM <- rbind(R1$S, R2$S, R3$S)
MOD <- rbind(R1$mod, R2$mod, R3$mod)

cat("KRITER KARSILASTIRMASI\n")
kt <- do.call(rbind, lapply(unique(TUM$etiket), function(e) {
  s <- TUM[TUM$etiket == e, ]
  data.frame(grid = e, uclu = nrow(s),
             sadece_boot = sum(s$boot_gecti & !s$binom_gecti),
             sadece_binom = sum(!s$boot_gecti & s$binom_gecti),
             ikisi_de = sum(s$boot_gecti & s$binom_gecti),
             hicbiri = sum(!s$boot_gecti & !s$binom_gecti),
             boot_orani = round(mean(s$boot_gecti), 3),
             binom_orani = round(mean(s$binom_gecti), 3),
             cift_orani = round(mean(s$boot_gecti & s$binom_gecti), 3)) }))
print(kt, row.names = FALSE, right = FALSE)
cat("\nTOPLAM: bootstrap gecen", sum(TUM$boot_gecti), "/", nrow(TUM),
    " binom gecen", sum(TUM$binom_gecti), "/", nrow(TUM),
    " cift kriter", sum(TUM$karar == "STABILIZATOR"), "/", nrow(TUM), "\n")

cat("\nISARET TESTI: sign(delta_l) ~ Bernoulli(0.5) mi\n")
cat("moderator duzeyi gozlem sayisi :", nrow(MOD), "\n")
cat("pozitif delta_l                :", sum(MOD$dl > 0),
    " (%", round(100*mean(MOD$dl > 0), 1), ")\n", sep = "")
bt <- binom.test(sum(MOD$dl > 0), nrow(MOD), 0.5)
cat("binom testi p                  :", round(bt$p.value, 4),
    "  %95 GA [", round(bt$conf.int[1], 3), ",", round(bt$conf.int[2], 3), "]\n")

cat("\nGRID BAZINDA ISARET ORANI\n")
print(do.call(rbind, lapply(unique(MOD$etiket), function(e) {
  m <- MOD[MOD$etiket == e, ]
  data.frame(grid = e, n = nrow(m), poz = sum(m$dl > 0),
             oran = round(mean(m$dl > 0), 3),
             p = round(binom.test(sum(m$dl > 0), nrow(m), 0.5)$p.value, 4)) })),
  row.names = FALSE, right = FALSE)

cat("\nMODERATOR BAZINDA ISARET ORANI\n")
print(do.call(rbind, lapply(sort(unique(MOD$moderator)), function(g) {
  m <- MOD[MOD$moderator == g, ]
  data.frame(moderator = g, n = nrow(m), poz = sum(m$dl > 0),
             oran = round(mean(m$dl > 0), 3),
             ort_dl = round(mean(m$dl), 3), medyan_dl = round(median(m$dl), 3),
             ort_cv0 = round(mean(m$cv0), 1)) })), row.names = FALSE, right = FALSE)

cat("\nPOZ_DL DAGILIMI vs BINOM(K, 0.5)\n")
for (kk in sort(unique(TUM$K))) {
  s <- TUM[TUM$K == kk, ]
  if (nrow(s) < 10) next
  gzl <- table(factor(s$poz_dl, levels = 0:kk))
  bek <- dbinom(0:kk, kk, 0.5) * nrow(s)
  cat("\nK =", kk, " (", nrow(s), " uclu)\n", sep = "")
  print(rbind(gozlenen = as.integer(gzl), beklenen = round(bek, 1)))
  ki <- suppressWarnings(chisq.test(as.integer(gzl), p = dbinom(0:kk, kk, 0.5)))
  cat("ki-kare =", round(ki$statistic, 2), " sd =", ki$parameter,
      " p =", round(ki$p.value, 4), "\n")
}

cat("\nDELTA_L DAGILIMININ KUYRUK YAPISI (Onerme 3, buyukluk ayagi)\n")
print(round(quantile(MOD$dl, c(.01,.05,.25,.5,.75,.95,.99)), 3))
cat("ortalama :", round(mean(MOD$dl), 3), "  medyan :", round(median(MOD$dl), 3),
    "  carpiklik :", round(psych::skew(MOD$dl), 2),
    "  basiklik :", round(psych::kurtosi(MOD$dl), 2), "\n")
cat("|dl| > 1 olan gozlem :", sum(abs(MOD$dl) > 1), " (%",
    round(100*mean(abs(MOD$dl) > 1), 1), ")\n", sep = "")
cat("bunlarin ortalama CV0'i :", round(mean(MOD$cv0[abs(MOD$dl) > 1]), 1),
    "  digerlerinin :", round(mean(MOD$cv0[abs(MOD$dl) <= 1]), 1), "\n")

cat("\nTABAN ETKISI: CV0 ile |delta_l| iliskisi\n")
print(do.call(rbind, lapply(list(c(0,1),c(1,5),c(5,10),c(10,25),c(25,100),c(100,1e6)), function(b) {
  i <- MOD$cv0 >= b[1] & MOD$cv0 < b[2]
  if (sum(i) < 5) return(NULL)
  data.frame(cv0_araligi = paste0("[", b[1], ",", b[2], ")"), n = sum(i),
             ort_abs_dl = round(mean(abs(MOD$dl[i])), 3),
             maks_abs_dl = round(max(abs(MOD$dl[i])), 3),
             poz_oran = round(mean(MOD$dl[i] > 0), 3)) })),
  row.names = FALSE, right = FALSE)

cizgi("BOLUM 10: KAYIT")
write.csv(TUM, paste0(kok, "svt_uclu_sonuclari.csv"), row.names = FALSE)
write.csv(MOD, paste0(kok, "svt_moderator_duzeyi.csv"), row.names = FALSE)
saveRDS(list(uclu = TUM, moderator = MOD, alt = alt, veri = d_tam),
        paste0(kok, "svt_tum_sonuclar.rds"))
cat("svt_uclu_sonuclari.csv  :", nrow(TUM), "satir\n")
cat("svt_moderator_duzeyi.csv:", nrow(MOD), "satir\n")
cat("svt_tum_sonuclar.rds yazildi\n")

cizgi("OZET")
cat("N (tam)          :", nrow(d_tam), "\n")
cat("N (dikkatli)     :", sum(alt[["Dikkatli"]]), "\n")
cat("toplam uclu      :", nrow(TUM), "\n")
cat("STABILIZATOR     :", sum(TUM$karar == "STABILIZATOR"), "\n")
cat("C2 + taban + karar:", sum(TUM$karar == "STABILIZATOR" & TUM$C2 == "gecti" & TUM$taban_ok), "\n")
print(table(TUM$etiket, TUM$karar))
print(sessionInfo()$otherPkgs$lavaan$Version)



kok     <- "data/derived/"
ham_kok <- Sys.getenv("SVT_RAW_DIR", unset = "raw-data/")
TUM <- read.csv(paste0(kok, "svt_uclu_sonuclari.csv"), stringsAsFactors = FALSE)
MOD <- read.csv(paste0(kok, "svt_moderator_duzeyi.csv"), stringsAsFactors = FALSE)
anah <- function(df) paste(df$etiket, df$X, df$Y, df$Z, sep = "|")
MOD$maxZ <- TUM$maxZ[match(anah(MOD), anah(TUM))]
MOD$C2 <- TUM$C2[match(anah(MOD), anah(TUM))]
MOD$poz <- as.integer(MOD$dl > 0)

cat("A) ISARET ORANI GRADYANI\n")
MOD$mz_bin <- cut(MOD$maxZ, c(-Inf, .1, .2, .3, .5, Inf),
                  labels = c("<.10", ".10-.20", ".20-.30", ".30-.50", ">=.50"))
print(do.call(rbind, lapply(split(MOD, MOD$mz_bin), function(m)
  data.frame(bin = as.character(m$mz_bin[1]), n = nrow(m),
             poz_oran = round(mean(m$poz), 3),
             medyan_dl = round(median(m$dl), 3),
             ort_dl = round(mean(m$dl), 3)))), row.names = FALSE, right = FALSE)
print(do.call(rbind, lapply(split(MOD, MOD$C2), function(m)
  data.frame(C2 = m$C2[1], n = nrow(m), poz_oran = round(mean(m$poz), 3),
             medyan_dl = round(median(m$dl), 3)))), row.names = FALSE, right = FALSE)
sp <- cor.test(MOD$maxZ, MOD$dl, method = "spearman", exact = FALSE)
cat("spearman(maxZ, dl) =", round(sp$estimate, 3), " p =", signif(sp$p.value, 3), "\n\n")

cat("B) MODERATORLER ARASI ISARET UYUMU\n")
W <- reshape(MOD[, c("etiket", "X", "Y", "Z", "moderator", "poz")],
             idvar = c("etiket", "X", "Y", "Z"), timevar = "moderator", direction = "wide")
S6 <- W[complete.cases(W[, grep("^poz\\.", names(W))]), grep("^poz\\.", names(W))]
names(S6) <- sub("poz\\.", "", names(S6))
cat("tam K=6 uclu :", nrow(S6), "\n")
p_m <- colMeans(S6)
print(round(p_m, 3))
mods <- names(S6)
uy_g <- uy_b <- matrix(NA, length(mods), length(mods), dimnames = list(mods, mods))
for (i in seq_along(mods)) for (j in seq_along(mods)) {
  uy_g[i, j] <- mean(S6[[i]] == S6[[j]])
  uy_b[i, j] <- p_m[i]*p_m[j] + (1-p_m[i])*(1-p_m[j])
}
cat("gozlenen uyum - bagimsizlik beklentisi:\n")
print(round(uy_g - uy_b, 3))
cat("ort fark (ust ucgen):", round(mean((uy_g - uy_b)[upper.tri(uy_g)]), 3), "\n")
S <- rowSums(S6)
for (e in unique(W$etiket)) {
  s <- S[W$etiket[complete.cases(W[, grep("^poz\\.", names(W))])] == e]
  if (length(s) < 20) next
  ph <- mean(s)/6
  cat(sprintf("%-10s n=%3d  p=%.3f  Var(S)=%.2f  binomVar=%.2f  DI=%.2f\n",
              e, length(s), ph, var(s), 6*ph*(1-ph), var(s)/(6*ph*(1-ph))))
}

cat("\nC) BOOTSTRAP TETIKLENEN 5 UCLU\n")
bg <- TUM[TUM$boot_gecti, c("etiket", "X", "Y", "Z")]
for (i in seq_len(nrow(bg))) {
  m <- MOD[MOD$etiket == bg$etiket[i] & MOD$X == bg$X[i] &
             MOD$Y == bg$Y[i] & MOD$Z == bg$Z[i], ]
  cat("\n", bg$etiket[i], bg$X[i], "->", bg$Y[i], "| Z =", bg$Z[i], "\n")
  print(m[, c("moderator", "dl", "ocr", "cv0", "cv1")], row.names = FALSE, digits = 3)
}

cat("\nD) DL DAGILIMI maxZ BINLERINDE\n")
print(do.call(rbind, lapply(split(MOD, MOD$mz_bin), function(m)
  data.frame(bin = as.character(m$mz_bin[1]),
             q05 = round(quantile(m$dl, .05), 2), q25 = round(quantile(m$dl, .25), 2),
             q75 = round(quantile(m$dl, .75), 2), q95 = round(quantile(m$dl, .95), 2),
             ort_abs = round(mean(abs(m$dl)), 2)))), row.names = FALSE, right = FALSE)

set.seed(9186)

mc_gruplar <- function(n, bag = 1) {
  u <- rnorm(n)
  kar <- function(w) bag*w*u + rnorm(n, 0, sqrt(max(0.02, 1 - (bag*w)^2)))
  data.frame(
    g1 = 1 + rbinom(n, 1, 0.37),
    g2 = 1 + as.integer(kar(0.7) > 0.15),
    g3 = as.integer(cut(kar(0.6), c(-Inf, -0.85, 0.75, Inf))),
    g4 = 1 + rbinom(n, 1, 0.35),
    g5 = as.integer(cut(rnorm(n), c(-Inf, -1.15, 0.95, Inf))),
    g6 = as.integer(cut(kar(0.8), c(-Inf, -1.25, -0.55, 0.25, 1.05, Inf))))
}

mc_veri <- function(n = 500, rho_m = 0, rho_t = 0.55, r_b = 0.65, sd_het = 0.15,
                    lam_zx = 0, stab = 0, s_a = 0.12, bag = 1) {
  G <- mc_gruplar(n, bag)
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het/sqrt(6)))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  lam <- min(0.95, 0.15*(stab > 0) + lam_zx)
  z <- lam*tx + sqrt(1 - lam^2)*rnorm(n)
  aj <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, s_a))
  ai <- rowSums(vapply(1:6, function(j) aj[[j]][G[[j]]], numeric(n)))
  ytrue <- bi*tx + stab*ai*z + sqrt(max(0.05, 1 - rho_t^2))*rnorm(n)
  ytrue <- as.numeric(scale(ytrue))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), G)
}

mc_stab <- function(dat, g, min_n = 25) {
  b0 <- b1 <- nn <- numeric(0)
  for (l in sort(unique(dat[[g]]))) {
    sub <- dat[dat[[g]] == l, ]
    if (nrow(sub) < min_n) next
    sx <- as.numeric(scale(sub$X)); sy <- as.numeric(scale(sub$Y)); sz <- as.numeric(scale(sub$Z))
    f0 <- coef(lm(sy ~ sx))[2]
    f1 <- coef(lm(sy ~ sx + sz))[2]
    if (any(is.na(c(f0, f1)))) next
    b0 <- c(b0, f0); b1 <- c(b1, f1); nn <- c(nn, nrow(sub))
  }
  if (length(b0) < 2) return(NULL)
  wm <- function(v) weighted.mean(v, nn)
  ws <- function(v) sqrt(sum(nn*(v - wm(v))^2)/sum(nn))
  mu0 <- wm(b0); mu1 <- wm(b1); sd0 <- ws(b0); sd1 <- ws(b1)
  cv0 <- 100*sd0/max(abs(mu0), 1e-4)
  cv1 <- 100*sd1/max(abs(mu1), 1e-4)
  list(dl = log(cv0 + 1e-8) - log(cv1 + 1e-8), cv0 = cv0, cv1 = cv1,
       ocr = mean(as.integer(abs(b1 - mu1) < abs(b0 - mu0))), n = sum(nn))
}

mc_svt <- function(dat, R = 500) {
  res <- lapply(paste0("g", 1:6), function(g) mc_stab(dat, g))
  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) < 3) return(NULL)
  dv <- vapply(res, function(r) r$dl, numeric(1))
  wv <- vapply(res, function(r) r$n, numeric(1))
  ov <- vapply(res, function(r) r$ocr, numeric(1))
  cv0v <- vapply(res, function(r) r$cv0, numeric(1))
  bo <- boot(dv, function(dd, idx) sum(dd[idx]*wv[idx]/sum(wv[idx])), R = R)
  md <- bo$t0; se <- sd(bo$t)
  pb <- if (is.finite(se) && se > 0) pnorm(md/se, lower.tail = FALSE) else 1
  kr2 <- (dv > 0) | (ov >= 0.5)
  data.frame(p_boot = pb, md = md,
             p_sign = binom.test(sum(dv > 0), length(dv), 0.5, "greater")$p.value,
             p_or = binom.test(sum(kr2), length(kr2), 0.5, "greater")$p.value,
             poz = sum(dv > 0), M = length(dv), med_dl = median(dv),
             cv0_ort = mean(cv0v), cv0_min = min(cv0v),
             isaret = paste(as.integer(dv > 0), collapse = ""),
             stringsAsFactors = FALSE)
}

mc_hucre <- function(rep = 30, R = 500, ...) {
  t0 <- Sys.time()
  out <- vector("list", rep)
  for (i in seq_len(rep)) {
    r <- NULL
    while (is.null(r)) r <- mc_svt(mc_veri(...), R = R)
    out[[i]] <- r
  }
  H <- do.call(rbind, out)
  H$boot_g <- H$p_boot < 0.05 & H$md > 0
  H$sign_g <- H$p_sign < 0.05
  H$or_g <- H$p_or < 0.05
  attr(H, "dk") <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
  H
}

ozet <- function(H, ad) cat(sprintf(
  "%-18s rep=%3d boot=%.3f sign=%.3f or=%.3f ciftS=%.3f ciftOR=%.3f p+=%.3f meddl=%+.3f cv0=%.1f %.2f dk\n",
  ad, nrow(H), mean(H$boot_g), mean(H$sign_g), mean(H$or_g),
  mean(H$boot_g & H$sign_g), mean(H$boot_g & H$or_g),
  mean(H$poz/H$M), median(H$med_dl), median(H$cv0_ort), attr(H, "dk")))

h1 <- mc_hucre(rho_m = 0.0, stab = 0)
ozet(h1, "rho0 null")
h2 <- mc_hucre(rho_m = 0.8, stab = 0)
ozet(h2, "rho.8 null")
h3 <- mc_hucre(rho_m = 0.0, stab = 1)
ozet(h3, "rho0 stab")
h4 <- mc_hucre(rho_m = 0.8, stab = 1)
ozet(h4, "rho.8 stab")
v0 <- mc_veri(rho_m = 0); v8 <- mc_veri(rho_m = 0.8)
cat("r(X,Y) rho0:", round(cor(v0$X, v0$Y), 2), " rho.8:", round(cor(v8$X, v8$Y), 2), "\n")


set.seed(9186)
kok     <- "data/derived/"
ham_kok <- Sys.getenv("SVT_RAW_DIR", unset = "raw-data/")

mc_gruplar <- function(n, bag = 1) {
  u <- rnorm(n)
  kar <- function(w) bag*w*u + rnorm(n, 0, sqrt(max(0.02, 1 - (bag*w)^2)))
  list(G = data.frame(
    g1 = 1 + rbinom(n, 1, 0.37),
    g2 = 1 + as.integer(kar(0.7) > 0.15),
    g3 = as.integer(cut(kar(0.6), c(-Inf, -0.85, 0.75, Inf))),
    g4 = 1 + rbinom(n, 1, 0.35),
    g5 = as.integer(cut(rnorm(n), c(-Inf, -1.15, 0.95, Inf))),
    g6 = as.integer(cut(kar(0.8), c(-Inf, -1.25, -0.55, 0.25, 1.05, Inf)))),
    u = u)
}

mc_veri <- function(n = 500, rho_m = 0, rho_t = 0.60, r_b = 0.55, sd_het = 0.18,
                    lam_zx = 0, stab = 0, s_a = 0.5, bag = 1) {
  gu <- mc_gruplar(n, bag)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  z <- lam_zx*tx + sqrt(1 - lam_zx^2)*rnorm(n)
  ytrue <- bi*tx + stab*s_a*gu$u*z + sqrt(1 - rho_t^2)*rnorm(n)
  ytrue <- as.numeric(scale(ytrue))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), G)
}

mc_stab <- function(dat, g, min_n = 25) {
  b0 <- b1 <- nn <- numeric(0)
  for (l in sort(unique(dat[[g]]))) {
    sub <- dat[dat[[g]] == l, ]
    if (nrow(sub) < min_n) next
    sx <- as.numeric(scale(sub$X)); sy <- as.numeric(scale(sub$Y)); sz <- as.numeric(scale(sub$Z))
    f0 <- coef(lm(sy ~ sx))[2]
    f1 <- coef(lm(sy ~ sx + sz))[2]
    if (any(is.na(c(f0, f1)))) next
    b0 <- c(b0, f0); b1 <- c(b1, f1); nn <- c(nn, nrow(sub))
  }
  if (length(b0) < 2) return(NULL)
  wm <- function(v) weighted.mean(v, nn)
  ws <- function(v) sqrt(sum(nn*(v - wm(v))^2)/sum(nn))
  mu0 <- wm(b0); mu1 <- wm(b1); sd0 <- ws(b0); sd1 <- ws(b1)
  cv0 <- 100*sd0/max(abs(mu0), 1e-4)
  cv1 <- 100*sd1/max(abs(mu1), 1e-4)
  list(dl = log(cv0 + 1e-8) - log(cv1 + 1e-8), cv0 = cv0, cv1 = cv1,
       ocr = mean(as.integer(abs(b1 - mu1) < abs(b0 - mu0))), n = sum(nn))
}

mc_svt <- function(dat, R = 500) {
  res <- lapply(paste0("g", 1:6), function(g) mc_stab(dat, g))
  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) < 3) return(NULL)
  dv <- vapply(res, function(r) r$dl, numeric(1))
  wv <- vapply(res, function(r) r$n, numeric(1))
  ov <- vapply(res, function(r) r$ocr, numeric(1))
  cv0v <- vapply(res, function(r) r$cv0, numeric(1))
  bo <- boot(dv, function(dd, idx) sum(dd[idx]*wv[idx]/sum(wv[idx])), R = R)
  md <- bo$t0; se <- sd(bo$t)
  pb <- if (is.finite(se) && se > 0) pnorm(md/se, lower.tail = FALSE) else 1
  kr2 <- (dv > 0) | (ov >= 0.5)
  data.frame(p_boot = pb, md = md,
             p_sign = binom.test(sum(dv > 0), length(dv), 0.5, "greater")$p.value,
             p_or = binom.test(sum(kr2), length(kr2), 0.5, "greater")$p.value,
             poz = sum(dv > 0), M = length(dv), med_dl = median(dv),
             cv0_ort = mean(cv0v), cv0_min = min(cv0v),
             isaret = paste(as.integer(dv > 0), collapse = ""),
             stringsAsFactors = FALSE)
}

mc_hucre <- function(rep = 30, R = 500, ...) {
  t0 <- Sys.time()
  out <- vector("list", rep)
  for (i in seq_len(rep)) {
    r <- NULL
    while (is.null(r)) r <- mc_svt(mc_veri(...), R = R)
    out[[i]] <- r
  }
  H <- do.call(rbind, out)
  H$boot_g <- H$p_boot < 0.05 & H$md > 0
  H$sign_g <- H$p_sign < 0.05
  H$or_g <- H$p_or < 0.05
  attr(H, "dk") <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
  H
}

mc_kalibre <- function(rho_vec = c(0, 0.3, 0.8), rep = 40) {
  do.call(rbind, lapply(rho_vec, function(r) {
    rk <- cvm <- numeric(rep)
    for (i in seq_len(rep)) {
      v <- mc_veri(rho_m = r)
      rk[i] <- cor(v$X, v$Y)
      cvv <- vapply(paste0("g", 1:6), function(g) {
        st <- mc_stab(v, g); if (is.null(st)) NA_real_ else st$cv0 }, numeric(1))
      cvm[i] <- mean(cvv, na.rm = TRUE)
    }
    data.frame(rho_m = r, r_xy = round(mean(rk), 3),
               cv0_med = round(median(cvm), 1),
               cv0_q25 = round(quantile(cvm, 0.25), 1),
               cv0_q75 = round(quantile(cvm, 0.75), 1))
  }))
}

cat("KALIBRASYON\n")
print(mc_kalibre(), row.names = FALSE, right = FALSE)

faz_a <- function(rho_vec = c(0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.85, 0.95),
                  rep = 1000, R = 500) {
  ozl <- list(); ham <- list()
  for (het in c(0.18, 0)) for (r in rho_vec) {
    H <- mc_hucre(rep = rep, R = R, rho_m = r, sd_het = het, stab = 0, lam_zx = 0)
    H$rho_m <- r; H$sd_het <- het
    et <- sprintf("het%.2f_rho%.2f", het, r)
    ham[[et]] <- H
    h6 <- H[H$M == 6, ]
    ph <- mean(h6$poz)/6
    di <- if (nrow(h6) > 30 && ph > 0 && ph < 1) var(h6$poz)/(6*ph*(1 - ph)) else NA
    ozl[[et]] <- data.frame(sd_het = het, rho_m = r,
                            boot = mean(H$boot_g), sign = mean(H$sign_g), or_ = mean(H$or_g),
                            cift_sign = mean(H$boot_g & H$sign_g), cift_or = mean(H$boot_g & H$or_g),
                            p_poz = round(ph, 3), med_dl = round(median(H$med_dl), 3),
                            cv0 = round(median(H$cv0_ort), 1), DI = round(di, 2), dk = attr(H, "dk"))
    cat(sprintf("%s  boot=%.3f sign=%.3f or=%.3f  p+=%.3f  cv0=%.1f  DI=%.2f  %.1f dk\n",
                et, mean(H$boot_g), mean(H$sign_g), mean(H$or_g), ph,
                median(H$cv0_ort), di, attr(H, "dk")))
  }
  A <- do.call(rbind, ozl)
  write.csv(A, paste0(kok, "mc_fazA_ozet.csv"), row.names = FALSE)
  write.csv(do.call(rbind, ham), paste0(kok, "mc_fazA_ham.csv"), row.names = FALSE)
  A
}

fazA <- faz_a()
cat("\nFAZ A OZET\n")
print(fazA, row.names = FALSE, right = FALSE)

faz_b <- function(lam_vec = c(0, 0.15, 0.30, 0.45, 0.60),
                  rho_vec = c(0.3, 0.8), rep = 1000, R = 500) {
  ozl <- list()
  for (r in rho_vec) for (lz in lam_vec) {
    H <- mc_hucre(rep = rep, R = R, rho_m = r, sd_het = 0.18, stab = 0, lam_zx = lz)
    et <- sprintf("rho%.1f_lam%.2f", r, lz)
    ozl[[et]] <- data.frame(rho_m = r, lam_zx = lz,
                            boot = mean(H$boot_g), sign = mean(H$sign_g), or_ = mean(H$or_g),
                            cift_or = mean(H$boot_g & H$or_g),
                            p_poz = round(mean(H$poz/H$M), 3), med_dl = round(median(H$med_dl), 3),
                            q25_dl = round(quantile(H$med_dl, 0.25), 3), dk = attr(H, "dk"))
    cat(sprintf("%s  boot=%.3f or=%.3f  p+=%.3f  meddl=%+.3f  %.1f dk\n",
                et, mean(H$boot_g), mean(H$or_g),
                mean(H$poz/H$M), median(H$med_dl), attr(H, "dk")))
  }
  B <- do.call(rbind, ozl)
  write.csv(B, paste0(kok, "mc_fazB.csv"), row.names = FALSE)
  B
}

c1_veri <- function(K = 8, ng = 60, rho_m = 0, s_a = 0.5, rho_t = 0.6,
                    r_b = 0.55, sd_het = 0.05) {
  n <- K*ng
  g <- rep(seq_len(K), each = ng)
  bk <- rho_t + rnorm(K, 0, sd_het)
  ak <- rnorm(K, 0, s_a)
  x <- rnorm(n)
  z <- 0.15*x + sqrt(1 - 0.15^2)*rnorm(n)
  y <- bk[g]*x + ak[g]*z + sqrt(1 - rho_t^2)*rnorm(n)
  y <- as.numeric(scale(y))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*x + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*y + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), g = g)
}

c1_svt <- function(dat, R = 500) {
  b0 <- b1 <- nn <- numeric(0)
  for (l in sort(unique(dat$g))) {
    sub <- dat[dat$g == l, ]
    sx <- as.numeric(scale(sub$X)); sy <- as.numeric(scale(sub$Y)); sz <- as.numeric(scale(sub$Z))
    b0 <- c(b0, coef(lm(sy ~ sx))[2])
    b1 <- c(b1, coef(lm(sy ~ sx + sz))[2])
    nn <- c(nn, nrow(sub))
  }
  wm <- function(v) weighted.mean(v, nn)
  mu0 <- wm(b0); mu1 <- wm(b1)
  dk <- abs(b0 - mu0) - abs(b1 - mu1)
  Ik <- as.integer(abs(b1 - mu1) < abs(b0 - mu0))
  bo <- boot(dk, function(dd, idx) sum(dd[idx]*nn[idx]/sum(nn[idx])), R = R)
  dw <- bo$t0; se <- sd(bo$t)
  pb <- if (is.finite(se) && se > 0) pnorm(dw/se, lower.tail = FALSE) else 1
  pbin <- binom.test(sum(Ik), length(Ik), 0.5, "greater")$p.value
  ws <- function(v) sqrt(sum(nn*(v - wm(v))^2)/sum(nn))
  cv0 <- 100*ws(b0)/max(abs(mu0), 1e-4)
  cv1 <- 100*ws(b1)/max(abs(mu1), 1e-4)
  data.frame(p_boot = pb, dw = dw, p_binom = pbin, S = sum(Ik),
             dl = log(cv0 + 1e-8) - log(cv1 + 1e-8), cv0 = cv0)
}

faz_c1 <- function(rho_vec = c(0, 0.3, 0.45, 0.6, 0.75, 0.9),
                   rep = 1000, R = 500, s_a = 0.5) {
  ozl <- list()
  for (r in rho_vec) {
    t0 <- Sys.time()
    H <- do.call(rbind, lapply(seq_len(rep), function(i)
      c1_svt(c1_veri(rho_m = r, s_a = s_a), R = R)))
    H$boot_g <- H$p_boot < 0.05 & H$dw > 0
    H$bin_g <- H$p_binom < 0.05
    dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
    ozl[[as.character(r)]] <- data.frame(rho_m = r,
                                         boot = mean(H$boot_g), binom = mean(H$bin_g), cift = mean(H$boot_g & H$bin_g),
                                         med_dl = round(median(H$dl), 3), med_S = median(H$S),
                                         cv0 = round(median(H$cv0), 1), dk = dk)
    cat(sprintf("C1 rho=%.2f  boot=%.3f binom=%.3f cift=%.3f  meddl=%+.3f  S=%.0f  cv0=%.1f  %.1f dk\n",
                r, mean(H$boot_g), mean(H$bin_g), mean(H$boot_g & H$bin_g),
                median(H$dl), median(H$S), median(H$cv0), dk))
  }
  C1 <- do.call(rbind, ozl)
  write.csv(C1, paste0(kok, "mc_fazC1.csv"), row.names = FALSE)
  C1
}

mc_veri_c2 <- function(n = 500, rho_m = 0, rho_t = 0.60, r_b = 0.55, sd_het = 0.18,
                       s_c = 0.4, art = 1, bag = 1) {
  gu <- mc_gruplar(n, bag)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  cj <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, s_c))
  ci <- rowSums(vapply(1:6, function(j) cj[[j]][G[[j]]], numeric(n)))
  z <- 0.15*tx + sqrt(1 - 0.15^2)*rnorm(n)
  ytrue <- bi*tx + art*ci*z + sqrt(1 - rho_t^2)*rnorm(n)
  ytrue <- as.numeric(scale(ytrue))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), G)
}

faz_c2 <- function(rho_vec = c(0, 0.3, 0.45, 0.6, 0.75, 0.9),
                   rep = 1000, R = 500, s_c = 0.4) {
  ozl <- list()
  for (r in rho_vec) {
    t0 <- Sys.time()
    out <- vector("list", rep)
    for (i in seq_len(rep)) {
      rr <- NULL
      while (is.null(rr)) rr <- mc_svt(mc_veri_c2(rho_m = r, s_c = s_c), R = R)
      out[[i]] <- rr
    }
    H <- do.call(rbind, out)
    H$boot_g <- H$p_boot < 0.05 & H$md > 0
    H$sign_g <- H$p_sign < 0.05
    H$or_g <- H$p_or < 0.05
    dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
    ozl[[as.character(r)]] <- data.frame(rho_m = r,
                                         boot = mean(H$boot_g), sign = mean(H$sign_g), or_ = mean(H$or_g),
                                         cift_sign = mean(H$boot_g & H$sign_g), cift_or = mean(H$boot_g & H$or_g),
                                         p_poz = round(mean(H$poz/H$M), 3), med_dl = round(median(H$med_dl), 3),
                                         cv0 = round(median(H$cv0_ort), 1), dk = dk)
    cat(sprintf("C2 rho=%.2f  boot=%.3f or=%.3f ciftOR=%.3f  p+=%.3f meddl=%+.3f  cv0=%.1f  %.1f dk\n",
                r, mean(H$boot_g), mean(H$or_g), mean(H$boot_g & H$or_g),
                mean(H$poz/H$M), median(H$med_dl), median(H$cv0_ort), dk))
  }
  C2 <- do.call(rbind, ozl)
  write.csv(C2, paste0(kok, "mc_fazC2.csv"), row.names = FALSE)
  C2
}

fazB <- faz_b()
fazC1 <- faz_c1()
fazC2 <- faz_c2()
cat("\nFAZ B\n"); print(fazB, row.names = FALSE, right = FALSE)
cat("\nFAZ C1\n"); print(fazC1, row.names = FALSE, right = FALSE)
cat("\nFAZ C2\n"); print(fazC2, row.names = FALSE, right = FALSE)


mc_veri_b2 <- function(n = 500, rho_m = 0.3, phi = 0.35, rho_t = 0.60,
                       r_b = 0.55, sd_het = 0.18) {
  gu <- mc_gruplar(n)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  ytrue <- as.numeric(scale(bi*tx + sqrt(1 - rho_t^2)*rnorm(n)))
  s <- rnorm(n)
  z <- phi*s + sqrt(1 - phi^2)*rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), G)
}

faz_b2 <- function(phi_vec = c(0.35, 0.60), rho_vec = c(0.3, 0.8),
                   rep = 1000, R = 500) {
  ozl <- list()
  for (r in rho_vec) for (ph in phi_vec) {
    t0 <- Sys.time()
    out <- vector("list", rep); rxz <- numeric(rep)
    for (i in seq_len(rep)) {
      rr <- NULL
      while (is.null(rr)) { v <- mc_veri_b2(rho_m = r, phi = ph); rr <- mc_svt(v, R = R) }
      out[[i]] <- rr; rxz[i] <- cor(v$X, v$Z)
    }
    H <- do.call(rbind, out)
    H$boot_g <- H$p_boot < 0.05 & H$md > 0
    H$or_g <- H$p_or < 0.05
    dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
    ozl[[paste(r, ph)]] <- data.frame(rho_m = r, phi = ph, r_xz = round(mean(rxz), 3),
                                      boot = mean(H$boot_g), or_ = mean(H$or_g), cift_or = mean(H$boot_g & H$or_g),
                                      p_poz = round(mean(H$poz/H$M), 3), med_dl = round(median(H$med_dl), 3), dk = dk)
    cat(sprintf("B2 rho=%.1f phi=%.2f  rxz=%+.2f  boot=%.3f or=%.3f cift=%.3f  p+=%.3f meddl=%+.3f  %.1f dk\n",
                r, ph, mean(rxz), mean(H$boot_g), mean(H$or_g),
                mean(H$boot_g & H$or_g), mean(H$poz/H$M), median(H$med_dl), dk))
  }
  B2 <- do.call(rbind, ozl)
  write.csv(B2, paste0(kok, "mc_fazB2.csv"), row.names = FALSE)
  B2
}

c1_veri <- function(K = 8, ng = 60, rho_m = 0, s_a = 1.0, lam = 0.25,
                    rho_t = 0.6, r_b = 0.55, sd_het = 0.05) {
  n <- K*ng
  g <- rep(seq_len(K), each = ng)
  bk <- rho_t + rnorm(K, 0, sd_het)
  ak <- rnorm(K, 0, s_a)
  x <- rnorm(n)
  z <- lam*x + sqrt(1 - lam^2)*rnorm(n)
  y <- bk[g]*x + ak[g]*z + sqrt(1 - rho_t^2)*rnorm(n)
  y <- as.numeric(scale(y))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*x + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*y + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), g = g)
}

faz_c1r <- function(rho_vec = c(0, 0.3, 0.45, 0.6, 0.75, 0.9),
                    rep = 1000, R = 500) {
  ozl <- list()
  for (r in rho_vec) {
    t0 <- Sys.time()
    out <- vector("list", rep); rxy <- numeric(rep)
    for (i in seq_len(rep)) {
      v <- c1_veri(rho_m = r)
      rxy[i] <- cor(v$X, v$Y)
      out[[i]] <- c1_svt(v, R = R)
    }
    H <- do.call(rbind, out)
    H$boot_g <- H$p_boot < 0.05 & H$dw > 0
    H$bin_g <- H$p_binom < 0.05
    dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
    ozl[[as.character(r)]] <- data.frame(rho_m = r, r_xy = round(mean(rxy), 3),
                                         boot = mean(H$boot_g), binom = mean(H$bin_g), cift = mean(H$boot_g & H$bin_g),
                                         med_dl = round(median(H$dl), 3), med_S = median(H$S),
                                         cv0 = round(median(H$cv0), 1), dk = dk)
    cat(sprintf("C1r rho=%.2f  rxy=%.2f  boot=%.3f binom=%.3f cift=%.3f  meddl=%+.3f  S=%.0f  cv0=%.1f  %.1f dk\n",
                r, mean(rxy), mean(H$boot_g), mean(H$bin_g),
                mean(H$boot_g & H$bin_g), median(H$dl), median(H$S), median(H$cv0), dk))
  }
  C1 <- do.call(rbind, ozl)
  write.csv(C1, paste0(kok, "mc_fazC1r.csv"), row.names = FALSE)
  C1
}

mc_veri_c2 <- function(n = 500, rho_m = 0, rho_t = 0.60, r_b = 0.55, sd_het = 0.18,
                       s_c = 0.4, lam = 0.25, bag = 1) {
  gu <- mc_gruplar(n, bag)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  cj <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, s_c))
  ci <- rowSums(vapply(1:6, function(j) cj[[j]][G[[j]]], numeric(n)))
  z <- lam*tx + sqrt(1 - lam^2)*rnorm(n)
  ytrue <- bi*tx + ci*z + sqrt(1 - rho_t^2)*rnorm(n)
  ytrue <- as.numeric(scale(ytrue))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx,
             Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
             Z = as.numeric(scale(z)), G)
}

faz_c2r <- function(rho_vec = c(0, 0.3, 0.45, 0.6, 0.75, 0.9),
                    rep = 1000, R = 500) {
  ozl <- list()
  for (r in rho_vec) {
    t0 <- Sys.time()
    out <- vector("list", rep); rxy <- numeric(rep)
    for (i in seq_len(rep)) {
      rr <- NULL
      while (is.null(rr)) { v <- mc_veri_c2(rho_m = r); rr <- mc_svt(v, R = R) }
      out[[i]] <- rr; rxy[i] <- cor(v$X, v$Y)
    }
    H <- do.call(rbind, out)
    H$boot_g <- H$p_boot < 0.05 & H$md > 0
    H$sign_g <- H$p_sign < 0.05
    H$or_g <- H$p_or < 0.05
    dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
    ozl[[as.character(r)]] <- data.frame(rho_m = r, r_xy = round(mean(rxy), 3),
                                         boot = mean(H$boot_g), sign = mean(H$sign_g), or_ = mean(H$or_g),
                                         cift_sign = mean(H$boot_g & H$sign_g), cift_or = mean(H$boot_g & H$or_g),
                                         p_poz = round(mean(H$poz/H$M), 3), med_dl = round(median(H$med_dl), 3),
                                         cv0 = round(median(H$cv0_ort), 1), dk = dk)
    cat(sprintf("C2r rho=%.2f  rxy=%.2f  boot=%.3f or=%.3f ciftOR=%.3f  p+=%.3f meddl=%+.3f  cv0=%.1f  %.1f dk\n",
                r, mean(rxy), mean(H$boot_g), mean(H$or_g),
                mean(H$boot_g & H$or_g), mean(H$poz/H$M),
                median(H$med_dl), median(H$cv0_ort), dk))
  }
  C2 <- do.call(rbind, ozl)
  write.csv(C2, paste0(kok, "mc_fazC2r.csv"), row.names = FALSE)
  C2
}

faz_e <- function(q_vec = c(0, 0.2, 0.35, 0.5), rep = 300, rho_m = 0.15) {
  E <- do.call(rbind, lapply(q_vec, function(q) {
    rxy <- cvm <- numeric(rep)
    for (i in seq_len(rep)) {
      v <- mc_veri(rho_m = rho_m)
      k <- which(runif(nrow(v)) < q)
      if (length(k)) {
        a <- rnorm(length(k))
        v$X[k] <- 0.9*a + sqrt(0.19)*rnorm(length(k))
        v$Y[k] <- 0.9*a + sqrt(0.19)*rnorm(length(k))
      }
      rxy[i] <- cor(v$X, v$Y)
      cvv <- vapply(paste0("g", 1:6), function(g) {
        st <- mc_stab(v, g); if (is.null(st)) NA_real_ else st$cv0 }, numeric(1))
      cvm[i] <- mean(cvv, na.rm = TRUE)
    }
    out <- data.frame(q = q, r_xy = round(mean(rxy), 3),
                      cv0_med = round(median(cvm), 1),
                      cv0_q25 = round(quantile(cvm, 0.25), 1),
                      cv0_q75 = round(quantile(cvm, 0.75), 1))
    cat(sprintf("E q=%.2f  rxy=%.2f  cv0=%.1f [%.1f, %.1f]\n",
                q, out$r_xy, out$cv0_med, out$cv0_q25, out$cv0_q75))
    out
  }))
  write.csv(E, paste0(kok, "mc_fazE.csv"), row.names = FALSE)
  E
}

fazB2 <- faz_b2()
fazC1r <- faz_c1r()
fazC2r <- faz_c2r()
fazE <- faz_e()
cat("\nFAZ B2\n"); print(fazB2, row.names = FALSE, right = FALSE)
cat("\nFAZ C1r\n"); print(fazC1r, row.names = FALSE, right = FALSE)
cat("\nFAZ C2r\n"); print(fazC2r, row.names = FALSE, right = FALSE)
cat("\nFAZ E\n"); print(fazE, row.names = FALSE, right = FALSE)



t2_veri <- function(n = 500, rho_m = 0.15, rho_t = 0.60, r_b = 0.55, sd_het = 0.18,
                    s_c = 0, lam = 0, yuk = 0.85) {
  gu <- mc_gruplar(n)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  if (s_c > 0) {
    cj <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, s_c))
    ci <- rowSums(vapply(1:6, function(j) cj[[j]][G[[j]]], numeric(n)))
  } else ci <- 0
  z <- lam*tx + sqrt(1 - lam^2)*rnorm(n)
  ytrue <- as.numeric(scale(bi*tx + ci*z + sqrt(1 - rho_t^2)*rnorm(n)))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  Tg <- sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx
  dat <- data.frame(Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
                    Z = as.numeric(scale(z)), G)
  for (j in 1:4) dat[[paste0("x", j)]] <- yuk*Tg + sqrt(1 - yuk^2)*rnorm(n)
  dat
}

t2_stab <- function(dat, g, min_n = 25) {
  m0 <- "F =~ x1 + x2 + x3 + x4\nY ~ c*F"
  m1 <- "F =~ x1 + x2 + x3 + x4\nZ ~ a*F\nY ~ cprime*F + b*Z"
  b0 <- b1 <- nn <- numeric(0)
  for (l in sort(unique(dat[[g]]))) {
    sub <- dat[dat[[g]] == l, ]
    if (nrow(sub) < min_n) next
    f0 <- tryCatch(suppressWarnings(sem(m0, sub, std.lv = TRUE)), error = function(e) NULL)
    f1 <- tryCatch(suppressWarnings(sem(m1, sub, std.lv = TRUE)), error = function(e) NULL)
    if (is.null(f0) || is.null(f1)) next
    q0 <- standardizedSolution(f0); q1 <- standardizedSolution(f1)
    v0 <- q0$est.std[q0$label == "c"][1]
    v1 <- q1$est.std[q1$label == "cprime"][1]
    if (any(is.na(c(v0, v1)))) next
    b0 <- c(b0, v0); b1 <- c(b1, v1); nn <- c(nn, nrow(sub))
  }
  if (length(b0) < 2) return(NULL)
  wm <- function(v) weighted.mean(v, nn)
  ws <- function(v) sqrt(sum(nn*(v - wm(v))^2)/sum(nn))
  mu0 <- wm(b0); mu1 <- wm(b1); sd0 <- ws(b0); sd1 <- ws(b1)
  cv0 <- 100*sd0/max(abs(mu0), 1e-4)
  cv1 <- 100*sd1/max(abs(mu1), 1e-4)
  list(dl = log(cv0 + 1e-8) - log(cv1 + 1e-8), cv0 = cv0, cv1 = cv1,
       ocr = mean(as.integer(abs(b1 - mu1) < abs(b0 - mu0))), n = sum(nn))
}

t2_svt <- function(dat, R = 500) {
  res <- lapply(paste0("g", 1:6), function(g) t2_stab(dat, g))
  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) < 3) return(NULL)
  dv <- vapply(res, function(r) r$dl, numeric(1))
  wv <- vapply(res, function(r) r$n, numeric(1))
  ov <- vapply(res, function(r) r$ocr, numeric(1))
  cv0v <- vapply(res, function(r) r$cv0, numeric(1))
  bo <- boot(dv, function(dd, idx) sum(dd[idx]*wv[idx]/sum(wv[idx])), R = R)
  md <- bo$t0; se <- sd(bo$t)
  pb <- if (is.finite(se) && se > 0) pnorm(md/se, lower.tail = FALSE) else 1
  kr2 <- (dv > 0) | (ov >= 0.5)
  data.frame(p_boot = pb, md = md,
             p_sign = binom.test(sum(dv > 0), length(dv), 0.5, "greater")$p.value,
             p_or = binom.test(sum(kr2), length(kr2), 0.5, "greater")$p.value,
             poz = sum(dv > 0), M = length(dv), med_dl = median(dv),
             cv0_ort = mean(cv0v), stringsAsFactors = FALSE)
}

t2_hucre <- function(ad, rep = 100, R = 500, ...) {
  t0 <- Sys.time()
  out <- vector("list", rep); atla <- 0
  for (i in seq_len(rep)) {
    r <- NULL
    while (is.null(r)) {
      r <- t2_svt(t2_veri(...), R = R)
      if (is.null(r)) atla <- atla + 1
    }
    out[[i]] <- r
  }
  H <- do.call(rbind, out)
  H$boot_g <- H$p_boot < 0.05 & H$md > 0
  H$sign_g <- H$p_sign < 0.05
  H$or_g <- H$p_or < 0.05
  dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
  oz <- data.frame(hucre = ad, boot = mean(H$boot_g), sign = mean(H$sign_g),
                   or_ = mean(H$or_g), cift_or = mean(H$boot_g & H$or_g),
                   p_poz = round(mean(H$poz/H$M), 3),
                   med_dl = round(median(H$med_dl), 3),
                   cv0 = round(median(H$cv0_ort), 1), atla = atla, dk = dk)
  cat(sprintf("T2 %-12s boot=%.3f or=%.3f cift=%.3f  p+=%.3f meddl=%+.3f  cv0=%.1f  atla=%d  %.1f dk\n",
              ad, oz$boot, oz$or_, oz$cift_or, oz$p_poz, oz$med_dl, oz$cv0, atla, dk))
  oz
}

T2 <- rbind(
  t2_hucre("null_r015", rho_m = 0.15, s_c = 0, lam = 0),
  t2_hucre("null_r080", rho_m = 0.80, s_c = 0, lam = 0),
  t2_hucre("guc_r030", rho_m = 0.30, s_c = 0.4, lam = 0.25))
write.csv(T2, paste0(kok, "mc_tier2.csv"), row.names = FALSE)
cat("\nTIER 2\n"); print(T2, row.names = FALSE, right = FALSE)

rds <- readRDS(paste0(kok, "svt_tum_sonuclar.rds"))
dt <- rds$veri; al <- rds$alt
ciftler <- list(c("ERISIM","MEMNUNIYET"), c("ERISIM","MTT"), c("ERISIM","FCV"),
                c("ERISIM","SOY_index"), c("MEMNUNIYET","SOY_index"))
fk <- do.call(rbind, lapply(ciftler, function(p) {
  v <- vapply(c("Tum orneklem", "Dikkatli", "Dikkatsiz"), function(a)
    cor(dt[[p[1]]][al[[a]]], dt[[p[2]]][al[[a]]], method = "spearman"), numeric(1))
  data.frame(yol = paste(p[1], "->", p[2]), tum = round(v[1], 3),
             dikkatli = round(v[2], 3), dikkatsiz = round(v[3], 3),
             fark = round(v[3] - v[2], 3))
}))
cat("\nSAHA KONTROLU: korelasyonlar altorneklemlere gore\n")
print(fk, row.names = FALSE, right = FALSE)


mi_profil <- function(etiket) {
  do.call(rbind, lapply(names(gost), function(x) do.call(rbind, lapply(tum_mod, function(g) {
    r <- mi_ham(x, g)
    if (is.null(r)) return(NULL)
    dl <- r$delta
    asima <- rownames(dl)[apply(dl > matrix(esik, 3, 4, byrow = TRUE), 1, any)]
    data.frame(orneklem = etiket, olcek = x, moderator = g,
               metric_cfi = round(dl["metric", "cfi"], 4),
               scalar_cfi = round(dl["scalar", "cfi"], 4),
               strict_cfi = round(dl["strict", "cfi"], 4),
               metric_rmsea = round(dl["metric", "rmsea"], 4),
               scalar_rmsea = round(dl["scalar", "rmsea"], 4),
               strict_rmsea = round(dl["strict", "rmsea"], 4),
               ihlal_asamalari = if (length(asima)) paste(asima, collapse = "+") else "yok",
               invaryant = r$invaryant, stringsAsFactors = FALSE)
  }))))
}

d <- d_tam
P_tam <- mi_profil("tam")
d <- d_tam[alt[["Dikkatli"]], ]
P_dik <- mi_profil("dikkatli")
d <- d_tam
PROF <- rbind(P_tam, P_dik)
write.csv(PROF, paste0(kok, "mi_asama_profili.csv"), row.names = FALSE)
print(PROF, row.names = FALSE, right = FALSE)
cat("\nasama bazinda ihlal sayisi:\n")
print(table(unlist(strsplit(PROF$ihlal_asamalari[PROF$ihlal_asamalari != "yok"], "\\+"))))

kok     <- "data/derived/"
ham_kok <- Sys.getenv("SVT_RAW_DIR", unset = "raw-data/")
ad <- c(uclu = "svt_uclu_sonuclari.csv", moderator = "svt_moderator_duzeyi.csv",
        mi_profil = "mi_asama_profili.csv", fazA_ozet = "mc_fazA_ozet.csv",
        fazA_ham = "mc_fazA_ham.csv", fazB = "mc_fazB.csv", fazB2 = "mc_fazB2.csv",
        fazC1 = "mc_fazC1r.csv", fazC2 = "mc_fazC2r.csv", fazE = "mc_fazE.csv",
        tier2 = "mc_tier2.csv")
oku <- function(f) if (file.exists(paste0(kok, f))) read.csv(paste0(kok, f), stringsAsFactors = FALSE) else NULL
paket <- lapply(ad, oku)
eksik <- names(ad)[vapply(paket, is.null, logical(1))]
if (length(eksik)) cat("eksik dosya:", paste(eksik, collapse = ", "), "\n")
paket$bilgi <- list(seed = 9186, boot_R_sim = 500, boot_R_ampirik = 1000,
                    rep = c(fazA = 1000, fazB = 1000, fazB2 = 1000, fazC = 1000,
                            fazE = 300, tier2 = 100),
                    tarih = format(Sys.time()), R_surum = R.version.string,
                    lavaan = as.character(packageVersion("lavaan")))
saveRDS(paket, paste0(kok, "mc_tum_sonuclar.rds"))
cat("mc_tum_sonuclar.rds yazildi:", length(paket) - 1, "tablo + bilgi\n")




cv_lat <- function(dd, xolcek, yad, min_n = 25) {
  md <- paste0("F =~ ", paste(gost[[xolcek]], collapse = " + "), "\n", yad, " ~ c*F")
  cvv <- vapply(tum_mod, function(g) {
    b <- nn <- numeric(0)
    for (l in sort(unique(dd[[g]]))) {
      sub <- dd[dd[[g]] == l, ]
      if (nrow(sub) < min_n) next
      f <- tryCatch(suppressWarnings(sem(md, sub, std.lv = TRUE)), error = function(e) NULL)
      if (is.null(f) || !lavInspect(f, "converged")) next
      q <- standardizedSolution(f)
      v <- q$est.std[q$label == "c"][1]
      if (is.na(v)) next
      b <- c(b, v); nn <- c(nn, nrow(sub))
    }
    if (length(b) < 2) return(NA_real_)
    wm <- weighted.mean(b, nn)
    ws <- sqrt(sum(nn*(b - wm)^2)/sum(nn))
    100*ws/max(abs(wm), 1e-4)
  }, numeric(1))
  mean(cvv, na.rm = TRUE)
}

yollar <- list(c("ERISIM", "MEMNUNIYET"), c("ERISIM", "MTT"), c("ERISIM", "FCV"),
               c("ERISIM", "SOY_index"), c("MEMNUNIYET", "SOY_index"))
altlar <- c("Tum orneklem", "Dikkatli", "Dikkatsiz")

cat("A) TABLO 1 PANEL B - LATENT PIPELINE\n")
T1B <- do.call(rbind, lapply(yollar, function(p) {
  v <- vapply(altlar, function(a) cv_lat(d_tam[alt[[a]], ], p[1], p[2]), numeric(1))
  data.frame(yol = paste(p[1], "->", p[2]), tum = round(v[1], 1),
             dikkatli = round(v[2], 1), dikkatsiz = round(v[3], 1))
}))
print(T1B, row.names = FALSE, right = FALSE)
write.csv(T1B, paste0(kok, "tablo1_panelB_latent.csv"), row.names = FALSE)


bilinen <- c(mad_e, mad_g, mad_s)
aday <- grep("^[A-Za-z]+[0-9]+$", names(d_tam), value = TRUE)
aday <- setdiff(aday, bilinen)
aday <- aday[vapply(d_tam[aday], is.numeric, logical(1))]
cat("aday kolonlar (", length(aday), "):\n")
print(aday)
r_m <- if (length(aday)) cor(rowMeans(d_tam[, aday]), d_tam$MEMNUNIYET) else NA
cat("rowMeans(aday) ~ MEMNUNIYET  r =", round(r_m, 4), "\n")
cat("rowMeans(GUV)  ~ GUVEN       r =", round(cor(rowMeans(d_tam[, mad_g]), d_tam$GUVEN), 4), "\n")
cat("rowMeans(SOY)  ~ SOY_index   r =", round(cor(rowMeans(d_tam[, mad_s]), d_tam$SOY_index), 4), "\n")
if (is.na(r_m) || r_m < 0.999 || length(aday) != 26) {
  cat("UYARI: aday seti 26 memnuniyet maddesiyle eslesmedi, listeden elle sec\n")
} else {
  mad_m <- aday
  mad78 <- c(mad_e, mad_m, mad_g, mad_s)
  cat("mad78 toplam:", length(mad78), "\n\n")
  
  cat("B) HARMAN - GERCEK TEK FAKTOR FA vs PCA (78 madde)\n")
  for (a in altlar) {
    dd <- d_tam[alt[[a]], mad78]
    pca <- eigen(cor(dd, use = "pairwise.complete.obs"))$values
    fa1 <- tryCatch(suppressWarnings(fa(dd, nfactors = 1, fm = "minres")),
                    error = function(e) NULL)
    fp <- if (is.null(fa1)) NA else 100*fa1$Vaccounted["Proportion Var", 1]
    cat(sprintf("%-14s PCA %%%.1f   FA %%%s\n", a, 100*pca[1]/length(pca),
                ifelse(is.na(fp), "hata", sprintf("%.1f", fp))))
  }
  
  cat("\nC) ESIK DUYARLILIK TARAMASI\n")
  DUY <- do.call(rbind, lapply(c(20, 25, 30, 35), function(ls_k) do.call(rbind, lapply(c(0, 1), function(db_k) {
    sec <- d_tam$Longstring < ls_k & d_tam$DuzBlok <= db_k & d_tam$FarkliSecenek >= 3
    dd <- d_tam[sec, ]
    ev <- eigen(cor(dd[, mad_m], use = "pairwise.complete.obs"))$values
    fa1 <- tryCatch(suppressWarnings(fa(dd[, mad78], nfactors = 1, fm = "minres")),
                    error = function(e) NULL)
    fp <- if (is.null(fa1)) NA else round(100*fa1$Vaccounted["Proportion Var", 1], 1)
    out <- data.frame(LS = ls_k, DB = db_k, n = sum(sec),
                      alpha_M = round(suppressWarnings(psych::alpha(dd[, mad_m])$total$raw_alpha), 3),
                      ev_oran = round(ev[1]/ev[2], 1), fa1_pct = fp,
                      cv_EM = round(cv_lat(dd, "ERISIM", "MEMNUNIYET"), 1))
    cat(sprintf("LS<%d DB<=%d  n=%3d  alphaM=%.3f  ev=%4.1f  fa1=%s  cvEM=%.1f\n",
                ls_k, db_k, out$n, out$alpha_M, out$ev_oran,
                ifelse(is.na(fp), "hata", paste0("%", fp)), out$cv_EM))
    out
  }))))
  write.csv(DUY, paste0(kok, "esik_duyarlilik.csv"), row.names = FALSE)
}



set.seed(9186)

if (!exists("mc_gruplar")) {
  mc_gruplar <- function(n, bag = 1) {
    u <- rnorm(n)
    kar <- function(w) bag*w*u + rnorm(n, 0, sqrt(max(0.02, 1 - (bag*w)^2)))
    list(G = data.frame(
      g1 = 1 + rbinom(n, 1, 0.37),
      g2 = 1 + as.integer(kar(0.7) > 0.15),
      g3 = as.integer(cut(kar(0.6), c(-Inf, -0.85, 0.75, Inf))),
      g4 = 1 + rbinom(n, 1, 0.35),
      g5 = as.integer(cut(rnorm(n), c(-Inf, -1.15, 0.95, Inf))),
      g6 = as.integer(cut(kar(0.8), c(-Inf, -1.25, -0.55, 0.25, 1.05, Inf)))),
      u = u)
  }
  mc_veri <- function(n = 500, rho_m = 0, rho_t = 0.60, r_b = 0.55, sd_het = 0.18,
                      lam_zx = 0, stab = 0, s_a = 0.5, bag = 1) {
    gu <- mc_gruplar(n, bag)
    G <- gu$G
    tx <- rnorm(n)
    del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
    bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
    z <- lam_zx*tx + sqrt(1 - lam_zx^2)*rnorm(n)
    ytrue <- as.numeric(scale(bi*tx + stab*s_a*gu$u*z + sqrt(1 - rho_t^2)*rnorm(n)))
    s <- rnorm(n)
    mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
    my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
    data.frame(X = sqrt(1 - rho_m)*tx + sqrt(rho_m)*mx,
               Y = sqrt(1 - rho_m)*ytrue + sqrt(rho_m)*my,
               Z = as.numeric(scale(z)), G)
  }
  mc_stab <- function(dat, g, min_n = 25) {
    b0 <- b1 <- nn <- numeric(0)
    for (l in sort(unique(dat[[g]]))) {
      sub <- dat[dat[[g]] == l, ]
      if (nrow(sub) < min_n) next
      sx <- as.numeric(scale(sub$X)); sy <- as.numeric(scale(sub$Y)); sz <- as.numeric(scale(sub$Z))
      f0 <- coef(lm(sy ~ sx))[2]
      f1 <- coef(lm(sy ~ sx + sz))[2]
      if (any(is.na(c(f0, f1)))) next
      b0 <- c(b0, f0); b1 <- c(b1, f1); nn <- c(nn, nrow(sub))
    }
    if (length(b0) < 2) return(NULL)
    wm <- function(v) weighted.mean(v, nn)
    ws <- function(v) sqrt(sum(nn*(v - wm(v))^2)/sum(nn))
    mu0 <- wm(b0); mu1 <- wm(b1); sd0 <- ws(b0); sd1 <- ws(b1)
    cv0 <- 100*sd0/max(abs(mu0), 1e-4)
    cv1 <- 100*sd1/max(abs(mu1), 1e-4)
    list(dl = log(cv0 + 1e-8) - log(cv1 + 1e-8), cv0 = cv0, cv1 = cv1,
         ocr = mean(as.integer(abs(b1 - mu1) < abs(b0 - mu0))), n = sum(nn))
  }
  mc_svt <- function(dat, R = 500) {
    res <- lapply(paste0("g", 1:6), function(g) mc_stab(dat, g))
    res <- res[!vapply(res, is.null, logical(1))]
    if (length(res) < 3) return(NULL)
    dv <- vapply(res, function(r) r$dl, numeric(1))
    wv <- vapply(res, function(r) r$n, numeric(1))
    ov <- vapply(res, function(r) r$ocr, numeric(1))
    cv0v <- vapply(res, function(r) r$cv0, numeric(1))
    bo <- boot(dv, function(dd, idx) sum(dd[idx]*wv[idx]/sum(wv[idx])), R = R)
    md <- bo$t0; se <- sd(bo$t)
    pb <- if (is.finite(se) && se > 0) pnorm(md/se, lower.tail = FALSE) else 1
    kr2 <- (dv > 0) | (ov >= 0.5)
    data.frame(p_boot = pb, md = md,
               p_sign = binom.test(sum(dv > 0), length(dv), 0.5, "greater")$p.value,
               p_or = binom.test(sum(kr2), length(kr2), 0.5, "greater")$p.value,
               poz = sum(dv > 0), M = length(dv), med_dl = median(dv),
               cv0_ort = mean(cv0v), cv0_min = min(cv0v),
               isaret = paste(as.integer(dv > 0), collapse = ""),
               stringsAsFactors = FALSE)
  }
  mc_hucre <- function(rep = 1000, R = 500, ...) {
    t0 <- Sys.time()
    out <- vector("list", rep)
    for (i in seq_len(rep)) {
      r <- NULL
      while (is.null(r)) r <- mc_svt(mc_veri(...), R = R)
      out[[i]] <- r
    }
    H <- do.call(rbind, out)
    H$boot_g <- H$p_boot < 0.05 & H$md > 0
    H$sign_g <- H$p_sign < 0.05
    H$or_g <- H$p_or < 0.05
    attr(H, "dk") <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
    H
  }
}

uyum_ozet <- function(H) {
  h6 <- H[H$M == 6, ]
  S <- do.call(rbind, lapply(strsplit(h6$isaret, ""), as.integer))
  p_m <- colMeans(S)
  fark <- numeric(0)
  for (i in 1:5) for (j in (i+1):6)
    fark <- c(fark, mean(S[, i] == S[, j]) - (p_m[i]*p_m[j] + (1-p_m[i])*(1-p_m[j])))
  ph <- mean(rowSums(S))/6
  c(DI = var(rowSums(S))/(6*ph*(1-ph)), ort_fark = mean(fark), maks_fark = max(fark))
}

cat("D1) BAGIMSIZ MODERATOR KOLU (bag=0)\n")
D1 <- do.call(rbind, lapply(c(0.15, 0.80), function(r) {
  H <- mc_hucre(rho_m = r, bag = 0)
  u <- uyum_ozet(H)
  out <- data.frame(rho_m = r, boot = mean(H$boot_g), sign = mean(H$sign_g),
                    cift_or = mean(H$boot_g & H$or_g), p_poz = round(mean(H$poz/H$M), 3),
                    DI = round(u["DI"], 3), ort_fark = round(u["ort_fark"], 4),
                    maks_fark = round(u["maks_fark"], 4), dk = attr(H, "dk"))
  cat(sprintf("bag0 rho=%.2f  boot=%.3f sign=%.3f cift=%.3f  p+=%.3f  DI=%.3f  fark(ort/maks)=%.4f/%.4f  %.1f dk\n",
              r, out$boot, out$sign, out$cift_or, out$p_poz, out$DI, out$ort_fark, out$maks_fark, attr(H, "dk")))
  out
}))
write.csv(D1, paste0(kok, "mc_fazD1_bag0.csv"), row.names = FALSE)

mc_veri_mask <- function(n = 500, rho_base = 0.45, d_rho = 0, rho_t = 0.60,
                         r_b = 0.55, sd_het = 0.18) {
  gu <- mc_gruplar(n)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  ytrue <- as.numeric(scale(bi*tx + sqrt(1 - rho_t^2)*rnorm(n)))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  rho_i <- pmin(pmax(rho_base + d_rho*c(-1, 0, 1)[G$g3], 0.05), 0.90)
  L <- abs(s)*sqrt(rho_i) + rnorm(n, 0, 0.3)
  data.frame(X = sqrt(1 - rho_i)*tx + sqrt(rho_i)*mx,
             Y = sqrt(1 - rho_i)*ytrue + sqrt(rho_i)*my,
             Z = as.numeric(scale(sqrt(rho_i)*s + sqrt(1 - rho_i)*rnorm(n))),
             G, L = L)
}

cat("\nD2) MASKELENME FAZI (grup-degisken rho_M)\n")
D2 <- do.call(rbind, lapply(c(0, 0.10, 0.20, 0.30), function(dr) {
  t0 <- Sys.time()
  out <- vector("list", 1000); e2 <- numeric(1000)
  for (i in seq_len(1000)) {
    v <- NULL
    while (is.null(attr(v, "ok"))) {
      dat <- mc_veri_mask(d_rho = dr)
      v <- mc_svt(dat[, 1:9])
      if (!is.null(v)) attr(v, "ok") <- TRUE
    }
    out[[i]] <- v
    a <- anova(lm(dat$L ~ factor(dat$g3)))
    e2[i] <- a[1, 2]/sum(a[, 2])
  }
  H <- do.call(rbind, out)
  H$boot_g <- H$p_boot < 0.05 & H$md > 0
  H$sign_g <- H$p_sign < 0.05
  H$or_g <- H$p_or < 0.05
  dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
  o <- data.frame(d_rho = dr, boot = mean(H$boot_g), sign = mean(H$sign_g),
                  or_ = mean(H$or_g), cift_sign = mean(H$boot_g & H$sign_g),
                  cift_or = mean(H$boot_g & H$or_g),
                  p_poz = round(mean(H$poz/H$M), 3), med_dl = round(median(H$med_dl), 3),
                  eta2_L = round(mean(e2), 3), dk = dk)
  cat(sprintf("d_rho=%.2f  boot=%.3f sign=%.3f ciftS=%.3f ciftOR=%.3f  p+=%.3f meddl=%+.3f  eta2L=%.3f  %.1f dk\n",
              dr, o$boot, o$sign, o$cift_sign, o$cift_or, o$p_poz, o$med_dl, o$eta2_L, dk))
  o
}))
write.csv(D2, paste0(kok, "mc_fazD2_mask.csv"), row.names = FALSE)

mc_veri_uneq <- function(n = 500, rho_x = 0.5, rho_y = 0.5, rho_t = 0.60,
                         r_b = 0.55, sd_het = 0.18) {
  gu <- mc_gruplar(n)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  ytrue <- as.numeric(scale(bi*tx + sqrt(1 - rho_t^2)*rnorm(n)))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_x)*tx + sqrt(rho_x)*mx,
             Y = sqrt(1 - rho_y)*ytrue + sqrt(rho_y)*my,
             Z = as.numeric(scale(rnorm(n))), G)
}

cat("\nD3) ESITSIZ PAY DOGRULAMASI (a = sqrt((1-rx)(1-ry)))\n")
D3 <- do.call(rbind, lapply(list(c(.2, .8), c(.8, .2), c(.6, .6), c(.5, .5)), function(p) {
  t0 <- Sys.time()
  out <- vector("list", 1000); rxy <- numeric(1000)
  for (i in seq_len(1000)) {
    r <- NULL
    while (is.null(r)) { v <- mc_veri_uneq(rho_x = p[1], rho_y = p[2]); r <- mc_svt(v) }
    out[[i]] <- r; rxy[i] <- cor(v$X, v$Y)
  }
  H <- do.call(rbind, out)
  H$boot_g <- H$p_boot < 0.05 & H$md > 0
  H$or_g <- H$p_or < 0.05
  dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
  o <- data.frame(rho_x = p[1], rho_y = p[2], a = round(sqrt((1-p[1])*(1-p[2])), 3),
                  r_xy = round(mean(rxy), 3), cv0 = round(median(H$cv0_ort), 1),
                  boot = mean(H$boot_g), cift_or = mean(H$boot_g & H$or_g), dk = dk)
  cat(sprintf("rx=%.1f ry=%.1f  a=%.3f  rxy=%.3f  cv0=%.1f  boot=%.3f cift=%.3f  %.1f dk\n",
              p[1], p[2], o$a, o$r_xy, o$cv0, o$boot, o$cift_or, dk))
  o
}))
write.csv(D3, paste0(kok, "mc_fazD3_uneq.csv"), row.names = FALSE)

cat("\nD4) TIER-2, 500 TEKRAR (uzun surer)\n")
if (exists("t2_hucre")) {
  T2b <- rbind(
    t2_hucre("null_r015", rep = 500, rho_m = 0.15, s_c = 0, lam = 0),
    t2_hucre("null_r080", rep = 500, rho_m = 0.80, s_c = 0, lam = 0),
    t2_hucre("guc_r030", rep = 500, rho_m = 0.30, s_c = 0.4, lam = 0.25))
  write.csv(T2b, paste0(kok, "mc_tier2_500.csv"), row.names = FALSE)
  print(T2b, row.names = FALSE, right = FALSE)
} else cat("t2 fonksiyonlari yuklu degil, tier-2 icin onceki oturum kodunu calistir\n")


set.seed(9186)

cat("D1r) BAGIMSIZ MODERATOR KOLU, rep=1000\n")
D1r <- do.call(rbind, lapply(c(0.15, 0.80), function(r) {
  H <- mc_hucre(rep = 1000, R = 500, rho_m = r, bag = 0)
  u <- uyum_ozet(H)
  out <- data.frame(rho_m = r, boot = mean(H$boot_g), sign = mean(H$sign_g),
                    or_ = mean(H$or_g), cift_or = mean(H$boot_g & H$or_g),
                    p_poz = round(mean(H$poz/H$M), 3), DI = round(u["DI"], 3),
                    ort_fark = round(u["ort_fark"], 4), maks_fark = round(u["maks_fark"], 4),
                    dk = attr(H, "dk"))
  cat(sprintf("bag0 rho=%.2f  boot=%.3f sign=%.3f or=%.3f cift=%.3f  p+=%.3f  DI=%.3f  fark=%.4f/%.4f  %.1f dk\n",
              r, out$boot, out$sign, out$or_, out$cift_or, out$p_poz, out$DI,
              out$ort_fark, out$maks_fark, out$dk))
  out
}))
write.csv(D1r, paste0(kok, "mc_fazD1_bag0_v2.csv"), row.names = FALSE)

mc_veri_mask2 <- function(n = 500, rho_base = 0.45, d_rho = 0, rho_t = 0.35,
                          r_b = 0.80, sd_het = 0.18, phi = 0.5) {
  gu <- mc_gruplar(n)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  ytrue <- as.numeric(scale(bi*tx + sqrt(1 - rho_t^2)*rnorm(n)))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  rho_i <- pmin(pmax(rho_base + d_rho*c(-1, 0, 1)[G$g3], 0.05), 0.90)
  L <- abs(s)*sqrt(rho_i) + rnorm(n, 0, 0.3)
  data.frame(X = sqrt(1 - rho_i)*tx + sqrt(rho_i)*mx,
             Y = sqrt(1 - rho_i)*ytrue + sqrt(rho_i)*my,
             Z = as.numeric(scale(phi*s + sqrt(1 - phi^2)*rnorm(n))),
             G, L = L)
}

cat("\nD2r) MASKELENME, plato-trait ayrismis (rho_t=.35, r_b=.80)\n")
hucre_mask <- function(dr, phi) {
  t0 <- Sys.time()
  out <- vector("list", 1000); e2 <- numeric(1000)
  for (i in seq_len(1000)) {
    v <- NULL
    while (is.null(v)) { dat <- mc_veri_mask2(d_rho = dr, phi = phi); v <- mc_svt(dat[, 1:9], R = 500) }
    out[[i]] <- v
    a <- anova(lm(dat$L ~ factor(dat$g3)))
    e2[i] <- a[1, 2]/sum(a[, 2])
  }
  H <- do.call(rbind, out)
  H$boot_g <- H$p_boot < 0.05 & H$md > 0
  H$sign_g <- H$p_sign < 0.05
  H$or_g <- H$p_or < 0.05
  dk <- round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2)
  o <- data.frame(d_rho = dr, phi = phi, boot = mean(H$boot_g), sign = mean(H$sign_g),
                  or_ = mean(H$or_g), cift_sign = mean(H$boot_g & H$sign_g),
                  cift_or = mean(H$boot_g & H$or_g), p_poz = round(mean(H$poz/H$M), 3),
                  med_dl = round(median(H$med_dl), 3), cv0 = round(median(H$cv0_ort), 1),
                  eta2_L = round(mean(e2), 3), dk = dk)
  cat(sprintf("d_rho=%.2f phi=%.1f  boot=%.3f ciftS=%.3f ciftOR=%.3f  p+=%.3f meddl=%+.3f  cv0=%.1f  eta2L=%.3f  %.1f dk\n",
              dr, phi, o$boot, o$cift_sign, o$cift_or, o$p_poz, o$med_dl, o$cv0, o$eta2_L, dk))
  o
}
D2r <- rbind(hucre_mask(0, 0.5), hucre_mask(0.10, 0.5), hucre_mask(0.20, 0.5),
             hucre_mask(0.30, 0.5), hucre_mask(0.30, 0))
write.csv(D2r, paste0(kok, "mc_fazD2_mask_v2.csv"), row.names = FALSE)

mc_veri_uneq2 <- function(n = 500, rho_x = 0.5, rho_y = 0.5, rho_t = 0.60,
                          r_b = 0.55, sd_het = 0.18) {
  gu <- mc_gruplar(n)
  G <- gu$G
  tx <- rnorm(n)
  del <- lapply(1:6, function(j) rnorm(max(G[[j]]), 0, sd_het))
  bi <- rho_t + rowSums(vapply(1:6, function(j) del[[j]][G[[j]]], numeric(n)))
  ytrue <- as.numeric(scale(bi*tx + sqrt(1 - rho_t^2)*rnorm(n)))
  s <- rnorm(n)
  mx <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  my <- sqrt(r_b)*s + sqrt(1 - r_b)*rnorm(n)
  data.frame(X = sqrt(1 - rho_x)*tx + sqrt(rho_x)*mx,
             Y = sqrt(1 - rho_y)*ytrue + sqrt(rho_y)*my, G)
}

d3_sd <- function(dat, min_n = 25) {
  vv <- vapply(paste0("g", 1:6), function(g) {
    b <- nn <- numeric(0)
    for (l in sort(unique(dat[[g]]))) {
      sub <- dat[dat[[g]] == l, ]
      if (nrow(sub) < min_n) next
      b <- c(b, coef(lm(scale(sub$Y) ~ scale(sub$X)))[2]); nn <- c(nn, nrow(sub))
    }
    if (length(b) < 2) return(c(NA, NA))
    wm <- weighted.mean(b, nn)
    c(sum(nn*(b - wm)^2)/sum(nn), wm)
  }, numeric(2))
  c(v = mean(vv[1, ], na.rm = TRUE), mu = mean(vv[2, ], na.rm = TRUE))
}

cat("\nD3r) SINYAL CIKARIMLI ESITSIZ PAY (het ikizleriyle)\n")
konf <- list(c(.2, .8), c(.8, .2), c(.6, .6), c(.5, .5))
D3r <- do.call(rbind, lapply(konf, function(p) {
  vh <- v0 <- mu <- numeric(600)
  for (i in seq_len(600)) {
    a1 <- d3_sd(mc_veri_uneq2(rho_x = p[1], rho_y = p[2], sd_het = 0.18))
    a0 <- d3_sd(mc_veri_uneq2(rho_x = p[1], rho_y = p[2], sd_het = 0))
    vh[i] <- a1["v"]; v0[i] <- a0["v"]; mu[i] <- a1["mu"]
  }
  sinyal <- sqrt(max(mean(vh, na.rm = TRUE) - mean(v0, na.rm = TRUE), 0))
  o <- data.frame(rho_x = p[1], rho_y = p[2], a = round(sqrt((1-p[1])*(1-p[2])), 3),
                  mu = round(mean(mu, na.rm = TRUE), 3),
                  sd_toplam = round(sqrt(mean(vh, na.rm = TRUE)), 4),
                  sd_taban = round(sqrt(mean(v0, na.rm = TRUE)), 4),
                  sinyal_sd = round(sinyal, 4),
                  oran = round(sinyal/sqrt((1-p[1])*(1-p[2])), 4))
  cat(sprintf("rx=%.1f ry=%.1f  a=%.3f  mu=%.3f  toplam=%.4f taban=%.4f  sinyal=%.4f  sinyal/a=%.4f\n",
              p[1], p[2], o$a, o$mu, o$sd_toplam, o$sd_taban, o$sinyal_sd, o$oran))
  o
}))
write.csv(D3r, paste0(kok, "mc_fazD3_uneq_v2.csv"), row.names = FALSE)


kok     <- "data/derived/"
ham_kok <- Sys.getenv("SVT_RAW_DIR", unset = "raw-data/")
## NOTE: this archive embeds participant-level objects and must never be
## committed to the public repository; it is written outside the repo tree.
hedef <- paste0(ham_kok, "svt_cmv_arsiv.rds")

oku <- function(f) { yol <- paste0(kok, f); if (file.exists(yol)) read.csv(yol, stringsAsFactors = FALSE) else NULL }
ad <- c(uclu = "svt_uclu_sonuclari.csv", moderator = "svt_moderator_duzeyi.csv",
        mi_profil = "mi_asama_profili.csv", t1_latent = "tablo1_panelB_latent.csv",
        esik = "esik_duyarlilik.csv", fazA_ozet = "mc_fazA_ozet.csv",
        fazA_ham = "mc_fazA_ham.csv", fazB = "mc_fazB.csv", fazB2 = "mc_fazB2.csv",
        fazC1 = "mc_fazC1r.csv", fazC2 = "mc_fazC2r.csv", fazE = "mc_fazE.csv",
        tier2 = "mc_tier2.csv", tier2_500 = "mc_tier2_500.csv",
        fazD1 = "mc_fazD1_bag0_v2.csv", fazD2 = "mc_fazD2_mask_v2.csv",
        fazD3 = "mc_fazD3_uneq_v2.csv")
tablolar <- lapply(ad, oku)
eksik <- names(ad)[vapply(tablolar, is.null, logical(1))]
if (length(eksik)) cat("eksik csv:", paste(eksik, collapse = ", "), "\n")
tablolar <- tablolar[!vapply(tablolar, is.null, logical(1))]

g <- globalenv()
isim <- ls(g)
fn_ad <- isim[vapply(isim, function(a) is.function(get(a, envir = g)), logical(1))]
fn_ad <- setdiff(fn_ad, c("oku"))
fonksiyonlar <- mget(fn_ad, envir = g)

al <- function(a) if (exists(a, envir = g)) get(a, envir = g) else NULL
veri <- list(d_tam = al("d_tam"), alt = al("alt"), gost = al("gost"),
             tum_mod = al("tum_mod"), mad_e = al("mad_e"), mad_m = al("mad_m"),
             mad_g = al("mad_g"), mad_s = al("mad_s"))
veri <- veri[!vapply(veri, is.null, logical(1))]

surum <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_
bilgi <- list(seed = 9186, tarih = format(Sys.time()), R = R.version.string,
              lavaan = surum("lavaan"), psych = surum("psych"), boot = surum("boot"),
              rep = c(fazA = 1000, fazB = 1000, fazB2 = 1000, fazC = 1000,
                      fazD1 = 1000, fazD2 = 1000, fazD3 = 600, fazE = 300, tier2 = 500))

paket <- list(tablolar = tablolar, fonksiyonlar = fonksiyonlar, veri = veri, bilgi = bilgi)
if (file.exists(hedef)) {
  eski <- readRDS(hedef)
  paket$tablolar <- modifyList(eski$tablolar, paket$tablolar)
  paket$fonksiyonlar <- modifyList(eski$fonksiyonlar, paket$fonksiyonlar)
  if (length(eski$veri)) paket$veri <- modifyList(eski$veri, paket$veri)
}
saveRDS(paket, hedef)
cat("svt_cmv_arsiv.rds:", length(paket$tablolar), "tablo,",
    length(paket$fonksiyonlar), "fonksiyon,", length(paket$veri), "veri nesnesi,",
    round(file.size(hedef)/1e6, 1), "MB\n")

