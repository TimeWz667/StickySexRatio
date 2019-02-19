burden <- read.csv("Input/WHO_TB_burden_2018-10-31.csv")
notification <- read.csv("Input/WHO_TB_notifications_2018-10-31.csv")

data.cdr <- burden[c("country", "year", "g_whoregion", "c_newinc_100k", 
                     "e_inc_100k", "e_inc_100k_lo", "e_inc_100k_hi",
                     "c_cdr", "c_cdr_lo", "c_cdr_hi")]

names(data.cdr) <- c("Country", "Year", "Region", "Notif", "Inc", "Inc_L", "Inc_H", "CDR", "CDR_L", "CDR_H")



# data.noti <- notification[c("country", "year", "c_newinc")]


data.cdr$CDR[data.cdr$CDR <= 0] <- NA
data.cdr$CDR_L[data.cdr$CDR_L <= 0] <- NA
data.cdr$CDR_H[data.cdr$CDR_H <= 0] <- NA

data.cdr$CDR[data.cdr$CDR > 100] <- NA
data.cdr$CDR_L[data.cdr$CDR_L > 100] <- NA
data.cdr$CDR_H[data.cdr$CDR_H > 100] <- NA



noti <- notification #notification[notification$country == "Argentina",]

gps <- c("514", "014", "1524", "2534", "3544", "4554", "5564", "65", "u")

noti = data.frame(Year=noti$year, Country=noti$country, notif=noti$c_newinc, 
                  m=rowSums(noti[, paste0("new_sp_", "m", gps)], na.rm = T) + 
                    rowSums(noti[, paste0("new_sn_", "m", gps)], na.rm = T) + 
                    rowSums(noti[, paste0("new_ep_", "m", gps)], na.rm = T) + 
                    rowSums(noti[, paste0("newrel_", "m", gps)], na.rm = T),
                  f=rowSums(noti[, paste0("new_sp_", "f", gps)], na.rm = T) + 
                    rowSums(noti[, paste0("new_sn_", "f", gps)], na.rm = T) + 
                    rowSums(noti[, paste0("new_ep_", "f", gps)], na.rm = T) + 
                    rowSums(noti[, paste0("newrel_", "f", gps)], na.rm = T))

noti$add <- noti$f + noti$m
noti <- subset(noti, add>0 & m>0 & f>0)
  
noti$m <- noti$m / noti$add
noti$f <- noti$f / noti$add

noti$lf <- log(noti$f / (noti$add - noti$f))
noti$lfm <- log(noti$f) - log(noti$m)
noti <- noti[abs(noti$add-noti$notif) < 0.1*noti$notif, c("Year", "Country", "m", "f", "lf", "lfm")]
noti <- merge(data.cdr, noti)

save(noti, file="Input/WHO_TB.rdata")
