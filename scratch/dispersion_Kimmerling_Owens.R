d1 <- readRDS("output/filtered_data_Kimmerling_threshold1.rds")
r1 <- d1$relative
g1 <- d1$groups

d2 <- readRDS("output/filtered_data_Owens_threshold1.rds")
r2 <- d2$relative
g2 <- d2$groups

temp <- data.frame(x = apply(r1[g1 == levels(g1)[1],], 2, function(x) 1/sd(x)),
                   y = "r1_A")
temp <- rbind(temp,
              data.frame(x = apply(r1[g1 == levels(g1)[2],], 2, function(x) 1/sd(x)),
                         y = "r1_B"))
temp <- rbind(temp,
              data.frame(x = apply(r2[g2 == levels(g2)[1],], 2, function(x) 1/sd(x)),
                         y = "r2_A"))
temp <- rbind(temp,
              data.frame(x = apply(r2[g2 == levels(g2)[2],], 2, function(x) 1/sd(x)),
                         y = "r2_B"))

ggplot(temp, aes(x = x, fill = y)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  xlim(c(0,2))
