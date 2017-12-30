library(car)
library(corrplot)

head(prestige)
dim(prestige)
str(prestige)

prestige.mod1 <- lm(prestige ~ education + log(income) + women, data = prestige)
summary(prestige.mod1)
anova(prestige.mod1)

prestige.mod2 <- lm(prestige ~ education + log(income), data = prestige)
summary(prestige.mod2)
anova(prestige.mod2)

anova(prestige.mod1, prestige.mod2)


prestige.mod.inter <- lm(prestige ~ education + log(income) + women + education:women, data = prestige)
summary(prestige.mod.inter)
anova(prestige.mod.inter)

anova(prestige.mod1, prestige.mod.inter)


####
income.mod1 = lm(income ~ education + prestige + women, data = prestige)
income.mod11 = lm(income ~ education + women + prestige, data = prestige)
income.mod111 = lm(income ~ prestige + women + education, data = prestige)
income.mod1111 = lm(income ~ women + prestige + education, data = prestige)

summary(income.mod1)
summary(income.mod11)
summary(income.mod111)
summary(income.mod1111)

anova(income.mod1)
anova(income.mod11)
anova(income.mod111)
anova(income.mod1111)

income.mod2 = lm(income ~ education + women, data = prestige)
summary(income.mod2)
anova(income.mod2)

anova(income.mod1, income.mod2)

datacor = cor(prestige[1:4])
corrplot(datacor, method = "number")
