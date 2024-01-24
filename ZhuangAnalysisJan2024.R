# Load the workspace with data:
setwd("~/Desktop/ZhuangRScript")
load("DZWorkspace.RData")

# Use parallel processing to speed up the GAMM models
require(parallel) 
nc <- 8   ## cluster size, set for example portability
    if (detectCores()>1) { ## no point otherwise
      cl <- makeCluster(nc) 
    } else cl <- NU

# Linear Mixed Model analysis of rhyme duration
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn) # needed for r.squaredGLMM function
library(car) # needed to run Type II ANOVA, rather than type III
library(ggplot2)
library(tictoc)

# Perform the analysis on all data together (nesting structures included)
lm1 <- lmer(Duration ~ checked + ReclassifiedVowelLength + ReclassifiedTone + (1|VowelQuality) + (1|Speaker) + (1|WordId), data = durdata)
summary(lm1)
Anova(lm1)

# Try dropping random effect of vowel quality (very small amount of variance explained)
lm2 <- lmer(Duration ~ checked + ReclassifiedVowelLength + ReclassifiedTone + (1|Speaker) + (1|WordId), data = durdata)
anova(lm2,lm1)
# Retain lm1

# Try an interaction effect between vowel length and tone
lm3 <- lmer(Duration ~ checked + ReclassifiedVowelLength * ReclassifiedTone + (1|VowelQuality) + (1|Speaker) + (1|WordId), data = durdata)
anova(lm1,lm3)
# Interaction not significant so adopt lm1

# R-squared
r.squaredGLMM(lm1)
# Note this returns R2m and R2c, where R2m is a marginal r-squared of the variance explained by the fixed factors only
# R2c is a conditional r-squared for the variance explained by the fixed and random factors (complete model)

# Post-hoc comparisons using emmeans:
syltype <- emmeans(lm1, specs = pairwise ~ checked)
vlength <- emmeans(lm1, specs = pairwise ~ ReclassifiedVowelLength)
tonesdur <- emmeans(lm1, specs = pairwise ~ ReclassifiedTone)

# EMMs for each tone:
tonesemm <- emmeans(lm1, ~ ReclassifiedTone, type="response")
tonesemmdf <- as.data.frame(tonesemm)
# Re-order the levels for ReclassifiedTone for plotting EMMs in Figure 2
tonesemmdf$ReclassifiedTone <- factor(tonesemmdf$ReclassifiedTone, levels(tonesemmdf$ReclassifiedTone)[c(10,9,8,7,6,5,4,3,2,1)])

# ==========================

# Linear mixed model for vowel duration (include Coda here since that affects vowel duration; syllable type is then redundant and is omitted from the model; coda stands in for it)
# Also exclude CV: syllables in durallvnoc to keep the symmetric 2 x 2 design (CVO, CVVO, CVN, CVVN)
# Need to exclude tone here also because coda x vowel length interaction would make it redundant also

# Model 1:
lmvdur1 <- lmer(Duration ~ ReclassifiedCoda * ReclassifiedVowelLength  + (1|Speaker) + (1|WordId) + (1|Repetition) + (1|VowelQuality), data = durallvnocv)
summary(lmvdur1)
anova(lmvdur1)

# Remove the interaction between coda & vowel length since it's insignificant
lmvdur2 <- lmer(Duration ~ ReclassifiedCoda + ReclassifiedVowelLength + (1|Speaker) + (1|WordId) + (1|Repetition) + (1|VowelQuality), data = durallvnocv)
Anova(lmvdur2)
anova(lmvdur1,lmvdur2)

# Add tone into the model now since it's no longer redundant with the interaction effect removed
lmvdur3 <- lmer(Duration ~ ReclassifiedCoda + ReclassifiedVowelLength + ReclassifiedTone + (1|Speaker) + (1|WordId) + (1|Repetition) + (1|VowelQuality), data = durallvnocv)
Anova(lmvdur3)
anova(lmvdur2,lmvdur3)
# Adopt lmvdur3 (lower AIC & significantly different from lmvdur2)

# Remove random effect of repetition
lmvdur4 <- lmer(Duration ~ ReclassifiedCoda + ReclassifiedVowelLength + ReclassifiedTone + (1|Speaker) + (1|WordId) + (1|VowelQuality), data = durallvnocv)
anova(lmvdur4,lmvdur3)
# Adopt lmvdur3

# Try adding Length x tone interaction
lmvdur5 <- update(lmvdur3, ~. + ReclassifiedVowelLength:ReclassifiedTone )
anova(lmvdur3,lmvdur5)
anova(lmvdur5)
# the interaction effect isn't significant in the Type III ANOVA so drop it and retain lmvdur3

r.squaredGLMM(lmvdur3)

# Post-hoc analysis
codavdur <- emmeans(lmvdur3, specs = pairwise ~ ReclassifiedCoda)
vlvdur <- emmeans(lmvdur3, specs = pairwise ~ ReclassifiedVowelLength)
tonevdur <- emmeans(lmvdur3, specs = pairwise ~ ReclassifiedTone)

# EMMs for each tone:
tonesvduremm <- emmeans(lmvdur3, ~ ReclassifiedTone , type="response")
tonesvduremmdf <- as.data.frame(tonesvduremm)
# Re-order the levels for ReclassifiedTone
tonesvduremmdf$ReclassifiedTone <- factor(tonesvduremmdf$ReclassifiedTone, levels(tonesvduremmdf$ReclassifiedTone)[c(8,7,10,9,6,5,4,3,2,1)])

# ==================

# Linear mixed model for Nasal C Duration in CV:N and CVN syllables
# Model 1
lmnas1 <- lmer(Duration ~ ReclassifiedTone*ReclassifiedVowelLength + (1|VowelQuality) + (1|Speaker) + (1|WordId), data = durnonly)

# Try with no interaction
lmnas2 <- lmer(Duration ~ ReclassifiedTone + ReclassifiedVowelLength + (1|VowelQuality) + (1|Speaker) + (1|WordId), data = durnonly)
# Likelihood Ratio Test
anova(lmnas1, lmnas2)
# significant difference and AIC is lower in lmnas1; adopt lmnas1
summary(lmnas1)

# Try without vowel quality:
lmnas3 <- lmer(Duration ~ ReclassifiedTone*ReclassifiedVowelLength + (1|Speaker) + (1|WordId), data = durnonly)
anova(lmnas1,lmnas3)
# Drop random effect of vowel quality - adopt lmnas3

# Type III ANOVA:
anova(lmnas3)
sumlmnas3 <- summary(lmnas3)

r.squaredGLMM(lmnas3)
# R2m is the marginal r-squared (excluding random effects)
# R2c is the conditional r-squared (including random effects)

# for tone:
nldz2 <- emmeans(lmnas3, specs = pairwise ~ ReclassifiedTone|ReclassifiedVowelLength)

# for vowel length:
nldz3 <- emmeans(lmnas3, specs = pairwise ~ ReclassifiedVowelLength|ReclassifiedTone)
# For tone :
nldz3b <- emmeans(lmnas3, specs = pairwise ~ ReclassifiedTone|ReclassifiedVowelLength)
# for vowel length (across tones):
nldz4 <- emmeans(lmnas3, specs = pairwise ~ ReclassifiedVowelLength)

# Post-hoc analysis for interaction effect
nldzuc<- emmeans(lmnas3, ~ReclassifiedTone* ReclassifiedVowelLength)
connl <- contrast(nldzuc, interaction = "pairwise")

# For plot for Fig. 5
nldzuc <- as.data.frame(nldzuc)
# Make a new factor for the interaction of Vowel Length & Tone
nldzuc$ToneVL <- interaction(nldzuc$ReclassifiedTone, nldzuc$ReclassifiedVowelLength, sep = " ")
# Re-order the levels for ToneVL
nldzuc$ToneVL <- factor(nldzuc$ToneVL, levels(nldzuc$ToneVL)[c(12,6,11,5,10,4,9,3,8,2,7,1)])

# ================================================================================================
# ========
# SS-ANOVA analysis of f0 by tone (not used anymore)
library(gss)
# Use X as a sequence that re-maps our results to a plot every 1% of normalized time
X = seq(0,1,by=0.01)

# Create random intercepts for Speaker, WordID:
ranf <- mkran(~1|WordId + 1|Speaker, pitchdatauc)

# SS-ANOVA model (unchecked)
ssafitf0uc <- ssanova(cents ~ ReclassifiedTone * normT, data = pitchdatauc, alpha = 1.4, seed = 1974, random = ranf)
sumfitf0uc <- summary(ssafitf0uc, diagnostics = T)
pif0uc <- sumfitf0uc$pi

# Check reduced models via a Kullback-Leibler projection
f0ucnoint <- project(ssafitf0uc, include = c("normT", "ReclassifiedTone"))
f0uctimeonly <- project(ssafitf0uc, include = "normT")
f0uctoneonly <- project(ssafitf0uc, include = "ReclassifiedTone")

# Display results
print(paste("SS-ANOVA for F0 - unchecked tones in Du'an Zhuang"))
print(sumfitf0uc)
print(paste("Percentage decomposition of explained variance into model terms"))
print(pif0uc)
print(paste("Feasibility of reduced model with no interaction via Kullback-Leibler projection"))
print(f0ucnoint)
print(paste("Feasibility of reduced model with only time via Kullback-Leibler projection"))
print(f0uctimeonly)
print(paste("Feasibility of reduced model with only tone via Kullback-Leibler projection"))
print(f0uctoneonly)

# Set grids for plotting model
gridf0uc <- expand.grid(normT = X, ReclassifiedTone = levels(pitchdatauc$ReclassifiedTone))
gridf0uc$ssafit <- predict(ssafitf0uc, gridf0uc, se = T)$fit
gridf0uc$ssaSE <- predict(ssafitf0uc, gridf0uc, se = T)$se.fit

# ========
# SS-ANOVA Model (checked)

# Create random intercepts for Speaker, Repetition, WordID:
ranfch <- mkran(~1|WordId + 1|as.ordered(Repetition), pitchdatach)

ssafitf0ch <- ssanova(cents ~ ReclassifiedTone * normT, data = pitchdatach, alpha = 1.4, seed = 1974, random = ranfch)
sumfitf0ch <- summary(ssafitf0ch, diagnostics = T)

# Get pi for each effect in the model, the percentage decomposition of explained variance into model terms
pif0ch <- sumfitf0ch$pi
# Check reduced models via a Kullback-Leibler projection
f0chnoint <- project(ssafitf0ch, include = c("normT", "ReclassifiedTone"))
f0chtimeonly <- project(ssafitf0ch, include = "normT")
f0chtoneonly <- project(ssafitf0ch, include = "ReclassifiedTone")

# Display results
print(paste("SS-ANOVA for F0 - checked tones in Du'an Zhuang"))
print(sumfitf0ch)
print(paste("Percentage decomposition of explained variance into model terms"))
print(pif0ch)
print(paste("Feasibility of reduced model with no interaction via Kullback-Leibler projection"))
print(f0chnoint)
print(paste("Feasibility of reduced model with only time via Kullback-Leibler projection"))
print(f0chtimeonly)
print(paste("Feasibility of reduced model with only tone via Kullback-Leibler projection"))
print(f0chtoneonly)

# Set grids for plotting model
gridf0ch <- expand.grid(normT = X, ReclassifiedTone = levels(pitchdatach$ReclassifiedTone))
gridf0ch$ssafit <- predict(ssafitf0ch, gridf0ch, se = T)$fit
gridf0ch$ssaSE <- predict(ssafitf0ch, gridf0ch, se = T)$se.fit

# ========

# GAMM analysis
library(mgcv)
library(itsadug)
    
# Identify first measurement in each trajectory:
pitchdatauc$start.event <- pitchdatauc$relTime == 0

# Create a model without auto-correlation, using the default thin-plate smooth class
# We use method = fREML when running the model from now on to save time
# Method = ML can be used when comparing models for significance but this is resource-intensive

tic()
gamucbase <- bam(cents ~ ReclassifiedTone +
	s(relTime, bs = "ad", k = 30) +
	s(relTime, by = ReclassifiedTone, k = 30) +
	s(Speaker, bs="re")+
	s(VowelQuality, bs="re")+
	s(WordId, bs="re")+
	s(OnsType, bs="re") +
	s(relTime, Speaker, by = ReclassifiedTone, bs = "fs", m = 1, k = 30) +
	s(relTime, VowelQuality, by = ReclassifiedTone, bs = "fs", m = 1, k=30)+
	s(relTime, WordId,  by = ReclassifiedTone, bs = "fs", m = 1, k=30)+
	s(relTime, OnsType,  by = ReclassifiedTone, bs = "fs", m = 1, k=30),
	data = pitchdatauc, method = "fREML", discrete = T, chunk.size = 5000, cluster = cl)
toc()
	
# Get a rough estimate of the correlation between adjacent errors:
ruc <- start_value_rho(gamucbase)

# Model for unchecked tones including residual autocorrelation
tic()
gamuc <- bam(cents ~ ReclassifiedTone +  
s(relTime, bs = "ad",  k=30)+
s(relTime, by= ReclassifiedTone,  k=30)+ 
s(Speaker, bs="re")+ #random intercept
s(VowelQuality, bs="re")+ #random intercept 
s(WordId, bs="re")+ #random intercept
s(OnsType, bs="re")+ #random intercept
s(relTime, Speaker,  by = ReclassifiedTone, bs = "fs", m = 1, k=30)+
s(relTime, VowelQuality, by = ReclassifiedTone, bs = "fs", m = 1, k=30)+
s(relTime, WordId,  by = ReclassifiedTone, bs = "fs", m = 1, k=30)+
s(relTime, OnsType,  by = ReclassifiedTone, bs = "fs", m = 1, k=30),
AR.start = pitchdatauc$start.event,  rho = ruc,
data= pitchdatauc, method="fREML",
discrete = T, chunk.size=5000,cluster=cl)
toc()

# For summary, re.test = F option saves time by removing the random effects. This is done because the summary doesn't finish within a few days even on a very fast processor.
tic()
gamucsummary <- summary(gamuc, re.test = F)
toc()
			
# Check the residual autocorrelation (should be below +/- 0.2 at lag 1 and beyond)
acf_resid(gamuc, split_pred="AR.start")

# =============

# Model for CVO tones only
# Identify first measurement in each trajectory:
pitchdatacvo$start.event <- pitchdatacvo$relTime == 0

# Create a model for CVO tones without auto-correlation and use it get an estimate of the amount of autocorrelation:
tic()
gamcvobase <- bam(cents ~ ReclassifiedTone +  
s(relTime, bs = "ad")+
s(relTime, by= ReclassifiedTone)+ 
s(Speaker, bs="re")+ #random intercept
s(VowelQuality, bs="re")+ #random intercept 
s(WordId, bs="re")+ #random intercept
s(OnsType, bs="re")+ #random intercept
s(relTime, Speaker,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, VowelQuality, by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, WordId,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, OnsType,  by = ReclassifiedTone, bs = "fs", m = 1),
data= pitchdatacvo, method="fREML",
discrete = T, chunk.size=5000,cluster=cl)
toc()
	
# Get a rough estimate of the correlation between adjacent errors:
rcvo <- start_value_rho(gamcvobase)

tic()
gamcvo <- bam(cents ~ ReclassifiedTone +  
s(relTime, bs = "ad")+
s(relTime, by= ReclassifiedTone)+ 
s(Speaker, bs="re")+ #random intercept
s(VowelQuality, bs="re")+ #random intercept 
s(WordId, bs="re")+ #random intercept
s(OnsType, bs="re")+ #random intercept
s(relTime, Speaker,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, VowelQuality, by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, WordId,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, OnsType,  by = ReclassifiedTone, bs = "fs", m = 1),
AR.start = pitchdatacvo$start.event,  rho = rcvo,
data= pitchdatacvo, method="fREML",
discrete = T, chunk.size=5000,cluster=cl)
toc()

gamcvosummary <- summary(gamcvo)

# Check the residual autocorrelation (should be below +/- 0.2 at lag 1 and beyond)
acf_resid(gamcvo, split_pred="AR.start")

# =============

# Model for CVVO tones only
# Identify first measurement in each trajectory:
pitchdatacvvo$start.event <- pitchdatacvvo$relTime == 0

# Create a model for CVVO tones without auto-correlation and use it get an estimate of the amount of autocorrelation:
tic()
gamcvvobase <- bam(cents ~ ReclassifiedTone +  
s(relTime, bs = "ad")+
s(relTime, by= ReclassifiedTone)+ 
s(Speaker, bs="re")+ #random intercept
s(VowelQuality, bs="re")+ #random intercept 
s(WordId, bs="re")+ #random intercept
s(OnsType, bs="re")+ #random intercept
s(relTime, Speaker,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, VowelQuality, by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, WordId,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, OnsType,  by = ReclassifiedTone, bs = "fs", m = 1),
data= pitchdatacvvo, method="fREML",
discrete = T, chunk.size=5000,cluster=cl)
toc()

# Get a rough estimate of the correlation between adjacent errors:
rcvvo <- start_value_rho(gamcvvobase)

# Model for CVVO tones including residual autocorrelation
tic()
gamcvvo <- bam(cents ~ ReclassifiedTone +  
s(relTime, bs = "ad")+
s(relTime, by= ReclassifiedTone)+ 
s(Speaker, bs="re")+ #random intercept
s(VowelQuality, bs="re")+ #random intercept 
s(WordId, bs="re")+ #random intercept
s(OnsType, bs="re")+ #random intercept
s(relTime, Speaker,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, VowelQuality, by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, WordId,  by = ReclassifiedTone, bs = "fs", m = 1)+
s(relTime, OnsType,  by = ReclassifiedTone, bs = "fs", m = 1),
AR.start = pitchdatacvvo$start.event,  rho = rcvvo,
data= pitchdatacvvo, method="fREML",
discrete = T, chunk.size=5000,cluster=cl)
toc()

gamcvvosummary <- summary(gamcvvo)

# Check the residual autocorrelation (should be below +/- 0.2 at lag 1 and beyond)
acf_resid(gamcvvo, split_pred="AR.start")


# ============================================================================

# This section creates plots and is in order of appearance by Figure in the associated JIPA article
library(ggplot2)
library(ggpubr)

# Set working directory to plots folder
setwd("~/Desktop/ZhuangRScript/Plots")

 # Create color scales used for plots
 # For unchecked 6-tone system (all tones actually, but the last 4 get dropped when we're only dealing with the unchecked ones)
 colorScale = c(
'#1f78b4',
'#a6cee3',
'#33a02c',
'#b2df8a',
'#e31a1c',
'#fb9a99',
'#ff7f00',
'#fdbf6f',
'#cab2d6',
'#6a3d9a') # created in http://colorbrewer2.org/

# For checked 4-tone system
colorScaleschecked = c(
'#fdbf6f',
'#ff7f00',
'#cab2d6',
'#6a3d9a')

colorScaleCVO = c(
'#fdbf6f',
'#ff7f00')

colorScaleCVVO = c(
'#cab2d6',
'#6a3d9a')

# Create a vector for black & white versions of the 6 unchecked tones (needed to make a legend in one row)
blackscale = rep(c("black"), times = 6)

# Create a greyscale color vector for B&W plots
greyscale = grey.colors(6, start = .2, end = .7)

# Set the subfolder for the by-speaker plots
plotpathspeakers = "~/Desktop/ZhuangRScript/Plots"

# Size & resolution for png plots
plotwidth = 9
plotheight = 9
plotunits = "cm"
plotres = 300

# Set the font size of the axis titles here
axtisize = 11

# Set the font size for axis text (ticks)
axticksize = 9

# Set the legend title font size here
legendtisize = 9

# Set the legend text font size here
legendtextsize = 9

# Set the font family to use here
fonttype = "Times New Roman"

# Set the thickness of the lines in geom_smooth plots
linethicknss = 0.6

# Set the theme to use in plots for by-speaker Figures here; Set the size for plot titles for speaker id; use hjust 0.5 to center it
themespeakers = theme(plot.title = element_text(size = 11, hjust = 0.5, family = fonttype),
	legend.position = "bottom",
	legend.box = "vertical",
	# Remove axis titles
	axis.title=element_blank(),
	legend.text=element_text(size= legendtextsize, family = fonttype),
	axis.text = element_text(size= axticksize, family = fonttype),
	legend.title=element_text(size= legendtisize, family = fonttype),
	legend.key.size = unit(0.5, "cm")
)

# Set the theme for the duration plots across speakers (no axis labels)
themedur = theme(plot.title = element_text(size = axtisize, hjust = 0.5, family = fonttype),
	legend.position = "bottom",
	# Remove axis titles
	axis.title=element_blank(),
	legend.text=element_text(size= legendtextsize, family = fonttype),
	axis.text = element_text(size= axticksize, family = fonttype),
	legend.title=element_text(size= legendtisize, family = fonttype))
	
# Set the theme for the duration plots across speakers (with axis labels)
themedurbasic = theme(plot.title = element_text(size = 14, hjust = 0.5, family = fonttype),
	axis.title = element_text(size = 12, family = fonttype),
	legend.title=element_text(size= 12, family = fonttype),
	axis.text = element_text(size= 12, family = fonttype)
)
	
# ===============================================================================================

# Figure 2. EMMs of rhyme duration (ms) by tone with 95% confidence intervals 
rdurbytone <- ggplot(tonesemmdf, aes(x = ReclassifiedTone, y = emmean, ymin = lower.CL, ymax = upper.CL))+
	geom_col(position = position_dodge(), color = 'black', alpha = .5)+
	geom_point(position = position_dodge(width = 0.9))+
	coord_flip()+
	geom_errorbar(position=position_dodge(width=0.9), width = 0.20)+
	labs(x = "Tone",
		y = "Rhyme duration (ms)")+
	themedurbasic
ggsave(filename = "Fig2.png", plot = rdurbytone, dpi = plotres)

# ===============================================================================================

# Figure 3. EMMs of vowel duration (ms) by tone with 95% confidence intervals 
vdurbytone <- ggplot(tonesvduremmdf, aes(x = ReclassifiedTone, y = emmean, ymin = lower.CL, ymax = upper.CL))+
	geom_col(position = position_dodge(), color = 'black', alpha = .5)+
	geom_point(position = position_dodge(width = 0.9))+
	coord_flip()+
	geom_errorbar(position=position_dodge(width=0.9), width = 0.20)+
	labs(x = "Tone",
		y = "Vowel duration (ms)")+
	themedurbasic
ggsave(filename = "Fig3.png", plot = vdurbytone, dpi = plotres)

# ===============================================================================================

# Figure 4 - Rhyme Duration split into vowel and nasal duration (ms) by syllable type, one speaker per panel 
# Plot duration by type for DZ speakers
# stacked barplot for vowel & coda duration shown together
for (i in 1:length(vdurspeakers)){
	spkr = vdurspeakers[i]

	# Get the subset for that speaker only
	spkrsub <- droplevels(subset(durall, durall$Speaker == spkr))
	
	# Plot vowel duration (you need to sort the data frame in reverse order to get vowel on the left in the plot...)
	# stat = "summary" is crucial to get mean Duration; otherwise, it counts data-points, giving a useless histogram
	assign(paste("vdurbytypegrey", spkr, sep = ""), ggplot(spkrsub, aes(fill = Label, y = Type, x = Duration)) +
	geom_bar(position = position_stack(reverse = T), stat = "summary") +
	theme_bw() +
	themespeakers +
	scale_fill_grey(start = 0.3, end = .6)+
	scale_colour_grey(start = 0.3, end = .6)+
	labs(fill = "Segment Type",
		title = spkr,
		x = "Duration (ms)",
		y = "Syllable Type"))
	}

# Plot Vowel Duration contrast by Syllable Type greyscale plot for DZ
figzgrey <- ggarrange(vdurbytypegreySpkr1, vdurbytypegreySpkr2, vdurbytypegreySpkr3, vdurbytypegreySpkr4, vdurbytypegreySpkr5, vdurbytypegreySpkr6, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
figzgrey <- annotate_figure(figzgrey, left = text_grob("Syllable Type", size = axtisize, rot = 90, family = fonttype), bottom = text_grob("Duration (ms)", size = axtisize, family = fonttype))
ggsave(filename = "VdurBySyllTypeHorizontal.png", figzgrey, path = plotpathspeakers, width = 2*plotwidth, height = 0.75*2*plotheight, units = plotunits, dpi = plotres)

# ===============================================================================================

# Figure 5. EMMs of nasal coda duration (ms) by tone and vowel length with 95% confidence intervals 
ndurbytone <- ggplot(nldzuc, aes(x = ToneVL, y = emmean, ymin = lower.CL, ymax = upper.CL))+
	geom_col(position = position_dodge(), color = 'black', alpha = .5)+
	geom_point(position = position_dodge(width = 0.9))+
	coord_flip()+
	geom_errorbar(position=position_dodge(width=0.9), width = 0.20)+
	labs(x = "Tone",
		y = "Nasal duration (ms)")+
	themedurbasic
ggsave(filename = "Fig5.png", plot = ndurbytone, dpi = plotres)

# =======================================================

# Figure 6: Smooth scatterplots for f0 for each tone plotted against normalized time, one panel per speaker. Note that the f0 ranges vary per speaker.
for (i in 1:length(dzspeakers)){
	spkr = dzspeakers[i]

	# Get subset for each speaker
	sppitchdatauc <- droplevels(subset(pitchdatauc, pitchdatauc$Speaker == spkr))
	sppitchdatach <- droplevels(subset(pitchdatach, pitchdatach$Speaker == spkr))
	# Color plot
	assign(paste("f0uc", spkr, "color", sep = ""), ggplot(sppitchdatauc, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents), linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	#coord_cartesian(ylim = c(f0min, f0max)) +
	scale_colour_manual(values=colorScale)+
	scale_fill_manual(values=colorScale)+
	guides(color = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		color = "Tone",
		# Use bquote() instead of expression() if you want to call values of constants/variables inside .()
		title = bquote(.(spkr))))
		
	# Greyscale plot 
	assign(paste("f0uc", spkr, "grey", sep = ""), ggplot(sppitchdatauc, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents, linetype = ReclassifiedTone), linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	#coord_cartesian(ylim = c(f0min, f0max)) +
	scale_colour_manual(values=blackscale)+
	scale_fill_manual(values=blackscale)+
	guides(color = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		color = "Tone",
		linetype = "Tone",
		title = bquote(.(spkr))))
	
	# Checked tones	

	# Color plot
	assign(paste("f0ch", spkr, "color", sep = ""), ggplot(sppitchdatach, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents), linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	#coord_cartesian(ylim = c(f0min, f0max)) +
	scale_colour_manual(values=colorScaleschecked)+
	scale_fill_manual(values=colorScaleschecked)+
	guides(color = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		color = "Tone",
		title = bquote(.(spkr))))
		
	# Greyscale plot 
	assign(paste("f0ch", spkr, "grey", sep = ""), ggplot(sppitchdatach, aes(x = normT, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents, linetype = ReclassifiedTone), color = "black", linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	#coord_cartesian(ylim = c(f0min, f0max)) +
	labs(fill = "Tone",
		color = "Tone",
		linetype = "Tone",
		title = bquote(.(spkr))))
}

# Color Plots for unchecked tones
f0DZuccolor <- ggarrange(f0ucSpkr1color, f0ucSpkr2color, f0ucSpkr3color, f0ucSpkr4color, f0ucSpkr5color, f0ucSpkr6color, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
f0DZuccolor <- annotate_figure(f0DZuccolor,
	top = text_grob("Unchecked Tones", size = axtisize, family = fonttype), fig.lab.face = c("bold"))
# Color plots for checked tones
f0DZchcolor <- ggarrange(f0chSpkr1color, f0chSpkr2color, f0chSpkr3color, f0chSpkr4color, f0chSpkr5color, f0chSpkr6color, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
f0DZchcolor <- annotate_figure(f0DZchcolor,
	 top = text_grob("Checked Tones", size = axtisize, family = fonttype), fig.lab.face = c("bold"))

# Join the unchecked and checked plots
f0DZcolor <- ggarrange(f0DZuccolor, f0DZchcolor, ncol = 1, nrow = 2)
f0DZcolor <- annotate_figure(f0DZcolor, left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype))
ggsave(filename = "f0DZ_color.png", f0DZcolor, path = plotpathspeakers, width = 1.5*plotwidth, height = 1.8*plotheight, units = plotunits, dpi = plotres)

# Greyscale Plots for unchecked tones
f0DZucgrey <- ggarrange(f0ucSpkr1grey, f0ucSpkr2grey, f0ucSpkr3grey, f0ucSpkr4grey, f0ucSpkr5grey, f0ucSpkr6grey, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
f0DZucgrey <- annotate_figure(f0DZucgrey,
	top = text_grob("Unchecked Tones", size = axtisize, family = fonttype), fig.lab.face = c("bold"))
# Greyscale plots for checked tones
f0DZchgrey <- ggarrange(f0chSpkr1grey, f0chSpkr2grey, f0chSpkr3grey, f0chSpkr4grey, f0chSpkr5grey, f0chSpkr6grey, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
f0DZchgrey <- annotate_figure(f0DZchgrey,
	top = text_grob("Checked Tones", size = axtisize, family = fonttype), fig.lab.face = c("bold"))

# Join the unchecked and checked plots
f0DZgrey <- ggarrange(f0DZucgrey, f0DZchgrey, ncol = 1, nrow = 2)
f0DZgrey <- annotate_figure(f0DZgrey, left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype))
ggsave(filename = "f0DZ_grey.png", f0DZgrey, path = plotpathspeakers, width = 1.5*plotwidth, height = 1.8*plotheight, units = plotunits, dpi = plotres)

# ===========================================================================================

# Plot for SS-ANOVA: F0 smooth splies for each tone in Du'an Zhuang plotted against normalized time with 95% credible intervals. Significant differences between contours correspond to non-overlapping periods

# Set the theme to use in ss-anova plots across speakers here
themessanova = theme(plot.title = element_text(size = 11, hjust = 0.5, family = fonttype),
	legend.position = "bottom",
	# Remove axis titles
	axis.title=element_blank(),
	legend.text=element_text(size= legendtextsize, family = fonttype),
	axis.text = element_text(size= axticksize, family = fonttype),
	legend.title=element_text(size= legendtisize, family = fonttype))

# Set min & max f0
dzf0minssa = -600
dzf0maxssa = 600

# f0 color plot - unchecked tones
ssaplotdzucf0color <- ggplot(gridf0uc, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_line(aes(y = ssafit), color = "black")+
	geom_ribbon(aes(ymin = ssafit - (1.96*ssaSE), ymax = ssafit + (1.96*ssaSE), fill = ReclassifiedTone), alpha = 0.5) +
	theme_bw()+
	themessanova +
	coord_cartesian(ylim = c(dzf0minssa, dzf0maxssa))+
	scale_colour_manual(values=colorScale)+
	scale_fill_manual(values=colorScale) +
	guides(colour = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		color = "Tone",
		title = "Unchecked Tones",
		x = "Normalized Time",
		y = "F0 (Cents re. median F0/Speaker)")
		
# f0 greyscale plot - unchecked tones
ssaplotdzucf0grey <- ggplot(gridf0uc, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_line(aes(y = ssafit, linetype = ReclassifiedTone), color = "black")+
	geom_ribbon(aes(ymin = ssafit - (1.96*ssaSE), ymax = ssafit + (1.96*ssaSE), fill = ReclassifiedTone), alpha = 0.5) +
	theme_bw()+
	themessanova +
	coord_cartesian(ylim = c(dzf0minssa, dzf0maxssa))+
	scale_fill_grey(start = 0.2, end = .7)+
	scale_colour_grey(start = 0.2, end = .7)+
	guides(colour = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		linetype = "Tone",
		color = "Tone",
		title = "Unchecked Tones",
		x = "Normalized Time",
		y = "F0 (Cents re. median F0/Speaker)")

# f0 color plot - checked tones
ssaplotdzcf0color <- ggplot(gridf0ch, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_line(aes(y = ssafit), color = "black")+
	geom_ribbon(aes(ymin = ssafit - (1.96*ssaSE), ymax = ssafit + (1.96*ssaSE), fill = ReclassifiedTone), alpha = 0.5) +
	theme_bw()+
	themessanova +
	coord_cartesian(ylim = c(dzf0minssa, dzf0maxssa))+
	scale_colour_manual(values=colorScaleschecked)+
	scale_fill_manual(values=colorScaleschecked) +
	guides(colour = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		color = "Tone",
		title = "Checked Tones",
		x = "Normalized Time",
		y = "F0 (Cents re. median F0/Speaker)")
		
# f0 greyscale plot - checked tones
ssaplotdzcf0grey <- ggplot(gridf0ch, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_line(aes(y = ssafit, linetype = ReclassifiedTone), color = "black")+
		geom_ribbon(aes(ymin = ssafit - (1.96*ssaSE), ymax = ssafit + (1.96*ssaSE), fill = ReclassifiedTone), alpha = 0.5) +
	theme_bw()+
	themessanova +
	coord_cartesian(ylim = c(dzf0minssa, dzf0maxssa))+
	scale_fill_grey(start = 0.2, end = .7)+
	scale_colour_grey(start = 0.2, end = .7)+
	guides(colour = guide_legend(nrow = 1)) +
	labs(fill = "Tone",
		linetype = "Tone",
		color = "Tone",
		title = "Checked Tones",
		x = "Normalized Time",
		y = "F0 (Cents re. median F0/Speaker)")
		
# Join the unchecked and checked plots
# Color
ssapf0DZcolor <- ggarrange(ssaplotdzucf0color, ssaplotdzcf0color, ncol = 2, nrow = 1)
ssapf0DZcolor <- annotate_figure(ssapf0DZcolor,
	left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype))
ggsave(filename = "ssanovaf0DZ_color_v3.png", ssapf0DZcolor, path = plotpathspeakers, width = 2*plotwidth, height = plotheight, units = plotunits, dpi = plotres)

# Greyscale
ssapf0DZgrey <- ggarrange(ssaplotdzucf0grey, ssaplotdzcf0grey, ncol = 2, nrow = 1)
ssapf0DZgrey <- annotate_figure(ssapf0DZgrey,
	left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype))
ggsave(filename = "ssanovaf0DZ_grey_v3.png", ssapf0DZgrey, path = plotpathspeakers, width = 2*plotwidth, height = plotheight, units = plotunits, dpi = plotres)

# =====================================================================================

# Figure 7 left side: Plot model for unchecked syllables with voiced stop onsets
png("GAMMUC_Jan19.png", width = 3.25, height = 3.25, units = "in", res = 300, pointsize = 8)
plot_smooth(gamuc, view = "relTime", plot_all="ReclassifiedTone", rug = F,  rm.ranef = T, col = colorScale, xlab = "Time (ms)", ylab = "F0 (cents)")
dev.off()

# Figure 7: model for CVO syllables: 
png("GAMMCVO_Jan19.png", width = 3.25, height = 3.25, units = "in", res = 300, pointsize = 8)
plot_smooth(gamcvo, view = "relTime", plot_all="ReclassifiedTone", rug = F,  rm.ranef = T, col = colorScaleCVO, xlab = "Time (ms)", ylab = "F0 (cents)")
dev.off()

# Figure 7: model for CVVO syllables:
png("GAMMCVVO_Jan19.png", width = 3.25, height = 3.25, units = "in", res = 300, pointsize = 8)
plot_smooth(gamcvvo, view = "relTime", plot_all="ReclassifiedTone", rug = F,  rm.ranef = T, col = colorScaleCVVO, xlab = "Time (ms)", ylab = "F0 (cents)")
dev.off()

# =====================================================================================

# Figures 8 to 12 - f0 contours for selections of tones, one speaker per panel

# Create the color scale matching the old shades used for each selection of tones
# Tones 5,6,9,10
colorscale56910 = c(
'#e31a1c',
'#fb9a99',
'#cab2d6',
'#6a3d9a'
)

# Tones 1,3,7
colorscale137 = c(
'#1f78b4',
'#33a02c',
'#ff7f00')

# Tone 2, 4, 6, 8
colorscale8 = c(
'#a6cee3',
'#b2df8a',
'#fb9a99',
'#fdbf6f',
'#6a3d9a'
)

# Create customized subsets of pitchdatadzdic for each set of tones we want to view together:

# Tones 5,6,9,10
dz56910 <- droplevels(subset(pitchdatadzdic, pitchdatadzdic$ReclassifiedTone == "5" | pitchdatadzdic$ReclassifiedTone == "6" | pitchdatadzdic$ReclassifiedTone == "9"| pitchdatadzdic$ReclassifiedTone == "10"))

# Tones 1,3,7
dz137 <- droplevels(subset(pitchdatadzdic, pitchdatadzdic$ReclassifiedTone == "1" | pitchdatadzdic$ReclassifiedTone == "3" | pitchdatadzdic$ReclassifiedTone == "7"))

dz246810 <- droplevels(subset(pitchdatadzdic, pitchdatadzdic$ReclassifiedTone == "2" | pitchdatadzdic$ReclassifiedTone == "4" | pitchdatadzdic$ReclassifiedTone == "6" | pitchdatadzdic$ReclassifiedTone == "8" | pitchdatadzdic$ReclassifiedTone == "10"))

# Set x-axis limit for time
xmax = 400

# Set y-axis limits for f0
ymin = -500
ymax = 250

# Choose which set of tones to display below and the appropriate color scale
tns <- dz246810
colsc <- colorscale8
pngname <- "DZAllotones8"

for (i in 1:length(dzspeakers)){
	spkr = dzspeakers[i]

	# Get subsets for each speaker
	sp <- droplevels(subset(tns, tns$Speaker == spkr))
			
	# F0 greyscale plot, absolute time
	assign(paste(spkr, "grey", sep = ""), ggplot(sp, aes(x = relTime, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents, linetype = ReclassifiedTone, fill = ReclassifiedTone), color = "black", linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	coord_cartesian(xlim = c(0, xmax)
	, ylim = c(ymin, ymax)
	) +
	scale_fill_grey(start = 0.2, end = .7)+
	scale_colour_grey(start = 0.2, end = .7)+
	labs(fill = "Tone",
		color = "Tone",
		linetype = "Tone",
		title = bquote(.(spkr)),
		x = "Time (ms)",
		y = "F0 (Cents re. median F0/Speaker)"))
		
	# F0 color plot, absolute time
	assign(paste(spkr, "color", sep = ""), ggplot(sp, aes(x = relTime, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents), linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	coord_cartesian(xlim = c(0, xmax)
	, ylim = c(ymin, ymax)
	)+
	scale_colour_manual(values = colsc)+
	scale_fill_manual(values = colsc)+
	labs(fill = "Tone",
		color = "Tone",
		title = bquote(.(spkr)),
		x = "Time (ms)",
		y = "F0 (Cents re. median F0/Speaker)"))
		
	# F0 greyscale plot, normalized time
	assign(paste(spkr, "greynorm", sep = ""), ggplot(sp, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents, linetype = ReclassifiedTone, fill = ReclassifiedTone), color = "black", linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	coord_cartesian(xlim = c(0, 1), ylim = c(ymin, ymax)) +
	scale_fill_grey(start = 0.2, end = .7)+
	scale_colour_grey(start = 0.2, end = .7)+
	labs(fill = "Tone",
		color = "Tone",
		linetype = "Tone",
		title = bquote(.(spkr)),
		x = "Normalized Time",
		y = "F0 (Cents re. median F0/Speaker)"))
		
	# F0 color plot, absolute time
	assign(paste(spkr, "colornorm", sep = ""), ggplot(sp, aes(x = normT, color = ReclassifiedTone, group = ReclassifiedTone))+
	geom_smooth(aes(y = cents), linewidth = linethicknss, method = "gam")+
	theme_bw()+
	themespeakers +
	coord_cartesian(xlim = c(0, 1), ylim = c(ymin, ymax)) +
	scale_colour_manual(values = colsc)+
	scale_fill_manual(values = colsc)+
	labs(fill = "Tone",
		color = "Tone",
		title = bquote(.(spkr)),
		x = "Normalized Time",
		y = "F0 (Cents re. median F0/Speaker)"))

}

# Color Plot for F0, absolute time
pcolor <- ggarrange(Spkr1color, Spkr2color, Spkr3color, Spkr4color, Spkr5color, Spkr6color, ncol = 3, nrow = 2, common.legend = T, legend = "right")
pcolor <- annotate_figure(pcolor,
	left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype), bottom = text_grob("Time (ms)", size = axtisize, family = fonttype))
ggsave(filename = paste(pngname,"Color.png"), pcolor, path = plotpathspeakers, width = 1.5*plotwidth, height = plotheight, units = plotunits, dpi = plotres)

# Greyscale Plot for F0, absolute time
pgrey <- ggarrange(Spkr1grey, Spkr2grey, Spkr3grey, Spkr4grey, Spkr5grey, Spkr6grey, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
pgrey <- annotate_figure(pgrey,
	left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype), bottom = text_grob("Time (ms)", size = axtisize, family = fonttype))
ggsave(filename = paste(pngname,"Grey.png"), pgrey, path = plotpathspeakers, width = 1.5*plotwidth, height = plotheight, units = plotunits, dpi = plotres)

# Color Plot for F0, normalized time
pcolornorm <- ggarrange(Spkr1colornorm, Spkr2colornorm, Spkr3colornorm, Spkr4colornorm, Spkr5colornorm, Spkr6colornorm, ncol = 3, nrow = 2, common.legend = T, legend = "right")
pcolornorm <- annotate_figure(pcolornorm,
	left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype), bottom = text_grob("Normalized Time", size = axtisize, family = fonttype))
ggsave(filename = paste(pngname,"Color_NormT.png"), pcolornorm, path = plotpathspeakers, width = 1.5*plotwidth, height = plotheight, units = plotunits, dpi = plotres)

# Greyscale Plot for F0, absolute time
pgreynorm <- ggarrange(Spkr1greynorm, Spkr2greynorm, Spkr3greynorm, Spkr4greynorm, Spkr5greynorm, Spkr6greynorm, ncol = 3, nrow = 2, common.legend = T, legend = "bottom")
pgreynorm <- annotate_figure(pgreynorm,
	left = text_grob("F0 (Cents re. median F0/Speaker)", size = axtisize, rot = 90, family = fonttype), bottom = text_grob("Normalized Time", size = axtisize, family = fonttype))
ggsave(filename = paste(pngname,"Grey_NormT.png"), pgreynorm, path = plotpathspeakers, width = 1.5*plotwidth, height = plotheight, units = plotunits, dpi = plotres)