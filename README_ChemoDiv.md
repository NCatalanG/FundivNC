Initially I will  add some functions on alpha and beta diversity based on vegan
======================================
```{r}
require (vegan)
require (BiodiversityR)
```
Vegan-based shannon
```{r}
data (BCI)#tree counts in Colorado Island
H <- diversity(BCI, index="shannon") #Shannon index
simp <- diversity(BCI, index="simpson") #Simpson index
isimp <- diversity(BCI, index= "inv") #inverse Simpson
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(BCI, 2) - 1

## Fisher alpha
alpha <- fisher.alpha(BCI)
## Plot all
pairs(cbind(H, simp, isimp, unbias.simp, alpha), pch="+", col="blue")
```
Richness and rarefaction:
Species richness increases with sample size, and differences in richness actually may be caused by differences in sample size.To solve this problem, we may try to rarefy species richness to the same number of individuals.
```{r}
quantile(rowSums(BCI))
S <- specnumber(BCI) #number of species
fS <- specnumber(BCI, margin=2) #frequency of species 

Srar <- rarefy(BCI, min(rowSums(BCI)))#richness for the same number of individuals
plot(S,Srar, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
(raremax <- min(rowSums(BCI)))
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

```

Species abundance and Eveness indices:
```{r}
J <- H/log(specnumber(BCI))#Pielou's Eveness (included in Eve)
eve1 <- Eve(A = BCI) #from fundiv. Ea= 
#eve2 <- Eve(A = dummy$abun, scales = c(0.25,0.5,2,4,8,Inf))
pairs(eve1)
#pairs(eve2)

#for a more specific look into Renyi diversities:
k <- sample(nrow(BCI), 6)
R <- renyi(BCI[k,])
plot(R)

#expected number of species with n individuals (Fisher)
k <- sample(nrow(BCI), 1)
fish <- fisherfit(BCI[k,])
fish
plot(prestondistr(BCI[k,]))#alternative to the lognormal model

#ranked abundances:
rad <- radfit(BCI[k,])
rad

```
Beta diversity and accumulation
Species accumulation models are similar to rarefaction: they study the accumulation of species when the number of sites increases.
Beta diversity defined as gamma/alpha - 1:
```{r}
beta<-ncol(BCI)/mean(specnumber(BCI)) - 1

#Sorensen index(distance)
beta <- vegdist(BCI, binary=TRUE)
mean(beta)
#other indices:
betadiver(help=TRUE)
#Function betadisper can be used to analyse beta diversities with respect to classes or factors
data(dune)
data(dune.env)
z <- betadiver(dune, "z")
mod <- with(dune.env, betadisper(z, Management))
mod

#accumulation
sac <- specaccum(BCI)
plot(sac, ci.type="polygon", ci.col="yellow")

```

Total species pool (unknown  number of species  in the community):
```{r}
specpool(BCI)#with a collection of sites

estimateR(BCI[k,])#number of unseen species for each single site.These
functions need counts of individuals, and species seen only once or twice, or other rare species, take the place of species with low frequencies.
veiledspec(BCI[k,])

#which unseen species might be where?it has been suggested as a method of estimating which of the missing species could occur in a site, or which sites are suitable for a species. The probability for each species at each site is assessed from other species occurring on the site.
smo <- beals(BCI)
j <- which(colnames(BCI) == "Ceiba.pentandra")
plot(beals(BCI, species=j, include=FALSE), BCI[,j],
ylab="Occurrence", main="Ceiba pentandra",
xlab="Probability of occurrence")
```
Taxonomic diversity using compound categories as e.g. genus.Maybe that is conceptually more correct?
````{r}
data(dune)
data(dune.taxon)
taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
mod <- taxondive(dune, taxdis)
plot(mod)
```
fundiv: analyzing functional trait diversity
========================================================
This we maybe can apply with compounds that we know whether degrade or no, whether are photoreactive or not, etc.
````{r}
tr <- hclust(taxdis, "aver")
mod <- treedive(dune, tr)#implements functional diversity defined as the total branch length in a trait dendrogram connecting all species, but excluding the unnecessary root segments of the tree

```

This is a wrapper to Petchey and Gaston FD indexes and Laliberte FD package. Hence most code is borrowed from them. It has some additions, most notably the dendogram measures can also be weighted by abundance as implemented in Gagic, Bartomeus et al. (2015 Proc B http://rspb.royalsocietypublishing.org/content/282/1801/20142620 ). There is also a function to calculate several Evennes indexes.

This version is functional, but in development, so please report bugs, etc...

To install the package run:

```{r}
install.packages("devtools")
require(devtools)
install_github("ibartomeus/fundiv")
require(fundiv)
```

To calculate dendogram indexes:

```{r}
FD_dendro <- FD_dendro(S = dummy$trait, A = dummy$abun, Cluster.method = "average", ord = "podani",
                    Weigthedby = "abundance")
FD_dendro
```

The dendrogram is ploted and returns a dataframe with the indexes.


You may also want indexes based in trait space in FD package. You can run 'dbFD {FD}',
or the wrapper below which also includes other diversity metrics.

```{r}
FD_all <- FDindexes(S = dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
                  Cluster.method= "average", corr= "cailliez", Weigthedby = "abundance")
FD_all
```

A nice addition if that you may want the indexes weigthed by biomass, not abundance. For bees and carabids, a function to convert body length to body mass is provided (see ?length_to_mass)

```{r}
FD_all_bm <- FDindexes(S = dummy$trait, A = dummy$abun, Distance.method= "gower", ord= "podani", 
                  Cluster.method= "average", corr= "cailliez", Weigthedby = "biomassValue",
                  biomassValue = c(1.2, 2.3, 0.6, 1.0, 3.4, 0.2, 1.6, 2.2))

FD_all_bm
```

We can see that both families of indexes are correlated

```plot(FD_all$FDpg ~ FD_all$Frich)```

Note that Functional richness indexes are highly correlated with richness, hence you may also want to know if those indexes are higher or lower than expected by its richness levels (See Rader et al. 2014 D&D for details http://onlinelibrary.wiley.com/doi/10.1111/ddi.12221/abstract)

```{r}
null <- null.FD(S = dummy$trait, A = dummy$abun, it = 100, w = NA)
null
```

No much significant results in this dataset, but is not surprising given:

```plot(FD_all$FDpg ~ FD_all$n_sp)```

You can calculate standardized FD indexes as in Rader et al. 2014 by 

```{r}
(null$FD - null$null_meanFD) / null$null_sdFD
```

Indexes implemented in Clark et al. 2012 (Plos One) are also available (`FD_Clark`), but in my opinion they perform very similar that the proposed FDw, but are computational time consuming because they need to 
create a dendogram for each community, instead of using only a general dendogram containing all species. This last approach is better according to the literature...

Lastly I added the calculation of some evenness indexes. See ?Eve for details.

```{r}
eve1 <- Eve(A = dummy$abun)
eve2 <- Eve(A = dummy$abun, scales = c(0.25,0.5,2,4,8,Inf))
pairs(eve1)
pairs(eve2)
```


