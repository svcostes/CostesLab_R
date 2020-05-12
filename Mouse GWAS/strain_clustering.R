# this is for clustering mouse strains based on FPG phenotypes 

#load packages ---------

library(tidyverse)
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # dendrogram visualization

#load data
setwd(getwd())
fpg_data <- read.csv("FPG_data.csv", header=TRUE, row.names=1)

fpg_scaled_data <- scale(fpg_data)    # centers and scales data

d <- dist(fpg_data, method = "euclidean")   #raw data
ds <- dist(fpg_scaled_data, method = "euclidean") #scaled data

# agglomerative hierarchical clustering  ------

# using agnes, compare various linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute agglomerative coefficients
ac <- function(x) {
  agnes(d, method = x)$ac
}

map_dbl(m, ac)  #display agglomerative coefficients

# agglomerative clustering using ward's method
hc_a_w <- agnes(d, method = "ward")     #raw data
hc_a_w <- as.dendrogram(hc_a_w)   # turn it into a dendrogram
colored_hc_a_w <- color_labels(hc_a_w, k=3)    # color the labels based on clusters
plot(colored_hc_a_w, main = "Agglomerative, ward, raw data")

grp_hc_a_w <- cutree(hc_a_w, k = 3)  # assign cluster group

s_hc_a_w <- agnes(ds, method = "ward")      #standardized data
s_hc_a_w <- as.dendrogram(s_hc_a_w)     # turn it into a dendrogram
colored_s_hc_a_w <- color_labels(s_hc_a_w, k=3)   # color the labels based on clusters
plot(colored_s_hc_a_w, main = "Agglomerative, ward, scaled data")

grp_s_hc_a_w <- cutree(s_hc_a_w, k = 3)    # assign cluster group

# agglomerative clustering using complete linkage 
hc_a_c <- agnes(d, method = "complete")     #raw data
hc_a_c <- as.dendrogram(hc_a_c)  # turn it into a dendrogram
colored_hc_a_c <- color_labels(hc_a_c, k=3)   # color the labels based on clusters
plot(colored_hc_a_c, main = "Agglomerative, complete, raw data") 

grp_hc_a_c <- cutree(hc_a_c, k = 3)  # assign cluster group

s_hc_a_c <- agnes(ds, method = "complete")      #standardized data
s_hc_a_c <- as.dendrogram(s_hc_a_c)  #turn it into a dendrogram
colored_s_hc_a_c <- color_labels(s_hc_a_c, k=3)  # color the labels based on clusters
plot(colored_s_hc_a_c, main = "Agglomerative, complete, scaled data") 

grp_s_hc_a_c <- cutree(s_hc_a_c, k = 3)  # assign cluster group

# agglomerative clustering using average linkage 
hc_a_a <- agnes(d, method = "average")     #raw data
hc_a_a <- as.dendrogram(hc_a_a)
colored_hc_a_a <- color_labels(hc_a_a, k=3)
plot(colored_hc_a_a, main = "Agglomerative, average, raw data") 

grp_hc_a_a <- cutree(hc_a_a, k = 3)  # assign cluster group

s_hc_a_a <- agnes(ds, method = "average")      #standardized data
s_hc_a_a <- as.dendrogram(s_hc_a_a)
colored_s_hc_a_a <- color_labels(s_hc_a_a, k=3)
plot(colored_s_hc_a_a, main = "Agglomerative, average, scaled data") 

grp_s_hc_a_a <- cutree(s_hc_a_a, k = 3)  # assign cluster group

# agglomerative clustering using single linkage 
hc_a_s <- agnes(d, method = "single")     #raw data
hc_a_s <- as.dendrogram(hc_a_s)
colored_hc_a_s <- color_labels(hc_a_s, k=3)
plot(colored_hc_a_s, main = "Agglomerative, single, raw data") 

grp_hc_a_s <- cutree(hc_a_s, k = 3)  # assign cluster group

s_hc_a_s <- agnes(ds, method = "single")      #standardized data
s_hc_a_s <- as.dendrogram(s_hc_a_s)
colored_s_hc_a_s <- color_labels(s_hc_a_s, k=3)
plot(colored_s_hc_a_s, main = "Agglomerative, single, scaled data") 

grp_s_hc_a_s <- cutree(s_hc_a_s, k = 3)  # assign cluster group

# Divisive hierarchical clustering -------------
# can't change linkage method, linkage is based on max avg dissimilarity, which is essentially "average"

hc_d <- diana(d)  # raw data
s_hc_d <- diana(ds)   #standardized data

# Divisive coefficient
hc_d$dc  # raw data
s_hc_d$dc  # standardized data

hc_d <- as.dendrogram(hc_d)
s_hc_d <- as.dendrogram(s_hc_d)

colored_hc_d <- color_labels(hc_d, k=3)
colored_s_hc_d <- color_labels(s_hc_d, k=3)

plot(colored_hc_d, main = "Divisive, raw data")  
plot(colored_s_hc_d, main = "Divisive, scaled data")  

grp_hc_d <- cutree(hc_d, k = 3)  # assign cluster group
grp_s_hc_d <- cutree(s_hc_d, k = 3)  # assign cluster group


# final table with all cluster group assignments
table1 <- mutate(fpg_data, A_W_R = grp_hc_a_w, A_W_S = grp_s_hc_a_w, 
                 A_C_R = grp_hc_a_c, A_C_S = grp_s_hc_a_c, 
                 A_A_R = grp_hc_a_a, A_A_S = grp_s_hc_a_a,
                 A_S_R = grp_hc_a_s, A_S_S = grp_s_hc_a_s,
                 D_R = grp_hc_d, D_S = grp_s_hc_d)

# comparing dendrograms --------
plot(dendlist(hc_a_w, s_hc_a_w))  
plot(dendlist(hc_a_c, s_hc_a_c))
plot(dendlist(hc_a_a, s_hc_a_a))
plot(dendlist(hc_a_s, s_hc_a_s))
plot(dendlist(hc_d, s_hc_d))
plot(dendlist(s_hc_a_w, s_hc_d))
plot(dendlist(s_hc_a_w, s_hc_a_c))
plot(dendlist(s_hc_a_w, s_hc_a_a))
plot(dendlist(s_hc_a_c, s_hc_a_a))
plot(dendlist(s_hc_a_a, s_hc_d))
