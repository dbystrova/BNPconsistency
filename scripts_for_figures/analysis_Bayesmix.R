

a = load("~/bayesmix/simdata/best_clustering.csv")
data = date_500$y

hist(data[,1])
best_clust <-read.table("~/bayesmix/simdata/best_clustering.csv", header = TRUE, sep = ",") 
data_in_clust_df = data.frame()
c = unique(best_clust)[,1]
for (i in 1:length(c)){
  data_in_clust =data[best_clust ==c[i],]
  data_in_clust_ = as.data.frame(data_in_clust)
  data_in_clust_$color = c[i]
  data_in_clust_df = rbind(data_in_clust_df,data_in_clust_)
}
names(data_in_clust_df) = c("x", "y", "color")
ggplot(data_in_clust_df, aes(x=x, y=y, color=color)) +
  geom_point()
data_in_clust =data[best_clust ==c,]
data_in_clust_ = as.data.frame(data_in_clust)

clust_chain <-read.table("~/bayesmix/simdata/clustering_chain.csv", header = TRUE, sep = ",") 

num_chain <-read.table("~/bayesmix/simdata/numclust_chain.csv", header = TRUE, sep = ",") 
hist(num_chain[,1])
#eval_dens <-read.table("~/bayesmix/simdata/eval_dens.csv", header = TRUE, sep = ",") 
#grid<-read.table("~/bayesmix/simdata/grid.csv", header = TRUE, sep = ",") 

