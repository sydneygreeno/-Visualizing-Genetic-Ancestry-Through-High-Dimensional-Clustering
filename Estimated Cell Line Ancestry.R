
## organize data set

cell_data <- cell_line_ancestry
cell_data <- cell_line_ancestry[, c(35:41)]
#take away cell data
cell_data2<-filter(cell_line_ancestry, p != "cell")
cell_data_data<-cell_data2[,c(35:41)]

library(Rtsne)
library(plotly)
library(dbscan)
library(RColorBrewer)
RColorBrewer::display.brewer.all()

## make a list for the color/marker assignments to super/sub population

color_assign <- c("#ce4278","#5bb645", "#46b9d1", "#646aaa", "#dc964d")
marker_assign <-c("circle-open", "square-open", "diamond-open", "cross-open", "triangle-up-open", 
                  "triangle-left-open", "triangle-right-open")
marker_assign <- setNames(marker_assign, c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"))
marker_assign2 <- setNames(marker_assign, c("CLM", "MXL", "PEL", "PUR"))
marker_assign3 <- setNames(marker_assign, c("CDX", "CHB", "CHS", "JPT", "KHV"))
marker_assign4 <- setNames(marker_assign, c("CEU", "FIN", "GBR", "IBS", "TSI"))
marker_assign5 <- setNames(marker_assign, c("BEB", "GIH", "ITU", "PJL", "STU"))
ma <- c(marker_assign,marker_assign2,marker_assign3,marker_assign4,marker_assign5)


## set the seed for best starting point
set.seed(60)
tsne <- Rtsne(cell_data_data, check_duplicates=FALSE, perplexity=30, max_iter=1000, theta=0.9)
pdb <- cbind(tsne$Y, cell_data2$p)
pdb <- data.frame(tsne$Y, p=cell_data2$p)
options(warn = -1)
fig <- plot_ly(data=pdb, x=~X1, y=~X2, type= "scatter", mode="markers", color= ~cell_data2$p, colors= "Set1")
fig <- fig %>%
  layout(title="R t-SNE NO CELL VER",legend=list(title=list(text="Ancestry Population")))
fig
## use DBSCAN for sub population symbols with clustering
dbscan_res <- dbscan(tsne$Y, eps=2, minPts=3)
fig2 <- plot_ly(data=pdb, x=~X1, y=~X2, type= "scatter", mode="markers", symbol=~cell_data2$sp, 
                symbols=ma, 
                color=~cell_data2$p, colors=color_assign, marker=list(size=5), opacity=0.7)
fig2 <- fig2 %>%
  layout(title="Estimated Cell Line Ancestry",legend=list(title=list(text="Ancestry")))
fig2


##  SQL in order to determine admixed populations  ##

# need row names to label
rownames(cell_data_data) <- cell_data2$id
names(dbscan_res$cluster) <- rownames(cell_data_data)
clusters <- tibble::enframe(dbscan_res$cluster)
## will assign each dbscan cluster(value) with a subpop(sp)
cluster_assign <- 
  ## add columns from y to x matching rows based on id
  dplyr::left_join(clusters, cell_data2, by=c("name"="id"))|>
  # isolates the var we want, name from celldata2 & value from clusters
  dplyr::select(name, value, sp) |>
  #  groups subpops together & labels
  dplyr::group_by(value, sp) |>
  ## counts the number of subpops in each cluster
  dplyr::summarise(n=dplyr::n()) |>
  dplyr::group_by(value) |>
  ## finds the max subpop in each cluster
  dplyr::mutate(max=max(n)) |>
  ##  filters out the non max subpops
  dplyr::filter(n==max)
assignments <- 
  # same as cluster_assign
  dplyr::left_join(clusters, cell_data2, by=c("name"="id"))|>
  dplyr::select(name, value, sp) |>
  dplyr::group_by(value, sp) |>
  # add columns from y to x matching rows based on cluster(value)
  ## sp.x is max subpop according to dbscan sp.y is the real subpop
  dplyr::left_join(cluster_assign, by=c("value"="value")) |>
  ## add new variable matches if subpop(sp) is equal to max cluster subpop
  dplyr::mutate(matches= sp.x==sp.y) |>
  ## filter the ones that are true
  dplyr::filter(matches==TRUE) |>
  ## leave data frame ungrouped
  dplyr::ungroup() |>
  ## get a numeric value
  dplyr::count()
assignments
## counts the non admixed groups


## will assign each dbscan cluster(value) with a subpop(sp)
define_if_admixed <- 
  ## add columns from y to x matching rows based on id
  dplyr::left_join(clusters, cell_data2, by=c("name"="id"))|>
  # isolates the var we want, name from celldata2 & value from clusters
  dplyr::select(name, value, sp) |>
  #  groups subpops together & labels
  dplyr::group_by(value, sp) |>
  ## counts the number of subpops in each cluster
  dplyr::summarise(n=dplyr::n()) |>
  dplyr::group_by(value) |>
  ## finds the max subpop in each cluster
  dplyr::mutate(max=max(n)) |>
  ##  filters out the non max subpops
  dplyr::filter(n==max)

isolate_admixed <- 
  # same as cluster_assign
  dplyr::left_join(clusters, cell_data2, by=c("name"="id"))|>
  dplyr::select(name, value, sp) |>
  dplyr::group_by(value, sp) |>
  # add columns from y to x matching rows based on cluster(value)
  ## sp.x is max subpop according to dbscan sp.y is the real subpop
  dplyr::left_join(define_if_admixed, by=c("value"="value")) |>
  ## add new variable matches if subpop(sp) is equal to max cluster subpop
  dplyr::mutate(matches= sp.x==sp.y)
isolate_admixed


## color non-admixed population in gray
## pink shows admixed populations

color_assign2= c("#ce4278" ,"#858a84")

set.seed(60)
tsne <- Rtsne(cell_data_data, check_duplicates=FALSE, perplexity=30, max_iter=1000, theta=0.9)
pdb <- cbind(tsne$Y, cell_data2$p)
pdb <- data.frame(tsne$Y, p=cell_data2$p)
options(warn = -1)
fig <- plot_ly(data=pdb, x=~X1, y=~X2, type= "scatter", mode="markers", color= ~cell_data2$p, colors= "Set1")
fig <- fig %>%
  layout(title="R t-SNE NO CELL VER",legend=list(title=list(text="Ancestry Population")))

dbscan_res <- dbscan(tsne$Y, eps=2, minPts=3)
fig3 <- plot_ly(data=pdb, x=~X1, y=~X2, type= "scatter", mode="markers", symbol=~cell_data2$sp, 
                symbols=ma, 
                color=~isolate_admixed$matches, colors=color_assign2, marker=list(size=5), opacity=0.7)
fig3 <- fig3 %>%
  layout(title="ADMIXED",legend=list(title=list(text="Ancestry")))
fig3

## color by dbscan clusters

set.seed(60)
tsne <- Rtsne(cell_data_data, check_duplicates=FALSE, perplexity=30, max_iter=1000, theta=0.9)
pdb <- cbind(tsne$Y, cell_data2$p)
pdb <- data.frame(tsne$Y, p=cell_data2$p)
options(warn = -1)
fig <- plot_ly(data=pdb, x=~X1, y=~X2, type= "scatter", mode="markers", color= ~cell_data2$p, colors= "Set1")
fig <- fig %>%
  layout(title="R t-SNE NO CELL VER",legend=list(title=list(text="Ancestry Population")))

dbscan_res <- dbscan(tsne$Y, eps=2, minPts=3)
fig4 <- plot_ly(data=pdb, x=~X1, y=~X2, type= "scatter", mode="markers", symbol=~cell_data2$sp, symbols=ma, 
                color=dbscan_res$cluster, colors="Paired",
                marker=list(size=5), opacity=0.7)
fig4 <- fig4 %>%
  layout(title="DBSCAN",legend=list(title=list(text="Ancestry")))
fig4

















