BiocManager::install("ggplot2")
BiocManager::install("ggdendro")
BiocManager::install("reshape2")

library(ggplot2)
library(ggdendro)
library(reshape2)

library("grid")

# Read in data
mouse <- read.csv(file = "trial_graphdata2.csv",
                  stringsAsFactors = TRUE)

head(mouse)

mouse.scaled <- mouse
mouse.scaled[, c(2:8)] <- scale(mouse.scaled[, 2:8])

mouse.matrix <- as.matrix(mouse.scaled[, -c(1)])
mouse.matrix <- as.matrix(mouse[, -c(1)])
rownames(mouse.matrix) <- mouse.scaled$Gene_ID
rownames(mouse.matrix) <- mouse$Gene_ID
mouse.dendro <- as.dendrogram(hclust(d = dist(x = mouse.matrix)))

dendro.plot <- ggdendrogram(data = mouse.dendro, rotate = TRUE)

print(dendro.plot)

dendro.plot <- dendro.plot + theme (axis.text.y = element_text(size = 3))

print(dendro.plot)

mouse.long <- melt(mouse.scaled, id = c("Gene_ID"))

mouse.long <- melt(mouse, id = c("Gene_ID"))

heatmap.plot <- ggplot(data = mouse.long, aes(x = variable, y = Gene_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 1))

print(heatmap.plot)

grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))

mouse.order <- order.dendrogram(mouse.dendro)

mouse.long$Gene_ID <- factor(x = mouse.long$Gene_ID,
                             levels = mouse.scaled$Gene_ID[mouse.order],
                             ordered = TRUE)

heatmap.plot <- ggplot(data = mouse.long, aes(x = variable, y = Gene_ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6))

print(heatmap.plot)
