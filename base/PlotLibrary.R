# Plot the frontier of the function
plot.frontier = function(data, derivative=TRUE, title=NULL) {
  df = data
  if (derivative) {
    dx = c(NA,df$OptimalCost[2:nrow(df)]-df$OptimalCost[1:(nrow(df)-1)])
    dx.r = range(dx, na.rm=T)
    df.r = range(df$OptimalCost)
    scale = (df.r[2]-df.r[1])/(dx.r[2]-dx.r[1])
    dx = dx*scale
    dx.r = range(dx, na.rm=T)
    offset = (df.r[1]-dx.r[1])
    dx = dx+offset
    df$Derivative = dx
  }
  p = ggplot(df, aes(x=Mutations, y=OptimalCost))
  if (derivative) {
    p = p + geom_line(aes(y=Derivative), color="red") + 
      scale_y_continuous(sec.axis = sec_axis(trans = ~(.-offset)/scale, name="Derivative"))
  }
  p = p + geom_line() + 
    xlab("Mutational Distance") + 
    ylab("Maximum Cost") + 
    ggtheme
  if (derivative) {
    p = p + theme(axis.line.y.right = element_line(color="red"), 
                  axis.ticks.y.right = element_line(color="red"), 
                  axis.text.y.right = element_text(color="red"), 
                  axis.title.y.right = element_text(color="red"))
  }
  if (!is.null(title)) {
    p = p + ggtitle(title) 
  }
  return(p)
}

# Plot mutational matrix
plot.mut.matrix = function(data, title=NULL) {
  input.seqs = data
  
  # Loop over all proteins and then all mutations
  sum.val = NULL
  base.sum.vals = NULL
  max.val = NULL
  base.max.vals = NULL
  for (curr.prot in 1:nrow(all.models[[1]])) {
    for (curr.seq in 1:nrow(input.seqs)) {
      # rescaled score needed
      scores = score.genome(genomicSequence = DNAString(input.seqs$OptimalSeq[curr.seq]), fits = all.models[[2]], 
                            index = curr.prot, mode = 1)/max.seq(all.models[[2]], index = curr.prot, mode=1)$MaxAffinity
      if (curr.seq==1) {
        base.sum = sum(scores)
        base.max = max(scores)
        base.sum.vals = c(base.sum.vals, base.sum)
        base.max.vals = c(base.max.vals, base.max)
        next
      } 
      sum.val = c(sum.val, sum(scores)/base.sum)
      max.val = c(max.val, max(scores)/base.max)
    }
  }
  sum.val = log2(sum.val)
  max.val = log2(max.val)
  
  # Create data.frame for plotting
  df = data.frame(Protein = rep(x = model.names, each=(nrow(input.seqs)-1)), 
                  Mutations = rep(x = 1:max(input.seqs$Mutations), nrow(all.models[[1]])),
                  Sum = sum.val, Max = max.val)
  
  # Plot
  p = ggplot(df, aes(Protein, Mutations)) + geom_tile(aes(fill=Sum)) + geom_text(aes(label = round(Sum, 2))) + scale_fill_gradient2(high="darkred", low="blue", midpoint=0, mid="red")
  print(p)
  p = ggplot(df, aes(Protein, Mutations)) + geom_tile(aes(fill=Max)) + geom_text(aes(label = round(Max, 2))) + scale_fill_gradient2(high="darkred", low="blue", midpoint=0, mid="red")
  print(p)
}