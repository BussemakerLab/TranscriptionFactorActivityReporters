library(ggplot2)
library(reshape2)
library(dplyr)

ggtheme = theme_bw() + theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
                             axis.title.x=element_text(margin=margin(20,0,0,0)), aspect.ratio=1, axis.line=element_line(color="black", size=1), 
                             axis.ticks=element_line(color="black", size=1), panel.border=element_blank(), legend.justification=c(1,0),
                             legend.title.align=0.5, legend.position=c(.99,0.01))

# Plot the frontier of the function
plot.frontier = function(data, derivative=TRUE, title=NULL) {
  if (class(data)=="list") { # Check for 'backwards compatability'
    df = data$MaximalSequences
  } else if (class(data)=="data.frame") { 
    df = data
  } else {
    stop(paste0("cannot deal with input datatype: ", class(data)))
  }
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

# Visualize the range of valid data
plot.cost.range = function(data, title=NULL) {
  if (class(data)=="list") { # Check for 'backwards compatability'
    df = data$SequenceExamples
  } else {
    stop(paste0("cannot deal with input datatype: ", class(data)))
  }
  df = (df!="NA")*1
  rownames(df) = as.character(0:(nrow(df)-1))
  cost.idx.0 = which(df[1,]==1)
  colnames(df) = as.character((0-cost.idx.0+1):(ncol(df)-cost.idx.0))
  df = melt(df)
  p = ggplot(df, aes(x=Var2, y=Var1)) + 
    geom_raster(aes(alpha=value), fill="orangered2") +
    xlab("Cost") +
    ylab("Mutational Distance") +
    scale_alpha(range=c(0,1)) +
    ggtheme +
    theme(legend.position = "none")
  if (!is.null(title)) {
    p = p + ggtitle(title) 
  }
  return(p)
}

# Create helper function to compute affinity scores for a given sequence (both max and sum)
score.seq.helper = function(sequences, base.seq, input.models) {
  # Create scoring object
  nModels = length(input.models)
  score.object = new(RapidScore)
  
  # Load models
  for (i in 1:nModels) {
    nuc = input.models[[i]]$NB
    dinuc = input.models[[i]]$DB
    score.object$add(nuc, dinuc)
  }
  
  # Score base
  base.scores = score.object$scoreBulkChar(base.seq)
  base.max = rep(0, nModels)
  base.sum = rep(0, nModels)
  for (i in 1:nModels) {
    base.max[i] = max(exp(base.scores[[i]]))
    base.sum[i] = sum(exp(base.scores[[i]]))
  }
  
  # Matrix output for sum/max
  sum.ratio = matrix(data = 0, ncol = nModels, nrow = length(sequences))
  max.ratio = sum.ratio
  for (curr.seq in 1:length(sequences)) {
    if (is.na(sequences[curr.seq])) {
      # Handle the "NA" sequence case
      sum.ratio[curr.seq,] = NA
      max.ratio[curr.seq,] = NA
    } else {
      curr.score = score.object$scoreBulkChar(sequences[curr.seq])
      for (j in 1:nModels) {
        sum.ratio[curr.seq, j] = sum(exp(curr.score[[j]]))/base.sum[j]
        max.ratio[curr.seq, j] = max(exp(curr.score[[j]]))/base.max[j]
      }
    }
  }
  # Name matrix for easy reading
  colnames(sum.ratio) = names(input.models)
  colnames(max.ratio) = names(input.models)
  
  return(list(Sum=sum.ratio, Max=max.ratio, BaseSum=base.sum, BaseMax=base.max))
}

# Plot mutational impact
plot.mut.impact = function(data, input.models, line.plot=FALSE, base.seq=NULL, mut.slice=NA, cost.slice=NA, 
                           title=NULL, reorder=NULL, matrix.plot.text=TRUE, matrix.plot.text.size=3, 
                           matrix.plot.gradient.sym=TRUE, matrix.plot.mut.divs=5, line.plot.legend.cols=1,
                           line.plot.legend.pos = c("side", "bottomleft", "bottomright", "topleft", "topright")) {
  leg.pos = match.arg(line.plot.legend.pos)
  if (class(data)=="list") { # Check for 'backwards compatability'
    if (!is.na(mut.slice)) {
      mut.slice = round(mut.slice) # Ensure integer value
      # Make a mutational slice. First, check to see if mut.slice value is within range
      if (mut.slice<1 || mut.slice >= nrow(data$SequenceExamples)) {
        stop("Mutational slice value is out of range!")
      }
      # Create slice and replace NAs
      seqs = data$SequenceExamples[mut.slice+1,]   # Mutational distance slice will be offset by 1
      seqs[seqs=="NA"] = NA
      cost.idx.0 = which(data$SequenceExamples[1,]!="NA")
      if (is.null(base.seq)) {
        # Find the base sequence just from the sequence examples themselves
        base.seq = data$SequenceExamples[1,cost.idx.0]
      }
      varying.label = "Cost"
      varying.ticks = (0-cost.idx.0+1):(ncol(data$SequenceExamples)-cost.idx.0)
      if (is.null(title)) {
        title = paste0("Mutational Slice at ", mut.slice)
      }
    } else if (!is.na(cost.slice)) {
      cost.slice = round(cost.slice) # Ensure integer value
      # Make a cost slice. First, check to see if cost.slice value is within range
      cost.idx.0 = which(cost.fun.out$SequenceExamples[1,]!="NA")
      cost.range = (0-cost.idx.0+1):(ncol(data$SequenceExamples)-cost.idx.0)
      if (cost.slice<range(cost.range)[1] || cost.slice>range(cost.range)[2]) {
        stop("Cost slice value is out of range!")
      }
      # Create slice and replace NAs
      seqs = data$SequenceExamples[,which(cost.range==cost.slice)]
      seqs = seqs[-1]   # Remove 0 mutational distance sequence
      seqs[seqs=="NA"] = NA
      if (is.null(base.seq)) {
        base.seq = cost.fun.out$MaximalSequences$OptimalSeq[1]
      }
      varying.label = "Mutations"
      varying.ticks = data$MaximalSequences$Mutations
      varying.ticks = varying.ticks[-1]
      if (is.null(title)) {
        title = paste0("Cost Slice at ", cost.slice)
      }
    } else {
      # Traditional slice
      seqs = data$MaximalSequences$OptimalSeq
      if (is.null(base.seq)) {
        # Find the base sequence just from the sequence examples themselves
        base.seq = seqs[1]
      }
      seqs = seqs[-1]
      varying.label = "Mutations"
      varying.ticks = data$MaximalSequences$Mutations
      varying.ticks = varying.ticks[-1]
      if (is.null(title)) {
        title = "Best Sequences"
      }
    }
  } else if (!is.na(mut.slice) || !is.na(cost.slice)) {
    stop("cannot create slice with given input!")
  } else {  # Only the maximal sequences are provided. Default to providing a sweep over the maximum
    seqs = data$OptimalSeq
    if (is.null(base.seq)) {
      base.seq = seqs[1]
    }
    seqs = seqs[-1]
    varying.label = "Mutations"
    varying.ticks = data$Mutations
    varying.ticks = varying.ticks[-1]
    if (is.null(title)) {
      # Find the base sequence just from the sequence examples themselves
      title = "Best Sequences"
    }
  }
  
  # Now score sequences and label them appropriately; also, log2-transform
  scored.output = score.seq.helper(sequences = seqs, base.seq = base.seq, input.models = input.models)
  sum.out = log2(scored.output$Sum)
  max.out = log2(scored.output$Max)
  rownames(sum.out) = varying.ticks
  sum.out = melt(sum.out, value.name = "Score", varnames = c(varying.label, "Protein"))
  rownames(max.out) = varying.ticks
  max.out = melt(max.out, value.name = "Score", varnames = c(varying.label, "Protein"))
  if (!is.null(reorder)) {
    sum.out$Protein = factor(sum.out$Protein, levels = reorder)
    max.out$Protein = factor(max.out$Protein, levels = reorder)
  }
  
  # Create line plot if requested
  
  # Things to implement: 5) plot area, legend, labels?
  
  if (line.plot) {
    # Need to compute shape/line type parameters
    nProts = length(input.models)
    if (nProts>36) {
      # Need to split amongst line color, shape, and type
      ltypes = rep(0, nProts)
      stypes = rep(0, nProts)
      ctypes = rep(0, nProts)
      nVariants = ceiling(nProts^(1/3))
      currIter = 1
      for (line.type in (linetype_pal()(nVariants))) {
        if (currIter > nProts) {
          break
        }
        for (color.type in (hue_pal()(nVariants))) {
          if (currIter > nProts) {
            break
          }
          for (shape.type in (shape_pal()(nVariants))) {
            if (currIter > nProts) {
              break
            }
            stypes[currIter] = shape.type
            ctypes[currIter] = color.type
            currIter = currIter + 1
          }
        }
      }
    } else {
      # Line type is fixed
      ltypes = rep(1, nProts)
      stypes = rep(0, nProts)
      ctypes = rep(0, nProts)
      nVariants = ceiling(sqrt(nProts))
      currIter = 1
      for (color.type in (hue_pal()(nVariants))) {
        if (currIter > nProts) {
          break
        }
        for (shape.type in (shape_pal()(nVariants))) {
          if (currIter > nProts) {
            break
          }
          stypes[currIter] = shape.type
          ctypes[currIter] = color.type
          currIter = currIter + 1
        }
      }
    }
    
    legend.theme = theme_bw() +
      theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
            axis.title.x=element_text(margin=margin(20,0,0,0)), aspect.ratio=1, axis.line=element_line(color="black", size=1), 
            axis.ticks=element_line(color="black", size=1), panel.border=element_blank(), legend.title = element_blank(), 
            legend.text = element_text(size=12))
    if (leg.pos=="bottomleft") {
      legend.theme = legend.theme + 
        theme(legend.box = "horizontal", legend.position = c(0.02, 0.02), legend.justification = c(0,0))
    } else if (leg.pos=="bottomright") {
      legend.theme = legend.theme +
        theme(legend.box = "horizontal", legend.position = c(1,.02), legend.justification = c(1,0))
    } else if (leg.pos=="topleft") {
      legend.theme = legend.theme +
        theme(legend.box = "horizontal", legend.position = c(.02, 1), legend.justification = c(0,1))
    } else if (leg.pos=="topright") {
      legend.theme = legend.theme +
        theme(legend.box = "horizontal", legend.position = c(1,1), legend.justification = c(1,1))
    }
    
    sum.plot = ggplot(sum.out, aes_string(x=varying.label, y="Score", color="Protein")) +
      geom_line(aes(linetype=Protein)) +
      scale_linetype_manual(values=ltypes) +
      scale_color_manual(values=ctypes) + 
      geom_point(aes(shape=Protein)) + 
      scale_shape_manual(values=stypes) + 
      ggtitle(paste0(title, " - Sum")) +
      guides(linetype = guide_legend(ncol=line.plot.legend.cols)) + 
      legend.theme
    max.plot = ggplot(max.out, aes_string(x=varying.label, y="Score", color="Protein")) +
      geom_line(aes(linetype=Protein)) +
      scale_linetype_manual(values=ltypes) +
      scale_color_manual(values=ctypes) + 
      geom_point(aes(shape=Protein)) + 
      scale_shape_manual(values=stypes) + 
      ggtitle(paste0(title, " - Max")) + 
      guides(linetype = guide_legend(ncol=line.plot.legend.cols)) + 
      legend.theme
  } else {
    # For rescaling the gradient 
    sum.score.range = range(sum.out$Score, na.rm = T)
    max.score.range = range(max.out$Score, na.rm = T)
    
    sum.plot = ggplot(sum.out, aes_string("Protein", varying.label)) + 
      geom_tile(aes(fill=Score)) + 
      ggtitle(paste0(title," - Sum")) + 
      theme_bw() +
      theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
            axis.title.x=element_text(margin=margin(20,0,0,0)), axis.line.x=element_line(color="black", size=1), 
            axis.ticks=element_line(color="black", size=1), panel.border=element_blank(),
            axis.text.x.bottom = element_text(angle = 45, hjust = 1), axis.text.x.top = element_text(hjust=.5, size=15),
            axis.title.x.top = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    max.plot = ggplot(max.out, aes_string("Protein", varying.label)) + 
      geom_tile(aes(fill=Score)) + 
      ggtitle(paste0(title," - Max")) + 
      theme_bw() +
      theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
            axis.title.x=element_text(margin=margin(20,0,0,0)), axis.line.x=element_line(color="black", size=1), 
            axis.ticks=element_line(color="black", size=1), panel.border=element_blank(),
            axis.text.x.bottom = element_text(angle = 45, hjust = 1), axis.text.x.top = element_text(hjust=.5, size=15),
            axis.title.x.top = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # Change color scale gradient preferences
    if (matrix.plot.gradient.sym) {
      sum.plot = sum.plot + scale_fill_gradient2(high="red", low="blue", midpoint=0, mid="white")
    } else {
      if (sum.score.range[1]>=0 || sum.score.range[2]<=0) {
        sum.plot = sum.plot + scale_fill_gradient2(high="red", low="blue", midpoint=0, mid="white")
      } else {
        sum.plot = sum.plot + scale_fill_gradientn(colors=c("blue", "white", "red"), values=rescale(c(sum.score.range[1],0,sum.score.range[2])))
      }
    }
    # Change color scale gradient preferences
    if (matrix.plot.gradient.sym) {
      max.plot = max.plot + scale_fill_gradient2(high="red", low="blue", midpoint=0, mid="white")
    } else {
      if (max.score.range[1]>=0 || max.score.range[2]<=0) {
        max.plot = max.plot + scale_fill_gradient2(high="red", low="blue", midpoint=0, mid="white")
      } else {
        max.plot = max.plot + scale_fill_gradientn(colors=c("blue", "white", "red"), values=rescale(c(max.score.range[1],0,max.score.range[2])))
      }
    }
    
    # Relabel and re-order mutation-axis
    if (varying.label=="Mutations") {
      divs = (1:floor(max(data$MaximalSequences$Mutations)/matrix.plot.mut.divs))*matrix.plot.mut.divs
      sum.plot = sum.plot + scale_y_reverse(breaks = c(0, divs), labels = c("WT", as.character(divs)))
      max.plot = max.plot + scale_y_reverse(breaks = c(0, divs), labels = c("WT", as.character(divs)))
      base.seq.y.idx = 0
    } else {
      base.seq.y.idx = max(sum.out[,1])+1
    }
    # Add explicit ratio values
    if (matrix.plot.text) {
      sum.plot = sum.plot + geom_text(aes(label=round(Score, 2)), size=matrix.plot.text.size)
      base.scores = data.frame(Protein=names(input.models), Score=as.character(signif(scored.output$BaseSum,2)))
      sum.plot = sum.plot + geom_text(data = base.scores, aes(label=Score, x=Protein, y=base.seq.y.idx), size=matrix.plot.text.size, fontface = "bold")
      max.plot = max.plot + geom_text(aes(label=round(Score, 2)), size=matrix.plot.text.size)
      base.scores = data.frame(Protein=names(input.models), Score=as.character(signif(scored.output$BaseMax,2)))
      max.plot = max.plot + geom_text(data = base.scores, aes(label=Score, x=Protein, y=base.seq.y.idx), size=matrix.plot.text.size, fontface = "bold")
    }
  }
  print(sum.plot)
  print(max.plot)
  return(list(Sum=sum.plot, Max=max.plot))
}