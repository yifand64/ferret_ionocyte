

# functions for 3D plotting in signature coords

#plot.3d(scores, clusters=clusters, color="Notch.Pathway")
#plot.3d(scores, clusters=clusters, color.gene="Kdm5b")

library(rgl)
#library(Rcmdr)
#library(plot3D)
library(squash)

plot.3d_OLD <- function(seurat.obj, 
					rotating=F, 
					point.size=0.5,
					color.by= "Cell.Cycle_mean",
					counts=NULL,
					save.movie=NULL,
					save.movie.fps=24,
					save.movie.duration=5,
					save.movie.delete.frames=T,
					label.cex=1.5,
					label.col="black",
					legend.title=NULL,
					legend.width=1,
					legend.title.cex=3,
					legend.axis.cex=3,
					legend.horizontal=T,
					legend.shrink=0.5,
					legend.line=3, ### expands the gap between legend title and colorbar
					ss2=T,
					res.x=640,
					res.y=480,
					#old.coords=T,
					label.cell.types=TRUE,
					use.normalised.sig=T,
					spheres=T,
					cols=material.heat(100),
					use.curated=F,
					box="wireframe", # can be 'grey' or 'none'
					disp.max=0,
					three.d=T)
{

	if(!toupper(box) %in% c("GREY", "WIREFRAME", "NONE")){
		stop("'box' argument must be 'wireframe', 'grey', or 'none'")
	}else{
		box=tolower(box)
	}
	info("Collecting scores")

	if(use.curated)
	{
		sig = "_curated_200_"
	}else
	{
		sig = "_"
	}
	
	if(use.normalised.sig){
		#s = seurat.obj@data.info$Stem.Perm_mean_normalised
		s = seurat.obj@data.info$Stem_mean_normalised
		p = seurat.obj@data.info[,paste0("Paneth", sig,"mean_normalised")]
		tuft = seurat.obj@data.info[,paste0("Tuft", sig,"mean_normalised")]
		e = seurat.obj@data.info[,paste0("Enterocyte", sig,"mean_normalised")]
		g = seurat.obj@data.info[,paste0("Goblet", sig,"mean_normalised")]
		#m = seurat.obj@data.info[,paste0("Microfold", sig,"mean_normalised")]
		end = seurat.obj@data.info[,paste0("Endocrine", sig,"mean_normalised")]
	}else{
		
		#s = seurat.obj@data.info$Stem.Perm_mean_normalised
		s = seurat.obj@data.info$Stem_mean
		p = seurat.obj@data.info[,paste0("Paneth", sig,"mean")]
		tuft = seurat.obj@data.info[,paste0("Tuft", sig,"mean")]
		e = seurat.obj@data.info[,paste0("Enterocyte", sig,"mean")]
		g = seurat.obj@data.info[,paste0("Goblet", sig,"mean")]
		#m = seurat.obj@data.info[,paste0("Microfold", sig,"mean_normalised")]
		end = seurat.obj@data.info[,paste0("Endocrine", sig,"mean")]
	
	}
	print("huh")
	
	# if(old.coords)
	# {
	if(ss2){
		z =  - s # + tuft
		y = p - end
		x = e - g
	}else{
		info("Using 10X coords (try and compensate for lack of Paneth cells)")
		z =  - s # + tuft
		y = tuft - end
		x = e - g
	}
		
	# }else{
	# 	z = -s
	# 	y = (e + p + 2*m) - (end + g + tuft)
	# 	x = (2*m + tuft) - (p + g)

	# }
	

	scores = data.frame(x,y,z)
	print(head(scores))


	# if(rotating)
	# {
	if(!is.null(color.by))
	{
		cat(sprintf("Coloring points by %s..\n", color.by))
		if(is.null(legend.title))
		{
			legend.title = color.by
		}
		# if(color.palette != "SamR")
		# {
		#  	if(length(color.palette) > 1)
		#  	{
		#  		cat(sprintf("Building heatmap colors from gradient using %i colors: ", length(color.palette)))
		#  		info(sprintf(color.palette))
		#  		col.function = colorRampPalette(color.palette)
	 # 		}else{
	 # 			cat(sprintf("Using brewer palette: %s", color.palette))
	 # 			col.function = function(x){rev(colorRampPalette(brewer.pal(11, color.palette)))(x)}
	 # 		}
 	# 	}else
 	# 	{
 	# 		col.function = material.heat #get.hmap.col
 	# 	}
 		#open3d()
 		
 		
		coi = as.numeric(unlist(fetch.data(seurat.obj, color.by)))
 		if(disp.max>0){coi[coi>disp.max] <- disp.max}
 		color.map = makecmap(coi, colFn = colorRampPalette(cols))
 		

		## clear scene:
		clear3d("all")

		## setup env:
		
		if(three.d)
		{
			#bg3d(color="#887777")
			#open3d(antialias=4)
			par3d(windowRect = c(0, 0, res.x, res.y)) # make the window large
			par3d(zoom = 1.1)


			light3d()	
			aspect3d(1,1,1)
			if(spheres){type="s"}else{type="p"}
			if(box=="wireframe"){b=T; a=T}else{b=F; a=F}
			plot3d(scores$x, scores$y, scores$z, col=cmap(coi, map=color.map),size=point.size, zlab="",xlab="",ylab="", type=type, box=b, axes=a)
			#movie3d(spin3d(axis = c(0,0,1), rpm = 10), duration=6,  type = "png")
			info("Rotating plot")

			if(label.cell.types)
			{
					if(ss2){
						rgl.texts(x=8, -3, -1.75,"Enterocyte", col=label.col, cex=label.cex)
						rgl.texts(x=0,15,-1.5,"Paneth", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=0,-8,-1,"Enteroendocrine", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=-8,4,-2,"Goblet", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=0,5,-3,"TA progenitor", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=-3,0,-5,"Lgr5+ Stem", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=-3,0,-1.5,"Tuft", col=label.col, cex=label.cex, useFreeType=T)
					}else{ #10X dataset needs d
						rgl.texts(x=2, -0.5, 0,"Enterocyte", col=label.col, cex=label.cex)
						#rgl.texts(x=0,15,-1.5,"Paneth", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=0.5,-2, 0,"Enteroendocrine", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=-2, 1, 0 ,"Goblet", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=0, -1,-0.3,"TA progenitor", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=-0.5,0, -1,"Stem", col=label.col, cex=label.cex, useFreeType=T)
						rgl.texts(x=0.5, 2, 0,"Tuft", col=label.col, cex=label.cex, useFreeType=T)
					}
					
				# }else{
				# 	text3d(x=0,8, 1.5,"Enterocyte", useFreeType=T)
				# 	text3d(x=-8,-4,1.5,"Goblet", useFreeType=T)
				# 	text3d(x=0,-8, 1.5,"Enteroendocrine", useFreeType=T, font="Arial")
				# 	text3d(x=-8,4, 1.5,"Paneth", useFreeType=T)
				# 	text3d(x=0, 5, 0,"TA Progenitor", useFreeType=T)
				# 	text3d(x=10,-6, 1.5,"Tuft",  useFreeType=T)
				# 	text3d(x=10, -2, 1.5,"Microfold",  useFreeType=T)
				# 	text3d(x=-3,-3,-3,"Lgr5+ Stem",  useFreeType=T)
				# }
				
			}
			## color legend
			library(fields)
			bgplot3d(suppressWarnings(
				image.plot( 
					axis.args = list(cex.axis = legend.axis.cex, line=1), ## set number size and move the ticks away from colour bar
					legend.args=list(text=legend.title, side=2, font=2, line=2.5, cex=legend.title.cex), # set legend title position, side (left of the colour bar) and distance from the colour bar
					legend.width=legend.width, 
					legend.only=TRUE, 
					zlim= c(min(coi), max(coi)), 
					nlevel=128, 
					col=cols, 
					lwd=0,
					smallplot=c(0.88, 0.9, .3,.7),
					bty = "n" ,
					#legend.lab=legend.title,
					#horizontal=legend.horizontal,
					legend.shrink=legend.shrink)
			))

			if(box=="grey"){rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('grey100'))}			

			if(!is.null(save.movie)){
				info(sprintf("Saving movie to --> %s", save.movie))
				movie3d(spin3d(axis = c(0,0,1), rpm = 10), 
					type = "gif", 
					movie=save.movie, dir=getwd(), 
					duration=save.movie.duration, 
					clean=save.movie.delete.frames,
					fps=save.movie.fps)
			}
		}else
		{
			plot3d(scores$x, scores$y, scores$z, col=cmap(coi, map=color.map),size=8, zlab="",xlab="",ylab="")
		}
		

		# print(coi)
		#plot3d(scores$x, scores$y, scores$z, col=cmap(coi, map=color.map),size=4)
		
		#legend3d("topright", legend = 'Proliferation', pch = 16, col = cmap(coi, map=color.map), cex=1, inset=c(0.02))

	}else{
		plot3d(scores$x, scores$y, scores$z)
	}
		
	# }else{
	# 	scatter3D(scores$x, scores$y, scores$z)
	# 		with(scores, scatter3D(x = x, y = y, z = z, colvar = Cell.Cycle,
	# 			pch = 16, cex = 1.5, xlab = "Entero/Goblet", ylab = "Paneth/Endo",
	# 			zlab = "Stemness", clab = c("Proliferation"),
	# 			main = "Gut Differentiation Traj.", ticktype = "detailed",
	# 			theta = 5, d = 2,
	# 			colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75))
	# 		)
	# }
	

	#return(scores)
}




