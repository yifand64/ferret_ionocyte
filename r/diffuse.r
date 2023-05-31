

# build a diffusion map and run spectral clustering
library(diffusionMap)
diffusion_map <- function(data=NULL, dist.m=NULL, run_clustering=T, eps.val=0.05, k=2)
{
	if(is.null(data))
	{
		if(is.null(dist.m))
		{
			print("ERROR: provide a non-null distance or data matrix!")
		}else
		{	
			d = dist.m
		}
	}else
	{
		d = dist(data)
	}
	info("Building diffusion map..")
	dmap = diffuse(d, eps.val=eps.val) # compute diffusion map
	print(dmap)
	plot(dmap)
	if(run_clustering)
	{
		info("Running spectral clustering..")
		dkmeans = diffusionKmeans(dmap, k) 
		c = dkmeans$part
	}else
	{
		c = NULL
	}
	return(list("dmap"=dmap, "clusters"=c))
}
