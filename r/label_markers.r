marker.status <- Vectorize(function(micro,rnaseq, t){

	r=NULL
	if(micro > t & rnaseq > t){
		r="Stem both"
	}else{
		if(micro < -t & rnaseq < -t){
			r="Paneth both"
		} else{
			if(micro > t){
				r="Stem Microarray"
			}else
	 		{
				if(rnaseq > t)
				{
					r="Stem RNA-seq"
				}else{
					if(micro < -t)
					{
						r="Paneth Microarray"
					}else
					{
						if(rnaseq < -t)
						{
							r="Paneth RNA-seq"
						}
					}
				}
			}
		}
	}
	if(is.null(r))
	{
		r = "None"
	}
	return (r)
})