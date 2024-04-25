  
 if (iepsout==1){
	    print(paste("postscipt driver",fileout))      		
	    setEPS()
	    par(family = "serif")
	    par(xpd=FALSE)
        postscript(file=fileout, width=8, height=9,horizontal=FALSE)        		 

} else if (iepsout==2) {
		par(family = "Arial")
		par(xpd=FALSE)
		cairo_ps(file=fileout,width=8,height=7,family="sans")
	    print(paste("cairo driver 2",fileout))   
	       		
}    else if (iepsout==3) {
		par(family = "Arial")
		par(xpd=FALSE)
		cairo_pdf(file=paste0(fileout,".pdf"),width=8,height=7)
	    print(paste("pdf driver",fileout))      		

} else if (iepsout==4) {
		par(family = "Arial")
		par(xpd=FALSE)
                par(mar = c(1, 1, 1, 1))
		cairo_pdf(file=paste0(fileout,".pdf"),width=7,height=7)
	    print(paste("pdf driver",fileout))      		

} else if (iepsout==5) {
		par(family = "Arial")
		par(xpd=FALSE)
                par(mar = c(1, 1, 1, 1))
                par(oma = c(0,0,0,0))
		pdf(file=paste0(fileout,".pdf"),width=7,height=7.5)
	         print(paste("pdf driver",fileout))      		

} else if (iepsout==6) {
		par(family = "Arial")
		par(xpd=FALSE)
                par(mar = c(0, 0, 0, 0))
                par(oma = c(0,0,0,0))
                par(plt=c(0.05,0.95,0.05,0.95))
		png(filename = fileout, width = 1000, height = 1000, units = "px", pointsize = 12,bg = "white",  res = NA)
	        print(paste("png",fileout))   
} else if (iepsout==7) {
		par(family = "Arial")
		par(xpd=FALSE)
		cairo_ps(file=fileout,width=7,height=7,family="sans")
	    print(paste("cairo driver 2",fileout))   
	       		
} else if (iepsout==8) {
		par(family = "Arial")
		par(xpd=FALSE)
		cairo_ps(file=fileout,width=6,height=9,family="sans")
	    print(paste("cairo driver 2",fileout))   
}  else if (iepsout==9) {
		par(family = "Arial")
		par(xpd=FALSE)
		cairo_ps(file=fileout,width=5,height=7,family="sans")
	    print(paste("cairo driver 2",fileout))   
}  else if (iepsout==10) {
		par(family = "Arial")
		par(xpd=FALSE)
		cairo_ps(file=fileout,width=6,height=7,family="sans")
	    print(paste("cairo driver 2",fileout))   
	       		
}  else if (iepsout==11) {
		par(family = "Arial")
		par(xpd=FALSE)
                par(mar = c(1, 1, 1, 1))
                par(oma = c(0,0,0,0))
		pdf(file=paste0(fileout,".pdf"),width=5,height=7.5)
	         print(paste("pdf driver",fileout))      		

}

	       		

