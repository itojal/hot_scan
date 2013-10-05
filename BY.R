args=(commandArgs(TRUE))

window = as.integer(args[1])
path   = args[2]
sample = args[3]
chr    = args[4]

file_in 	= paste(path,'/',sample,'/',window,'/',chr,'.txt',sep='');

if( (file.exists(file_in)) && (file.info(file_in)$size) ){

	data 		= read.table(file_in,header=TRUE)
	tot_rows  	= length(data[,1]) 

	p.unadjusted  = c()
	p.adjusted    = c()

	p.unadjusted  = data[,6]

	l = 1
	p.adjusted = p.adjust(p.unadjusted,"BY")
	capture.output(
		while ( l <= tot_rows ){ 
    		if ( p.adjusted[l] >= 0 & p.adjusted[l] < 0.05 ){
	      		x = as.numeric(data[l,])
    	  		print(paste(x[2],x[3],x[4],x[5],p.adjusted[l],' ',sep=" "))
    		} 
	    	l = l + 1
	  	} 
	    ,file = paste(path,'/',sample,'/',window,'/BY/',chr,'.txt',sep='')
	)
}
