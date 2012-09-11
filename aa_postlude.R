

DrawRefLine <- function(yval){
print('\n');
num <- -log10(yval)
#downViewport('pvalsClipped');
panel.abline(h=transformation(num),col="red3");
}


MakeData <- function(yank_loc){

x <- read.delim(yank_loc)
dataFrame <- as.data.frame(x)
print(dataFrame);
return(dataFrame);

}

LabelRefSNP <- function(ref){
specialref <- c('rs2476601','rs9273363','rs3842727')
reflabel = ref;
if (ref %in% specialref){
	reflabel <- paste(ref,"*",sep='')
	}
return(reflabel);
}

LowerLabel <- function(label_pos, ymax){
if (label_pos > ymax){
new_pos = label_pos - .2;
label_pos <- LowerLabel(new_pos, ymax)
}
else{
return(label_pos);
}
}

LabelSNPs <- function(dataFrame){
#mycols <- c("chr","snp_name","snp_position","p.value",'pch')
mycols <- c("chr","snp_name","snp_position","p.value","locuszoom_snp")
indices <- match(mycols, colnames(dataFrame))
print(indices);
text_col = "black";


for (i in 1:nrow(dataFrame)){
row <- dataFrame[i,]

if(length(row)!=0) {
	print(row);

	gold_ref = row[indices[5]]
	print(gold_ref);
	lz_ref = args[['refsnp']];
	print(lz_ref);
	pos = row[[indices[3]]];
	chr = row[[indices[1]]];
	gold_chrpos <- paste("chr",chr,":",pos,sep='')
	gold_lz = row[[indices[5]]];
	print(gold_lz);
	
	if (gold_chrpos == lz_ref){
		print(row[indices[2]]);
		orig_p =as.numeric(row[indices[4]])
		reflabel <- LabelRefSNP(gold_lz)
		if ( orig_p == 0){
		new_p = 100;
		label_y = 110;
		}
		else {
			new_p <- -log10(orig_p)
			label_y = new_p + 1.5;
			}
		ymax = yRange[2];
		ydiff = ymax - new_p;
		
		if (ymax < 8){
		label_y = new_p + .8;
		label_y <- LowerLabel(label_y, ymax);
		}
		else if (ymax < 15){
		label_y = new_p + 1.2;
		label_y <- LowerLabel(label_y, ymax);
		}
		else if (ymax < 30){
		label_y = new_p + 2;
		label_y <- LowerLabel(label_y, ymax);
		}
		else if (ymax < 50){
		label_y = new_p + 3;
		label_y <- LowerLabel(label_y, ymax);
		}
		else if (ymax < 70){
		label_y = new_p + 5;
		label_y <- LowerLabel(label_y, ymax);
		}
		panel.text(x=(as.numeric(row[indices[3]])*1e-6), y = label_y,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = NULL);
    	}
	else{next;}
	}
else {return;}
     }

}


yank_loc = '/home/jkb4y/ubs/work/results/July2012/Logistic_pc8/RegionYank/AA_8_yank_postlude.tbl'
downViewport('pvalsClipped');

DrawRefLine(0.05)
dataFrame <- MakeData(yank_loc)
#LabelSNPs(dataFrame)
