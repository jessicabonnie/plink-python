
DrawRefLine <- function(x){
print('\n');
num <- -log10(3.5e-7)
downViewport('pvalsClipped');
panel.abline(h=transformation(num),col="red3");
}

MakeData <- function(yank_loc){

x <- read.delim(yank_loc)
dataFrame <- as.data.frame(x)
print(dataFrame);
return(dataFrame);

}

LabelRefSNP <- function(ref){
specialref <- c('rs2476601','rs3087243','rs9273363','rs61839660',
				'rs12416116','rs3842727','rs1701704','rs653178',
				'rs12927355','rs3826110')
reflabel = ref;
if (ref %in% specialref){
	reflabel <- paste(ref,"*",sep='')
	}
return(reflabel);
}

ColorSNP <- function(pos, ypos, shape){
xpos = pos * 1e-6;

grid.points(x=xpos,y=ypos, 
				gp=gpar(col=args[['refsnpColor']],fill=args[['refsnpColor']],
				cex= if (args[['bigDiamond']] & args[['showRefsnpAnnot']]) 1.6*args[['refDot']] else args[['refDot']]),
				pch=21,
				default.units='native'
				);

}


LabelSNPs <- function(dataFrame){
#mycols <- c("chr","snp_name","snp_position","p.value",'pch')
mycols <- c("chr","snp_name","snp_position","p.value")
indices <- match(mycols, colnames(dataFrame))
print(indices);
text_col = "darkblue";

if (as.numeric(args[['chr']])== 8){
	panel.text(x=56554275*1e-6, y = -log10(6.11e-6)+1.7,
	labels = 'rs7839768', adj = NULL, col=text_col,vfont = NULL,
    font = 2);
    return;
    }
for (i in 1:nrow(dataFrame)){
row <- dataFrame[i,]
print(args[['chr']]);
if(length(row)!=0) {
	gold_chr = as.numeric(row[indices[1]]);
	lz_chr = as.numeric(args[['chr']])
	if (gold_chr == lz_chr){
	print(row[indices[2]]);
	
	new_p <- -log10(as.numeric(row[indices[4]]))
	reflabel <- LabelRefSNP(row[[indices[2]]])
		if ( new_p > 16.2){
		new_p = 16.1
		}
		if (lz_chr== 9){
		panel.text(x=(as.numeric(row[indices[3]])*1e-6+5), y = new_p+1.5,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else if (lz_chr== 20){
		panel.text(x=(as.numeric(row[indices[3]])*1e-6)+3, y = new_p+1.5,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else if (lz_chr== 11){
		panel.text(x=10, y = new_p+1.5,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else if (lz_chr == 1){
		xloc = as.numeric(row[indices[3]])*1e-6;
		if (row[[indices[2]]] =='rs3024493'){
			xloc = xloc + 20;
			}
		panel.text(x=xloc, y = new_p+1,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else if (lz_chr== 23){
		panel.text(x=145, y = new_p+1.5,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else if (row[[indices[2]]] == 'rs61839660'){
		panel.text(x=(as.numeric(row[indices[3]])*1e-6)+5, y = new_p+1.5,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else if (row[[indices[2]]] == 'rs12922409'){
		panel.text(x=(as.numeric(row[indices[3]])*1e-6)+.5, y = new_p+1.3,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
		}
		else{
		panel.text(x=(as.numeric(row[indices[3]])*1e-6), y = new_p+1.5,
		labels = reflabel, adj = NULL, col=text_col,vfont = NULL,
    	font = 2);
    	}
    }
	else{next;}
	}
else {return;}
     }

}

HideGenes<- function(x){
upViewport(4);
#upViewport('refFlatInner');
downViewport("refFlatInner");
panel.rect(0, 0, 100, 100, density = NULL, angle = 45,
     col = "red", border = NULL, lty = NULL, lwd = par("lwd"),
     xpd = NULL);
}

yank_loc = '/home/jkb4y/work/results/eurmeta/hg19/RegionYank_05242012/eurmeta_yank_chrpostlude.tbl'

DrawRefLine()
dataFrame <- MakeData(yank_loc)
LabelSNPs(dataFrame)

#HideGenes()
#DrawAxis()

