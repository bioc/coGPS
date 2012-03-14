PCOPA<-function(exprslist,alpha,side,type)
{
	##################################################################################
	Emperical_P_up<-function(baseline,tumor)
	{
		result<-sapply(1:length(tumor),function(i) length(which(baseline>tumor[i])))
		result<-result/length(baseline)
	}
	Emperical_P_down<-function(baseline,tumor)
	{
		result<-sapply(1:length(tumor),function(i) length(which(baseline<tumor[i])))
		result<-result/length(baseline)
	}
	Emperical_P_twosided<-function(baseline,tumor)
	{
		result<-sapply(1:length(tumor),function(i) 
			length(which(baseline<tumor[i]))+length(which(baseline<tumor[i])))
		result<-result/length(baseline)
	}
	#################################################################################
	NG<-nrow(exprslist[[1]]$exprs)
	genename<-rownames(exprslist[[1]]$exprs)
	Pmatrices<-list()
	for(j in 1:length(exprslist))
	{
		if(side[j]=="up")
		{
			Pmatrices[[j]]<-sapply(1:NG,function(i) 
				Emperical_P_up(exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==0)],
				exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==1)]))
		}
		if(side[j]=="down")
		{
			Pmatrices[[j]]<-sapply(1:NG,function(i) 
				Emperical_P_down(exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==0)],
				exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==1)]))
		}
		if(side[j]=="twosided")
		{
			Pmatrices[[j]]<-sapply(1:NG,function(i) 
				Emperical_P_twosided(exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==0)],
				exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==1)]))
		}
		Pmatrices[[j]]<-t(Pmatrices[[j]])
	}
	#############################################################################
	#pvalues to 0-1, summed up
	##############################################################################
	Cutoff_p<-function(vector,cut)
	{
		result<-length(which(vector<cut))
	}
	if(type=="uniform")
	{
		int.Pmatrices<-Pmatrices[[1]]
		if(length(exprslist)>=2)
		{
			for(j in 2:length(exprslist))
			{
				int.Pmatrices<-cbind(int.Pmatrices,Pmatrices[[j]])
			}
		}
		totaltummor<-sum(exprslist[[1]]$classlab)
		if(length(exprslist)>=2)
		{
			for(j in 2:length(exprslist))
			{
				totaltummor<-totaltummor+sum(exprslist[[j]]$classlab)
			}
		}
		result<-sapply(1:NG,function(i) Cutoff_p(int.Pmatrices[i,],alpha/totaltummor))
	}
	
	if(type=="subtype")
	{
		score.matrices<-list()
		for(j in 1:length(exprslist))
		{
			score.matrices[[j]]<-sapply(1:NG,function(i) 
							Cutoff_p(Pmatrices[[j]][i,],alpha/(sum(exprslist[[j]]$classlab))))

		}
		int.Pmatrices<-score.matrices[[1]]
		if(length(exprslist)>=2)
		{
			for(j in 2:length(exprslist))
			{
				int.Pmatrices<-cbind(int.Pmatrices,score.matrices[[j]])
			}
		}
		result<-int.Pmatrices
		if(length(exprslist)>=2)
		{
			result<-sapply(1:NG,function(i) max(int.Pmatrices[i,]))
		}
	}
	names(result)<-genename
	result
}
PatientSpecificGeneList<-function(exprslist,alpha,side,type,TopGeneNum)
{
	##################################################################################
	Emperical_P_up<-function(baseline,tumor)
	{
		result<-sapply(1:length(tumor),function(i) length(which(baseline>tumor[i])))
		result<-result/length(baseline)
	}
	Emperical_P_down<-function(baseline,tumor)
	{
		result<-sapply(1:length(tumor),function(i) length(which(baseline<tumor[i])))
		result<-result/length(baseline)
	}
	Emperical_P_twosided<-function(baseline,tumor)
	{
		result<-sapply(1:length(tumor),function(i) 
			length(which(baseline<tumor[i]))+length(which(baseline<tumor[i])))
		result<-result/length(baseline)
	}
	#################################################################################
	NG<-nrow(exprslist[[1]]$exprs)
	genename<-rownames(exprslist[[1]]$exprs)
	Pmatrices<-list()
	for(j in 1:length(exprslist))
	{
		if(side[j]=="up")
		{
			Pmatrices[[j]]<-sapply(1:NG,function(i) 
				Emperical_P_up(exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==0)],
				exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==1)]))
		}
		if(side[j]=="down")
		{
			Pmatrices[[j]]<-sapply(1:NG,function(i) 
				Emperical_P_down(exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==0)],
				exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==1)]))
		}
		if(side[j]=="twosided")
		{
			Pmatrices[[j]]<-sapply(1:NG,function(i) 
				Emperical_P_twosided(exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==0)],
				exprslist[[j]]$exprs[i,which(exprslist[[j]]$classlab==1)]))
		}
		Pmatrices[[j]]<-t(Pmatrices[[j]])
	}
	#############################################################################
	#pvalues to 0-1, summed up
	##############################################################################
	Cutoff_p<-function(vector,cut)
	{
		result<-length(which(vector<cut))
	}
	if(type=="uniform")
	{
		int.Pmatrices<-Pmatrices[[1]]
		if(length(exprslist)>=2)
		{
			for(j in 2:length(exprslist))
			{
				int.Pmatrices<-cbind(int.Pmatrices,Pmatrices[[j]])
			}
		}
		totaltummor<-sum(exprslist[[1]]$classlab)
		if(length(exprslist)>=2)
		{
			for(j in 2:length(exprslist))
			{
				totaltummor<-totaltummor+sum(exprslist[[j]]$classlab)
			}
		}
		result<-sapply(1:NG,function(i) Cutoff_p(int.Pmatrices[i,],alpha/totaltummor))
 		######################################################################################
		#outlier gene list for each patient
		######################################################################################
		outliergene_bypatient<-list()
		for(k in 1:ncol(Pmatrices[[1]]))
		{
			outliergene_bypatient[[k]]<-list()
			for(j in 1:length(exprslist))
			{
				#outliergene_bypatient[[k]][[j]]<-list()
				outliergene_bypatient[[k]][[j]]<-genename[which(Pmatrices[[j]][,k]<alpha/totaltummor & order(result)<=TopGeneNum)]		
			}
		}
	}
	
	if(type=="subtype")
	{
		score.matrices<-list()
		for(j in 1:length(exprslist))
		{
			score.matrices[[j]]<-sapply(1:NG,function(i) 
							Cutoff_p(Pmatrices[[j]][i,],alpha/(sum(exprslist[[j]]$classlab))))

		}
		int.Pmatrices<-score.matrices[[1]]
		if(length(exprslist)>=2)
		{
			for(j in 2:length(exprslist))
			{
				int.Pmatrices<-cbind(int.Pmatrices,score.matrices[[j]])
			}
		}
		result<-sapply(1:NG,function(i) max(int.Pmatrices[i,]))
		outlier_datatype<-sapply(1:NG,function(i) which.max(int.Pmatrices[i,]))
		for(i in 1:NG)
		{
			for(j in 1:length(exprslist))
			{
				if(j!=outlier_datatype[i])
				{
					Pmatrices[[j]][i,]<-1
				}		
			}
		}	
		outliergene_bypatient<-list()
		for(k in 1:ncol(Pmatrices[[1]]))
		{
			outliergene_bypatient[[k]]<-list()
			for(j in 1:length(exprslist))
			{
				#outliergene_bypatient[[k]][[j]]<-list()
				outliergene_bypatient[[k]][[j]]<-genename[which(Pmatrices[[j]][,k]<alpha/(sum(exprslist[[j]]$classlab))  & order(result)<=TopGeneNum)]			
			}
		}

	}
	names(outliergene_bypatient)<-paste("Patient",
						colnames(exprslist[[1]]$exprs)[which(exprslist[[1]]$classlab==1)])
	outliergene_bypatient
}

plotCOPA<- function(expData, index,labt,namelist,type) {
      plotData <- expData[index,]
      tumors <- plotData[which(labt==1)]
      normals <- plotData[which(labt==0)]
	  colorlab<-c(rep("blue",length(which(labt==0))),
			rep("red",length(which(labt==1))))

      orderT <- order(tumors)
      tumors <- tumors[orderT]
      orderN <- order(normals)
      normals <- normals[orderN]

      plotData <- c(normals,tumors)
	  temp<-plotData
      barplot(temp,las=2,cex.names=0.5,main=paste(namelist[index],type),col=colorlab)
}

PlotTopPCOPA<-function(exprslist,PCOPAresult,topcut,typelist)
{
	layout(matrix(c(1:4),byrow=TRUE,ncol=2))
	toplist<-order(PCOPAresult,decreasing=TRUE)
	for(i in 1:topcut)
	{
		for(j in 1:length(exprslist))
		{
			plotCOPA(exprslist[[j]]$exprs,toplist[i],exprslist[[j]]$classlab,rownames(exprslist[[j]]$exprs),type=typelist[j])
		}
		plot(PCOPAresult[toplist[i]],main="score",ylim=range(PCOPAresult,na.rm=TRUE),ylab="score")
		if((length(exprslist)+1)%%4!=0)
		{
			for(k in 1:(4-length(exprslist)-1%%4))
			{
				plot.new()
			}
		}
	}
}

permCOPA<-function(exprslist, alpha=0.05, side, type, perms=100){

       nType <- length(exprslist)

       tempList <- list()
       classLab <- exprslist[[1]]$'classlab'
       nSample <- length(classLab)
       nGene <- nrow(exprslist[[1]]$'exprs')

       permResult <- matrix(nrow=nGene, ncol=perms)

       for (i in 1:perms) {
               temp <- sample(1:nSample,nSample)
               classLab2 <- classLab[temp]
               for (j in 1:nType) {
                       matExpr <- exprslist[[j]]$exprs
                       temp <- list(exprs=matExpr, classlab = classLab2)
                       tempList[[j]] <- temp
               }
               permResult[1:nGene,i] <- PCOPA(tempList, alpha, side, type)
       }
       matExpr <- exprslist[[1]]$exprs
       rownames(permResult) <- rownames(matExpr)
       return(permResult)
}
