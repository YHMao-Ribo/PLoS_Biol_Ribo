library(poolr)

Args<-commandArgs()
total_gene<-as.numeric(Args[6])

result=list()
p_value=list()

count=1
myfile<-file("./potential_uTIS.gz","r")
while(count<total_gene){
	mylines=readLines(myfile,n=1)
	
	dataLine<-strsplit(mylines,split="\t")
	data_len<-length(dataLine[[1]])
	
	data<-data.frame('pos'=seq(0,data_len-4),'value'=as.numeric(dataLine[[1]][4:data_len]))
	
	data1<-subset(data,pos%%3==0)
	data2<-subset(data,pos%%3==1)
	data3<-subset(data,pos%%3==2)
	
	data2$mean_outframe<-(data2$value+data3$value)/2
	pval1<-wilcox.test(data1$value,data2$value,alternative='g')
	pval2<-wilcox.test(data1$value,data3$value,alternative='g')
	pcomb<-stouffer(c(pval1$p.value,pval2$p.value))
	#pval<-wilcox.test(data1$value,data2$mean_outframe,alternative='g')
	
	inframe<-sum(data1$value)/sum(data$value)
	output<-paste(dataLine[[1]][1],dataLine[[1]][2],dataLine[[1]][3],inframe,sep="\t")
	
	result[[count]]<-output
	p_value[[count]]<-pcomb$p
	#p_value[[count]]<-pval$p.value
	count<-count+1
}
close(myfile)

#p_adjust<-p.adjust(p_value,method='fdr')
p_adjust<-p_value
output<-paste(result,p_adjust,sep="\t")
write(output,file="potential_uTIS_pvalue")



