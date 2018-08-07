#
# ClueGO cyREST KEGG Batch Workflow
#
library(RJSONIO)
library(httr)

text.to.data.frame <- function(table.text) {
  table <- NULL
  rows <- unlist(strsplit(result.table.text, split="\n"))
  header <- t(unlist(strsplit(rows[1], split="\t")))
  for(i in 2:length(rows)) {
    if(is.null(table)) {
      table <- t(unlist(strsplit(rows[i], split="\t")))
    } else {
      table <- rbind(table,t(unlist(strsplit(rows[i], split="\t"))))
    }
  }
  table <- as.data.frame(table)
  names(table) <- header
  return(table)
}



if(TRUE){
      #Basic settings for cyREST
      home.folder = "F:/min-labs_paper/work/XEJ2/XEJ2/重排znt1高低表达的样本/kegg" # Linux / MacOSX
      cluego.home.folder = paste(home.folder,"",sep="") # Linux / MacOSX
      port.number = 1234
      host.address <- "localhost" # "192.168.0.20" #
      
      cytoscape.base.url = paste("http://",host.address,":", toString(port.number), "/v1", sep="")
      cluego.base.url = paste(cytoscape.base.url,"apps/cluego/cluego-manager", sep="/")
      
      # list all available organisms
      response <- GET(paste(cluego.base.url,"organisms","get-all-installed-organisms",sep="/"))
      print(content(response))
      
      
      #开启cluego
      
      response<-POST(url = "http://localhost:1234/v1/apps/cluego/start-up-cluego")
      response<-PUT(url = "http://localhost:1234/v1/apps/cluego/cluego-manager/ontologies/TRUE/1.0")
      print(response)
      # set organism to 'Homo Sapiens'
      organism.name = "Homo Sapiens"
      # organism.name = "Mus Musculus"
      response <- PUT(url=paste(cluego.base.url,"organisms","set-organism",URLencode(organism.name),sep="/"), encode = "json")
      print(response)
      #set ontology to Kegg
      # list all available ontologies
      response <- GET(paste(cluego.base.url,"ontologies","get-ontology-info",sep="/"))
      print(content(response))
      
      # set ontologies
      selected.ontologies <- toJSON(c("8;Rectangle"))
      response <- PUT(url=paste(cluego.base.url,"ontologies","set-ontologies",sep="/"), body=selected.ontologies, encode = "json", content_type_json())
}



######开始读取差异表达CSV文件并且进行KEGG富集分析
if(TRUE){
      # set gene list for cluster 1
      dir_infile="F:/min-labs_paper/work/XEJ2/XEJ2/重排znt1高低表达的样本/output"
      (files_csv=list.files(path = dir_infile,pattern = "*.csv",full.names = T,recursive = T))
      setwd(dir_infile)
      data_thecsv=data.frame()
      for(onecsv in files_csv){
        currentrowindex=nrow(data_thecsv)+1
        data_thecsv[currentrowindex,"full_filename"]=onecsv
        currentfile_dir=stringi::stri_extract(str = onecsv,regex =".*\\/(?=[A-Za-z0-9_\\.]+$)")
        data_thecsv[currentrowindex,"dir"]=currentfile_dir
        filename=stringi::stri_extract(str = onecsv,regex ="(?<=\\/)[A-Za-z0-9_\\.]+$")
        data_thecsv[currentrowindex,"filename"]=filename
        message(onecsv)
      }
      
      
      ##开始进行KEGG富集分析
      apply(X =data_thecsv,MARGIN = 1,FUN = function(x){
        all_out_file=paste(x[["dir"]],"Kegg_",x[["filename"]],sep = "")
        message(all_out_file)
        data_thecsv_filtered=read.csv(file = x[["full_filename"]],stringsAsFactors = FALSE)
        ##存在则跳过
        if(file.exists(all_out_file)){
          return(paste("存在",all_out_file))
        }
        ###分列
        for(rowindex in c(1:nrow(data_thecsv_filtered))){
          theGeneGeneName=data_thecsv_filtered[rowindex,"X"]
          theGeneGeneNameS=stringi::stri_split(str = theGeneGeneName,regex = "_",simplify = TRUE)
          data_thecsv_filtered[rowindex,"GeneName"]=theGeneGeneNameS[1]
          data_thecsv_filtered[rowindex,"EntrezId"]=theGeneGeneNameS[2]
        }
        
        ##筛选
        selectedrowindex_pvalue=which(as.numeric(data_thecsv_filtered[,"PValue"])<=0.05)
        selectedrowindex_logFC=which(abs(as.numeric(data_thecsv_filtered[,"logFC"]))>1)
        selectedrowindex=intersect(selectedrowindex_pvalue,selectedrowindex_logFC)
        data_thecsv_filtered=data_thecsv_filtered[selectedrowindex,]
        ###################
        cluster1 = "1"
        gene.list <- toJSON(c(data_thecsv_filtered$EntrezId))
        response <- PUT(url=paste(cluego.base.url,"cluster","upload-ids-list",URLencode(cluster1),sep="/"), body=gene.list, encode = "json", content_type_json())
        print(response) 
        # run the analysis
        analysis.name <- onecsv
        response <- GET(paste(cluego.base.url,URLencode(analysis.name),sep="/"))
        print(content(response, encode = "text",encoding = "UTF-8"))
        
        # get network id
        current.network.suid <- as.numeric(stringi::stri_extract(str = content(response, encode = "text",encoding = "UTF-8"),regex = "(?<=SUID: )\\d+"))
        print(current.network.suid)
        
        # get analysis results
        response <- GET(paste(cluego.base.url,"analysis-results","get-cluego-table",current.network.suid,sep="/"))
        result.table.text <- content(response, encode = "text", encoding = "UTF-8")
        result.table.data.frame <- text.to.data.frame(result.table.text)
        table.file.name = paste(,paste("Kegg",onecsv,sep = "_"),sep="/") # Linux / MacOSX
        # table.file.name = paste(home.folder,"ClueGOExampleResultTable.txt",sep="\\") # Windows
        write.table(result.table.data.frame,file=table.file.name,row.names=FALSE, na="",col.names=TRUE, sep=",")
        # print(result.table.data.frame)
        ######delete_all_the_network
        # response <- GET(url =paste(cytoscape.base.url,"/networks",sep = ""))
        response<-DELETE(url =paste(cytoscape.base.url,"/session",sep = ""))
        content(response, encode = "text",encoding = "UTF-8")
        
      })
}










