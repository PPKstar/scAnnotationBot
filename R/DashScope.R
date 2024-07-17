#' @title Annotates cell clusters in single-cell sequencing analysis.

#'

#' @description This package annotates cell clusters in single-cell sequencing analysis by using ERNIEBot or DashScope models.

#' @param api_key character. Obtained from the website https://dashscope.console.aliyun.com/apiKey.
#' @param model character. Check from website https://dashscope.console.aliyun.com/billing, the default value is "qwen-max-longcontext".
#' @param inputTable data.frame returned by FindAllMarkers.
#' @param tissueName The tissue origin of single-cell sequencing data, such as "mouse brain".
#' @param useGeneNumber numeric. The number of genes used for annotating cell types, with a default value of 10.

#' @return character vector

#' @export DashScopeCellType
#' @name DashScopeCellType

library(httr)
library(jsonlite)
DashScopeChat = function(api_key, contentText, model = "qwen-max-longcontext") {
  url = 'https://dashscope.aliyuncs.com/api/v1/services/aigc/text-generation/generation'
  headers = add_headers("Content-Type" = "application/json",
                        "Authorization" = paste("Bearer", api_key))
  json_input = list(
    model = model,
    input = list(
      messages = list(
        list(role = "system", content = "You are a tool used for annotating cell types in single-cell sequencing data, strictly executing annotation tasks according to the requirements!"),
        list(role = "user", content = contentText)
      )
    ),
    parameters = list(result_format = "message")
  )
  response = POST(url, headers, body = json_input, encode = "json")
  response = content(response)
  prompt_tokens = response$usage$input_tokens
  completion_tokens = response$usage$output_tokens
  total_tokens = response$usage$total_tokens
  tokenMessage = paste0('Tokens usage report: ',
                        '\n\t-prompt tokens: ',prompt_tokens,
                        '\n\t-completion tokens: ',completion_tokens,
                        '\n\t-total tokens: ',total_tokens,'\n')
  cat(tokenMessage)
  return(response$output$choices[[1]]$message$content)
}

DashScopeCellType = function(api_key = NULL, model = "qwen-max-longcontext", inputTable, tissueName = NULL, useGeneNumber = 5) {
  if(is.null(api_key)) {
    print('API key is not available!')
    return()
  }
  clusterNum = length(levels(factor(inputTable$cluster)))
  minCluster = as.numeric(levels(factor(inputTable$cluster))[1])
  sectionSize = 30
  sliceNum = ceiling(clusterNum/sectionSize)
  responseList = list()
  for (i in 1:sliceNum) {
    cat(paste0('----------section ',i,' start----------\n'))
    section = inputTable[inputTable$cluster<=sectionSize*i & inputTable$cluster>=minCluster+sectionSize*(i-1),]
    markedClusters = tapply(section$gene, list(section$cluster), function(i) paste0(i[1:useGeneNumber], collapse = ","))
    message = paste0("Use the following markers to identify the cell types of ", tissueName, " tissue cells. Provide only the cell type names in English, formatted as number followed by a colon and the cell type name in English, with no additional description after the name and the end of the response.",
                     "\n", paste0(names(markedClusters), ":", unlist(markedClusters), collapse = "\n"))
    response = DashScopeChat(api_key = api_key, content = message, model = model)
    responseList = append(responseList, response)
  }
  combinedRes = paste(unlist(responseList), collapse = "\n")
  clearedRes = gsub(pattern = "\n", replacement = "", x = combinedRes)
  clearedRes = strsplit(clearedRes, "\\d+:")[[1]]
  clearedRes = clearedRes[-1]
  cellTypes = setNames(clearedRes, levels(factor(inputTable$cluster)))
  return(cellTypes)
}
