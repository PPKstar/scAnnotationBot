#' @title Annotates cell clusters in single-cell sequencing analysis.

#'

#' @description This package annotates cell clusters in single-cell sequencing analysis by using ERNIEBot or DashScope models.

#' @param account list. Obtained from the website https://console.bce.baidu.com/qianfan/ais/console/applicationConsole/application?_=1720524456151. client_id is the API Key, and client_secret is the Secret Key.
#' @param model character. Check from website https://console.bce.baidu.com/qianfan/ais/console/onlineService. the default value is "completions_pro".
#' @param inputTable data.frame returned by FindAllMarkers.
#' @param tissueName The tissue origin of single-cell sequencing data, such as "mouse brain".
#' @param useGeneNumber numeric. The number of genes used for annotating cell types, with a default value of 10.

#' @return character vector

#' @export ERNIEBotCellType
#' @name ERNIEBotCellType

library(httr)
library(jsonlite)
getAccessToken_ERNIEBot = function(client_id, client_secret) {
  url = 'https://aip.baidubce.com/oauth/2.0/token'
  headers = c('Content-Type' = 'application/json')
  body = list(grant_type = 'client_credentials',
              client_id = client_id,
              client_secret = client_secret)
  response = POST(url, add_headers(headers), body = body)
  access_token = content(response)$access_token
  return(access_token)
}

ERNIEBotChat = function(access_token, contentText, model = "completions_pro") {
  url = paste0('https://aip.baidubce.com/rpc/2.0/ai_custom/v1/wenxinworkshop/chat/',model)
  full_url = paste0(url, "?access_token=", access_token)
  headers = c('Content-Type' = 'application/json')
  data = list(
    messages = list(
      list(role = "user", content = contentText)
    ),
    system = 'You are a tool used for annotating cell types in single-cell sequencing data, strictly executing annotation tasks according to the requirements!'
  )
  json_data = toJSON(data, auto_unbox = TRUE)
  response = POST(full_url, add_headers(headers), body = json_data)
  response = content(response)
  prompt_tokens = response$usage$prompt_tokens
  completion_tokens = response$usage$completion_tokens
  total_tokens = response$usage$total_tokens
  tokenMessage = paste0('Tokens usage report: ',
                        '\n\t-prompt tokens: ',prompt_tokens,
                        '\n\t-completion tokens: ',completion_tokens,
                        '\n\t-total tokens: ',total_tokens, '\n')
  cat(tokenMessage)
  return(response$result)
}

ERNIEBotCellType = function(account = list(client_id = NULL, client_secret = NULL), model = "completions_pro", inputTable, tissueName = NULL, useGeneNumber = 5) {
  if(is.null(account$client_id) || is.null(account$client_secret)) {
    print('client_id or client_secret is not available!')
    return()
  }
  access_token = getAccessToken_ERNIEBot(client_id = account$client_id, client_secret = account$client_secret)
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
    response = ERNIEBotChat(access_token = access_token, content = message, model = model)
    responseList = append(responseList, response)
  }
  combinedRes = paste(unlist(responseList), collapse = "\n")
  clearedRes = gsub(pattern = "\n", replacement = "", x = combinedRes)
  clearedRes = strsplit(clearedRes, "\\d+:")[[1]]
  clearedRes = clearedRes[-1]
  cellTypes = setNames(clearedRes, levels(factor(inputTable$cluster)))
  return(cellTypes)
}
