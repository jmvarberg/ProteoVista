#UI output for Help Sections

#HTML document for DEA help
output$dea_md <- renderUI({
    includeHTML("markdown/dea_algorithms.html")
})

#HTML document for FragPipe help
output$fragpipe_md <- renderUI({
    includeHTML("markdown/proteovista_dda_help.html")
})

