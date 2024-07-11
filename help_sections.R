#UI output for Help Sections

#HTML document for main panel
output$dea_md <- renderUI({
    includeHTML("markdown/dea_algorithms.html")
})
