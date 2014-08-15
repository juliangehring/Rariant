rariantInspect <- function(x) {
  
    ui = pageWithSidebar(
        headerPanel("Rariant Explorer"),
        sidebarPanel(
            #selectInput("chr", "Seqname:", seqlevels(globalVariantSet)),
            textInput("chr", "Seqname:", "1"), ## TODO
            numericInput("start", "Start position:", 1),
            numericInput("end", "End position:", 1e9),
            checkboxInput("show_somatic", "Show somatic events", TRUE),
            checkboxInput("show_hetero", "Show hetero events", TRUE),
            checkboxInput("show_undecided", "Show undecided events", TRUE),
            checkboxInput("show_powerless", "Show powerless events", TRUE),
            checkboxInput("select", "Show only calls", TRUE),
            numericInput("alpha", "Threshold on adjusted p-values", 0.1),
            submitButton("Call Variants"),
            actionButton("close_app", "Return to R")
            ),
        mainPanel(
            tabsetPanel(
                tabPanel("CI Plot", plotOutput("ci_plot")),
                tabPanel("Single Plot", plotOutput("single_plot")),
                tabPanel("Shift Plot", plotOutput("shift_plot")),
                tabPanel("Support Plot", plotOutput("support_plot")),
                tabPanel("Table", dataTableOutput("summary"))
                )
            )
        )
       
    server = function(input, output) {
        data = reactive({
            if(input$close_app > 0)
                stopApp()
            res = globalVariantSet
            idx_show = logical(length(res))
            if(input$show_somatic)
                idx_show[res$eventType %in% "somatic"] = TRUE
            if(input$show_hetero)
                idx_show[res$eventType %in% "hetero"] = TRUE
            if(input$show_undecided)
                idx_show[res$eventType %in% "undecided"] = TRUE
            if(input$show_powerless)
                idx_show[res$eventType %in% "powerless"] = TRUE
            res = res[idx_show]
            if(input$select)
                res = res[res$padj < input$alpha]
            roi = GRanges(input$chr, IRanges(input$start, input$end))
            res = subsetByOverlaps(res, roi)
            res
        })       
        output$ci_plot = renderPlot({
            print(plotConfidenceIntervals(data(), color = "eventType")) # + scale_color_hue(drop = FALSE) ## import
        })
        output$single_plot = renderPlot({
            d = data()
            p1 = print(plotConfidenceIntervals(d[d$eventType %in% "somatic", ])) + ylab("")
            p2 = print(plotConfidenceIntervals(d[d$eventType %in% "hetero", ])) + ylab("")
            p3 = print(plotConfidenceIntervals(d[d$eventType %in% "undecided", ])) + ylab("")
            p4 = print(plotConfidenceIntervals(d[d$eventType %in% "powerless", ])) + ylab("")
            t = ggbio::tracks(p1, p2, p3, p4)
            print(t)
        })
        output$shift_plot = renderPlot({
            print(plotAbundanceShift(data()))
        })
        output$support_plot = renderPlot({
            print(plotShiftSupport(data()))
        })
        output$summary = renderDataTable({
            as(data(), "data.frame")
        })
    }
    
    globalVariantSet = x
    runApp(list(ui = ui, server = server))

}
