
library(shiny)
library(shinythemes)
library(dplyr) 
library(DT)
library(png)
library(shinyWidgets)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  #theme = shinytheme("united"),
  #shinythemes::themeSelector(),
  #theme = shinytheme("slate"),
  # shinythemes::themeSelector(),
  # setBackgroundImage(src = "www/hope.png") ,
  navbarPage(
    #themeSelector(),
    "Home",
    theme = shinytheme("united"),
    # theme = shinytheme("cerulean"),
    tabPanel('Load Data',
             sidebarLayout(

                sidebarPanel(
                       h3(strong("Uploading file:"))
                       ,fileInput("upload.count.df", "Choose table File(csv)",multiple = FALSE)
                       ,fileInput("upload.meta.df", "Choose meta File(csv)",multiple = FALSE )    
                       ,h6( 'Note:The first column about SampleID should be set for the table. ')
                       ,selectInput("SelectP1_main_group", "Select column of meta", choices  = c("Group", "SampleID","Times.at.enrollment") ,selected = "Group")
                       ,selectizeInput("SelectP1_covar", "Select columns about covariance  in meta",c("Item A", "Item B"),multiple = TRUE,options = list(maxItems = 2L) )
                       ,selectizeInput("SelectP1_main_colname_contain", "Select content of column in meta",c("Item A", "Item B"),multiple = TRUE,options = list(maxItems = 2L) )
                       ,h5(strong("Note: The order of groups need be set."))
                       # ,sliderInput("prevalence_values", "Prevalence values",min = 0.0, max = 1, value = 0.25)
                       # ,actionButton("SelectP1_run1_uploading", "Submit")
                       ,br() 
                       ,actionButton("run_cleandata", "run_cleandata")
                       ,br() 
                       ,h6( 'You can download a example here. ')
                       ,tabPanel("downloadexample",downloadButton('downloadexample', 'Download the example'))
                       ,h6( 'You can watch the video,and it tell u how to use! ')
                       ,tabPanel("downloadvideo",downloadButton('downloadvideo', 'Download the video'))
                       ),
              
                mainPanel(
                      titlePanel("Analysis Platform-GV")

                      ,tabsetPanel(type = "tabs"
                                   ,tabPanel("main", img(src = "Analysis_Platform_v2.png", height = 300, width = 800) )  
                                   # ,tabPanel("main", img(src = "hope.png", height = 500, width = 800)  )
                                   ,tabPanel("meta", DT::dataTableOutput("meta_table")  )
                                   ,tabPanel("EDA"  
                                             ,selectInput("pca_group", "Select pca group", choices  =c("Group") ,selected = "Group")
                                             ,selectInput("EDA_picture_selected", "Select one picture", choices  =c("vis","vis_miss","screeplot","pca") ,selected = "vis")
                                             ,sliderInput("eda_width", "eda_width:", min = 0, max = 10000, value = 750 )
                                             ,sliderInput("eda_height", "eda_height:", min = 0, max = 10000, value = 700 )
                                             # ,imageOutput("EDA")
                                             ,h5("Note:better size")
                                             ,h5("vis.height:1300;vis.width:3000")
                                             ,h5("vis_miss.height:1300;vis_miss.width:3000")
                                             ,h5("pca_height:600;pca_width:600")
                                             ,plotOutput("EDA_picture_selected")
                                              )  
                                   ,tabPanel("lm"
                                             ,tableOutput("lm_table")
                                             ,selectInput("lm_group", "Select group(var)", choices  =c("Group") ,selected = "Group")
                                             ,selectInput("Select_lm_x", "Select x(var) of lm", choices  = c("x1", "x2","x3") ,selected = "x1")
                                             ,selectInput("Select_lm_y", "Select y(var) of lm", choices  =  c("y1", "y2","y3")  ,selected = "x2")
                                             ,selectInput("lm_picture_selected", "Select one picture", choices  =c("ggplot","ggscatterstats") ,selected = "ggplot")
                                             ,sliderInput("fig_lm_fit_width", "fig_lm_fit_width:", min = 0, max = 5000, value = 1300 )
                                             ,sliderInput("fig_lm_fit_height", "fig_lm_fit_height:", min = 0, max = 5000, value = 700 )

                                             ,h5("Note:better size")
                                             ,h5("ggscatterstats.height:3000;ggscatterstats.width:1300")
                                                                   
                                             ,imageOutput("AD_plot_lm_fit")
                                             # ,imageOutput("plot_rare_fig",width = 100,height = 100) 
                                             ) 
                                   ,tabPanel("heatmap" 
                                             ,sliderInput("heatmap_width", "heatmap_width:", min = 0, max = 3000, value = 900 )
                                             ,sliderInput("heatmap_height", "heatmap_height:", min = 0, max = 1000, value = 500 )
                                             ,imageOutput("AD_plot_heatmap",width = 800,height = '800')
                                             ) 
                                   ,tabPanel("boxplot"
                                             ,h6( 'Note:You also can download the figures by button from More(download) .')
                                             ,plotOutput("AD_plot_box",width = 600,height = '600'))

                                   ,tabPanel("hcpc"
                                             ,sliderInput("fig_hcpc_map_width", "fig_hcpc_map_width:", min = 0, max = 5000, value = 1200 )
                                             ,sliderInput("fig_hcpc_map_height", "fig_hcpc_map_height:", min = 0, max = 5000, value = 1000 )
                                             ,imageOutput("AD_plot_hcpc")
                                             
                                   ) 
                                   ,tabPanel("oplsda"
                                             ,imageOutput("AD_plot_oplsda_model")
                                             ,br()
                                             ,br()
                                             ,br()
                                             ,br()
                                             ,br()
                                             ,br()
                                             ,sliderInput("fig_oplsda_width", "fig_oplsda_width:", min = 0, max = 2000, value = 600 )
                                             ,sliderInput("fig_oplsda_height", "fig_oplsda_height:", min = 0, max = 2000, value = 400 )
                                             ,imageOutput("AD_plot_oplsda")
                                             # ,imageOutput("AD_plot_oplsda_model")
                                   )


                                     
                      
                                  ,tabPanel("ggbetweenstats"
                                            ,br()
                                            ,br()
                                            ,selectInput("gts.plot.type", "Select gts.plot.type", choices  = c("boxviolin", "box","violin") ,selected = "boxviolin")
                                            ,selectInput("gts.type", "Select gts.type", choices  = c("noparametric", "parametric","robust","bayes") ,selected = "noparametric")
                                            # ,selectInput("gts.data.type", "Select raw or log data", choices  = c("rawdata", "log") ,selected = "rawdata")
                                            # ,selectInput("gts.centrality.plotting", "Select TRUE or FALSE in pairwise.comparisons", choices  = c("TRUE", "FALSE") ,selected = "TRUE")
                                            ,selectInput("pairwise.comparisons", "Select TRUE or FALSE in centrality.plotting", choices  = c("TRUE", "FALSE")  ,selected = "TRUE")
                                            ,selectInput("results.subtitle", "Select TRUE or FALSE in subtitle", choices  = c("TRUE", "FALSE")  ,selected = "TRUE")
                                            ,selectInput("asterisk_label", "Select TRUE or FALSE in asterisk_label", choices  = c("TRUE", "FALSE")  ,selected = "TRUE")
                                            # ,sliderInput("gts.type", "fig_oplsda_width:", min = 0, max = 2000, value = 600 )
                                            # ,sliderInput("fig_oplsda_height", "fig_oplsda_height:", min = 0, max = 2000, value = 400 )
                                            ,h5("If you have some questions,you can refer https://github.com/IndrajeetPatil/ggstatsplot")
                                            ,imageOutput("gts_box",width = 300,height = '300')
                                  )
                      )

               ) ) )
    

    ,navbarMenu("More",
                tabPanel("Refresh",actionButton('switchtab',"Click this") )
                ,tabPanel("downloadData",downloadButton('downloadData', 'Download the results'))
                ,tabPanel("R",verbatimTextOutput("sessionInfo")  )      
    )        



))







