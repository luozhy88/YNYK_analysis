library(shiny)
library(shinythemes)
library(plyr) 
library(dplyr) 
library(DT)
# library(metacoder)
library(purrr)
library(ggplot2)
library(phyloseq)
library(mixOmics)
library(metagMisc)
library(png)
# install.packages("scales")
library("scales")
library(future)
library(promises)
library(vegan)
library(optparse)
library(tibble)
library(microbiome)
library(parameters)

`%+%` <- function(a,b) {paste0(a,b)}

system("pwd")
system("rm -rf output/*")

# options(shiny.error = browser)
# options(shiny.error = recover)

message( paste0("###############################################################Start data:",Sys.time()) )
message( paste0("###############################################################Start dir:",getwd() ) )

server <- function(input, output, session) {
    rt <- reactiveValues()
    system("pwd")
    system("rm -rf output/*")
    
    observe({
    print("#######################################update meta observe#############################################")
        meta.df<-rt$meta.df 
       
       print("第一个参数："%+% input$SelectP1_main_group)
       filter_str_lie_name <- unique(meta.df[,input$SelectP1_main_group]) 
       print(class(filter_str_lie_name))
       print(paste0("ddd 第2个参数:",filter_str_lie_name ))
       updateSelectizeInput(session, "SelectP1_main_colname_contain",label = paste("Select two classficiations(example: F,M )", length( filter_str_lie_name )),choices = c("NULL",filter_str_lie_name ) ,options = list(maxItems = 1000L) ) 
       
       infer_var<- c(input$SelectP1_main_group,input$SelectP1_covar)
       updateSelectizeInput(session, "pca_group",label = paste("Select pca group", length( infer_var )),choices = c("All",infer_var ) ,options = list(maxItems = 1L) ) 
       
    })
    
    observe({
       print("#######################################update count observe#############################################")
       
       count.df <- rt$count.df 
       
       colname_<-colnames(count.df)
       updateSelectizeInput(session, "Select_lm_x",label = paste("Select x(var) of lm", length( colname_ )),choices = c(colname_), options = list(maxItems = 1L))
       updateSelectizeInput(session, "Select_lm_y",label = paste("Select y(var) of lm", length( colname_ )),choices = c(colname_), options = list(maxItems = 1L))
       
    })
    
    observeEvent(input$upload.meta.df,{
    print("########################update data observeEvent##########################################################")
            input.upload.meta.df.datapath<-input$upload.meta.df$datapath
            
            meta.df<-read.csv(input.upload.meta.df.datapath)
            colnames(meta.df) <- make.names(colnames(meta.df))
            meta.df <- column_to_rownames(meta.df, var = "SampleID")
            rt$meta.df <- meta.df

            colname_<-colnames(meta.df)
            updateSelectizeInput(session, "SelectP1_main_group",label = paste("Select one group(example: Group )", length( colname_ )),choices = c(colname_), options = list(maxItems = 1L))

            updateSelectizeInput(session, "SelectP1_covar",label = paste("Select two covariant(example: age,bmi )", length( colnames(meta.df) )),choices = c("NULL",colnames(meta.df) ) ,options = list(maxItems = 5L) ) 
            print(paste0("ddd 第3个参数:",input$SelectP1_covar ))
            updateSelectizeInput(session, "lm_group",label = paste("Select one group(example: Group)", length( colnames(meta.df) )),choices = c("NULL",colnames(meta.df) ) ,options = list(maxItems = 1L) ) 
          })
    observeEvent(input$upload.count.df,{
       print("########################update count observeEvent##########################################################")
       input.upload.count.df.datapath<-input$upload.count.df$datapath
       
       count.df<-read.csv(input.upload.count.df.datapath)
       colnames(count.df) <- make.names(colnames(count.df))
       count.df <- column_to_rownames(count.df, var = "SampleID")
       rt$count.df <- count.df
       
       colname_<-colnames(count.df)
       updateSelectizeInput(session, "Select_lm_x",label = paste("Select x(var) of lm", length( colname_ )),choices = c(colname_), options = list(maxItems = 1L))
       updateSelectizeInput(session, "Select_lm_y",label = paste("Select y(var) of lm", length( colname_ )),choices = c(colname_), options = list(maxItems = 1L))
       print(paste0("ddd 第3个参数:",input$SelectP1_covar ))
       
    })

# Clean data --------------------------------------------------------------
   print("################################################ dataframe clean ########################################")
   df_clean_list <- eventReactive(input$run_cleandata, {
            input.upload.count.df.datapath<-input$upload.count.df$datapath
            input.upload.meta.df.datapath<-input$upload.meta.df$datapath
            # browser()
            # input.upload.count.df.datapath<-"input/SHJW_cytokine_HC_DP_20210923_filter30_filterimpute_count.csv"
            # input.upload.meta.df.datapath<-"input/SHJW_cytokine_HC_DP_meta.csv"
            cytokine <- read.csv(input.upload.count.df.datapath)
            cytokine <- column_to_rownames(cytokine, var = "SampleID")
            colnames(cytokine) <- make.names(colnames(cytokine))
            
            meta.df<-read.csv(input.upload.meta.df.datapath)
            meta.df <- column_to_rownames(meta.df, var = "SampleID")
            colnames(meta.df) <- make.names(colnames(meta.df))
            
            rt$meta.df <- meta.df
            print("###################filter ######################")
            # browser()
            # selected_var<-c("Times.at.enrollment","Age","Sex") 
            selected_var<-c(input$SelectP1_main_group,input$SelectP1_covar) 
            main_colname<-input$SelectP1_main_group
            # main_colname_contain<- c("1", "2")
            
            main_colname_contain<- input$SelectP1_main_colname_contain
            meta.df <- dplyr::filter(meta.df, grepl(paste0(main_colname_contain,collapse = "|"),  meta.df[,main_colname] )) #筛选分组
            
            cytokine<-cytokine[rownames(meta.df),]
            
            # groups <- c("1", "2")
            # cytokine$Times.at.enrollment <- factor(cytokine$Times.at.enrollment, levels = groups)
            
            # df <- cytokine[, -c(1:4)]
            df_raw<-merge(meta.df[,selected_var] , cytokine ,by=0)
            df_raw<- column_to_rownames(df_raw,"Row.names")
            print("################### trans ######################")
            cytokine_tans <- cytokine %>% dplyr::mutate(across(where(is.numeric),function(x)log(x + 0.00001))) 
            df_trans<-merge(meta.df[,selected_var] , cytokine_tans ,by=0)
            df_trans<- column_to_rownames(df_trans,"Row.names")
            df_clean_list<-list()
            df_clean_list[["raw"]]<-df_raw
            df_clean_list[["trans"]]<-df_trans
            df_clean_list
          })
   output$meta_table <- DT::renderDataTable({  rt$meta.df  })
   
   


# lm ----------------------------------------------------------------------
   print("##################################################################lm ################################################")
   output$lm_table <- renderTable({
                    print("run..")
                    # browser()
                    withProgress(message='Program running:', detail="", {
                      incProgress(0.2, detail = paste("To become the richest man in the cemetery doesn't matter to me... 
                                                      Said the night to go to bed we have done a great thing... It is important for me."))
                      
                      df_clean_list<- df_clean_list()
                      df<-df_clean_list[["trans"]]
                      selected_var<-c(input$SelectP1_main_group,input$SelectP1_covar) 
                      main_colname<-input$SelectP1_mai
                      # selected_var<-c("Times.at.enrollment","Age","Sex") # the first of selected_var is main group.
                      # main_colname<-selected_var[1]
                      main_colname_contain<- input$SelectP1_main_colname_contain
                      df[,main_colname] <- factor(df[,main_colname])
                      ## linear regression (adjust sex age)
                      output.dir<-"output/" 
                      dir.create(  output.dir %+% paste0(main_colname_contain, collapse = "_vs_"),recursive = T  )
                      lm.fit.1 <- list()
                      ano.res <- list()
                      # browser()
                      for (i in colnames(df)[-c(1:length(selected_var))]   ) {
                        print(i)
                        # i<-"r5"
                        lm.fit.1[[i]] <- lm(as.formula( paste0(i, " ~ " %+% paste0(selected_var,collapse = "+"))  ) , data = df) 
                        ano.res[[i]] <- parameters::model_parameters(lm.fit.1[[i]])
                      }
                      # browser()
                      ano_dat <- do.call(rbind, ano.res)
                      
                      # Group
                      ano_df_list<-list()
                      for (var_name in selected_var){
                        print(var_name)
                        
                        ano_df <- ano_dat[grepl(var_name, ano_dat$Parameter),]
                        # ano_df <- ano_dat
                        ano_df <- add_column(ano_df, features = colnames(df)[-c(1:length(selected_var))] , .before = colnames(ano_df)[1] )
                        ano_df$p.adj <- p.adjust(ano_df$p, method = "fdr")
                        openxlsx::write.xlsx(ano_df, output.dir %+% paste0(paste0(main_colname_contain, collapse = "_vs_"),"/linear_regression_adjust_" %+% var_name %+%  ".xlsx"))
                        
                        
                        ano_df_list[[var_name]]<-ano_df
                      incProgress(0.9, detail = paste("To become the richest man in the cemetery doesn't matter to me... 
                                                      Said the night to go to bed we have done a great thing... It is important for me."))
                      message( paste0("#############################END lm##################################") )
                      
                      }
                      # browser()
                      rt$diff_table<-ano_df_list[[1]]
                      ano_df_list[[1]]
                        })           
                    })
   print("##################################################################lm fit################################################")

# lm fit ------------------------------------------------------------------
   output$AD_plot_lm_fit <- renderImage({         
      # 加上进度条显示
      withProgress(message='Program running:', detail="", {
         incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
         ######处理数据
         # browser()
         df_clean_list<- df_clean_list()
         data<-df_clean_list[["trans"]]
         var_x<-input$Select_lm_x
         var_y<-input$Select_lm_y
         print("######################################lm fit################################################")
         group_fit_line<-input$lm_group
         print("group_fit_line:" %+% group_fit_line)
         main_colname_contain<-data[,group_fit_line] %>% unique()
         data[,group_fit_line]<- factor(data[,group_fit_line])
         
         output.dir<-"output/" 
         dir.create(  output.dir %+% paste0(main_colname_contain, collapse = "_vs_"),recursive = T  )
         
         fig_lm_fit_height<- 700
         fig_lm_fit_width<-1300
         
         fig_lm_fit_height2<- 2000
         fig_lm_fit_width2<-2000
         
         if(input$lm_picture_selected =='ggplot'){
            fig_lm_fit_height<- input$fig_lm_fit_height
            fig_lm_fit_width<-input$fig_lm_fit_width
            print('neiceng:ggplot')
            print(fig_lm_fit_height)
            print(fig_lm_fit_width)
         }
         if(input$lm_picture_selected =='ggscatterstats'){
            fig_lm_fit_height2<- input$fig_lm_fit_height
            fig_lm_fit_width2<-input$fig_lm_fit_width
            print('neiceng:ggscatterstats')
            print(fig_lm_fit_height2)
            print(fig_lm_fit_width2)
         }
         #fig1
         lm_fit<-ggplot(data, aes_string(x = var_x , y = var_y, color = group_fit_line))  +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE)
         out.fig.name.fitline<-output.dir %+% paste0(main_colname_contain, collapse = "_vs_")  %+% '/' %+% var_x %+% '_' %+% var_y  %+% ".png"
         ggsave(lm_fit, device="png",
                filename = out.fig.name.fitline  , 
                 height = fig_lm_fit_height, width = fig_lm_fit_width , units = "px")
         
         #fig2
         library(ggstatsplot)
         out.fig.name.fitline2<-output.dir %+% paste0(main_colname_contain, collapse = "_vs_")  %+% '/' %+% var_x %+% '_' %+% var_y  %+% "_ggscatterstats.png"
         grouped.ggscatterstats<-grouped_ggscatterstats(data,x = !!var_x , y = !!var_y, grouping.var = !!group_fit_line
            ,type= "nonparametric"
            # label.var        = title,
            # label.expression = length > 200,
            # ,xlab             = "IMDB rating"
            ,ggtheme          = ggplot2::theme_grey()
            # ggplot.component = list(ggplot2::scale_x_continuous(breaks = seq(2, 9, 1), limits = (c(2, 9)))),
            # ,plotgrid.args    = list(nrow = 2,ncol=2),
            # annotation.args  = list(title = "Relationship between movie length and IMDB ratings")
                     )
         ggsave(grouped.ggscatterstats, device="png",
                filename = out.fig.name.fitline2  , 
                height = fig_lm_fit_height2, width = fig_lm_fit_width2 , units = "px")
         
         filename.out<-if_else(input$lm_picture_selected =='ggscatterstats',out.fig.name.fitline2,out.fig.name.fitline)
         ######处理数据
         incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
         list(src = filename.out,contentType = "image/png")
      })
   })

   
   
# boxplot -----------------------------------------------------------------
   print("################################################ boxplot ########################################")
   output$AD_plot_box <- renderPlot({         
            # 加上进度条显示
            withProgress(message='Program running:', detail="", {
              incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
              ######处理数据
              library(ggstatsplot)
              library(ggpubr)
              library(readxl)
              library(dplyr)
              # my_colors <- c("#647687","#1ba1e2","#008a00","#FFCE30", "#e51400")
              # my_colors <- c("#E0E0E0", "#FFB6C1")
              my_colors <- c("#00AFBB",  "#FC4E07")
              plot.title.size<-15
              point_size <-0.5
              x_lable_size<-15
              y_lable_size<-15
              
              library(readr)
              df_clean_list<- df_clean_list()
              data_<-df_clean_list[["trans"]]
              
              # data_ <- read.csv("input/SHJW_cytokine_HC_DP_20210923_filter30_filterimpute.csv") %>% data_frame()
              colnames(data_) <- make.names(colnames(data_))
              
              colnames(data_)[1] <- 'Groups'
              selected_var<-c(input$SelectP1_main_group,input$SelectP1_covar) 
              main_colname<-input$SelectP1_main_group
              
              main_colname_contain<- input$SelectP1_main_colname_contain
              
              ## dif import
              # diff1 <- read_excel("output/1_vs_2/linear_regression_adjust_Times.at.enrollment.xlsx")
              diff1 <- rt$diff_table
              diff1 <- subset(diff1, p < 0.05)
              
              
              data_ <- data_[, colnames(data_) %in% c("Groups",diff1$features) ]
              data_$Groups <- factor(data_$Groups)
              
              
              ## boxplot_original
              dir.create("output/boxplot")
              
              
              library(ggplot2)
              library(ggprism)
              
              for (i in names(data_[, -1])) {
                print(i)
                #get sig
                df_p_val <- rstatix::t_test(data_, as.formula(paste0(i,' ~ Groups')) , ref.group = main_colname_contain[1] ) %>% rstatix::add_xy_position() %>%   rstatix::adjust_pvalue(p.col = "p", method = "holm") %>% rstatix::add_significance(p.col = "p.adj")
                
                #plot
                base <- ggplot( data_, aes(x = Groups, y = get(i) ))  + 
                  geom_violin(aes( colour = Groups,fill = Groups), alpha=.5,trim = FALSE) +   
                  scale_colour_manual(values = my_colors)+
                  scale_fill_manual(values =  my_colors) +
                  geom_boxplot(aes(fill = Groups), alpha=.5,width = 0.2, colour = "black")  +
                  add_pvalue(df_p_val, label = "p = {p.adj}",label.size=5) +  
                  theme(legend.position = "none",
                        plot.title = element_text(size = plot.title.size, face="bold") ,
                        axis.text.y = element_text(size=y_lable_size) ,
                        axis.text.x = element_text(angle = 0,size=x_lable_size),#plot.margin = unit(c(0,0,0,0), "cm")
                        axis.title.y = element_blank(),
                        axis.title.x = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(color = 'black'),
                        panel.background = element_blank()
                  ) +   
                  # theme_bw() +
                  ggtitle(i)
                base
                ggsave(paste0("output/boxplot/",i, "_boxplot.pdf"))
              }   
              
              

              ######处理数据
              incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
            })
            base
          })

# OPLSDA ------------------------------------------------------------------
   print("################################################ OPLSDA ########################################")
   output$AD_plot_oplsda <- renderImage({         
      # 加上进度条显示
      withProgress(message='Program running:', detail="", {
         incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
         ######处理数据
         df_clean_list<- df_clean_list()
         data_<-df_clean_list[["trans"]]
         colnames(data_) <- make.names(colnames(data_))
         
         library(ropls)
         library(ggpubr)
         # browser()
         selected_var<-c(input$SelectP1_main_group,input$SelectP1_covar) 
         mainvar.name<-input$SelectP1_main_group
         df_count<-dplyr::select(data_,!one_of(selected_var))

         
            
         method_use<-"oplsda"
         group <- input$SelectP1_main_colname_contain
         X = df_count
         Y = factor(data_[,c(mainvar.name)], levels = group)
         
         dir.create(paste0("output/",method_use,"/figure/"),recursive = T)
         dir.create(paste0("output/",method_use,"/table/"),recursive = T)
         out.name.AD_plot_oplsda_model<-paste("output/",method_use,"/figure/",method_use,"_model",".png", sep = "")
         # png(out.name.AD_plot_oplsda_model, width = 1000, height = 1000)
         png(out.name.AD_plot_oplsda_model)
         data.oplsda <- opls(X, Y, orthoI = if_else(method_use=="oplsda",1,0) )
         dev.off()
         rt$out.name.AD_plot_oplsda_model <- out.name.AD_plot_oplsda_model
         
         vip <- data.frame(data.oplsda@vipVn)
         loading <- data.frame(data.oplsda@loadingMN)
         all(rownames(vip) == rownames(loading))
         
         # get important features
         res <- cbind(vip, loading)
         colnames(res)[1:2] <- c("vip", "loading")
         res$VIP <- ifelse(res$loading > 0, res$vip, -res$vip)
         res$VIP_direction <- ifelse(res$loading > 0, group[2], group[1])
         res <- rownames_to_column(res, var = "feature")
         res$VIP_direction <- factor(res$VIP_direction)
         
         #top30
         final_table_plsda <- dplyr::select(res, feature, vip)
         colnames(final_table_plsda) <- c("var", "value")
         # final_table_plsda <- dplyr::arrange(final_table_plsda, desc(value))[c(1:30),]
         biomarker_select_df<- dplyr::select(data_, final_table_plsda$var, mainvar.name)
         
         
         p3 <- ggbarplot(res, x = "feature", y = "VIP",
                         fill = "VIP_direction",           # change fill color by mpg_level
                         color = "white",            # Set bar border colors to white
                         palette = "lancet",            # jco journal color palett. see ?ggpar
                         sort.val = "desc",          # Sort the value in descending order
                         sort.by.groups = FALSE,     # Don't sort inside each group
                         x.text.angle = 90,          # Rotate vertically x axis texts
                         ylab = "VIP",
                         xlab = "Features",
                         legend.title = "VIP",
                         rotate = TRUE,
                         ggtheme = theme_minimal()
         )
         
         otu.name.oplsda.p3<-paste("output/",method_use,"/figure/","variable_importance_barplot_",method_use,".png", sep = "")
         png(otu.name.oplsda.p3, width = input$fig_oplsda_width, height = input$fig_oplsda_height)
         print(p3)
         dev.off()
         
         write.csv(biomarker_select_df, paste0("output/",method_use,"/table/","biomarker_select_df",  "_",method_use, ".csv")   )
         write.csv(res, paste0("output/",method_use,"/table/","res",  "_",method_use,  ".csv")   )
         
         
         
         
         
         
         ######处理数据
         incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
      })
      list(src = otu.name.oplsda.p3,contentType = "image/png") 
   })
   output$AD_plot_oplsda_model <- renderImage({    
      print(input$fig_oplsda_width)
      # browser()
      # print(input$pca_width)
      list(src = rt$out.name.AD_plot_oplsda_model,contentType = "image/png")   })

# HCPC --------------------------------------------------------------------
   print("################################################ HCPC ########################################")
   output$AD_plot_hcpc <- renderImage({         
      # 加上进度条显示
      withProgress(message='Program running:', detail="", {
         incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
         ######处理数据
         library(factoextra)
         library(FactoMineR)
         library(impute)
         library(ggpubr)
         
         df_clean_list<- df_clean_list()
         data<-df_clean_list[["trans"]]

         # browser()
         ## diff input
         diff_input<- rt$diff_table
         diff <- subset(diff_input, p < 0.05)
         df_sel <- data[, c(input$SelectP1_main_group, diff$features)]
         
         # dfs <- t(scale(data[, -c(1:3)]))
         df_sels <- scale(df_sel[, -1])
         # colnames(dfs) <- colnames(df_sels) <- rownames(data)
         
         group <- input$SelectP1_main_group
         
         dir.create("output/HCPC/result",recursive = T)
         dir.create("output/HCPC/heatmap",recursive = T)
         dir.create("output/HCPC/clusterplot",recursive = T)
         
         df <- df_sels
         # Compute PCA with ncp = 3
         res.pca <- PCA(df, ncp = 20, graph = FALSE)
         # Compute hierarchical clustering on principal components
         res.hcpc <- HCPC(res.pca, graph = FALSE, nb.clust = 2)
         
         tmp2 <- res.hcpc$desc.var$quanti
         for (i in 1:length(tmp2)) {
            tmp2[[i]] <- data.frame(tmp2[[i]])
            tmp2[[i]]$clust <- names(tmp2)[i]
         }
         # browser()
         tmp3 <- do.call(rbind.data.frame, tmp2[1:2])
         tmp3$feature <- stringr::str_split(rownames(tmp3), "\\.", simplify = T)[,2]
         
         write.csv(tmp3, paste0("output/HCPC/result/", group, "_2cluster_hcpc.csv"), row.names = F)
         
         write.csv(tmp2[[1]], paste0("output/HCPC/result/", group, "_cluster1_hcpc.csv"))
         write.csv(tmp2[[2]], paste0("output/HCPC/result/", group, "_cluster2_hcpc.csv"))
         #write.csv(tmp2[[3]], paste0("output/HCPC/result/", group, "_cluster3_hcpc.csv"))
         
         
         hcpc_data_table <- res.hcpc$data.clust
         write.csv(hcpc_data_table, file = paste0("output/HCPC/result/", group, "_cluster_detail.csv"))
         
         
         ## Plot
         fviz_cluster(res.hcpc,
                      repel = TRUE, # Avoid label overlapping
                      axes = c(1, 2),
                      geom = "point",
                      show.clust.cent = FALSE, # Show cluster centers
                      palette = c("black", as.character(see::metro_colors(c("blue", "red", "purple")))),
                      ggtheme = theme_minimal(),
                      main = "HCPC on Omics "
         )
         ggsave(paste0("output/HCPC/clusterplot/", group, "_factor_cluster_map.pdf"), width = 6, height = 5)
         ##
         # browser()
         
         out.fig.name.hcpc.cluster_dend_map <- paste0("output/HCPC/clusterplot/", group, "_factor_cluster_dend_map.png")
         fviz_dend(res.hcpc,show_labels = T,cex = 0.4)
         ggsave(filename = out.fig.name.hcpc.cluster_dend_map  , height = input$fig_hcpc_map_height, width = input$fig_hcpc_map_width , units = "px")
         
         ## vtest barplot
         desc.var <- res.hcpc$desc.var
         quanti.var <- desc.var$quanti.var
         quanti.1 <- data.frame(desc.var$quanti$`1`)
         quanti.2 <- data.frame(desc.var$quanti$`2`)
         #quanti.3 <- data.frame(desc.var$quanti$`3`)
         
         quanti.1 <- rownames_to_column(quanti.1, var = "feature")
         quanti.2 <- rownames_to_column(quanti.2, var = "feature")
         #quanti.3 <- rownames_to_column(quanti.3, var = "feature")
         
         quanti.1$vtest <- ifelse(quanti.1$v.test > 0, "positive", "negative")
         quanti.2$vtest <- ifelse(quanti.2$v.test > 0, "positive", "negative")
         #quanti.3$vtest <- ifelse(quanti.3$v.test > 0, "positive", "negative")
         
         p1 <- ggbarplot(quanti.1, x = "feature", y = "v.test",
                         fill = "vtest",               # change fill color by cyl
                         color = "white",            # Set bar border colors to white
                         palette = as.character(see::metro_colors(c("blue", "red"))),
                         sort.val = "desc",          # Sort the value in dscending order
                         sort.by.groups = FALSE,     # Don't sort inside each group
                         x.text.angle = 90,           # Rotate vertically x axis texts
                         rotate = TRUE,
                         ggtheme = theme_minimal(),
                         title = "cluster1"
         ) +
            theme(axis.text.y = element_text(size = 6),
                  axis.text.x = element_text(size = 8))
         
         p2 <- ggbarplot(quanti.2, x = "feature", y = "v.test",
                         fill = "vtest",               # change fill color by cyl
                         color = "white",            # Set bar border colors to white
                         palette = as.character(see::metro_colors(c("blue", "red"))),
                         sort.val = "desc",          # Sort the value in dscending order
                         sort.by.groups = FALSE,     # Don't sort inside each group
                         x.text.angle = 90,           # Rotate vertically x axis texts
                         rotate = TRUE,
                         ggtheme = theme_minimal(),
                         title = "cluster2"
         ) +
            theme(axis.text.y = element_text(size = 6),
                  axis.text.x = element_text(size = 8))
      
         ggsave(p1, device="png", filename = paste0("output/HCPC/clusterplot/", group, "quanti.1.png"), dpi = 300, height = input$fig_hcpc_map_height, width = input$fig_hcpc_map_width , units = "px")
         ggsave(p2, device="png", filename = paste0("output/HCPC/clusterplot/", group, "quanti.2.png"), dpi = 300, height = input$fig_hcpc_map_height, width = input$fig_hcpc_map_width , units = "px")
         
         print("################################################ END HCPC ########################################")
         ######处理数据
         incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	 
         list(src = out.fig.name.hcpc.cluster_dend_map,contentType = "image/png")
      })
   })

# EDA ---------------------------------------------------------------------
   print("################################################ EDA ########################################")
   # output$EDA <- renderImage({        observeEvent
   # EDA <- eventReactive(input$pca_group, {
   # EDA <- observeEvent(input$pca_group, {
   EDA <- eventReactive( list(input$pca_group,input$eda_width ,input$eda_height ), {
      # 加上进度条显示
      withProgress(message='Program running:', detail="", {
         incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
         ######处理数据
         library(rasterpdf)
         library(visdat)
         df_clean_list<- df_clean_list()
         df<-df_clean_list[["trans"]]
         df<-df_clean_list[["trans"]]
         df.raw<-df_clean_list[["raw"]]
         # browser()
         df.raw.vis<-vis_dat(df.raw)
         df.raw.vis.miss<-vis_miss(df.raw)
         dir.create("output/EDA",recursive = T)
         
         out.fig.name.raw.vis<-"output/EDA/df.raw_vis" %+% ".png"
         print('##############run out.fig.name.raw.vis###########')
         # png(file=out.fig.name.raw.vis,width = input$pca_width,height = input$pca_height, units = "px")
         # print(df.raw.vis)
         # dev.off()
         vis.height<-1300
         vis.width<-3000
         vis_miss.height<-1300
         vis_miss.width<-3000
         pca_width<-600
         pca_height<-600

         if(input$EDA_picture_selected =='vis'){
            vis.height<- input$eda_height
            vis.width<-input$eda_width
            print('neiceng:vis')
            print(vis.height)
            print(vis.width)
         }
         if(input$EDA_picture_selected =='vis_miss'){
            vis_miss.height<- input$eda_height
            vis_miss.width<-input$eda_width
            print('neiceng:vis_miss')
            print(vis_miss.height)
            print(vis_miss.width)
         }
         if(input$EDA_picture_selected =='pca'){
            pca_height<- input$eda_height
            pca_width<-input$eda_width
         }

         ggsave(df.raw.vis, device="png", filename =  out.fig.name.raw.vis, 
                height = vis.height, width = vis.width, units = "px")
         
         out.fig.name.raw.vis.miss<-"output/EDA/df.raw.vis.miss" %+% ".png"
         # png(file=out.fig.name.raw.vis.miss,width = input$pca_width,height = input$pca_height, units = "px")
         # print(df.raw.vis.miss)
         # dev.off()
         ggsave(df.raw.vis.miss, device="png", filename =  out.fig.name.raw.vis.miss, dpi=300,
                height = vis_miss.height, width = vis_miss.width, units = "px")
         
         
         pca_group<-input$pca_group
         pca_group<-"Group"
         df[,pca_group] <- factor(df[,pca_group])
         immune_t <- dplyr::select( df,  -contains( c(input$SelectP1_main_group,input$SelectP1_covar ))   )  %>% base::t()

         meta<- df[,c(input$SelectP1_main_group,input$SelectP1_covar)]


         p <- PCAtools::pca(immune_t, metadata = meta, removeVar = 0.1,scale = T)
         lab_show<-p$metadata
         lab_show<-lab_show[,pca_group]
         
         out.fig.name.screeplot<-"output/EDA/screeplot_label_group_" %+% pca_group %+% ".png"
         png(file=out.fig.name.screeplot)
         screeplot(p, axisLabSize = 18, titleLabSize = 22)
         dev.off()
         

         # browser()
         out.fig.name.pca<-"output/EDA/pca_label_group_" %+% pca_group %+% ".png"
         png(file=out.fig.name.pca,width = pca_width,height = pca_height, units = "px")
            plot( PCAtools::biplot(p,
                   lab = lab_show,
                   colby = pca_group,
                   hline = 0, vline = 0,
                   legendPosition = 'right') )
         dev.off()
         rt$out.fig.name.pca<-out.fig.name.pca
         rt$out.fig.name.raw.vis<-out.fig.name.raw.vis
         rt$out.fig.name.raw.vis.miss<-out.fig.name.raw.vis.miss
         rt$out.fig.name.screeplot<-out.fig.name.screeplot
         ######处理数据
         incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
         # list(src = out.fig.name.screeplot,contentType = "image/png")
      })
   })
   

   
   output$EDA_picture_selected <- renderImage({    
      print("waiceng:")
      print(input$pca_group)
      print(input$EDA_picture_selected)
      print(input$eda_height)
      print(input$eda_width)
      EDA()
      source_file<-list()
      source_file[["vis"]]=rt$out.fig.name.raw.vis
      source_file[["vis_miss"]]=rt$out.fig.name.raw.vis.miss
      source_file[["screeplot"]]=rt$out.fig.name.screeplot
      source_file[["pca"]]=rt$out.fig.name.pca
      
      list(src = source_file[[input$EDA_picture_selected ]],contentType = "image/png")
      })
   


# Gts ---------------------------------------------------------------------

   output$gts_box <- renderImage({         
      # 加上进度条显示
      withProgress(message='Program running:', detail="", {
         incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
         print("################################################ Gts ########################################")
         ######处理数据
         library(ggstatsplot)
         library(dplyr)
         library(ggplot2)
         library(ggsignif)
         
         df_clean_list<- df_clean_list()
         # df_trans<-df_clean_list[["trans"]]
         # gts.data.type<- if_else(input$gts.data.type == "rawdata","raw","trans")
         # df__<-df_clean_list[[gts.data.type]]
         df__raw<-df_clean_list[["raw"]]
         df__<-df_clean_list[["trans"]]
         # browser()
         # pca_group<-input$pca_group
         # df[,pca_group] <- factor(df[,pca_group])
         # immune_t <- dplyr::select( df,  -contains( c(input$SelectP1_main_group,input$SelectP1_covar ))   )  %>% base::t()
         # meta<- df[,c(input$SelectP1_main_group,input$SelectP1_covar)]
         out.dir<-"output/ggbetweenstats/"
         dir.create(out.dir,recursive = T)
         # plot
         # df__<-"iris"
         # group_<-"Species"
         # df__<-df_trans
         group_<-input$SelectP1_main_group
         cor_var<-input$SelectP1_covar
         var_start_num<-1+length(cor_var)
         # plot.type<-"boxviolin" 
         # type <- "noparametric"
         plot.type<-input$gts.plot.type
         type <- input$gts.type

         
         df__[,group_]<-factor(df__[,group_],levels = input$SelectP1_main_colname_contain)
         df__raw[,group_]<-factor(df__raw[,group_],levels = input$SelectP1_main_colname_contain)
         print("开始进入for 循环！")
         df_all<-NULL
         for(i in colnames(df__)[var_start_num:ncol(df__) ] ){
            # i<-colnames(df__)[5]
            # i<-"X.2.methoxyethoxy.propanoic.acid.isomer"
            df_<-NULL
            print("遍历每个特征")
            print(i)
            filename=paste0(out.dir, i, "_",plot.type, "_",type, ".png")
            # if(input$pairwise.comparisons== FALSE){pairwise.comparisons=FALSE}
            # if(input$pairwise.comparisons== TRUE &input$asterisk_label == TRUE){pairwise.comparisons=FALSE}
            # if(input$pairwise.comparisons== TRUE &input$asterisk_label == FALSE){pairwise.comparisons=TRUE}

            
            # pairwise.comparisons<- if_else(input$pairwise.comparisons== TRUE  &  input$asterisk_label == "TRUE","FALSE","TRUE")
            # pairwise.comparisons<- if_else(input$pairwise.comparisons== FALSE ,"FALSE","TRUE")
            # print("pairwise.comparisons:")
            # print(pairwise.comparisons)
            #get mpc_df
            p_ <-ggbetweenstats( data  = df__,
                                x = !!group_, 
                                y = !!i ,
                                plot.type = plot.type,
                                type = type,
                                centrality.plotting=input$gts.centrality.plotting,#可以删去图中median的标注
                                pairwise.comparisons = TRUE,
                                results.subtitle = input$results.subtitle,
                                title = " ")
            df_add <- p_$plot_env$mpc_df
            
             #get plot
             p <-ggbetweenstats( data  = df__raw,
                                 x = !!group_, 
                                 y = !!i ,
                                 plot.type = plot.type,
                                 type = type,
                                 centrality.plotting=input$gts.centrality.plotting,#可以删去图中median的标注
                                 pairwise.comparisons = FALSE,
                                 results.subtitle = input$results.subtitle,
                                 title = " ")
            # ggsave(filename = filename, width = 6, height = 5)
            
            
            ## using `pairwise_comparisons()` function to create a dataframe with results
            set.seed(123)
            # df_add2 <-pairwise_comparisons(df__,  x = !!group_, y = !!i ) %>%
            #    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
            #    dplyr::arrange(group1)
            # 
            
            df_add<- df_add  %>% dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
               dplyr::arrange(group1)
            df_add$asterisk_label <-  dplyr::case_when(df_add$p.value > 0.05 ~ "",
                                                       (df_add$p.value > 0.01)&(df_add$p.value <= 0.05)  ~ "*",
                                                       (df_add$p.value > 0.001)&(df_add$p.value < 0.01)  ~ "**",
                                                       df_add$p.value <=0.001 ~ "***"
                                                        )
            df_add<-df_add %>% dplyr::filter(asterisk_label != "")
            max.point<- df__raw[,i] %>% max() 
            
            if(max.point>100){seq_step<-max.point/10
            }else if(max.point>1000){seq_step<-max.point/10 
            }else if(max.point>10000){seq_step<-max.point/10 
            }else{seq_step<-0.2}
            
            y_position<- seq(from=0,by=seq_step,length.out = nrow(df_add) ) +max.point
            
            if( input$asterisk_label==TRUE & nrow(df_add) != 0 ){
               print("asterisk_label:")
               print(input$asterisk_label)
               ## adding pairwise comparisons using `ggsignif`
               p<-p +
                  ggsignif::geom_signif(
                     comparisons = df_add$groups,
                     map_signif_level = TRUE,
                     annotations = df_add$asterisk_label,
                     y_position = y_position,
                     test = NULL,
                     na.rm = TRUE
                     )
                  }
               p
               ggsave(filename = filename, width = 6, height = 5)
               
               df_add$variable<-i
   
               df_all<-rbind(df_all,df_add)
         }
         
         # df_all<- dplyr::select(df_all, variable, everything())
         df_all<- df_all %>% data.frame()
         df_all$groups <-paste0( df_all$groups %>% unlist(),collapse = "_" )
         
         # df_add <- apply(df_add,2,as.character)
         # fwrite(df_all, file ="myDT.csv")
         # df_all = as.matrix(df_all)
         write.csv(  df_all,paste0(out.dir,"ggbetweenstats_compar_",plot.type,type,".all.csv")  )
         

         ######处理数据
         incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
         list(src = filename,contentType = "image/png")
      })
   })
   print("################################################End Gts #####################################")

# Heatmap -----------------------------------------------------------------
   print("################################################ heatmap ########################################")
   output$AD_plot_heatmap <- renderImage({         
     # 加上进度条显示
     withProgress(message='Program running:', detail="", {
       incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
       ######处理数据
       library(readxl)
       library(dplyr)
       library(purrr)
       library(tibble)
       library(MASS)
       library(broomExtra)
       library(rstatix)
       library(ComplexHeatmap)
       library(readxl)
       library(see)
       library(circlize)
       library(readr)

       df_clean_list<- df_clean_list()
       data_<-df_clean_list[["trans"]]
       
       # data_ <- read.csv("input/SHJW_cytokine_HC_DP_20210923_filter30_filterimpute.csv") %>% data_frame()
       colnames(data_) <- make.names(colnames(data_))
       
       colnames(data_)[1] <- 'Groups'
       selected_var<-c(input$SelectP1_main_group,input$SelectP1_covar) 
       main_colname<-input$SelectP1_main_group
       
       main_colname_contain<- input$SelectP1_main_colname_contain
       
       ## dif import
       # diff1 <- read_excel("output/1_vs_2/linear_regression_adjust_Times.at.enrollment.xlsx")
       diff1 <- rt$diff_table
       diff1 <- subset(diff1, p < 0.05)
       
       
       data_ <- data_[, colnames(data_) %in% c("Groups",diff1$features) ]
       data_$Groups <- factor(data_$Groups)
       
       
       # browser()
       dfs1 <- t(apply(data_[,-1], 2, scale))
       colnames(dfs1) <- data_ %>% rownames()
       all_group1 <- data_$Groups
       
       
       groups<-main_colname_contain
       mycol <- c("#90C2D3", "#0B318F", "#C4780A", "#7E318E", "#9B0A33")
       col_fun <- circlize::colorRamp2(seq(-4,4),c("#0D8CFF","#43A6FF","#5EB3FF","#93CCFF","#FFFFFF","#FFBBAA","#FF7755","#FF6039","#FF3300"))
       
       ht_list <- NULL
       for (i in 1:length(groups)) {
         ht_list <- ht_list + Heatmap(dfs1[,all_group1 == groups[i]],
                                      cluster_columns = T, cluster_rows = T,
                                      col = col_fun, 
                                      row_names_gp = gpar(fontsize =10),
                                      row_title_gp = gpar(fontsize = 10),
                                      row_title_rot = 0,
                                      column_title = groups[i],
                                      column_title_gp = gpar(fontsize=11),
                                      column_names_gp = gpar(fontsize = 7),
                                      show_column_names = T)
       }
       dir.create("output/heatmap/",recursive = T)
       # browser()
       
       
       png(paste0("output/heatmap/differential_lm_diff_heatmap_adjusted.png"),width = input$heatmap_width,height = input$heatmap_height,units = "px")
       draw(ht_list, show_heatmap_legend = F)
       dev.off()
       

        filename <-"output/heatmap/differential_lm_diff_heatmap_adjusted.png"
        # list(src = filename,contentType = "image/pdf",width = 900,height = 800)

       ######处理数据
       incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
       # list(src = filename,contentType = "image/png",width = 1500,height = 300)
       list(src = filename,contentType = "image/png")
     })
   })

# Download ----------------------------------------------------------------
   ###downloadData and send the zip of result
   output$downloadexample <- downloadHandler(
      filename = 'example.zip',
      content = function(fname) {
         file.copy("www/example.zip", fname)
      }
   )
   
   output$downloadvideo <- downloadHandler(
      filename = 'video.mp4.zip',
      content = function(fname) {
         file.copy("www/video.mp4.zip", fname)
      }
   )
   
   print("################################################ Download all the result ########################################")
   output$downloadData <- downloadHandler(
      
      filename = 'files.zip',
      content = function(fname) {
         withProgress(message='Program running:', detail="Still water run deep", {
            incProgress(0.2, detail = paste("Still water run deep"))  
            #tmpdir <- tempdir()
            # tempdir_ <- "/data/zhiyu/software/16s/rarecure_16s/output"
            print("uuu")
            
            tempdir_ <- getwd()
            # setwd(tempdir_)
            print("tempdir:" %+% tempdir_)
            print("getwd():" %+% getwd())
            fs<- dir(tempdir_,recursive = T)
            
            print("dir_contains:" %+% fs)
            zip(zipfile=fname, files=fs)
            
            incProgress(0.3, detail = paste("Still water run deep"))   
         })
      },
      contentType = "application/zip"
   )

# Refreshing --------------------------------------------------------------
   print("Initializing")
   print("################################################ refreshing ########################################")
   ###refreshing
   observeEvent(input$switchtab,{
      aggg_result = -1
      if(aggg_result == -1)
      {
         session$reload()
         return()
         print("session reload not working")
      }
      print("Code running this line")
      output$code_ran <- renderText("code Ran this line without refreshing")
   })

# Session -----------------------------------------------------------------
   print("################################################ sessionInfo ########################################")
   ###refreshing sessionInfo
   output$sessionInfo <- renderPrint({ capture.output(sessionInfo())})
   
   
   
   
   
}







