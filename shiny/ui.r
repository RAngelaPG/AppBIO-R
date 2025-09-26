library(adegenet)
library(ape)
library(cluster)
library(car)
library(chromoMap)
library(corehunter)
library(data.table)
library(dendextend)
library(dplyr)
library(DT)
library(ff)
library(factoextra)
library(ggtree)
library(ggplot2)
library(Hmisc)
library(htmlwidgets)
library(IntAssoPlot)
library(naturalsort)
library(pedigreemm)
library(plotly)
library(plyr)
library(rJava)
library(RColorBrewer)
library(reshape)
library(reshape2)
library(rtracklayer)
library(shiny)
library(shinyBS)
library(shinyFiles)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyalert)
library(statgenGWAS)
library(stringr)
library(vcfR)
library(vegan)



#BIOR
source("local/GeneraXML.R")
source("local/Biodv2.R")
source("local/gdiv.R")
source("local/CountFreq.R")
source("local/CountFreqNoImpute2.R")
source("local/functCORE.R")
source("local/Functions.R")
source("local/utils.R")
source("local/CompChrom.R")

radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){

  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

body<-dashboardBody(
		tags$head(tags$style(HTML(' /* main sidebar */
                                           .skin-green .main-sidebar {
                                           background-color: #0b6b21;
                                           font-family: "Georgia",Times,"Times New Roman", serif;
                                           font-weight: bold;
                                           font-size: 15px;
                                           }
                                          /* active selected tab in the sidebarmenu */
                                           .skin-green .main-sidebar .sidebar .sidebar-menu .active a{
                                          background-color: #400080;
                                          } 
                                          
                                          /* other links in the sidebarmenu when hovered */
                                          .skin-green .main-sidebar .sidebar .sidebar-menu a:hover{
                                          background-color: #6600cd;
                                          }
                                         
										  .content-wrapper {background-color:#E4E7E4;}
										  .skin-green .left-side, .skin-green .wrapper {
                                           background-color: #E4E7E4;
                                           }										  
									  ')
                                     )
                          ),
  ###################################################################################################################################################  
  ##################################################################################################################################################  
  tabItems( 
    ###################################################################################################################################################  
    ##################################################################################################################################################  
    tabItem(tabName="datos",
            fluidRow(
               tags$style(".nav-tabs {background-color: #FFF;}

               .nav-tabs-custom .nav-tabs li.active:hover a, .nav-tabs-custom .nav-tabs li.active a {
               background-color: transparent;
               border-color: transparent;
               }

               .nav-tabs-custom .nav-tabs li.active {
               border-top-color: #006747;
               background-color: #77BC1F;
               }
			   "),
              column(width = 4,
                     box(title = "Read Genetic Data File",  background="olive", width = NULL,  collapsible=F, collapsed=F,                         
						 fluidRow(
						 column(4,radioButtons("startAna", "Type analysis",choices = c(Chromosome="StarChrom",PAVs="pavs",Biodiversity = "StarBio"))),
						 column(12,
						 column(4,shinyFilesButton('filegen', 'Choose geno file', 'Select File', FALSE)),
						 column(4,shinyFilesButton('fileRdata', 'Choose RData file', 'Select File', FALSE)),
						 column(4,actionButton("readgenofile", "Upload Data",icon("cloud-arrow-up"), style = "color: white; background-color: #400080;border-color: #400080;"))
						 )
						),
                         # Horizontal line ----
                         tags$hr(),
                         fluidRow(
                           column(4,radioButtons("typedata", "Code Data",choices = c( VCF="vcfile", AlleleFrequency = "FREQ", SNP = "SNP", Counts="CUENTA", DistanceMatrix="DistMat"))),
                           column(8,
                                  column(4,textInput("ht1","AA","0")),
                                  column(4,textInput("ht2","Aa","2")),
                                  column(4,textInput("ht3","aa","1"))								  
                           )
                         ),
                         tags$hr(),
                         h4("Filters"),
                         fluidRow(
                           column(4,textInput("missval","%NA Genotype","0")),
                           column(4,textInput("mayorque",">polymorphism","0.95")),
                           column(4,textInput("menorque","<polymorphism","0.05"))
                         )
                     ),
					 box(title="Read second VCF file to compare", background="green", width = NULL,  collapsible=F, collapsed=F,                         
                      shinyFilesButton('fileVCF', 'Add VCF file', 'Select VCF file', FALSE)                  
					 ),
					 box(title="About", background="purple", width = NULL,  collapsible=F, collapsed=F, align="center",                        
                        h3("BIO-R (Copyright 2016 CIMMYT) Version 5.0 (2024-Nov)"),
						h3("maintainer: r.a.pacheco@cgiar.org"),
						h4("Authors:"),
						h4("Angela Pacheco"),
						h4("Guadalupe Valdez"),
						h4("Juan Burgueno"),
						h4("Jose Crossa"),
						h4("Keith Gardner"),
						h4("Fernando Toledo"),
						h4("Gregorio Alvarado"),
						h4("Francisco Rodriguez"),
						tags$hr(),
						h5("This program is based in R. Any R component of this program as well as the program as"),
						h5("a whole developed by CIMMYT are hereby licensed as per the terms of the GNU General"),
						h5("Public License version 3(available at http://www.gnu.org/licenses/gpl-3.0.html), as"),
						h5("specified below."),
						h5("This program is free software; you can redistribute it and/or modify it under the terms"),
						h5("of the GNU General Public License as published by the Free Software Foundation; either"),
						h5("version 3 of the License, or any later."),
						h5("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;"),
						h5("without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."),
						h5("See the GNU General Public License for more details."),
						h5("You should have received a copy of the GNU General Public License along with this program; if not,"),
						h5("to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA."),
						h5("For further information about the license, please contact CIMMYT at CIMMYT-Knowledge-Center@cgiar.org"),
						h5("or at Km. 45 Carretera Mexico-Veracruz, El Batan, Texcoco, Estado de Mexico, Mexico, C.P. 56237."),
						h5("We have invested a lot of time and effort in creating this platform, please cite it when using it for data analysis."),
						img(src="Bank_logo.png",width="40%")
					 )
             ),           
            
            column(width=8,
              tabBox(width =NULL,title = tagList(shiny::icon("database"), "Data"),side = "left",
                     tabPanel(h5(style="font-weight: bold","Chromo map"), chromoMap::chromoMapOutput("seeChromPlot")), 
					 tabPanel(h5(style="font-weight: bold","Chromo map comparison blast"), chromoMap::chromoMapOutput("BChromPlot")),
					 tabPanel(h5(style="font-weight: bold","PAVs"), plotlyOutput("seePAVS",height = "950px",width = "1250px")),	
					 tabPanel(h5(style="font-weight: bold","Data for Bio"), DT::dataTableOutput("seeDataGen"))				 
              )
            )
          )
    ),
	
  ###################################################################################################################################################  
  ###################################################################################################################################################  
  tabItem(tabName="bio",
          fluidRow(
            column(width = 8,
            box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"MDSplot2D"), solidHeader=T, closable=F, collapsible=T, collapsed=T, status="success", 
			sidebar=boxSidebar(id="boxMDS2D",width = 25,
                      background = "#8f8d8d",
                      useShinyalert(force=TRUE),         
                      selectInput('xcol', 'X Variable',choices =""),
                      selectInput('ycol', 'Y Variable',choices =""),
                      selectInput('zcol', 'Z Variable',choices =""),
                      tags$hr(),
                      selectInput('catv','Group',choices = '',selected=NULL),					  
                      #Creacion de color para cambiar el color de los puntos del grafico de acuerdo a una lista de colores generada
                      selectInput('color','Choose a color',choices = '',selected="",multiple=T),
                      tags$hr(),
                      textInput('tp','Plot Title',value='MDS Plot'),
                      selectInput('pnc','Title Color',choices=colors()[-c(1:25,27:30,37:46,48,57:67,69,80,82,83,85:89,92:93,97:98,101:106,108:113,123:127,
                                                                          130:131,138,140:141,152:253,255:256,260:366,377:392,394:447,449,451,463:464,468:469,
                                                                          478:489,492,499:500,504:505,509:511,513:534,536:546,548:550,553:554,557:561,563:564,
                                                                          569:570,579:583,585:586,589:609,611:612,617,620:629,631:632,636:637,642:644,646:651)],
                                  selected = 'blue'),
                      sliderInput('ts','Title Size',min=5,max=20,step = 1,value = 12),
                      tags$hr(),
                      textInput('tx','X Axis Label',value = 'Factor 1 ( 40.51 %)'),
                      textInput('ty','Y Axis Label',value = 'Factor 2 ( 28.09 %)'),
                      textInput('tz','Z Axis Label',value = 'Factor 3 ( 8.91 %)'),
                      selectInput('ac','Axes color',choices=colors()[-c(1:25,27:30,37:46,48,57:67,69,80,82,83,85:89,92:93,97:98,101:106,108:113,123:127,
                                                                        130:131,138,140:141,152:253,255:256,260:366,377:392,394:447,449,451,463:464,468:469,
                                                                        478:489,492,499:500,504:505,509:511,513:534,536:546,548:550,553:554,557:561,563:564,
                                                                        569:570,579:583,585:586,589:609,611:612,617,620:629,631:632,636:637,642:644,646:651)],
                                  selected = 'red'),
                      sliderInput('szl','Axes Label Size',min=5,max=20,step=1,value=12),
                      tags$hr(),
                      #Creacion de bkgp para cambiar el color de fondo en el grafico
                      selectInput('bkgp','Plot Background color',
                                  choices=c("white","antique","azure","beige",               
                                            "black","blanchedalmond","burlywood",           
                                            "cornsilk","darkgray","dimgray","floralwhite",         
                                            "gainsboro","ghostwhite","gray","honeydew",            
                                            "ivory","lavender","lavenderblush","lemonchiffon",        
                                            "lightcyan","lightgoldenrodyellow",
                                            "lightgray","lightslategray","lightsteelblue","linen","mintcream","mistyrose",
                                            "moccasin","oldlace","palegoldenrod","papayawhip","peachpuff","seashell","snow",
                                            "tan","thistle","whitesmoke"),selected=NULL),
                      tags$hr(),
                      selectInput('eti','Label points',choices=''),           
                      sliderInput('size','Points size',min=5,max=25,value=7)
                    ),
                    #cuadro de texto que se mostrara por default con una direccion en la que se encuentra el grafico 2d
                    verbatimTextOutput("default1",placeholder = TRUE),
                    #grafico 2d
                    div(plotlyOutput("try",height = "750px",width = "950px"),align="center")
            ),
      
       box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"MDSplot3D"), solidHeader=T, closable=F, collapsible=T, collapsed=T, status="success", 
	   sidebar=boxSidebar(id="boxMDS3D",width = 25,
                          background = "#8f8d8d",
                          selectInput('xcol3D', 'X Variable',choices =""),
                          selectInput('ycol3D', 'Y Variable',choices =""),
                          selectInput('zcol3D', 'Z Variable',choices =""),
                          tags$hr(),
                          selectInput('catv3D','Group',choices = '',selected=NULL),					  
                          #Creacion de color para cambiar el color de los puntos del grafico de acuerdo a una lista de colores generada
                          selectInput('color3D','Choose a color',choices = '',selected="",multiple=T),
                          tags$hr(),
                          textInput('tp3D','Plot Title',value='MDS Plot'),
                          selectInput('pnc3D','Title Color',choices=colors()[-c(1:25,27:30,37:46,48,57:67,69,80,82,83,85:89,92:93,97:98,101:106,108:113,123:127,
                                                                              130:131,138,140:141,152:253,255:256,260:366,377:392,394:447,449,451,463:464,468:469,
                                                                              478:489,492,499:500,504:505,509:511,513:534,536:546,548:550,553:554,557:561,563:564,
                                                                              569:570,579:583,585:586,589:609,611:612,617,620:629,631:632,636:637,642:644,646:651)],
                                      selected = 'blue'),
                          sliderInput('ts3D','Title Size',min=5,max=20,step = 1,value = 12),
                          tags$hr(),
                          textInput('tx3D','X Axis Label',value = 'Factor 1 ( 40.51 %)'),
                          textInput('ty3D','Y Axis Label',value = 'Factor 2 ( 28.09 %)'),
                          textInput('tz3D','Z Axis Label',value = 'Factor 3 ( 8.91 %)'),
                          selectInput('ac3D','Axes color',choices=colors()[-c(1:25,27:30,37:46,48,57:67,69,80,82,83,85:89,92:93,97:98,101:106,108:113,123:127,
                                                                            130:131,138,140:141,152:253,255:256,260:366,377:392,394:447,449,451,463:464,468:469,
                                                                            478:489,492,499:500,504:505,509:511,513:534,536:546,548:550,553:554,557:561,563:564,
                                                                            569:570,579:583,585:586,589:609,611:612,617,620:629,631:632,636:637,642:644,646:651)],
                                      selected = 'red'),
                          sliderInput('szl3D','Axes Label Size',min=5,max=20,step=1,value=12),
                          tags$hr(),
                          #Creacion de bkgp para cambiar el color de fondo en el grafico
                          selectInput('bkgp3D','Plot Background color',
                                      choices=c("white","antique","azure","beige",               
                                                "black","blanchedalmond","burlywood",           
                                                "cornsilk","darkgray","dimgray","floralwhite",         
                                                "gainsboro","ghostwhite","gray","honeydew",            
                                                "ivory","lavender","lavenderblush","lemonchiffon",        
                                                "lightcyan","lightgoldenrodyellow",
                                                "lightgray","lightslategray","lightsteelblue","linen","mintcream","mistyrose",
                                                "moccasin","oldlace","palegoldenrod","papayawhip","peachpuff","seashell","snow",
                                                "tan","thistle","whitesmoke"),selected=NULL),
                          tags$hr(),
                          selectInput('eti3D','Label points',choices=''),           
                          sliderInput('size3D','Points size',min=5,max=25,value=7)
                        ),
                        #cuadro de texto que se mostrara por default con una direccion en la que se encuentra el grafico 3d
                        verbatimTextOutput("default3d",placeholder = TRUE),
                        #grafico 3d
                        div(plotlyOutput("try3d",height = "750px",width = "950px"),align="center")
       ),
       box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"Heatmap for Distance Matrix"), solidHeader=T, closable=F, collapsible=T, collapsed=T,status="success", 
	       sidebar=boxSidebar(id="boxHeat",width = 25,
                 background = "#8f8d8d",
                 selectInput('colorheat', 'Scale color',choices =c("Jet","Reds","Blues","Viridis"), selected="Jet")
               ),
               verbatimTextOutput("defaultheat",placeholder = TRUE),
               #grafico heatmap
               div(plotlyOutput("heat",height = "750px",width = "950px"),align="center")
       ),
       box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"Dendogram"), solidHeader=T, closable=F, collapsible=T, status="success", 
	   sidebar=boxSidebar(id="boxdend",width = 25,
                 background = "#8f8d8d", 
                 sliderInput('sizeline','Size cluster line',min=0.1,max=2,value=0.9),
                 sliderInput('sizelab','Size labels',min=0.5,max=5,value=3),
                 sliderInput('space','Spaces',min=0.1,max=2,value=0.2),
                 tags$hr(),
                 selectInput('poslen','Position legend',choices=c('left','top','bottom','right'),selected='left'),
                 tags$hr(),
                 selectInput('typeclust','Type',choices=c('circular','rectangular'), selected='rectangular'),
                 selectInput('catvdend','Group',choices = '',selected=NULL),					  
                 #Creacion de color para cambiar el color de los puntos del grafico de acuerdo a una lista de colores generada
                 selectInput('colordend','Choose a color',choices = '',selected="",multiple=T)
               ),
               #radioButtons('desdend','File Type',choices = list('png','pdf'),selected = 'png'),
               verbatimTextOutput("defaultdend",placeholder = TRUE),
               fluidRow(
                 column(11, align="center",plotOutput("dend",height = "950px",width = "750px")))
       )
       #close column
       ),
       column(width=4,
              box(width=12,title=tagList(img(src="Biopng.png",height="25"),"More options"), solidHeader=T, closable=F, collapsible=T, status="success", 
                      shinyFilesButton('fileenvbio', 'Add external group', 'Select csv file', FALSE),
                      checkboxInput("quitomono","Remove monomorphic markers",value=FALSE),
				  	  checkboxInput("gapS","Optimize number of cluster",value=FALSE),
                      textInput('nclust','No. Clusters',value='3'),
                      radioButtons("distk", "Genetic Distance",choices = c(Rogers = "Rogers", Nei = "Nei"))
              ),
              tabBox(width=12,side="right",title = tagList(img(src="Biopng.png",height="25"),"Results"), 
                tabPanel(h5(style="font-weight: bold","Summary Diversity"),DT::dataTableOutput("seeDataDiver")),
                tabPanel(h5(style="font-weight: bold","Population structure"),DT::dataTableOutput("seeDataGDiver")),
				 tabPanel(h5(style="font-weight: bold","AMOVA"),DT::dataTableOutput("seeDataGAmova"))
                )
       #close column   
       )
       
      )    
    ),
  ###################################################################################################################################################  
  ###################################################################################################################################################  
  tabItem(tabName="core",
          box(solidHeader = TRUE,title = tagList(img(src="Biopng.png",height="25"),"Select Options"),width = 6,status = "success",
            fluidRow(
              tags$head(
                tags$style(type="text/css", "#inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                #inline .form-group { display: table-row;}")
              ),
              column(width=6,
                    useShinyjs(),
                    boxPad(color="gray",
                       checkboxGroupInput("datause", "Information to use:",c("Phenotypic Data" = "phendat",
                       "Genetic Data" = "gendat","Matrix Distance" = "distdat")),
                       shinyFilesButton('filedistbio', 'Choose distance matrix (csv file)', 'Select File', FALSE),
					   shinyFilesButton('filephendatbio', 'Choose phenotypic data (csv file)', 'Select File', FALSE)
                           )
                ),
              column(width=6,
                    boxPad(color="gray",tags$div(id = "inline", textInput('score','Size of Core Subset:',value = '0.2',width = '50%')))
              )
            ),
            tags$br(),
            fluidRow(
              column(width=12,
              boxPad(color="gray",
                     tags$b("OPTIMIZATION: Fill in the blanks using numbers between [0,1], this is the weight assigned to the objective
                     when maximizing a weighted index, such that their sum is equal to one.")
                     ),
              tags$br()
              ),
              tags$br(),
              column(width=6,
                     boxPad(color="gray",
                         tags$b("EN: Maximizes the average distances between each selected individual and the closest other selected item in the core."),
                         tags$hr(),
                         tags$div(id = "inline", textInput('mrdEN','Modified Rogers Distance',value = '0'),
                         textInput('csedEN','Cavalli-Sforza and Edwards Distance',value = '0'),
                         textInput('gdEN','Gower Distance',value = '0'),
                         textInput('pcdEN','Precomputed Distance',value = '0'))
                        ),
                     tags$br(),
                     boxPad(color="gray",
                            tags$b("EE: Maximizes the average distances between each pair of selected individuals in the core."),
                            tags$hr(),
                            tags$div(id = "inline",textInput('mrdEE','Modified Rogers Distance',value = '0'),
                            textInput('csedEE','Cavalli-Sforza and Edwards Distance',value = '0'),
                            textInput('gdEE','Gower Distance',value = '0'),
                            textInput('pcdEE','Precomputed Distance',value = '0'))
                        )
                    ),
              column(width=6,
                     boxPad(color="gray",
                            tags$b("AN: Minimizes the average distances between eachi individual closest selected item in the core."),
                            tags$hr(),
                            tags$div(id = "inline",textInput('mrdAN','Modified Rogers Distance',value = '0'),
                            textInput('csedAN','Cavalli-Sforza and Edwards Distance',value = '0'),
                            textInput('gdAN','Gower Distance',value = '0'),
                            textInput('pcdAN','Precomputed Distance',value = '0'))
                           ),
                     tags$br(),
                     boxPad(color="gray",
                            tags$div(id = "inline",textInput('SH','SH: Maximizes the entropy, as used in information theory, of the selected core.',value = '0')),
                            tags$hr(),
                            tags$div(id = "inline",textInput('HE','HE: Maximizes the expected proportion of heterozygous loci in offspring produced from random crossing within the selected core.',value = '0')),
                            tags$hr(),
                            tags$div(id = "inline",textInput('CV','CV: Maximizes the proportions of alleles observed in the full dataset that are retained in the selected core.',value = '0'))
                     )
                     )
            )
          ),
          box(solidHeader = T, title=tagList(img(src="Biopng.png",height="25"),"Warnings"),status = "success",  
              uiOutput(outputId = "defaultcore")
          )
  ),
  #close tab core
  ###################################################################################################################################################  
  ###################################################################################################################################################  
  tabItem(tabName="gwas",
		fluidRow(
		column(width=4,
          box(solidHeader = TRUE,title = tagList(img(src="Biopng.png",height="25"),"Select options"),width = 12,status = "success",     
              tags$head(
                tags$style(type="text/css", "#inline label{ display: table-cell; text-align: left; vertical-align: middle; } 
                #inline .form-group { display: table-row;}")
              ),
              column(width=12,
				boxPad(color="gray", width=12,
                    useShinyjs(),
						h3("Add filters for markers information"),
					   textInput('gwastraits', 'Write traits to analyze separated by commas',value=""),
					   textInput('nmissgeno','Genotypes with high proportion of missing values in the marker matrix can be removed',value='0.2'),
					   textInput('nmissind','SNPs with a high proportion of missing values can be removed',value='0.2'),
					   textInput('nmf','Can be used to remove SNPs with a lower MAF',value='0.05'),					   
					   tags$hr(),					   
					   h3("Load phenotypic data"),
                       shinyFilesButton('filephen', 'Choose phenotypic info (csv file)', 'Select File', FALSE),
					   tags$br(),
					   tags$hr(),
					   h3("Load genome reference data"),
					   textInput('RegGenPlotraits', 'Write trait to associate',value=""),
					   textInput('RegGenPlotChr', 'Chromosome',value="1"),
					   textInput('RegGenPlotChrThr','LOD-threshold',value='3.5'), 
					   textInput('RegGenPlotrigth','Rigth-threshold',value='Search suggestion in ChromDim file'),
					   textInput('RegGenPlotleft','Left-threshold',value='Search suggestion in ChromDim file'),
					   tags$br(),
					   shinyFilesButton('fileplants', 'Choose genome ref (gtf.gz file)', 'Select File', FALSE),
					   tags$br(),
					   tags$br()
					  ) 
					)
			)
		),
		column(width=8,
			box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"Regional Gen Plot"), solidHeader=T, closable=F, collapsible=T, collapsed=T, status="success", 
				#cuadro de texto que se mostrara por default con una direccion en la que se encuentra el grafico 2d
              verbatimTextOutput("defaultRegGenPlotChr",placeholder = TRUE),  
                
              fluidRow(
                 column(11, align="center",plotOutput("RegGenPlotChrPlot",height = "650px",width = "1050px"))
				)			  
            ),
			box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"Manhattan Plot"), solidHeader=T, closable=F, collapsible=T, collapsed=T, status="success", 
			sidebar=boxSidebar(id="boxManhattan",width = 25,
                      background = "#8f8d8d",
                      useShinyalert(force=TRUE),         
                      selectInput('ManTrait', 'Trait',choices =""),
					  textInput('ManyThr','LOD-threshold',value='3.5'),
					  textInput('Manchr','Chromosomes to be plotted ',value='All'),
					  selectInput('MancolPalette', 'Color',choices =c("Blue","Gray","Green","Colors"), selected="Blue")
					),
				#cuadro de texto que se mostrara por default con una direccion en la que se encuentra el grafico 2d
              verbatimTextOutput("defaultMan",placeholder = TRUE),
                fluidRow(
                 column(11, align="center",plotOutput("ManPlot",height = "550px",width = "950px"))
				)  
            ),
			box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"QTL Plot"), solidHeader=T, closable=F, collapsible=T, collapsed=T, status="success", 
			sidebar=boxSidebar(id="boxQTL",width = 25,
                      background = "#8f8d8d",
                      useShinyalert(force=TRUE),         
                      textInput('QTLyThr','LOD-threshold',value='3.5'),
					  textInput('QTLchr','Chromosomes to be plotted ',value='All')
					),
				#cuadro de texto que se mostrara por default con una direccion en la que se encuentra el grafico 2d
              verbatimTextOutput("defaultQTL",placeholder = TRUE),
                fluidRow(
                 column(11, align="center",plotOutput("QTLPlot",height = "750px",width = "950px"))
				)
            ),
			box(loadingState(),width=12,title=tagList(img(src="Biopng.png",height="25"),"QQ Plot"), solidHeader=T, closable=F, collapsible=T, collapsed=F, status="success", 
			sidebar=boxSidebar(id="boxQQ",width = 25,
                      background = "#8f8d8d",
                      useShinyalert(force=TRUE),         
                      selectInput('QQTrait', 'Trait',choices ="")  
					),
				#cuadro de texto que se mostrara por default con una direccion en la que se encuentra el grafico 2d
              verbatimTextOutput("defaultQQ",placeholder = TRUE),
                fluidRow(
                 column(11, align="center",plotOutput("QQPlot",height = "750px",width = "750px"))
				)
            )
		)
		)
	)
	#close tab gwas

)
)

notificationItemWithAttr <- function(text, icon = shiny::icon("warning"), status = "success", href = NULL, ...) {
  if (is.null(href)) 
    href <- "#"
    icon <- tagAppendAttributes(icon, class = paste0("text-", status))
    tags$li(a(href = href, icon, text, ...))
}
#class="logo-lg"
ui <- dashboardPage(skin = "green",title="BSU-CIMMYT",
                        dashboardHeader(
                           title =tagList(
                                       span(class="logo-lg","BSU&QG"), 
                                       img(src = "Allicon.png",style= 'margin-left:-10px', height="40")
                                  ),
                          tags$li(class = "dropdown",  
                                  img(src='logocimmytA.png',style= 'margin-right:50px',width="350",height="45")
                          ),
                          dropdownMenu(
                            type = "notifications", 
                            icon = icon("question-circle"),
                            badgeStatus = NULL,
                            headerText = "See Manuals:",
                            notificationItemWithAttr("BIO-R", icon = icon("file"),
                                             href = "BIO-R User Manual.pdf",target = "_blank")
                          )
                        ),
                        dashboardSidebar(
                          sidebarMenuOutput("menu")
                        ),
                        body						
                          
)

