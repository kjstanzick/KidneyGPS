
source("packages.r")
load("data/workspace.RData")

##########################
#
### Design User Interface
#
##########################
ui <- 
  fluidPage(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script('
                                        $(document).on("keyup", function(e) {
                                        if(e.keyCode == 13){
                                        Shiny.onInputChange("keyPressed", Math.random());
                                        }
                                        });
                          ')
    ),
  #useShinyjs(),
  
  tags$h1(list(tags$img(id="kidney_image",
                   height="50px",
                   width="50px",
                   src= "kidney-g8a69f0709_1920.png"),
          tags$em("Kidney "), tags$b("G"),"ene",tags$b("P"),"rioriti",tags$b("S"),"ation - ",tags$b("KidneyGPS"))),
  tags$p(list("This platform summarizes information on each of 5906 genes overlapping the 424 loci identified for ",
              tags$abbr(title="estimated glomerular filtration rate from serum creatinine",tags$b("eGFRcrea")), " based on a ",tags$abbr(title="genome-wide association studies","GWAS")," meta-analysis of ", tags$a(href="https://www.ukbiobank.ac.uk/", target="_blank", "UK Biobank data"), "and ", tags$a(href="https://ckdgen.imbi.uni-freiburg.de", target="_blank", "CKDGen consortium"), "data", tags$b("(n=1,201,909)"), 
         tags$strong(tags$a(href="https://www.nature.com/articles/s41467-021-24491-0#Sec31", target="_blank", "[Stanzick et al. Nat. Commun. 2021]")),
         ". Focussed on European ancestry ",tags$b("(n=1,004,040)")," , 634 independent signals across the 424 loci were identified using approximate conditioned analysis. "
         , "For each variant in these 634 signals, the posterior probability of association (PPA) was computed and, for each signal, a 99% credible set of variants was derived "
         ,"(i.e. smallest set of variants with >99% cumultative PPA). ",
         "Credible (set) variants are considered the most likely variants to drive the association signal (particularly those with high PPA).")),
  
   
  
  #tags$br(),
  
  navbarPage(tags$p(tags$b("G"),"ene",tags$b("P"),"rioriti",tags$b("S"),"ation"),
             
             
             
             # tabPanel("Home",
             #          tags$p("This platform summarizes information on genes and variants of ",tags$b("eGFRcrea"), " associated loci. Origin of this association is a GWAS meta-analysis of ", tags$a(href="https://www.ukbiobank.ac.uk/", "UK Biobank data"), "and ", tags$a(href="https://ckdgen.imbi.uni-freiburg.de", "CKDGen consortium"), "data", tags$b("(n=1,201,909)"), 
             #                 ". Loci were selected based on an association p-value treshhold of P<5E-8. In total 424 eGFRcrea asscoiated loci and 634 independent association signals were found."),
             #          tags$p("The displayed data is based on published work, that can be found ", tags$strong(tags$a(href="https://www.nature.com/articles/s41467-021-24491-0#Sec31", "here."))),
             #          
             #          ),
             
             tabPanel("Genes",
                      
                      tags$h3("Gene Search"),
                      p("Search for genes overlapping any of the 424 loci:"),
                      
                      sidebarLayout(

                       # Sidebar panel for inputs ----
                       sidebarPanel(
                        textInput(inputId="Genename", label="Gene Name", value = ""),
                        # tags$h4(actionButton(inputId="question_gene_batch",label="?",class="question"),"Upload Gene list:"),
                        tags$br(),
                        
                        div(class="gene-batch",textAreaInput("gene_batch",label="Paste a list of genes",value=NULL)),

                        tags$hr(),
                        # Input: Select a file ----
                        fluidRow(
                            column(12,
                                   div(id="gene-upload-section",
                                       fileInput("file2", "Choose txt file with a list of genes (max. 2000)",
                                       multiple = FALSE,
                                       accept = c("text/txt",
                                                  "text/comma-separated-values,text/plain",
                                                  ".txt")),
                                   actionButton(inputId="reset_gene", label="reset", class="go"))
                                   ),
                         ),

                         # Horizontal line ----
                         tags$hr(),
                         shinyjs::hidden(
                           div(id="gene_batch_upload_options",
                           # Input: Checkbox if file has header ---- not reqired if just a list should be uploaded
                           #checkboxInput("header_gene", "Header", TRUE),

                           # Input: Select separator ----
                               radioButtons("sep_gene", "Separator",
                                            choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ","),

                           # # Input: Select quotes ----
                           #     radioButtons("quote_gene", "Quote",
                           #                  choices = c(None = "",
                           #                              "Double Quote" = '"',
                           #                              "Single Quote" = "'"),
                           #                  selected = '"')

                           )),
                        fluidRow(
                          column(2, offset=10,
                                 shinyjs::hidden(div(id="gene_go_div",actionButton(inputId="gene_go", label="go", class="go"))),
                          )
                        )
                       ),
                       
                       mainPanel(
                         wellPanel(
                           
                           
                           
                           
                           tags$h4("Specification of functional evidence for the searched gene (s) displayed in GPS table:",actionButton(inputId="question_GPS_gene",label="?",class="question")),
                           
                           
                           div(class="CADD",checkboxGroupInput(inputId="columns1", label=span("Gene contains credible variant that is protein-relevant with high predicted deleteriousness (CADD-Phred ",HTML("<span>&#8805;</span>"),"15)"),  selected = c(11:13),
                                                               inline = FALSE, width = NULL, choiceNames = c("stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other functional consequence"), choiceValues = c(11:13))),
                           
                           div(class="eqtl",checkboxGroupInput(inputId="columns2", label="Gene maps to credible variant that modulates gene expression (eQTLs) or splicing (sQTLs), FDR < 5%",  selected = c(14:19,22),
                                                               inline = FALSE, width = NULL, choiceNames = c("eQTL in glomerular tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL in tubulo-interstitial tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL in kidney cortex tissue (GTEx)","eQTL in other tissue (GTEx)","sQTL in kidney cortex tissue (GTEx)","sQTL in other tissue (GTEx)", "Colocalization of gene-expression signal and eGFRcrea signal (NEPTUNE)"), choiceValues = c(14:19,22))),
                           
                           div(class="pheno",checkboxGroupInput(inputId="columns3", label="Gene has known kidney phenotype in mouse or human:",  selected = c(20,21),
                                                                inline = FALSE, width = NULL, choiceNames = c("mouse (Mouse Genome Informatics, MGI)","human (online mendelian inheritance in man, OMIM, or Groopman et al. [N.Engl.J.Med,2019])"), choiceValues = c(20,21))),
                           tags$br(),
                           tags$h4("Specification of additional information",actionButton(inputId="question_details_gene",label="?",class="question")),
                           
                           
                           checkboxGroupInput(inputId="detail_evidence_locus_gene", label=NULL,  selected = NULL,
                                              inline = FALSE, width = NULL, choiceNames = c("Locus information for loci containing the searched genes", "List of all credible variants in the locus containing the searched gene"), choiceValues = c("summary","cred_var")),
                           
                           div(class="genelocuszoom",checkboxGroupInput(inputId="gene.locuszoom", label=NULL,  selected = NULL,
                                                                          inline = FALSE, width = NULL, choiceNames = c("Regional association plot of the eGFRcrea locus containing the searched gene (LocusZoom)"), choiceValues = c("locus_zoom")),
                                 
                           ),
                           
                           tags$br(),
                           tags$h5(tags$b("Restrict results to credible variants with a minimum PPA* of:")),
                           numericInput(inputId="ppa", label=NULL, value =0, min=0, max=1, width="20%"),
                           tags$p("*PPA: posterior probability of the eGFRcrea association  (max. PPA 1.0)"),
                           
                           
                           
                         )            
                      )
                      ),
                      
                      
                      
                      useShinyjs(),
                      verbatimTextOutput("testforme"),
                      shinyjs::hidden(        
                        div(id="Gene_div",        
                            tags$p(textOutput("namecheck_gene")),
                            
                            
                            shinyjs::hidden(
                              div( id ="GPS_gene",
                                   tags$h4("GPS Table", class="Output_header"),
                                   DT::dataTableOutput("GPSrows"),
                                   br(),
                                  # downloadButton("downloadGPS", "Download"),
                              )
                            ),
                            hr(),
                            
                            #textOutput("extras_gene"),
                            
                            
                            shinyjs::hidden(div(id="div_cadd_gene",tags$h4(span("Protein-relevant variants with predicted deleteriousness (CADD-Phred",HTML("<span>&#8805;</span>") ,"15):"), class="Output_header"),
                                                textOutput("CADD_gene_text"),
                                                br(),
                                                DT::dataTableOutput("CADD_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_glo_gene",tags$h4("eQTLs in glomerular tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021]):", class="Output_header"),
                                                textOutput("glo_gene_text"),
                                                br(),
                                                DT::dataTableOutput("glo_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_tub_gene",tags$h4("eQTLs in tubulo-interstitial tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021]):", class="Output_header"),
                                                textOutput("tub_gene_text"),
                                                br(),
                                                DT::dataTableOutput("tub_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_GTEx_eqtl_gene",tags$h4("eQTLs in kidney-cortex tissue (GTEx):", class="Output_header"),
                                                textOutput("eqtl_gene_text"),
                                                br(),
                                                DT::dataTableOutput("eqtl_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_GTEx_sqtl_gene",tags$h4("sQTLs in kidney-cortex tissue (GTEx):", class="Output_header"),
                                                textOutput("sqtl_gene_text"),
                                                br(),
                                                DT::dataTableOutput("sqtl_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_GTEx_wo_kidney_eqtl_gene",tags$h4("eQTLs in other tissues (GTEx):", class="Output_header"),
                                                textOutput("eqtl_wo_kidney_gene_text"),
                                                br(),
                                                DT::dataTableOutput("eqtl_wo_kidney_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_GTEx_wo_kidney_sqtl_gene",tags$h4("sQTLs in other tissues (GTEx):", class="Output_header"),
                                                textOutput("sqtl_wo_kidney_gene_text"),
                                                br(),
                                                DT::dataTableOutput("sqtl_wo_kidney_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_mgi_gene",tags$h4("Kidney phenotypes in mouse:", class="Output_header"),
                                                textOutput("mgi_gene_text"),
                                                br(),
                                                DT::dataTableOutput("mgi_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_omim_gene",tags$h4("Kidney phenotypes in human:", class="Output_header"),
                                                textOutput("omim_gene_text"),
                                                br(),
                                                DT::dataTableOutput("omim_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_summary",tags$h4("Summary of loci and signals containing the searched gene(s):", class="Output_header"),
                                                textOutput("summary_text"),
                                                br(),
                                                DT::dataTableOutput("summary"),
                                                hr())),
                            shinyjs::hidden(div(id="div_cred_var_gene",tags$h4("List of all credible variants in the locus containing the searched gene:", class="Output_header"),
                                                p("Shown are results from GWAS meta-analyses in European ancestry."),
                                                br(),
                                                DT::dataTableOutput("cred_var_gene"),
                                                hr())),
                            shinyjs::hidden(div(id="div_locus_zoom_gene",tags$h4("Locus Zoom Plot:", class="Output_header"),
                                                textOutput("LocusZoom_gene_text"),
                                                uiOutput("plot_download_links_gene"),
                                                
                                                )),
                            
                            
                        )),
             ),
             tabPanel("Variants",
                      p(tags$h3("SNP Search",actionButton(inputId="question_SNP",label="?",class="question"),style="display:inline")),
                      p("Search for SNPs which are associated with log(eGFRcrea) at P<5E-8 (all ancestries, unconditioned, n=1,201,909) or for SNPs being a credible variant for any of the 634 signals (European ancestry, conditioned, n=1,004,040):"),
                      useShinyjs(),
                      sidebarLayout(
                      
                      sidebarPanel(
                        
                        textInput(inputId="snp", label="Single SNP search", value=NULL, placeholder = "rs12345"),
                        # tags$h4(actionButton(inputId="question_SNP_batch",label="?",class="question"),"Upload SNP list:"),
                        # tags$br(),
                        tags$br(),
                        
                        div(class="SNP-batch", textAreaInput("SNP_batch",label="Paste a list of RSIDs",value=NULL, placeholder= "rs1234, rs4567")),
                        
                        tags$hr(),
                        # Input: Select a file ----
                        fluidRow(
                          
                          column(12,
                                 div(id="snp-upload-section",
                                 fileInput("file1", "Choose txt file with a list of RSIDs (max. 2000)",
                                           multiple = FALSE,
                                           accept = c("text/txt",
                                                      "text/comma-separated-values,text/plain",
                                                      ".txt")),
                                 actionButton(inputId="reset_variant", label="reset", class="go")
                          ))
                          ,
                        ),
                        # Horizontal line ----
                        tags$hr(),
                        shinyjs::hidden(
                          div(id="SNP_batch_upload_options",
                              # Input: Checkbox if file has header ----not reqired if just a list should be uploaded
                              #checkboxInput("header", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep_snp", "Separator",
                                           choices = c(Comma = ",",Semicolon = ";",Tab = "\t",Space=" "),selected = ","),
                              
                              # Input: Select quotes ----
                              # radioButtons("quote_snp", "Quote",
                              #              choices = c(None = "",
                              #                          "Double Quote" = '"',
                              #                          "Single Quote" = "'"),
                              #              selected = '"'),
                              
                              
                              
                          )
                        ),
                        fluidRow(
                          column(12,shinyjs::hidden(div(id="snp_go_div",actionButton(inputId="snp_go", label="go", class="go"))))
                        ),
                        
                      
                      ),
                      
                      mainPanel(
                      
                        wellPanel(
                          div(class="snp-search-options",
                          h4("Search options:"),
                        
                          div(class="gwscredinputbox",checkboxGroupInput(inputId="variant-search-options", label=NULL,  selected = c("gws","credvar"),
                                                                   inline = FALSE, width = NULL, choiceNames = c("search for SNPs at P<5E-8", "search for credible variants"), choiceValues = c("gws","credvar")),
                          ),
                          uiOutput("detail.evidence.snp"),
                          
                          
                          
                          div(class="snplocuszoom",checkboxGroupInput(inputId="snp.locuszoom", label="Specification of locus-based information",  selected = NULL,
                                                                      inline = FALSE, width = NULL, choiceNames = c("Regional association plot of the eGFRcrea locus containing the searched SNP"), choiceValues = c("locus_zoom")),
                              
                          ),
                          
                        ))
                      )
                      ),
                      
                      shinyjs::hidden(
                        div(id="SNP_div",
                            textOutput("snpinputcheck"),
                            div(id="gws_div",tags$h4("eGFRcrea Association (P<5E-8 in all-ancestry GWAS meta-analyis, unconditioned):", class="Output_header"),
                            textOutput("namecheck"),
                            DT::dataTableOutput("all_var"),
                            tags$hr()),
                            div(id="cred_var_div",tags$h4("Credible variant(s):", class="Output_header"),   
                            textOutput("namecheck2"),
                            DT::dataTableOutput("cred_var"),
                            tags$hr()),
                            shinyjs::hidden(div(id="div_cadd_snp",tags$h4(span("Protein-relevance and predicted deleteriousness (CADD-Phred",HTML("<span>&#8805;</span>") ,"15):"), class="Output_header"),
                                                #textOutput("extras"),
                                                textOutput("CADD_text"),
                                                DT::dataTableOutput("CADD"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_glo_snp",tags$h4("eQTL in glomerular tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021]):", class="Output_header"),
                                                textOutput("glo_text"),
                                                DT::dataTableOutput("glo"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_tub_snp",tags$h4("eQTL in tubulo-interstitial tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021]):", class="Output_header"),
                                                textOutput("tub_text"),
                                                DT::dataTableOutput("tub"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_GTEx_eqtl_snp",tags$h4("eQTL in kidney-cortex tissue (GTEx):", class="Output_header"),
                                                textOutput("eqtl_text"),
                                                DT::dataTableOutput("eqtl"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_GTEx_sqtl_snp",tags$h4("sQTL in kidney-cortex tissue (GTEx):", class="Output_header"),
                                                textOutput("sqtl_text"),
                                                DT::dataTableOutput("sqtl"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_GTEx_eqtl_wo_kidney_snp",tags$h4("eQTL in other tissues (GTEx):", class="Output_header"),
                                                textOutput("eqtl_wo_kidney_text"),
                                                DT::dataTableOutput("eqtl_wo_kidney"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_GTEx_sqtl_wo_kidney_snp",tags$h4("sQTL in other tissues (GTEx):", class="Output_header"),
                                                textOutput("sqtl_wo_kidney_text"),
                                                DT::dataTableOutput("sqtl_wo_kidney"),
                                                tags$hr())),
                            shinyjs::hidden(div(id="div_locus_zoom_snp",tags$h4("Locus Zoom Plot:", class="Output_header"),
                                                textOutput("LocusZoom_SNP_text"),
                                                uiOutput("plot_download_links_snp"),
                                            ))
                            
                        ))
             ),
             tabPanel("GPS",
                      
                      p(tags$h3("Gene Prioritisation - Overview",actionButton(inputId="question_region",label="?",class="question"))),
                      
                      tabsetPanel(
                        tabPanel("no PPA-filter",
                                 div(class="overflow",DT::dataTableOutput("GPS"))
                                 ),
                        tabPanel("PPA > 10%",
                                 div(class="overflow",DT::dataTableOutput("GPS_10"))
                        ),
                        tabPanel("PPA > 50%",
                                 div(class="overflow",DT::dataTableOutput("GPS_50"))
                        )
                      )
                      #tags$br(),
                      
                      
             ),
             
             tabPanel("Region",
                      useShinyjs(),
                      tags$h3("Region search",actionButton(inputId="question_region_search",label="?",class="question")),
                      wellPanel(
                        fluidRow(
                          column(4, 
                                 selectInput("Chromosome", label = "Chromosome", 
                                             choices = list("Chromosome 1" = 1, "Chromosome 2" = 2, "Chromosome 3" = 3,"Chromosome 4" = 4,"Chromosome 5" = 5,"Chromosome 6" = 6,"Chromosome 7" = 7,"Chromosome 8" = 8,"Chromosome 9" = 9,"Chromosome 10" = 10,"Chromosome 11" = 11,"Chromosome 12" = 12,
                                                            "Chromosome 13" = 13, "Chromosome 14" = 14,"Chromosome 15" = 15,"Chromosome 16" = 16,"Chromosome 17" = 17,"Chromosome 18" = 18,"Chromosome 19" = 19,"Chromosome 20" = 20,"Chromosome 21" = 21,"Chromosome 22" = 22), 
                                             selected = 1)
                          ),
                          column(8,
                                 numericInput(inputId="pos_start", label="Start of region [bp]", value = 0),
                                 numericInput(inputId="pos_end", label="End of region [bp]", value = 0),
                                 shinyjs::hidden(div(id="region_go_div",p(actionButton(inputId="region_go", label="go", class="go"),class="go"))),
                          ),
                        )
                      ),
                      #verbatimTextOutput("value"),
                      shinyjs::hidden(
                        div(id="region_div",
                            span(textOutput("valid_gene_test")),
                            bsCollapse(id="regionCollapse",multiple=TRUE,
                                 bsCollapsePanel("Overlapping genetic eGFRcrea signals",span(textOutput("ov_loci_test"),style="color:red"),
                                                 DT::dataTableOutput("ov_loci")
                                                 ),
                                 bsCollapsePanel("GPS entries for genes in overlapping eGFRcrea loci", 
                                                 shinyjs::hidden(
                                                   div(id="ov_genes_down",
                                                       DT::dataTableOutput("GPS_region_search"),
                                                       )
                                                 )
                                 )
                                 # ,
                                 # bsCollapsePanel("Locus Zoom","hier mitunter", uiOutput("LocusZoom2"), textOutput("hier"))
                      ))),
             ),
             
						tabPanel("Documentation & Help",
						         h3("Data Sources:"),
						         br(),
						         h4("Association with eGFRcrea:"),
						         p("Selection of genetic loci and therefore the genes within those loci, which are presented here, is based on a ",tags$abbr(title="genome-wide association studies","GWAS")," meta-analysis for eGFRcrea of ", tags$a(href="https://www.ukbiobank.ac.uk/", target="_blank", "UK Biobank data"), "and ", tags$a(href="https://ckdgen.imbi.uni-freiburg.de", target="_blank", "CKDGen consortium"), "data", tags$b("(n=1,201,909)."),"
						           Detailed information on the selection process can be found ",tags$a(href="https://www.nature.com/articles/s41467-021-24491-0", target="_blank", "here.")," A GWAS-meta-analysis of the same data restricted to European-ancestry (n=1,004,040) was used to identify independent
						           association signals and to calculate posterior probabilities of association (PPA) for all variants for each signal. The 99% credible variant sets of each signal contain the variants with the highest PPA with the sum of the PPA of all credible variants for a signal being larger than 99%. 
						           These are the variants with the highest probability to be causal for the association signal under the assumption that there is only one causal variant per association signal and that this variant is included in the analysis."),
						         br(),
						         h4("Association with other phenotypes:"),
						         p("GWAS meta-analyses were also performed for the other kidney function biomarkers eGFR estimated from serum cystatin C (eGFRcys, n=460,826) and blood urea nitrogen (BUN, n=852,678). The information if the locus lead variant (variant with the smallest association p-value in a locus) is also nominal significantly associated with these additional phenotypes can be found in each GPS table.
						           Summary statistics of these analyses can be downloaded ", tags$a(href="https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/gwas-summary-statistics/index.html", target="_blank", "here.")),
						         br(),
						         h4("CADD:"),
						         p("The combined annotation dependend depletion (CADD) score is a measurement of the deleteriousness of a genetic variant. By integrating multiple annotations, it contrasts variants that survived natural selection with simulated mutations. 
						           CADD evaluated ~8.6 billion SNPs and the CADD-Phred Score used on this website represents the rank of variant compared to all annotated variants. Therefore, the CADD-Phred Score", HTML("<span>&#8805;</span>"),"15 cut-off restricts our analysis to the 3.2% most deleterious variants.
						           Further, the analysis is restricted to variants within the affected gene as overlap with eQTLs and sQTLs should be minimized to avoid overscoring particular genes and variants. For additional information regarding CADD, please vistit the ", tags$a(href="https://cadd.gs.washington.edu/info", target="_blank", "CADD website.")),
						         p("Used version: v1.6 [2020-03-23]"),
						         br(),
						         h4("eQTL and sQTL data:"),
						         p("All credible variants were searched in expression quantitative trait loci (eQTL) databases. Three sources for eQTL data were used:"),
						         h5("NEPTUNE",style="font-weight:bold"),
						         p("eQTL data from the NEPTUNE study includes ",tags$em("cis-"),"eQTLs, which are variants that influence expression of genes within a 1Mb region centred around the variant. The association between a variant and the expression of a gene was deemed to be significant if the false dicovery rate (FDR) was <0.05.
						           This eQTL data was obtained from glomerular and tubulo-interstitial tissue. Further information about the NEPTUNE study can be found on the webpage of the ", tags$a(href="https://www.neptune-study.org/about", target="_blank", "study"), " and on the "
						           , tags$a(href="http://nephqtl.org/about", target="_blank", "NephQTL browser.")),
						         p("Version from [2017-09-25]"),
						         h5("Susztaklab (Sheng et al)", style="font-weight:bold"),
						         p("The Susztaklab also provides comprehensive kidney omics data. We integrated the eQTL data from glomerular und tubulo-interstitial tissue published by Sheng et al.", tags$a(href="https://www.nature.com/articles/s41588-021-00909-9", target="_blank", " (Sheng, X. et al., Nature Genetics, 2021).")),
						         h5("GTEx", style="font-weight:bold"),
						         p("In contrast to the other two eQTL sources, the GTEx project is not restricted to kidney tissue. Furthermore, additional splicing altering variants (sQTLs) were investigated. Thus, the here integrated GTEx data includes ",tags$em("cis-"),"eQTL and -sQTL information from 48 different tissues with a mapping window of 1Mb up- and downstream of the transcription start site.
						            Further information about GTEx can be found ",tags$a(href="https://www.gtexportal.org/home/", target="_blank", "here.")),
						         p("Used version: GTEx Release v7 [2017-09-05]"),
						         br(),
						         h4("Mouse phenotypes:"),
						         p("Information on genes with kidney-relevant phenotypes in mice origin from the Mouse Genome Informatics database (MGI). This includes all phenotypes subordinate to \"abnormal kidney morphology\" (MP:0002135) and \"abnormal kidney physiology\" (MP:0002136).
						            Further information how this data was collected can be found on the ", tags$a(href="http://www.informatics.jax.org/", target="_blank", "MGI webpage.")),
						         p("Version from [2020-06-03]"),
						         br(),
						         h4("Human phenotypes:"),
						         p("We used two sources to identify genes causing genetic disorders with kidney phenotype in human:"),
						         h5("OMIM", style="font-weight:bold"),
						         p("The Online Mendelian Inheritance in Man (OMIM) database was queried for phenotype entries subordinate to the clinical synopsis class \"kidney\". Diseases with \"kidney\"-phenotype entries being: \"normal kidneys\", \"normal renal ultrasound at ages 4 and 7 (in two family)\", \"no kidney disease\", \"no renal disease; normal renal function\", \"normal renal function; no kidney disease\" and \"no renal findings\" were manually excluded.
						           Be aware that OMIM entries missing a clinical synopsis entry are not included in kidneyGPS regardless of a potential kidney involvement. Further information on the diseases can be found at the ",tags$a(href="https://www.omim.org/", target="_blank", "OMIM webpage.")),
						         p("Version from [2020-08-07]"),
						         h5("Groopman et al.", style="font-weight:bold"),
						         p(tags$a(href="(http://www.columbiamedicine.org/divisions/gharavi/files/Kidney_Gene_List_625.xlsx", target="_blank","A list of 625 genes "), "associated with Mendelian forms of kidney and genitourinary disease was published by Groopman et al. in 2019 in the New England Journal of Medicine. The original article \"Diagnostic Utility of Exome Sequencing for Kidney Disease\" can be found ",tags$a(href="https://www.nejm.org/doi/10.1056/NEJMoa1806891", target="_blank", "here.")," 
						           Please notice that not all 625 genes are included in any eGFRcrea locus and thus cannot be found in kidneyGPS."),
						         tags$hr(),
						         h3("Citation"),
						         p("We encourage the use and publication of results generated from these data. We request that if you use data from KidneyGPS browser cite that you cite our paper, \"KidneyGPS: an easily accessible web application to prioritize kidney function genes and variants based on evidence from genome-wide association studies\" currently published in preprint and and the original data sources listed on this page."),
						         tags$hr(),
						         h3("Contact"),
						         p("If you have any questions not answered by this page, please contact: Kira Stanzick ", tags$a(href="mailto:kira-julia.stanzick@klinik.uni-regensburg.de", "(kira-julia.stanzick@klinik.uni-regensburg.de)") ),
						         tags$hr(),
						         h3("Release history"),
						         p("KidneyGPS 1.1.2 [2022-06] - small layout changes"),
						         p("KidneyGPS 1.1.1 [2022-04] - minor fixes"),
						         p("KidneyGPS 1.1 [2022-03]- integration of eQTL data from Susztaklab"),
						         p("first publication KidneyGPS 1.0 [2022-03]")
						         )
	)
							
  
)



###############################################
#     server part
#
#
###############################################
server <- function(input, output, session) {
  
#overall:
  sketchcadd = htmltools::withTags(table(
    class = 'display',
    
    thead(
      tr(
        th('RSID'),
        th(class="PPA hover", 'PPA', p(class="PPAcomment", "Posterior probability of association of the SNP (probability of driving the association of the respective signal)")),
        th('Affected gene'),
        th('Locus ID'),
        th(class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
        th('Functional consequence'),
        th(class="REF hover",'Reference allele', span(class="REFcomment", "Allele for reference aminoacid")),
        th(class="ALT hover",'Alternative allele', span(class="ALTcomment", "Allele for altered aminoacid")),
        th('Aminoacid position'),
        th('Reference aminoacid'),
        th('Alternative aminoacid'),
        th(class="CADDscore hover",'CADD Phred-Score',span(class="CADDscorecomment", "Scaled score for the deleteriousness of this allelic exchange"))
      )
    )
  ))
  
  sketchcredvar = htmltools::withTags(table(
    class = 'display',
    
    thead(
      tr(
        th( rowspan=2, 'RSID*'),
        th(rowspan=2,'Locus ID'),
        th(rowspan=2,class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
        th(rowspan=2,class="PPA hover", 'PPA', p(class="PPAcomment", "Posterior probability of association of the SNP (probability of driving the association of the respective signal)")),
        th(rowspan=2,'N'),
        th(rowspan=2,'Effect allele'),
        th(rowspan=2,'Other allele'),
        th(rowspan=2,'EAF'),
        th(colspan=3, 'Association with log(eGFRcrea) conditioned on other signal index variants'),
        th(colspan=3, 'Association with log(eGFRcrea) unconditioned'),
        th(rowspan=2,'Chr'),
        th(rowspan=2, class="pos hover",'Pos',span(class="poscomment","SNP position from GRCh37")),
        ),
      tr(
        lapply(c('beta', 'StdErr', 'P-value','beta', 'StdErr', 'P-value'), th)
      )
    )
  ))
  
  sketchcredvargene = htmltools::withTags(table(
    class = 'display',
    
    thead(
      tr(
        th( rowspan=2, 'RSID*'),
        th(rowspan=2,'Locus ID'),
        th(rowspan=2,class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
        th(rowspan=2,class="PPA hover", 'PPA', p(class="PPAcomment", "Posterior probability of association of the SNP (probability of driving the association of the respective signal)")),
        th(rowspan=2, 'Functional impact on the searched gene'),
        th(rowspan=2,'N'),
        th(rowspan=2,'Effect allele'),
        th(rowspan=2,'Other allele'),
        th(rowspan=2,'EAF'),
        th(colspan=3, 'Association with log(eGFRcrea) conditioned on other signal index variants'),
        th(colspan=3, 'Association with log(eGFRcrea) unconditioned'),
        th(rowspan=2,'Chr'),
        th(rowspan=2, class="pos hover",'Pos',span(class="poscomment","SNP position from GRCh37")),
      ),
      tr(
        lapply(c('beta', 'StdErr', 'P-value','beta', 'StdErr', 'P-value'), th)
      )
    )
  ))
  
  sketchallsig = htmltools::withTags(table(
    class = 'display',
    
    thead(
      tr(
        th( rowspan=2, 'RSID*'),
        th(rowspan=2,'Chr'),
        th(rowspan=2,'Pos'),
        th(rowspan=2,'Effect allele'),
        th(rowspan=2,'Other allele'),
        th(rowspan=2,'N'),
        th(rowspan=2,'EAF'),
        th(colspan=3,'Association with log(eGFRcrea) unconditioned**'),
        th(rowspan=2,'Nearest gene'),
        th(rowspan=2,'Distance to nearest gene'),
        th(rowspan=2,'Locus ID')
      ),
      tr(
        lapply(c('beta', 'StdErr', 'P-value'), th)
      )
    )
  ))

# GPS ---------------------------------------------------------------------

  
  sketch = htmltools::withTags(table(
    class = 'display',
    
    thead(
      tr(
        th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
        th(rowspan = 2, 'Locus name**'),
        th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
        th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
        th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
        th(rowspan = 2, '# credible variants in signal'),
        th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
        th(rowspan = 2, 'overall gene score'),
        th(colspan = 3, class="CADD", span('# protein-relevant credible variants in the gene ')),
        th(colspan = 6, class="eqtl", '# credible variants that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
        th(rowspan = 2, class="eqtl", 'Coloc in NEPTUNE tissue'),
        th(rowspan = 2, class="pheno",'# kidney phenotypes in mouse'),
        th(rowspan = 2, class="pheno", '# kidney phenotypes in human'),
        th(rowspan = 2, class="dm", 'Diabetes specific effect'),
        th(rowspan = 2, class="decline", 'eGFRcrea decline effect'),
        th(rowspan = 2, class="drug", 'Drugable'),
        th(rowspan = 2, class="pathway", 'In enriched pathway'),
        th(rowspan = 2, 'Distance to locus lead variant'),
        th(rowspan = 2, 'Chromosome'),
        th(rowspan = 2, 'Position gene start'),
        th(rowspan = 2, 'Position gene end')
      ),
      tr(
        lapply(c('stop-gained, stop-lost, non-synonymous','canonical splice, noncoding change, synonymous, splice site'),th),
               th('other deleterious variant (CADD-Phred',HTML("<span>&#8805;</span>"),'15)'),
        lapply(c('eQTL glomerulus (NEPTUNE, or Sheng et al [Nat.Genet. 2021])','eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat.Genet. 2021])','eQTL kidney cortex (GTEx)','eQTL other tissue (GTEx)','sQTL kidney cortex (GTEx)','sQTL other tissue (GTEx)'), th)
      )
    )
  ))
 
  headerCallback <- "function( thead, data, start, end, display ) {
  $(thead).closest('thead').find('th').eq(8).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(9).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(10).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(11).css('background-color', '#b3ffb3');
  $(thead).closest('thead').find('th').eq(12).css('background-color', '#b3ffb3');
  $(thead).closest('thead').find('th').eq(21).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(22).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(23).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(24).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(25).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(26).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(27).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(28).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(29).css('background-color', '#ffdb99');
  }"
  
  output$GPS <- renderDataTable({
    
    
    datatable(
      GPS_show[,c(4,1:3,9,10,24,23,11:19,22,20,21,25:28,5:8)], rownames=FALSE,

        options = list(
                    columnDefs = list(list(className = 'dt-left', targets = "_all"))
                    , scrollX=TRUE
                    ,headerCallback = JS(headerCallback),
                     dom = 'rtBip',
                    buttons= list(c('copy', 'excel'))

        ),
      extensions= 'Buttons',
      class ='display',
      filter = list(
        position = 'top', clear = TRUE, plain = TRUE
      ),
      escape = FALSE,

      caption = htmltools::tags$caption(
        style = 'caption-side: bottom; text-align: left; font-weight: normal;',
        '* Link to GeneCards, **nearest gene to locus lead-variant'),
      container = sketch
    )
  })
  output$GPS_10 <- renderDataTable({
    
    
    datatable(
      GPS_10_show[,c(4,1:3,9,10,24,23,11:19,22,20,21,25:28,5:8)], rownames=FALSE,
      
      options = list(
        columnDefs = list(list(className = 'dt-left', targets = "_all"))
        , scrollX=TRUE
        ,headerCallback = JS(headerCallback),
        dom = 'rtBip',
        buttons= list(c('copy', 'excel'))
        
      ),
      extensions= 'Buttons',
      class ='display',
      filter = list(
        position = 'top', clear = TRUE, plain = TRUE
      ),
      escape = FALSE,
      
      caption = htmltools::tags$caption(
        style = 'caption-side: bottom; text-align: left; font-weight: normal;',
        '* Link to GeneCards, **nearest gene to locus lead-variant'),
      container = sketch
    )
  })
  output$GPS_50 <- renderDataTable({
    
    
    datatable(
      GPS_50_show[,c(4,1:3,9,10,24,23,11:19,22,20,21,25:28,5:8)], rownames=FALSE,
      
      options = list(
        columnDefs = list(list(className = 'dt-left', targets = "_all"))
        , scrollX=TRUE
        ,headerCallback = JS(headerCallback),
        dom = 'rtBip',
        buttons= list(c('copy', 'excel'))
        
      ),
      extensions= 'Buttons',
      class ='display',
      filter = list(
        position = 'top', clear = TRUE, plain = TRUE
      ),
      escape = FALSE,
      
      caption = htmltools::tags$caption(
        style = 'caption-side: bottom; text-align: left; font-weight: normal;',
        '* Link to GeneCards, **nearest gene to locus lead-variant'),
      container = sketch
    )
  })
                                

# region page -------------------------------------------------------------
  
  sketch2 = htmltools::withTags(table(
    class = 'display',
    
    thead(
      tr(
        th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
        th(rowspan = 2, 'Locus name**'),
        th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
        th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
        th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
        th(rowspan = 2, '# credible variants in signal'),
        th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
        th(colspan = 3, class="CADD", span('# protein-relevant credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
        th(colspan = 6, class="eqtl", '# credible variants that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
        th(rowspan = 2, class="eqtl", 'Coloc in NEPTUNE tissue'),
        th(rowspan = 2, class="pheno",'# kidney phenotypes in mouse'),
        th(rowspan = 2, class="pheno", '# kidney phenotypes in human'),
        th(rowspan = 2, class="dm", 'Diabetes specific effect'),
        th(rowspan = 2, class="decline", 'eGFRcrea decline effect'),
        th(rowspan = 2, class="drug", 'Drugable'),
        th(rowspan = 2, class="pathway", 'In enriched pathway'),
        th(rowspan = 2, 'Distance to locus lead variant'),
        th(rowspan = 2, 'Chromosome'),
        th(rowspan = 2, 'Position gene start'),
        th(rowspan = 2, 'Position gene end')
      ),
      tr(
        lapply(c('stop-gained, stop-lost, non-synonymous','canonical splice, noncoding change, synonymous, splice site','other deleterious variant','eQTL glomerulus (NEPTUNE, or Sheng et al [Nat.Genet. 2021])','eQTL tubulo-interstitium (NEPTUNE, or Sheng et al [Nat.Genet. 2021])','eQTL kidney cortex (GTEx)','eQTL other tissue (GTEx)','sQTL kidney cortex (GTEx)','sQTL other tissue (GTEx)'), th)
      )
    )
  ))
  

  headerCallback2 <- "function( thead, data, start, end, display ) {
  $(thead).closest('thead').find('th').eq(7).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(8).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(9).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(10).css('background-color', '#b3ffb3');
  $(thead).closest('thead').find('th').eq(11).css('background-color', '#b3ffb3');
  $(thead).closest('thead').find('th').eq(20).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(21).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(22).css('background-color', '#ccffff');
  $(thead).closest('thead').find('th').eq(23).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(24).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(25).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(26).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(27).css('background-color', '#ffdb99');
  $(thead).closest('thead').find('th').eq(28).css('background-color', '#ffdb99');
  }"
  
  
  #region page
  observe(
    if(input$pos_start!=0&input$pos_end!=0){
      shinyjs::showElement(id="region_go_div")
    }
  )
  
  #combine ENTER and "go" click to one reactive value only reacting if there is region input
  go <- reactiveValues('region'=c())
  
  observeEvent(input$region_go, {
    if(input$pos_start!=0&input$pos_end!=0){
      if(is.null(go$region)){go$region=1}
      else{go$region <- go$region+1}
      
    }
  })
  observeEvent(input$keyPressed, {
    if(input$pos_start!=0&input$pos_end!=0){
      if(is.null(go$region)){go$region=1}
      else{go$region <- go$region+3}
    }
  })
  
  region_start <- eventReactive(go$region,{input$pos_start})
  region_end <- eventReactive(go$region,{input$pos_end})
  
  valid_reg <- eventReactive( go$region,{(
    if(!is.na(region_start())& !is.na(region_end())){(is.integer(region_start()) & is.integer(region_end())& input$pos_start>0 & input$pos_end>0 & input$pos_end>input$pos_start)}
    else{FALSE}
    
  )}
  )
  
  chr <- eventReactive(go$region,{input$Chromosome})
  
  output$valid_gene_test <- renderText({
    if(!valid_reg()){
      if(is.na(region_start())){"Error: invalid region start position"}
      else{
        if(!is.integer(region_start())|region_start()<0) {"Error: invalid region start position"}
        else{
          if(is.na(region_end())){"Error: invalid region end position"}
          else{
            if(!is.integer(region_end())| region_end()<0) {"Error: invalid region end position"}
            else{
              if(region_end()<region_start()) {"Error: region start must be smaller than region end position"}
              else{"undefined error"}
            }
          }
        }
      }
    }
  })
  
  
  ov_loci<- reactive({
    if(valid_reg()){
      
      region_table_show[which(chr()==as.numeric(region_table$chr) & ((region_start()<=region_table$region_start_500kb & region_end()>=region_table$region_start_500kb) | (region_start()>=region_table$region_start_500kb & region_start()<=region_table$region_end_500kb))),]
    }else{0}
  })
  
  overlap_loci <- eventReactive(go$region,{nrow(ov_loci())!=0})
  
  output$ov_loci <- DT::renderDataTable({
    if(valid_reg() & overlap_loci()){
      datatable(ov_loci(),rownames=FALSE,
                options = list(
                  columnDefs = list(list(className = 'dt-left', targets = "_all"))
                  , scrollX=TRUE,
                  dom = 'lfrtBip',
                  buttons= list(c('copy', 'excel'))
                ),
                extensions = 'Buttons',
                caption = htmltools::tags$caption(
                  style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                  '* SNPs, which are associated with log(eGFRcrea) in the all-ancestry GWAS meta-analysis (unconditioned) with P<5E-8')
                
                )
    }
  })
  
  output$ov_loci_test <- renderText({
    
    if(valid_reg()& !overlap_loci()) {
      "no eGFRcrea associated loci found in this region"
    }else{
      if(valid_reg() & overlap_loci()){
        paste(length(unique(ov_loci()[,1]))," loci (including ",nrow(ov_loci())," signals) have been found in your searched region.", sep="")
      }
    }
  })
  
  GPS_region_search <- reactive({
    if(valid_reg()&overlap_loci()){
      x = GPS_show[which(GPS$locus_id%in%ov_loci()[,1]),c(4,1:3,9,10,24,11:19,22,20,21,25:28,5:8)]
      y = GPS[which(GPS$locus_id%in%ov_loci()[,1]),]
      v=c()
      for(i in 1:nrow(y)){
      v[i] = 11-length(which(y[i,c(9:19)]==0))
      }
      # x['Genescore']=v
      # x=x[,c(1:6,23,7:22)]
      x[order(v,decreasing=T),]
      
    }else{NULL}
  })
  
  output$GPS_region_search <- 
    DT::renderDataTable({
      if(valid_reg()& !is.null(GPS_region_search())){
        datatable(GPS_region_search(),rownames=FALSE,
                  
                  
                  options = list(
                    columnDefs = list(list(className = 'dt-left', targets = "_all"))
                    , scrollX=TRUE
                    ,headerCallback = JS(headerCallback2),
                    dom = 'lfrtBip',
                    buttons= list(c('copy', 'excel'))
                    
                  ),
                  extensions= 'Buttons',
                  class ='display',
                  filter = list(
                    position = 'top', clear = TRUE, plain = TRUE
                  ),
                  escape = FALSE,
                  
                  caption = htmltools::tags$caption(
                    style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                    '* Link to GeneCards, **nearest gene to locus lead-variant'),
                  container = sketch2
        )%>% formatPercentage('max PPA') 
      } 
    })
  
  observeEvent(go$region,{shinyjs::showElement(id="region_div")})
  observeEvent(go$region, ({
    
    if(valid_reg()){
      updateCollapse(session, "regionCollapse", open = c("Overlapping genetic eGFRcrea signals"))
      if(!is.null(GPS_region_search())){
          updateCollapse(session, "regionCollapse", open = c("Overlapping genetic eGFRcrea signals", "GPS entries for genes in overlapping eGFRcrea loci"))
        
      }
    }
  }))
  
  
  observe({
    
    shinyjs::toggleElement(id="ov_genes_down", condition = (!is.null(GPS_region_search())) )
    
  })
  
  
  
  observeEvent(input$question_region_search,{
    showModal(modalDialog(
      titel="Help region search",
      tags$h4("Region search - Help:"),
      p("Select a chromosome and enter a start/end position [bp] of a region you are interested in. The results show you a list of overlapping eGFRcrea loci and the GPS entries for the respective genes if there are any. Attation: The region end must be larger than region start and region start can't be 0. For further information check the Documentation page."),
      easyClose = TRUE,
      footer=NULL,
      size = "m"
    ))
  })

# snp page ----------------------------------------------------------------
  
  #define input type based on click for further analyses
  inputtype <- reactiveValues('snp'=c())
  onclick("snp", (inputtype$snp=c("single")), add=T)
  onclick("SNP_batch", (inputtype$snp=c("list")), add=T)
  onclick("snp-upload-section", (inputtype$snp=c("upload")), add=T)
  
  #show search button:
  
  observe(
    if(isTruthy(input$file1)|isTruthy(input$snp)|isTruthy(input$SNP_batch)){
      shinyjs::showElement(id="snp_go_div")
    }
  )
  
  
  #combine ENTER and "go" click to one reactive value only reacting if there is snp input
  go <- reactiveValues('snp'=c())
  
  observeEvent(input$snp_go, {
    if(isTruthy(input$file1)|isTruthy(input$snp)|isTruthy(input$SNP_batch)){
      if(is.null(go$snp)){go$snp=1}
      else{go$snp <- go$snp+1}
      
    }
  })
  observeEvent(input$keyPressed, {
    if(isTruthy(input$file1)|isTruthy(input$snp)|isTruthy(input$SNP_batch)){
      if(is.null(go$snp)){go$snp=1}
      else{go$snp <- go$snp+3}
    }
  })
  
  #upload options
  observe({
    shinyjs::toggleElement(id="SNP_batch_upload_options", condition = (isTruthy(input$file1)) )
  })
  
  #disable credvar only functions, when not searched for cred vars
  observe({
    
    if(any(input$`variant-search-options`=="credvar")){
    output$detail.evidence.snp <- renderUI(
    div(class="snpinputbox",checkboxGroupInput(inputId="detail_evidence_variant", label=span(class="snpsearch","Specification of additional functional SNP information",span(class="snpsearchcomment","functional information only available for credible variants")),  selected = c("CADD","NEPTUNE_glo","NEPTUNE_tub","GTEx_eQTL","GTEx_sQTL"),
                                               inline = FALSE, width = NULL, choiceNames = c("deleteriousness and functional annotation (CADD)","eQTL in glomerulus (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL in tubulo-interstitium (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL in kidney cortex (GTEx)","sQTL in kidney cortex (GTEx)","eQTL in other tissues (GTEx)","sQTL in other tissues (GTEx)"), choiceValues = c("CADD","NEPTUNE_glo","NEPTUNE_tub","GTEx_eQTL","GTEx_sQTL","GTEx_wo_kidney_eQTL","GTEx_wo_kidney_sQTL")),
    ))}else{
      output$detail.evidence.snp <- renderUI(
        div(class="snpinputbox",disabled(checkboxGroupInput(inputId="detail_evidence_variant", label=span(class="snpsearch","Specification of additional functional SNP information or locus based information",span(class="snpsearchcomment","functional information only available for credible variants")),
                                                   inline = FALSE, width = NULL, choiceNames = c("deleteriousness and functional annotation (CADD)","eQTL in glomerulus (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL in tubulo-interstitium (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL in kidney cortex (GTEx)","sQTL in kidney cortex (GTEx)","eQTL in other tissues (GTEx)","sQTL in other tissues (GTEx)"), choiceValues = c("CADD","NEPTUNE_glo","NEPTUNE_tub","GTEx_eQTL","GTEx_sQTL","GTEx_wo_kidney_eQTL","GTEx_wo_kidney_sQTL"))),
        ))
    }
  })
  
  ###first snp input processing: 
  
  single_snp <- eventReactive(go$snp, {input$snp})
  
  #build vector of snps from text input
  input_snp_batch_list <- eventReactive(go$snp,{
    
    req(input$SNP_batch)
    
    ivbl=unlist(strsplit(input$SNP_batch,"[^[:alnum:]^-]",perl=TRUE)) ## excludes all non-word characteres except "-" that might be used as separators
    if(any(ivbl=="")){
      ivbl[-which(ivbl=="")]
    }else{ivbl}
  })
  
  
  #build vector of snps from file input
  
  input_snp_batch_upload <- reactiveValues('input'= c())
  
  observe({
    req(input$file1)
    
    # df <- read.table(input$file1$datapath,
    #                  sep = input$sep_snp,
    #                  quote = input$quote_snp)
    # 
    # input_snp_batch_upload$input = paste(df[1,])
    input_snp_batch_upload$input <- scan(file=input$file1$datapath,
                                         sep = input$sep_snp,
                                         what="character")
    
  })
  
  
  #reset file input, when reset button is clicked
  observeEvent(input$reset_variant,{
    reset("file1")
    input_snp_batch_upload$input <- NULL
  })
  
  #final selection:
  snp <- eventReactive(go$snp, {
    
    if(inputtype$snp=="single"){
      single_snp()
    }else{
      if(inputtype$snp=="list"){
        input_snp_batch_list()
      }else{
        if(inputtype$snp=="upload"){
          input_snp_batch_upload$input
        }else{NULL}
      }
    }
    })
  
  inputtype <- reactiveValues('snp_final'=c()) #use inputtype$snp_final as it only reacts on "search"-click
  observeEvent(go$snp,{
    inputtype$snp_final = inputtype$snp 
  })
  
  #snp page	
  cred_vars <- eventReactive(go$snp, {
    if(!is.null(snp())){
      if(any(snp()%in%cred_var$rsid)){
        snp_in_cred=(snp()%in%cred_var$rsid)
        snp()[which(snp_in_cred)]
      }else{NULL}
      }else{NULL}
    })
  
  valid <- eventReactive(go$snp, {!is.null(cred_vars())}) #any of the searched variants is a credible variant
  
  detail_evidence_variant <- eventReactive(go$snp, {input$detail_evidence_variant})
  
  gws_snps <- eventReactive(go$snp,{
    if(!is.null(snp())){
      if(any(snp()%in%all_sig$RSID)){
      snp_in_gws = (snp()%in%all_sig$RSID)
      snp()[which(snp_in_gws)]
      }else{NULL}
      }else{NULL}
  })
  
  gws<-eventReactive(go$snp, {!is.null(gws_snps())}) # any of the searched variants is gws
  
  # open all relevant divs
  observeEvent(go$snp,{shinyjs::showElement(id="SNP_div")})
  observeEvent(go$snp,{runjs('
      document.getElementById("SNP_div").scrollIntoView({ left: 0, block: "end", behavior: "smooth" });
    ')})
  observeEvent(go$snp,{shinyjs::toggleElement(id="gws_div", condition = any(input$`variant-search-options`=="gws"))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="cred_var_div", condition = any(input$`variant-search-options`=="credvar"))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_cadd_snp", condition= (any(detail_evidence_variant()=="CADD")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_glo_snp", condition= (any(detail_evidence_variant()=="NEPTUNE_glo")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_tub_snp", condition= (any(detail_evidence_variant()=="NEPTUNE_tub")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_GTEx_eqtl_snp", condition= (any(detail_evidence_variant()=="GTEx_eQTL")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_GTEx_sqtl_snp", condition= (any(detail_evidence_variant()=="GTEx_sQTL")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_GTEx_eqtl_wo_kidney_snp", condition= (any(detail_evidence_variant()=="GTEx_wo_kidney_eQTL")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_GTEx_sqtl_wo_kidney_snp", condition= (any(detail_evidence_variant()=="GTEx_wo_kidney_sQTL")&valid()))})
  observeEvent(go$snp,{shinyjs::toggleElement(id="div_locus_zoom_snp", condition= (any(input$snp.locuszoom=="locus_zoom")&valid()))})
  
  
  #gws output
  output$namecheck <- renderText({
    if(!is.null(snp())){
      if(inputtype$snp_final=="single"){
        if(gws()){
          paste(snp(), "is associated with log(eGFRcrea) in at least one of the 424 loci at P<5E-8 (all-ancestry GWAS meta-analysis, unconditioned):", sep=" ")
        }else{
          paste(snp(), "is not included or not associated with log(eGFRcrea) in one of the 424 loci at P<5E-8 (all-ancestry GWAS meta-analysis, unconditioned).", sep=" ")
        }
      }else{
          if(gws()){
            if(length(gws_snps())>1){
              paste("Of your ",length(snp()), " queried SNPs, ", length(gws_snps()), " SNPs are associated with log(eGFRcrea) in at least one of the 424 loci at P<5E-8 (all-ancestry GWAS meta-analysis, unconditioned):", sep="")
            }else{
              paste("Of your ",length(snp()), " queried SNPs, ", length(gws_snps()), " SNP (", gws_snps(),") is associated with log(eGFRcrea) in at least one of the 424 loci at P<5E-8 (all-ancestry GWAS meta-analysis, unconditioned):", sep="") 
              }
          }else{
            paste("None of your ",length(snp()), " queried SNPs is associated with log(eGFRcrea) in one of the 424 loci at P<5E-8 (all-ancestry GWAS meta-analysis, unconditioned).", sep="") 
          }
        }
      
    }else{
        "Please check if you correctly entered or uploaded your SNPs and if the right input field is selected."
      }
  })
  
  output$all_var <- DT::renderDataTable({
    if(gws()){
      datatable(all_sig_show[which(all_sig$RSID%in%gws_snps()),c(15,11,12,2:8,22,23,24)], rownames=FALSE,
                options = list(
                  columnDefs = list(list(className = 'dt-left', targets = "_all"))
                  , scrollX=TRUE,
                  dom = 'lfrtBip',
                  buttons= list(c('copy', 'excel'))
                ),
                extensions = 'Buttons',
                escape=FALSE,
                caption = htmltools::tags$caption(
                  style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                  htmltools::tags$p('* Link to dbSNP'),
                  htmltools::tags$p('** adjusted for age, sex and other study-specific covariates'),
                  htmltools::tags$p('Effect allele: eGFRcrea lowering allele, N: Sample size of the all-ancestry GWAS meta-analysis for this SNP, EAF: Effect allele frequency'),
                  htmltools::tags$p('beta: Effectsize of the effect allele on log(eGFRcrea)')
                  ),
                class='display',
                container = sketchallsig
                )
    }
  })
  
  #cred_var ouput
  output$namecheck2 <- renderText({
    if(!is.null(snp())){
      if(inputtype$snp_final=="single"){
        if(valid()){
          paste(snp(), "is a 99% credible variant in at least one of the 634 eGFRcrea signals. Shown are results from GWAS meta-analyses in European ancestry:", sep=" ")
        }else{
          paste(snp(), "is not a 99% credible variant in any of the 634 eGFRcrea signals.", sep=" ")
        }
      }else{
        if(gws()){
          if(length(cred_vars())>1){
            paste("Of your ",length(snp()), " queried SNPs, ", length(cred_vars()), " SNPs are 99% credible variants in at least one of the 634 eGFRcrea signals. Shown are results from GWAS meta-analyses in European ancestry:", sep="")
          }else{
            paste("Of your ",length(snp()), " queried SNPs, ", length(cred_vars()), " SNP (", cred_vars(),") is a 99% credible variant in at least one of the 634 eGFRcrea signals. Shown are results from GWAS meta-analyses in European ancestry:", sep="") 
          }
        }else{
          paste("None of your ",length(snp()), " queried SNPs is a 99% credible variant in any of the 634 eGFRcrea signals.", sep="") 
        }
      }
      
    }
  })
  
  output$cred_var <- DT::renderDataTable({
    if(valid()){
      datatable(cred_var_show[which(cred_var$rsid%in%cred_vars()),], rownames=FALSE,
                options = list(
                  columnDefs = list(list(className = 'dt-left', targets = "_all"))
                  , scrollX=TRUE,
                  dom = 'lfrtBip',
                  buttons= list(c('copy', 'excel'))
                ),
                extensions = 'Buttons',
                escape=FALSE,
                caption = htmltools::tags$caption(
                  style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                  htmltools::tags$p('* Link to dbSNP'),
                  htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal)'),
                  htmltools::tags$p('N: Sample size of the European-ancestry GWAS meta-analysis for this SNP, Effect allele: eGFRcrea lowering allele, EAF: Effect allele frequency'),
                  htmltools::tags$p('Unconditioned association statistics origin from a european ancestry only GWAS meta-analysis for eGFRcrea. For loci with multiple independent signals, the conditioned association statistics origin from conditioning on the signal index variants of the other signals in the locus.')
                  ),
                class='display',
                container = sketchcredvar
                )%>% formatPercentage('PPA')
    }
  })
  
  output$extras <-  renderText({
    if(valid()){
      if(length(detail_evidence_variant())!=0){
        
          paste("The following functional evidence exists for the ",length(cred_vars())," 99% credible variant(s) included in your search:",sep="")
        
      }
    }else{
      "Your request didn't include any 99% credible variant so no search for functional evidence can be performed."
    }
  })
  
  
  #output CADD
  output$CADD_text <- renderText({
    
      if(any(detail_evidence_variant()=="CADD")){
        if(any(CADD$RSID%in%cred_vars())){
          if(length(cred_vars())==1){
            paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%CADD$RSID))," SNP resides in a gene and is predicted deleterious:",sep="")
          }else{
            paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%CADD$RSID))," SNP(s) reside(s) in a gene and is/are predicted deleterious:",sep="")  
          }
        }else{
          if(length(cred_vars())==1){
            paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not reside in a gene or is not predicted deleterious.",sep="")
          }else{
            paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, resides in a gene and is predicted deleterious.",sep="")  
          }
        }
      }
  })
  
  output$CADD <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="CADD")){
        if(any(CADD$RSID%in%cred_vars())){
          x = CADD_show[which(CADD$RSID%in%cred_vars()),]
          if(any(duplicated(x))) x = x[-which(duplicated(x)),]
          datatable(x, rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal)'),
                      htmltools::tags$p(' Reference and alternative allele from human reference genome (GRCh37)')
                    ),
                    class ='display',
                    container = sketchcadd
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  # Neptune + Susztak glomerulus
  output$glo_text <- renderText({
    
      if(any(detail_evidence_variant()=="NEPTUNE_glo")){
        if(any(glo_show$RSID%in%cred_vars())){
          if(length(cred_vars())==1){
            paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%glo_show$RSID))," SNP modulates gene expression in glomerular tissue:",sep="")
          }else{
            paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%glo_show$RSID))," SNP(s) modulate(s) gene expression in glomerular tissue:",sep="")  
          }
        }else{
          if(length(cred_vars())==1){
            paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not modulate gene expression in glomerular tissue.",sep="")
          }else{
            paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, modulates gene expression in glomerular tissue.",sep="")  
          }
        }
      }
    
  })
  
  output$glo <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="NEPTUNE_glo")){
        if(any(glo_show$RSID%in%cred_vars())){
          datatable(glo_show[which(glo_show$RSID%in%cred_vars()),], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                      htmltools::tags$p('** NEPTUNE, or Sheng et al, Nat.Genet. 2021')
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  #Neptune tubulo-interstitium
  output$tub_text <- renderText({
    
      if(any(detail_evidence_variant()=="NEPTUNE_tub")){
        if(any(tub_show$RSID%in%cred_vars())){
          if(length(cred_vars())==1){
            paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%tub_show$RSID))," SNP modulates gene expression in tubulo-interstitial tissue:",sep="")
          }else{
            paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%tub_show$RSID))," SNP(s) modulate(s) gene expression in tubulo-interstitial tissue:",sep="")  
          }
        }else{
          if(length(cred_vars())==1){
            paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not modulate gene expression in tubulo-interstitial tissue.",sep="")
          }else{
            paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, modulates gene expression in tubulo-interstitial tissue.",sep="")  
          }
        }
      }
    
  })
  
  output$tub <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="NEPTUNE_tub")){
        if(any(tub_show$RSID%in%cred_vars())){
          datatable(tub_show[which(tub_show$RSID%in%cred_vars()),], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                      htmltools::tags$p('** NEPTUNE, or Sheng et al, Nat.Genet. 2021')
                      )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  #GTEx eQTL kidney
  output$eqtl_text <- renderText({
    
      if(any(detail_evidence_variant()=="GTEx_eQTL")){
        if(any(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
          if(length(cred_vars())==1){
            paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%GTEx_eQTL$rs_id_dbSNP151_GRCh38p7))," SNP modulates gene expression in kidney cortex tissue:",sep="")
          }else{
            paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%GTEx_eQTL$rs_id_dbSNP151_GRCh38p7))," SNP(s) modulate(s) gene expression in kidney cortex tissue:",sep="")  
          }
        }else{
          if(length(cred_vars())==1){
            paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not modulate gene expression in kidney cortex tissue.",sep="")
          }else{
            paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, modulates gene expression in kidney cortex tissue.",sep="")  
          }
        }
      }
    
  })
  
  output$eqtl <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="GTEx_eQTL")){
        if(any(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
          datatable(GTEx_eQTL_show[which(GTEx_eQTL$rs_id_dbSNP151_GRCh38p7%in%cred_vars()),], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  #GTEx sQTL kidney
  output$sqtl_text <- renderText({
    
      if(any(detail_evidence_variant()=="GTEx_sQTL")){
        if(any(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
          if(length(cred_vars())==1){
            paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%GTEx_sQTL$rs_id_dbSNP151_GRCh38p7))," SNP modulates gene splicing in kidney cortex tissue:",sep="")
          }else{
            paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%GTEx_sQTL$rs_id_dbSNP151_GRCh38p7))," SNP(s) modulate(s) gene splicing in kidney cortex tissue:",sep="")  
          }
        }else{
          if(length(cred_vars())==1){
            paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not modulate gene splicing in kidney cortex tissue.",sep="")
          }else{
            paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, modulates gene splicing in kidney cortex tissue.",sep="")  
          }
         }
      }
    
  })
  
  output$sqtl <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="GTEx_sQTL")){
        if(any(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
          datatable(GTEx_sQTL_show[which(GTEx_sQTL$rs_id_dbSNP151_GRCh38p7%in%cred_vars()),], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* sQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  #GTEx eQTL other tissue
  output$eqtl_wo_kidney_text <- renderText({
    
    if(any(detail_evidence_variant()=="GTEx_wo_kidney_eQTL")){
      if(any(gtex_without_kidney$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
        if(length(cred_vars())==1){
          paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%gtex_without_kidney$rs_id_dbSNP151_GRCh38p7))," SNP modulates gene expression in other tissues:",sep="")
        }else{
          paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%gtex_without_kidney$rs_id_dbSNP151_GRCh38p7))," SNP(s) modulate(s) gene expression in other tissues:",sep="")  
        }
      }else{
        if(length(cred_vars())==1){
          paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not modulate gene expression in other tissues.",sep="")
        }else{
          paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, modulates gene expression in other tissues.",sep="")  
        }
      }
    }
    
  })
  
  output$eqtl_wo_kidney <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="GTEx_wo_kidney_eQTL")){
        if(any(gtex_without_kidney$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
          datatable(gtex_without_kidney_show[which(gtex_without_kidney$rs_id_dbSNP151_GRCh38p7%in%cred_vars()),c(1:13)], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
          )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  #GTEx sQTL other tissue
  output$sqtl_wo_kidney_text <- renderText({
    
    if(any(detail_evidence_variant()=="GTEx_wo_kidney_sQTL")){
      if(any(sqtl_without_kidney$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
        if(length(cred_vars())==1){
          paste("Of your ",length(cred_vars())," searched SNP that is a credible variant, ",length(which(cred_vars()%in%sqtl_without_kidney$rs_id_dbSNP151_GRCh38p7))," SNP modulates gene splicing in other tissues:",sep="")
        }else{
          paste("Of your",length(cred_vars()) ," searched SNPs that are credible variants, ",length(which(cred_vars()%in%sqtl_without_kidney$rs_id_dbSNP151_GRCh38p7))," SNP(s) modulate(s) gene splicing in other tissues:",sep="")  
        }
      }else{
        if(length(cred_vars())==1){
          paste("Your ",length(cred_vars())," searched SNP, that is a credible variant, does not modulate gene splicing in other tissues.",sep="")
        }else{
          paste("None of your ",length(cred_vars())," searched SNPs, that are credible variants, modulates gene splicing in other tissues.",sep="")  
        }
      }
    }
    
  })
  
  output$sqtl_wo_kidney <- DT::renderDataTable({
    if(valid()){
      if(any(detail_evidence_variant()=="GTEx_wo_kidney_sQTL")){
        if(any(sqtl_without_kidney$rs_id_dbSNP151_GRCh38p7%in%cred_vars())){
          datatable(sqtl_without_kidney_show[which(sqtl_without_kidney$rs_id_dbSNP151_GRCh38p7%in%cred_vars()),c(1:12)], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* sQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
          )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  observeEvent(input$question_SNP,{
    showModal(modalDialog(
      titel="Help SNP search",
      tags$h4("Variant search - Help:"),
      p("Enter a RSid of a SNP you are interested in. Other formats than RSids are not supported, yet. If you search without any additional options the association statistics for unconditioned eGFRcrea (if the SNP has an association p-value <5E-8) are displayed. If the SNP is a 99% credible variant for any eGFRcrea signal also conditioned results (if there are other independent signals in the respective locus) and the posterior propabitliy of association (PPA) is shown in a second table. For further information check the Documentation & Help page."),
      easyClose = TRUE,
      footer=NULL,
      size = "m"
    ))
  })
  
  #Locus Zoom
  Locus_ID_snp<-reactive({
    if(valid())
    {cred_var$`Locus Id`[which(cred_var$rsid==snp())]}
  })
  
  #create HTML syntax for link to plots
  snp_plots <- reactive({
    if(valid()){
      ids=which(links_id%in%Locus_ID_snp())[1] ## preliminary only one plot per locus, should be overall plot
      paste("<p><a href=\"",links[ids],"\"target=\"_blank\">Regional association plot of locus",links_id[ids],"</a></p>")
    }
  })
  
  #render links
  output$plot_download_links_snp <- renderUI({
    list(lapply(snp_plots(),HTML))
  })
  
  
  # filename1_snp <- reactive({
  #   if(valid()&inputtype$snp_final=="single"){
  #     if(any(input$snp.locuszoom=="locus_zoom")){
  #       as.character(lz_files[which(locus_id_files==as.character(Locus_ID_snp()))][1])
  #       
  #     }
  #   }
  # })
  # 
  # 
  # output$LocusZoom_SNP <- renderUI({
  #   req(filename1_snp())
  #   
  #   tags$img(height="100%",
  #            width="100%",
  #            src=filename1_snp())
  # })
  
  
  

# gene page ---------------------------------------------------------------
  
  
  #define input type based on click for further analyses
  inputtype <- reactiveValues('gene'=c())
  onclick("Genename", (inputtype$gene=c("single")), add=T)
  onclick("gene_batch", (inputtype$gene=c("list")), add=T)
  onclick("gene-upload-section", (inputtype$gene=c("upload")), add=T)
  
  
  observe({
    
    shinyjs::toggleElement(id="gene_batch_upload_options", condition = (isTruthy(input$file2)) )
    
  })
  
  #show search button:
  
  observe(
    if(isTruthy(input$file2)|isTruthy(input$Genename)|isTruthy(input$gene_batch)){
      shinyjs::showElement(id="gene_go_div")
    }
  )
  
  
  #combine ENTER and "go" click to one reactive value only reacting if there is gene input
  go <- reactiveValues('gene'=c(c()))
  
  observeEvent(input$gene_go, {
    if(isTruthy(input$file2)|isTruthy(input$Genename)|isTruthy(input$gene_batch)){
      if(is.null(go$gene)){go$gene=1}
      else{go$gene <- go$gene+1}
        
    }
  })
  observeEvent(input$keyPressed, {
    if(isTruthy(input$file2)|isTruthy(input$Genename)|isTruthy(input$gene_batch)){
      if(is.null(go$gene)){go$gene=1}
      else{go$gene <- go$gene+3}
    }
    })
  
   
  #output$testforme <- renderText({as.integer(go$gene)}) <- test option
  ###first gene input processing: 
  single_gene <- eventReactive(go$gene, {input$Genename})
  
  ###build vector of genes from text input
  input_gene_batch_list <- eventReactive(go$gene,{
    
    req(input$gene_batch)
    
    igbl=unlist(strsplit(input$gene_batch,"[^[:alnum:]^-]",perl=TRUE)) ## excludes all non-word characteres except "-" that might be used as separators
    if(any(igbl=="")){
      igbl[-which(igbl=="")]
    }else{igbl}
  })
  
  
  ###build vector of genes from file input
  
  input_gene_batch_upload <- reactiveValues('input'= c())
    
  observe({
    req(input$file2)
    
    # df <- read.table(input$file2$datapath,
    #                  sep = input$sep_gene,
    #                  quote = input$quote_gene)
    # 
    # input_gene_batch_upload$input = paste(df[1,])
    input_gene_batch_upload$input <- scan(file=input$file2$datapath,
                                          sep = input$sep_gene,
                                          what="character")
    
  }) 
  
  #reset file input, when reset button is clicked
  observeEvent(input$reset_gene,{
    reset("file2")
    input_gene_batch_upload$input <- NULL
  })
  
  #final selection:
  gene <- eventReactive(go$gene, {
    
    if(inputtype$gene=="single"){
      toupper(single_gene())
    }else{
      if(inputtype$gene=="list"){
        toupper(input_gene_batch_list())
      }else{
        if(inputtype$gene=="upload"){
          toupper(input_gene_batch_upload$input)
        }else{NULL}
      }
    }
  })
  
  inputtype <- reactiveValues('gene_final'=c()) #use inputtype$gene_final as it only reacts on "search"-click
  observeEvent(go$gene,{
    inputtype$gene_final = inputtype$gene 
  })
  
  
  
  
  # skript for gene page
  
  correct_genenames <- eventReactive(go$gene,{
    if(!is.null(gene())){
      if(any(gene()%in%toupper(genes_sorted$Gene))){
        gene()[which(gene()%in%toupper(genes_sorted$Gene))]
      }else{NULL}
    }else{NULL}
  })
  
  
  genes_in_gps <- eventReactive(go$gene,{
    if(!is.null(gene())){
      if(any(gene()%in%toupper(GPS$Gene))){
        gene()[which(gene()%in%toupper(GPS$Gene))]
      }else{NULL}
    }else{NULL}
  })
  
  valid_gene <- eventReactive(go$gene, {!is.null(genes_in_gps())})
  
  columns <- eventReactive(go$gene, {if(isTruthy(input$columns1)|isTruthy(input$columns2)|isTruthy(input$columns3)){ as.numeric(c(input$columns1,input$columns2,input$columns3))}else{NULL}})
  #detail_evidence_gene <- eventReactive(go$gene, {c(input$detail_evidence_gene_1,input$detail_evidence_gene_2,input$detail_evidence_gene_3)})
  
  detail_evidence_locus_gene <- eventReactive(go$gene, {input$detail_evidence_locus_gene})
  
  ppa <- eventReactive(go$gene, {input$ppa})
  
  incorrect_genenames <- eventReactive(go$gene,{
    length(gene())-length(correct_genenames())
  })
  
  output$namecheck_gene <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.null(correct_genenames())){
        paste(gene()," is no valid/official gene name from HGNC, please check your input.",sep="")
      }else{
        if(is.null(genes_in_gps())){
          paste(gene()," is not included in any eGFRcrea locus.", sep="")
        }else{
          paste("Following information has been found for ", genes_in_gps(),":",sep="")
        }
      }
    }else{
      if(is.null(correct_genenames)){
        if(incorrect_genenames()==1){
          paste("You queried ",length(gene())," gene name. ",gene(),"  is no official/valid gene name from HGNC. No further searches were performed. Please check your input.",sep="")
        }else{
          paste("You queried ",length(gene())," gene names. All of these are no official/valid gene names from HGNC. No further searches were performed. Please check your input.",sep="")
        }
      }else{
        if(is.null(genes_in_gps())){
          if(incorrect_genenames()==1){
            paste("You queried ",length(gene())," genes. ",incorrect_genenames()," of these is no official/valid gene name from HGNC. However, none of your ", length(correct_genenames()), " queried official genes is not included in any of the 424 eGFRcrea loci. Thus, no GPS entries are available.",sep="")
          }else{
            if(incorrect_genenames()==0){
              paste("You queried ",length(gene())," genes. All of these are official/valid gene names from HGNC. However, none of your ", length(correct_genenames()), " queried official genes is included in any of the 424 eGFRcrea loci. Thus, no GPS entries are available.",sep="")
            }else{
              paste("You queried ",length(gene())," genes. ",incorrect_genenames()," of these are no official/valid gene names from HGNC. However, none of your ", length(correct_genenames()), " queried official genes is included in any of the 424 eGFRcrea loci. Thus, no GPS entries are available.",sep="")
            }
          }
        }else{
          if(length(correct_genenames())==1){
            if(incorrect_genenames()==1){
              paste("You queried ",length(gene())," genes. ",incorrect_genenames()," of these is no official/valid gene names from HGNC. ",length(genes_in_gps()), " (",genes_in_gps(),")"," of your ", length(correct_genenames()), " queried official genes is included in any of the 424 eGFRcrea loci and thus included in the GPS.",sep="")
            }else{
              if(incorrect_genenames()==0){
                paste("You queried ",length(gene())," genes. All of these are official/valid gene names from HGNC. ",length(genes_in_gps())," (",genes_in_gps(),")"," of your ", length(correct_genenames()), " queried official genes is included in any of the 424 eGFRcrea loci and thus included in the GPS.",sep="")
              }else{
                paste("You queried ",length(gene())," genes. ",incorrect_genenames()," of these are no official/valid gene names from HGNC. ",length(genes_in_gps())," (",genes_in_gps(),")"," of your ", length(correct_genenames()), " queried official genes is included in any of the 424 eGFRcrea loci and thus included in the GPS.",sep="")
              }
            }
          }else{
            if(incorrect_genenames()==1){
              paste("You queried ",length(gene())," genes. ",incorrect_genenames()," of these is no official/valid gene names from HGNC. ",length(genes_in_gps())," of your ", length(correct_genenames()), " queried official genes are included in any of the 424 eGFRcrea loci and thus included in the GPS.",sep="") 
            }else{
              if(incorrect_genenames()==0){
                paste("You queried ",length(gene())," genes. All of these are official/valid gene names from HGNC. ",length(genes_in_gps())," of your ", length(correct_genenames()), " queried official genes are included in any of the 424 eGFRcrea loci and thus included in the GPS.",sep="") 
              }else{
                paste("You queried ",length(gene())," genes. ",incorrect_genenames()," of these are no official/valid gene names from HGNC. ",length(genes_in_gps())," of your ", length(correct_genenames()), " queried official genes are included in any of the 424 eGFRcrea loci and thus included in the GPS.",sep="") 
              }
            }
          }
        }
      }
    }
  })
  
 
  #CADD
  output$CADD_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()%in%c(11:13))){
        if(any(toupper(CADD$GeneName)%in%genes_in_gps())){
          x = CADD_show[which(toupper(CADD$GeneName)%in%genes_in_gps()&CADD$PPA>=ppa()),]
          if(any(duplicated(x))) x = x[-which(duplicated(x)),]
          datatable(x, rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal)'),
                      htmltools::tags$p('Reference and alternative allele from human reference genome (GRCh37)')
                      
                    )
                    ,
                    class ='display',
                    container = sketchcadd
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  cadd_gene <- reactive({
    if(valid_gene()){
      if(any(columns()%in%c(11:13))){
        if(any(toupper(CADD$GeneName)%in%genes_in_gps())){
          x = CADD[which(toupper(CADD$GeneName)%in%genes_in_gps()&CADD$PPA>=ppa()),]
          if(any(duplicated(x[,-c(19,20)]))){
            l = which(duplicated(x[,-c(19,20)]))
            x = x[-l,]
          } 
          x
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$CADD_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(cadd_gene())){
          paste(genes_in_gps()," contains at least one credible variant that is predicted deleterious and has a functional consequence for ",genes_in_gps() ,":",sep="")
      }else{
          paste(genes_in_gps()," does not contain a credible variant that is predicted deleterious and has a functional consequence for ",genes_in_gps(),".",sep="")
      }
    }else{
      if(is.data.frame(cadd_gene())){
        if(length(genes_in_gps())==1){
          paste(genes_in_gps()," contains at least one credible variant that is predicted deleterious and has a functional consequence for ",genes_in_gps() ,":",sep="")
        }else{
          if(length(unique(cadd_gene()[,15])==1)){
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, ",length(unique(cadd_gene()[,15]))," gene contains at least one credible variant that is predicted deleterious and has a functional consequence for this gene:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, ",length(unique(cadd_gene()[,15]))," genes contain at least one credible variant that is predicted deleterious and has a functional consequence for the respective gene:",sep="")
          }
        }
      }else{
        paste("None of your ",length(genes_in_gps()), " searched genes included in the GPS contains at least one credible variant that is predicted deleterious and has a functional consequence for the respective gene.",sep="")
      }
    }
  })
  
  #Neptune glomerulus
  output$glo_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==14)){
        if(any(toupper(glo_show$`Affected gene`)%in%genes_in_gps())){
          datatable(glo_show[which(toupper(glo_show$`Affected gene`)%in%genes_in_gps() & glo_show$PPA>=ppa()),], rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                      htmltools::tags$p('** NEPTUNE, or Sheng et al, Nat.Genet. 2021')
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  glo_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==14)){
        if(any(toupper(glo$Approved.symbol)%in%genes_in_gps())){
          glo[which(toupper(glo$Approved.symbol)%in%genes_in_gps() & glo$PPA>=ppa()),]
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$glo_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(glo_gene())){
        paste("Expression of ",genes_in_gps()," is modulated by a credible variant in glomerular tissue:",sep="")
      }else{
        paste("No credible variant modulates expression of ",genes_in_gps()," in glomerular tissue.",sep="")
      }
    }else{
      if(is.data.frame(glo_gene())){
        if(length(genes_in_gps())==1){
          paste("Expression of ",genes_in_gps()," is modulated by a credible variant in glomerular tissue:",sep="")
        }else{
          if(length(unique(glo_gene()[,5]))==1){
           paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(glo_gene()[,5]))," gene is modulated by a credible variant in glomerular tissue:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(glo_gene()[,5]))," genes is modulated by a credible variant in glomerular tissue:",sep="")
          }
        }
        
      }else{
        paste("Expression of none of your ",length(genes_in_gps()), " searched genes included in the GPS is modulated by a credible variant in glomerular tissue.",sep="")
      }
    }
  })
  
  #Neptune tubulo-interstitium
  output$tub_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==15)){
        if(any(toupper(tub_show$`Affected gene`)%in%genes_in_gps())){
          datatable(tub_show[which(toupper(tub_show$`Affected gene`)%in%genes_in_gps()&tub_show$PPA>=ppa()),],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                      htmltools::tags$p('** NEPTUNE, or Sheng et al, Nat.Genet. 2021')
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  tub_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==15)){
        if(any(toupper(tub$Approved.symbol)%in%genes_in_gps())){
          tub[which(toupper(tub$Approved.symbol)%in%genes_in_gps()&tub$PPA>=ppa()),]
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$tub_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(tub_gene())){
        paste("Expression of ",genes_in_gps()," is modulated by a credible variant in tubulo-interstitial tissue:",sep="")
      }else{
        paste("No credible variant modulates expression of ",genes_in_gps()," in tubulo-interstitial tissue.",sep="")
      }
    }else{
      if(is.data.frame(tub_gene())){
        if(length(genes_in_gps())==1){
          paste("Expression of ",genes_in_gps()," is modulated by a credible variant in tubulo-interstitial tissue:",sep="")
        }else{
          if(length(unique(tub_gene()[,5]))==1){
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(tub_gene()[,5]))," gene is modulated by a credible variant in tubulo-interstitial tissue:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(tub_gene()[,5]))," genes is modulated by a credible variant in tubulo-interstitial tissue:",sep="")
          }
          
        }
        
      }else{
        paste("Expression of none of your ",length(genes_in_gps()), " searched genes included in the GPS is modulated by a credible variant in tubulo-interstitial tissue.",sep="")
      }
    }  
  })
  
  #GTEx eQTL kidney
  output$eqtl_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==16)){
        if(any(toupper(GTEx_eQTL$gene)%in%genes_in_gps())){
          datatable(GTEx_eQTL_show[which(toupper(GTEx_eQTL$gene)%in%genes_in_gps()&GTEx_eQTL$PPA>=ppa()),],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  eqtl_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==16)){
        if(any(toupper(GTEx_eQTL$gene)%in%genes_in_gps())){
          GTEx_eQTL[which(toupper(GTEx_eQTL$gene)%in%genes_in_gps()&GTEx_eQTL$PPA>=ppa()),]
          
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$eqtl_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(eqtl_gene())){
        paste("Expression of ",genes_in_gps()," is modulated by a credible variant in kidney cortex tissue:",sep="")
      }else{
        paste("No credible variant modulates expression of ",genes_in_gps()," in kidney cortex tissue.",sep="")
      }
    }else{
      if(is.data.frame(eqtl_gene())){
        if(length(genes_in_gps())==1){
          paste("Expression of ",genes_in_gps()," is modulated by a credible variant in kidney cortex tissue:",sep="")
        }else{
          if(length(unique(eqtl_gene()[,17]))==1){
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(eqtl_gene()[,17]))," gene is modulated by a credible variant in kidney cortex tissue:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(eqtl_gene()[,17]))," genes is modulated by a credible variant in kidney cortex tissue:",sep="")
          }
          
        }
        
      }else{
        paste("Expression of none of your ",length(genes_in_gps()), " searched genes included in the GPS is modulated by a credible variant in kidney cortex tissue.",sep="")
      }
    }
  })
  
  #GTEx other tissue eQTL
  
  output$eqtl_wo_kidney_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==19)){
        if(any(toupper(gtex_without_kidney_show[,3])%in%genes_in_gps())){
          datatable(gtex_without_kidney_show[which(toupper(gtex_without_kidney_show[,3])%in%genes_in_gps()& gtex_without_kidney_show$PPA >= ppa()),c(1:13)],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions='Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* eQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
          )%>% formatPercentage('PPA')
        }
      }
    }
  })
  eqtl_wo_kidney_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==19)){
        if(any(toupper(gtex_without_kidney_show[,3])%in%genes_in_gps())){
          gtex_without_kidney_show[which(toupper(gtex_without_kidney_show[,3])%in%genes_in_gps()& gtex_without_kidney_show$PPA >= ppa()),]
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$eqtl_wo_kidney_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(eqtl_wo_kidney_gene())){
        paste("Expression of ",genes_in_gps()," is modulated by a credible variant in other tissues:",sep="")
      }else{
        paste("No credible variant modulates expression of ",genes_in_gps()," in other tissues.",sep="")
      }
    }else{
      if(is.data.frame(eqtl_wo_kidney_gene())){
        if(length(genes_in_gps())==1){
          paste("Expression of ",genes_in_gps()," is modulated by a credible variant in other tissues:",sep="")
        }else{
          if(length(unique(eqtl_wo_kidney_gene()[,3]))==1){
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(eqtl_wo_kidney_gene()[,3]))," gene is modulated by a credible variant in other tissues:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," expression of ",length(unique(eqtl_wo_kidney_gene()[,3]))," genes is modulated by a credible variant in other tissues:",sep="")
          }
          
        }
        
      }else{
        paste("Expression of none of your ",length(genes_in_gps()), " searched genes included in the GPS is modulated by a credible variant in other tissues.",sep="")
      }
    }
  })
  
  
  #sQTL kidney
  output$sqtl_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==18)){
        if(any(toupper(GTEx_sQTL$gene)%in%genes_in_gps())){
          datatable(GTEx_sQTL_show[which(toupper(GTEx_sQTL$gene)%in%genes_in_gps()& GTEx_sQTL$PPA >= ppa()),],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions='Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* sQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
                    )%>% formatPercentage('PPA')
        }
      }
    }
  })
  
  sqtl_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==18)){
        if(any(toupper(GTEx_sQTL$gene)%in%genes_in_gps())){
          GTEx_sQTL[which(toupper(GTEx_sQTL$gene)%in%genes_in_gps()& GTEx_sQTL$PPA >= ppa()),]
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$sqtl_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(sqtl_gene())){
        paste("Splicing of ",genes_in_gps()," is modulated by a credible variant in kidney cortex tissue:",sep="")
      }else{
        paste("No credible variant modulates splicing of ",genes_in_gps()," in kidney cortex tissue.",sep="")
      }
    }else{
      if(is.data.frame(sqtl_gene())){
        if(length(genes_in_gps())==1){
          paste("Splicing of ",genes_in_gps()," is modulated by a credible variant in kidney cortex tissue:",sep="")
        }else{
          if(length(unique(sqtl_gene()[,18]))==1){
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," splicing of ",length(unique(sqtl_gene()[,18]))," gene is modulated by a credible variant in kidney cortex tissue:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," splicing of ",length(unique(sqtl_gene()[,18]))," genes is modulated by a credible variant in kidney cortex tissue:",sep="")
          }
        }
        
      }else{
        paste("Splicing of none of your ",length(genes_in_gps()), " searched genes included in the GPS is modulated by a credible variant in kidney cortex tissue.",sep="")
      }
    }
  })
  
  
  
  #sQTL without kidney
  output$sqtl_wo_kidney_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==19)){
        if(any(toupper(sqtl_without_kidney_show[,3])%in%genes_in_gps())){
          datatable(sqtl_without_kidney_show[which(toupper(sqtl_without_kidney_show[,3])%in%genes_in_gps()& sqtl_without_kidney_show$PPA >= ppa()),c(1:13)],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions='Buttons',
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* sQTL effect allele is the eGFRcrea lowering allele'),
                      htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal), EAF: Effect allele frequency'),
                    )
          )%>% formatPercentage('PPA')
        }
      }
    }
  })
  sqtl_wo_kidney_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==19)){
        if(any(toupper(sqtl_without_kidney_show[,3])%in%genes_in_gps())){
          sqtl_without_kidney_show[which(toupper(sqtl_without_kidney_show[,3])%in%genes_in_gps()& sqtl_without_kidney_show$PPA >= ppa()),]
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$sqtl_wo_kidney_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(sqtl_wo_kidney_gene())){
        paste("Splicing of ",genes_in_gps()," is modulated by a credible variant in other tissues:",sep="")
      }else{
        paste("No credible variant modulates splicing of ",genes_in_gps()," in other tissues.",sep="")
      }
    }else{
      if(is.data.frame(sqtl_wo_kidney_gene())){
        if(length(genes_in_gps())==1){
          paste("Splicing of ",genes_in_gps()," is modulated by a credible variant in other tissues:",sep="")
        }else{
          if(length(unique(sqtl_wo_kidney_gene()[,3]))==1){
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," splicing of ",length(unique(sqtl_wo_kidney_gene()[,3]))," gene is modulated by a credible variant in other tissues:",sep="")
          }else{
            paste("Of your ",length(genes_in_gps()), " searched genes included in the GPS, "," splicing of ",length(unique(sqtl_wo_kidney_gene()[,3]))," genes is modulated by a credible variant in other tissues:",sep="")
          }
          
        }
        
      }else{
        paste("Splicing of none of your ",length(genes_in_gps()), " searched genes included in the GPS is modulated by a credible variant in other tissues.",sep="")
      }
    }
  })
  
  #MGI
  output$mgi_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==20)){
        if(any(toupper(MGI$human_symbol)%in%genes_in_gps())){
          datatable(MGI_show[which(toupper(MGI$human_symbol)%in%genes_in_gps()),],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons'
                    )
        }
      }
    }
  })
  
  mgi_gene<-reactive({
    if(valid_gene()){
      if(any(columns()==20)){
        if(any(toupper(MGI$human_symbol)%in%genes_in_gps())){
          (MGI_show[which(toupper(MGI$human_symbol)%in%genes_in_gps()),]
          )
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$mgi_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(mgi_gene())){
        paste("Following kidney phenotypes in mice are described for ",genes_in_gps()," :",sep="")
      }else{
        paste("No kidney phenotypes in mice were described for ",genes_in_gps()," .",sep="")
      }
    }else{
      if(is.data.frame(mgi_gene())){
        if(length(genes_in_gps())==1){
          paste("Following kidney phenotypes in mice are described for ",genes_in_gps()," :",sep="")
        }else{
          if(length(unique(mgi_gene()[,1]))==1){
            paste("Of your ", length(genes_in_gps())," searched genes included in the GPS, ",length(unique(mgi_gene()[,1]))," gene is described with kidney phenotypes in mice:",sep="")
          }else{
            paste("Of your ", length(genes_in_gps())," searched genes included in the GPS, ",length(unique(mgi_gene()[,1]))," genes are described with kidney phenotypes in mice:",sep="")
          }
        }
        
      }else{
        paste("No kidney phenotypes in mice were described for any of your ",length(genes_in_gps()), " searched genes included in the GPS.",sep="")
      }
    }
    
  })
  
  #OMIM
  output$omim_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(columns()==21)){
          if(any(toupper(OMIM$Gene)%in%genes_in_gps())){
          datatable(OMIM_MIM[which(toupper(OMIM$Gene)%in%genes_in_gps()),c(1,4,2,3)],rownames=FALSE,
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = "_all"))
                      , scrollX=TRUE,
                      dom = 'lfrtBip',
                      buttons= list(c('copy', 'excel'))
                    ),
                    extensions = 'Buttons',
                    escape = FALSE
                    ,caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                      htmltools::tags$p('* Link to OMIM entry'),
                      htmltools::tags$p('** Online Mendelian Inheritance in Man, OMIM or Groopman et al, N.Engl.J.Med. 2019')
                    )
                    )
        }
      }
    } 
  })
  
  
  
  omim_gene <- reactive({
    if(valid_gene()){
      if(any(columns()==21)){
        if(any(toupper(OMIM$Gene)%in%genes_in_gps())){
          OMIM_MIM[which(toupper(OMIM$Gene)%in%genes_in_gps()),c(1,4,2,3)]
        }else{0}
      }else{0}
    }else{0}
  })
  
  output$omim_gene_text <- renderText({
    if(inputtype$gene_final=="single"){
      if(is.data.frame(omim_gene())){
        paste("Following genetic disorders with kidney phenotypes in human are described for ",genes_in_gps()," :",sep="")
      }else{
        paste("No genetic disorders with kidney phenotypes in human were described for ",genes_in_gps()," .",sep="")
      }
    }else{
      if(is.data.frame(omim_gene())){
        if(length(genes_in_gps())==1){
          paste("Following genetic disorders with kidney phenotypes in human are described for ",genes_in_gps()," :",sep="")
        }else{
          if(length(unique(omim_gene()[,1]))==1){
            paste("Of your ", length(genes_in_gps())," searched genes included in the GPS, ",length(unique(omim_gene()[,1]))," gene is described with genetic disorders with kidney phenotypes in human:",sep="")
          }else{
            paste("Of your ", length(genes_in_gps())," searched genes included in the GPS, ",length(unique(omim_gene()[,1]))," genes are described with genetic disorders with kidney phenotypes in human:",sep="")
          }
          
        }
        
      }else{
        paste("No genetic disorders with kidney phenotypes in human were described for any of your ",length(genes_in_gps()), " searched genes included in the GPS.",sep="")
      }
    }
    
  })
  
  #credible variants
  output$cred_var_gene <- DT::renderDataTable({
    if(valid_gene()){
      if(any(detail_evidence_locus_gene()=="cred_var")){
        
        a=unique(GPS$locus_id[which(toupper(GPS$Gene)%in%genes_in_gps())])
        
        z=cred_var_show[which(cred_var$`Locus Id`%in%a &cred_var$ppa >= ppa()),]
        
        snps=cred_var$rsid[which(cred_var$`Locus Id`%in%a &cred_var$ppa >= ppa())]
        eqtl=rep("-",length(snps))
        if(is.data.frame(tub_gene())){eqtl[which(snps%in%tub_gene()[,1])]="eQTL"}
        if(is.data.frame(glo_gene())){eqtl[which(snps%in%glo_gene()[,1])]="eQTL"}
        if(is.data.frame(eqtl_gene())){eqtl[which(snps%in%eqtl_gene()[,1])]="eQTL"}
        sqtl=rep("-",length(snps))
        if(is.data.frame(sqtl_gene())){sqtl[which(snps%in%sqtl_gene()[,1])]="sQTL"}
        cadd=rep("-",length(snps))
        if(is.data.frame(cadd_gene())){cadd[which(snps%in%cadd_gene()[,1])]="protein-relevant SNP"}
        
        z['functional.consequence'] = ifelse(eqtl=="-", ifelse(sqtl=="-", ifelse(cadd=="-", "-", cadd),ifelse(cadd=="-", sqtl, paste(sqtl,cadd,sep=", "))),ifelse(sqtl=="-", ifelse(cadd=="-", eqtl, cadd),ifelse(cadd=="-", paste(eqtl,sqtlsep=", "), paste(eqtl,sqtl,cadd, sep=", "))))
        z=z[,c(1:4,17,5:16)]
        
        
        
        datatable(z, rownames=FALSE,
                  options = list(
                    columnDefs = list(list(className = 'dt-left', targets = "_all"))
                    , scrollX=TRUE,
                    dom = 'lfrtBip',
                    buttons= list(c('copy', 'excel'))
                  ),
                  extensions = 'Buttons',
                  escape=FALSE,
                  caption = htmltools::tags$caption(
                    style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                    htmltools::tags$p('* Link to dbSNP'),
                    htmltools::tags$p('PPA: posterior probability of association of the SNP (probability of driving the association of the respective signal)'),
                    htmltools::tags$p('N: Sample size of European-ancestry GWAS meta-analysis for this SNP, Effect allele: eGFRcrea lowering allele, EAF: Effect allele frequency'),
                    htmltools::tags$p('For loci with multiple independent signals, the association results are provided conditioned and unconditioned on the signal index variants of the other signals in the locus.')
                  ),
                  class='display',
                  container = sketchcredvargene
        )%>% formatPercentage('PPA')
      }
    }
  })
  
  #GPS
  
  GPS_gene <- reactive({
    if(valid_gene()){
      
        
        x=GPS_show[which(toupper(GPS$Gene)%in%genes_in_gps()),c(1:10,24,20,21,22,25:28)]

        y=GPS[which(toupper(GPS$Gene)%in%genes_in_gps()),]
        rows=nrow(x)
        cadd1=rep(0,rows)
        cadd2=rep(0,rows)
        cadd3=rep(0,rows)
        tub_g=rep(0,rows)
        glo_g=rep(0,rows)
        eqtl_kid=rep(0,rows)
        sqtl_kid=rep(0,rows)
        eqtl_all=rep(0,rows)
        sqtl_all=rep(0,rows)
        for(i in 1:rows){
          
          if(is.data.frame(cadd_gene())){
            cadd1[i]=nrow(cadd_gene()[which(cadd_gene()['signal_id']==y[i,26] & cadd_gene()['ConsScore']>6 & cadd_gene()['GeneName']==y[i,3] & cadd_gene()['region_id']==y[i,2]),])
            
            cadd2[i]=nrow(cadd_gene()[which(cadd_gene()['signal_id']==y[i,26] & cadd_gene()['ConsScore']>4& cadd_gene()['ConsScore']<=6& cadd_gene()['GeneName']==y[i,3] & cadd_gene()['region_id']==y[i,2]),])
            
            cadd3[i]=nrow(cadd_gene()[which(cadd_gene()['signal_id']==y[i,26] & cadd_gene()['ConsScore']<=4& cadd_gene()['GeneName']==y[i,3] & cadd_gene()['region_id']==y[i,2]),])
          }
          if(is.data.frame(tub_gene())){
            
            tub_g[i]=length(unique(tub_gene()[which(tub_gene()['signal_id']==y[i,26]& tub_gene()['Approved.symbol']==y[i,3] & tub_gene()['region_id']==y[i,2]),1]))
          }
          if(is.data.frame(glo_gene())){
            
            glo_g[i]=length(unique(glo_gene()[which(glo_gene()['signal_id']==y[i,26]& glo_gene()['Approved.symbol']==y[i,3] & glo_gene()['region_id']==y[i,2]),1]))
          }
          if(is.data.frame(eqtl_gene())){
            eqtl_kid[i]=nrow(eqtl_gene()[which(eqtl_gene()['signal_id']==y[i,26]& eqtl_gene()['gene']==y[i,3] & eqtl_gene()['region_id']==y[i,2]),])
          }
          if(is.data.frame(eqtl_wo_kidney_gene())){
            eqtl_all[i]=length(unique(eqtl_wo_kidney_gene()[which(eqtl_wo_kidney_gene()['signal_id']==y[i,26]& eqtl_wo_kidney_gene()['Affected gene']==y[i,3] & eqtl_wo_kidney_gene()['region_id']==y[i,2]),1]))
          }
          if(is.data.frame(sqtl_gene())){
            sqtl_kid[i]=nrow(sqtl_gene()[which(sqtl_gene()['signal_id']==y[i,26]& sqtl_gene()['gene']==y[i,3] & sqtl_gene()['region_id']==y[i,2]),])
          }
          if(is.data.frame(sqtl_wo_kidney_gene())){
            sqtl_all[i]=length(unique(sqtl_wo_kidney_gene()[which(sqtl_wo_kidney_gene()['signal_id']==y[i,26]& sqtl_wo_kidney_gene()['Affected gene']==y[i,3] & sqtl_wo_kidney_gene()['region_id']==y[i,2]),1]))
          }
        }
        x[,19]=cadd1
        x[,20]=cadd2
        x[,21]=cadd3
        x[,22]=glo_g
        x[,23]=tub_g
        x[,24]=eqtl_kid
        x[,25]=eqtl_all
        x[,26]=sqtl_kid
        x[,27]=sqtl_all
        x=x[,c(1:11,19:27,12,13,14,15:18)]
        names(x)=c("Locus name**","Locus ID", "Signal ID", "Gene*", "Distance to locus lead variant","Chromosome","Start of gene","End of gene", "eGFRcys or BUN validation", "# credible variants in signal","max PPA","stop-gained, stop-lost, non-synonymus","canonical splice, noncoding change, synonymous, splice site","other deleterious variant","eQTL glomerular tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL tubulo-interstitial tissue (NEPTUNE, or Sheng et al [Nat.Genet. 2021])","eQTL kidney cortex tissue (GTEx)","eQTL other tissue (GTEx)","sQTL kidney cortex tisse (GTEx)","sQTL other tissue (GTEx)","mouse phenotype","human phenotype", "Coloc in NEPTUNE tissue","Diabetes specific effect","eGFRcrea decline effect","Drugable","In enriched pathway")
        x[,c(4,1:3,9,10,11,(columns()+1),24:27,5:8)] #max PPA column is new -> alternativly change column values everywhere
        
      
    }else{0}
  })
  
  
  
  rowsGPS <- reactive({
    if(valid_gene()){nrow(GPS_gene())}
    else{0}
  })
  
  #header GPS
  sketch_gene <- reactive({
    if(any(columns()%in%c(11:13))&any(columns()%in%c(14:19))&any(columns()%in%c(20:22))){
       htmltools::withTags(table(
              class = 'display',

              thead(
                tr(
                  th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                  th(rowspan = 2, 'Locus name**'),
                  th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                  th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                  th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                  th(rowspan = 2, '# credible variants in signal'),
                  th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                  # lapply(c('Gene*','Locus name**','Locus ID','Signal ID','eGFRcys or BUN validation','# credible variants in signal'),th, rowspan=2),
                  th(colspan = length(which(columns()%in%c(11:13))), class="CADD", span('# protein-relevant credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                  th(colspan = length(which(columns()%in%c(14:19))), class="eqtl", '# credible variants that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                  lapply(names(GPS_show[columns()[which(columns()==(22))]]), class="eqtl", th, rowspan=2),
                  lapply(names(GPS_show[columns()[which(columns()%in%c(20:21))]]), class="pheno", th, rowspan=2),
                  lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)
                ),
                tr(
                  lapply(names(GPS_show[columns()[which(columns()%in%c(11:13))]]), th, class="CADD"),
                  lapply(names(GPS_show[columns()[which(columns()%in%c(14:19))]]), th, class="eqtl")
                )
              )
            ))
          }else{

            if(!any(columns()%in%c(11:13))&any(columns()%in%c(14:19))&any(columns()%in%c(20:22))){
              htmltools::withTags(table(
                class = 'display',

                thead(
                  tr(
                    th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                    th(rowspan = 2, 'Locus name**'),
                    th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                    th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                    th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                    th(rowspan = 2, '# credible variants in signal'),
                    th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                    #th(colspan = length(which(columns()%in%c(11:13))), span('# deleterious credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                    th(colspan = length(which(columns()%in%c(14:19))), class="eqtl", '# credible variants, that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                    lapply(names(GPS_show[columns()[which(columns()==(22))]]), class="eqtl", th, rowspan=2),
                    lapply(names(GPS_show[columns()[which(columns()%in%c(20:21))]]), class="pheno", th, rowspan=2),
                    lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                  ),
                  tr(
                    lapply(names(GPS_show[columns()[which(columns()%in%c(11:13))]]), th, class="CADD"),
                    lapply(names(GPS_show[columns()[which(columns()%in%c(14:19))]]), th, class="eqtl")
                  )
                )
              ))
            }else{

              if(any(columns()%in%c(11:13))&!any(columns()%in%c(14:19))&any(columns()%in%c(20:22))){
                htmltools::withTags(table(
                  class = 'display',

                  thead(
                    tr(
                      th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                      th(rowspan = 2, 'Locus name**'),
                      th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                      th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                      th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                      th(rowspan = 2, '# credible variants in signal'),
                      th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                      th(colspan = length(which(columns()%in%c(11:13))), class="CADD", span('# protein-relevant credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                      #th(colspan = length(which(columns()%in%c(14:19))), '# credible variants that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                      lapply(names(GPS_show[columns()[which(columns()==(22))]]), class="eqtl", th, rowspan=2),
                      lapply(names(GPS_show[columns()[which(columns()%in%c(20:21))]]), class="pheno", th, rowspan=2),
                      lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                    ),
                    tr(
                      lapply(names(GPS_show[columns()[which(columns()%in%c(11:13))]]), th, class="CADD"),
                      lapply(names(GPS_show[columns()[which(columns()%in%c(14:19))]]), th, class="eqtl")
                    )
                  )
                ))
              }else{
                if(any(columns()%in%c(11:13))&any(columns()%in%c(14:19))&!any(columns()%in%c(20:22))){
                  htmltools::withTags(table(
                    class = 'display',

                    thead(
                      tr(
                        th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                        th(rowspan = 2, 'Locus name**'),
                        th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                        th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                        th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                        th(rowspan = 2, '# credible variants in signal'),
                        th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                        th(colspan = length(which(columns()%in%c(11:13))), class="CADD", span('# protein-relevant credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                        th(colspan = length(which(columns()%in%c(14:19))), class="eqtl", '# credible variants that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                        #lapply(names(GPS_show[columns()[which(columns()%in%c(20:22))]]), th, rowspan=2),
                        lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                      ),
                      tr(
                        lapply(names(GPS_show[columns()[which(columns()%in%c(11:13))]]), th, class="CADD"),
                        lapply(names(GPS_show[columns()[which(columns()%in%c(14:19))]]), th, class="eqtl")
                      )
                    )
                  ))
                }else{
                  if(!any(columns()%in%c(11:13))&!any(columns()%in%c(14:19))&any(columns()%in%c(20:22))){
                    htmltools::withTags(table(
                      class = 'display',

                      thead(
                        tr(
                          th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                          th(rowspan = 2, 'Locus name**'),
                          th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                          th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                          th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                          th(rowspan = 2, '# credible variants in signal'),
                          th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                          #th(colspan = length(which(columns()%in%c(11:13))), span('# deleterious credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                          #th(colspan = length(which(columns()%in%c(14:19))), '# credible variants, that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                          lapply(names(GPS_show[columns()[which(columns()==(22))]]), class="eqtl", th, rowspan=2),
                          lapply(names(GPS_show[columns()[which(columns()%in%c(20:21))]]), class="pheno", th, rowspan=2),
                          lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                        )
                      )
                    ))
                  }else{
                    if(!any(columns()%in%c(11:13))&any(columns()%in%c(14:19))&!any(columns()%in%c(20:22))){
                      htmltools::withTags(table(
                        class = 'display',

                        thead(
                          tr(
                            th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                            th(rowspan = 2, 'Locus name**'),
                            th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                            th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                            th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                            th(rowspan = 2, '# credible variants in signal'),
                            th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                            #th(colspan = length(which(columns()%in%c(11:13))), span('# deleterious credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                            th(colspan = length(which(columns()%in%c(14:19))), class="eqtl", '# credible variants, that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                            #lapply(names(GPS_show[columns()[which(columns()%in%c(20:22))]]), th, rowspan=2),
                            lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                          ),
                          tr(
                            lapply(names(GPS_show[columns()[which(columns()%in%c(11:13,14:19))]]), th)
                          )
                        )
                      ))
                    }else{
                      if(any(columns()%in%c(11:13))&!any(columns()%in%c(14:19))&!any(columns()%in%c(20:22))){
                        htmltools::withTags(table(
                          class = 'display',

                          thead(
                            tr(
                              th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                              th(rowspan = 2, 'Locus name**'),
                              th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                              th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                              th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                              th(rowspan = 2, '# credible variants in signal'),
                              th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                              th(colspan = length(which(columns()%in%c(11:13))), class="CADD", span('# protein-relevant credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                              #th(colspan = length(which(columns()%in%c(14:19))), '# credible variants that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                              #lapply(names(GPS_show[columns()[which(columns()%in%c(20:22))]]), th, rowspan=2),
                              lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                            ),
                            tr(
                              lapply(names(GPS_show[columns()[which(columns()%in%c(11:13))]]), th, class="CADD"),
                              lapply(names(GPS_show[columns()[which(columns()%in%c(14:19))]]), th, class="eqtl")
                            )
                          )
                        ))
                      }else{
                        if(!is.numeric(columns())){
                          htmltools::withTags(table(
                            class = 'display',

                            thead(
                              tr(
                                th(rowspan = 2, class="Geneofficial hover", 'Gene*', span(class="Genecomment", "Official gene name from HGNC")),
                                th(rowspan = 2, 'Locus name**'),
                                th(rowspan = 2, class="Locusidwithcomment hover",'Locus ID', span(class="Locusidcomment", "k: known locus from Wuttke et al. Nat.commun.2019, n: novel locus in Stanzick et al. Nat.commun. 2021")),
                                th(rowspan = 2, class="Signalidwithcomment hover",'Signal ID', span(class="Signalidcomment", "ID of independent signals within a locus")),
                                th(rowspan = 2, class="validationwithcomment hover",'eGFRcys or BUN validation', span(class="validationcomment",list( "Information whether the lead variant of the respective is nominal significant (P<0.05) associated with concordent effect direction with", tags$abbr(title="glomerular filtration rate estimated from serum cystatin C","eGFRcys"),"or",tags$abbr(title="blood urea nitrogen","BUN")))),
                                th(rowspan = 2, '# credible variants in signal'),
                                th(rowspan = 2, class="maxPPA hover",'max PPA', span(class="maxPPAcomment", "max. probability of a credible variant in this signal to be causal")),
                                #th(colspan = length(which(columns()%in%c(11:13))), span('# deleterious credible variants in the gene (CADD-Phred',HTML("<span>&#8805;</span>"),'15)')),
                                #th(colspan = length(which(columns()%in%c(14:19))), '# credible variants, that modulate gene expression (eQTLs) or splicing (sQTLs), FDR < 5%'),
                                #lapply(names(GPS_show[columns()[which(columns()%in%c(20:22))]]), th, rowspan=2),
                                lapply(c('Diabetes specific effect','eGFRcrea decline effect','Drugable','In enriched pathway','Distance to locus lead variant','Chromosome','Position gene start','Position gene end'),th, rowspan=2)

                              )
                            )
                          )
                          )
                        }
                      }
                    }
                  }
                }
              }
            }
          }
   })
  
  
  
  output$GPSrows <-  DT::renderDataTable({
    
    
    datatable(GPS_gene(), rownames=FALSE,
              options = list(
                columnDefs = list(list(className = 'dt-left', targets = "_all"))
                , scrollX=TRUE,
                lengthChange=TRUE,
                dom = 'lfrtBip',
                buttons= list(c('copy', 'excel'))
              ),
              extensions = 'Buttons',
              escape=FALSE,
              caption = htmltools::tags$caption(
                style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                '* Link to GeneCards, **nearest gene to locus lead-variant'),
               container=sketch_gene()
    )%>% formatPercentage('max PPA')
  })
  
  #Locus summary
  output$summary <- renderDataTable({
        if(valid_gene()){
          if(any(detail_evidence_locus_gene()=="summary")){
            a=unique(GPS$locus_id[which(toupper(GPS$Gene)%in%genes_in_gps())])
            b=region_table_show[which(region_table$locus_id%in%a),]
            datatable(b[order(b[,1],b[,2]),], rownames=F,
                      options = list(
                        columnDefs = list(list(className = 'dt-left', targets = "_all"))
                        , scrollX=TRUE,
                        dom = 'lfrtBip',
                        buttons= list(c('copy', 'excel'))),
                      extensions = 'Buttons',
                      caption = htmltools::tags$caption(
                        style = 'caption-side: bottom; text-align: left; font-weight: normal;',
                        '* SNPs, which are associated with log(eGFRcrea) in the all-ancestry GWAS meta-analysis (unconditioned) with P<5E-8')
                      )
          }
        }
   })
  
  #Locus Zoom
  Locus_ID<-reactive({
    # if(inputtype$gene_final=="single")
    if(valid_gene())
    {GPS$locus_id[which(GPS$Gene%in%genes_in_gps())]}
  })
  filename1 <- reactive({
    if(valid_gene()&inputtype$gene_final=="single"){
      if(any(input$gene.locuszoom=="locus_zoom")){
        as.character(lz_files[which(locus_id_files==as.character(Locus_ID()))][1])
        
      }
    }
  })
  
  output$LocusZoom <- renderUI({
    req(filename1())
    
    tags$img(id="GPS_image",
               height="100%",
             width="100%",
             src=filename1())
  })
  
  #create HTML syntax for link to plots
  gene_plots <- reactive({
    if(valid_gene()){
      ids=which(links_id%in%Locus_ID())[1] ## preliminary only one plot
      paste("<p><a href=\"",links[ids],"\"target=\"_blank\">Regional association plot of locus",links_id[ids],"</a></p>")
    }
  })
  
  #render links
  output$plot_download_links_gene <- renderUI({
    list(lapply(gene_plots(),HTML))
  })
  
  
  observeEvent(input$question_GPS_gene,{
    showModal(modalDialog(
      titel="GPS table",
      tags$h4("GPS table help:"),
      p("Enter a gene name and choose every evidence you want to see in the GPS table. Numbers within the CADD-, eQTL- and sQTL-columns display the number of credible variants in the association signal with the respective characteristics. Numbers in the MGI and OMIM column represent the number of different phenotypes/diseases found for the respective gene. Fur further information check the Documentation page."),
      easyClose = TRUE,
      footer=NULL,
      size = "m"
    ))
  })
  
  observeEvent(input$question_details_gene,{
    showModal(modalDialog(
      titel="Details GPS Table",
      tags$h4("Details on GPS table - Help:"),
      p("Click here for details on the entries you find in the GPS table for a gene. If a gene has no entry in the respective column the selected option will not give you any result. Details can be shown also without selection the respective column in the GPS table.Fur further information check the Documentation page."),
      easyClose = TRUE,
      footer=NULL,
      size = "m"
    ))
  })
  
  #show relevant sections:
  observe({
    
    shinyjs::toggleElement(id="GPS_gene", condition = (rowsGPS()!=0) )
    
  })
  observeEvent(go$gene,{shinyjs::showElement(id="Gene_div")})
  observeEvent(go$gene,{runjs('
      document.getElementById("Gene_div").scrollIntoView({ left: 0, block: "end", behavior: "smooth" });
    ')})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_cadd_gene", condition= (any(columns()%in%c(11:13))&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_glo_gene", condition= (any(columns()==14)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_tub_gene", condition= (any(columns()==15)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_GTEx_eqtl_gene", condition= (any(columns()==16)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_GTEx_sqtl_gene", condition= (any(columns()==18)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_GTEx_wo_kidney_eqtl_gene", condition= (any(columns()==17)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_GTEx_wo_kidney_sqtl_gene", condition= (any(columns()==19)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_mgi_gene", condition= (any(columns()==20)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_omim_gene", condition= (any(columns()==21)&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_summary", condition= any(detail_evidence_locus_gene()=="summary"&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_cred_var_gene", condition= any(detail_evidence_locus_gene()=="cred_var"&valid_gene()))})
  observeEvent(go$gene,{shinyjs::toggleElement(id="div_locus_zoom_gene", condition= any(input$gene.locuszoom=="locus_zoom"&valid_gene()))})
  
  
}


shinyApp(ui = ui, server = server)