#
# Documentation here
#
# Attribution starts after this:
#
#    https://github.com/repository-address-once-it-goes-public
#


# load libraries (install if necessary)
packages = c("shiny", "shinyWidgets", "shinycustomloader", "bslib", "tidyr",
             "dplyr", "stringr", "BiocManager", "trackViewer", "RColorBrewer")

for (pack in packages) {
  if (pack %in% rownames(installed.packages())) {
    do.call("library", list(pack))
  } else {
    install.packages(pack)
    do.call("library", list(pack))
  }
}

# necessary for deployment
options(repos=BiocManager::repositories())

# Client-side UI
ui = fluidPage(

    # suppress warning messages
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"),

    # change to appropriate title
    titlePanel("Mute_Drug_Spot"),

    # for UI customization
    # (for a list of themes, check https://bootswatch.com)
    theme = bs_theme(version=4, bootswatch="lux"),

    # UI elements
    fluidRow(

      # custom UI for 3 column layout
      # left panel for data input
      column(3,
             wellPanel(
               fluidRow(
                 column(6,
                        numericInput("seq_length",
                                     h6("Length of Sequence"),
                                        value=200,
                                        min=2),
                        ),
                 column(6,
                        numericInput("num_domains",
                                     h6("Number of Domains"),
                                     value=1,
                                     min=1),
                        ),
               ),
               uiOutput("pos_inputs"),
               fileInput("mutation_list",
                         h5("List of Mutations"),
                         accept=c(".txt", ".csv", ".tsv")),
               selectInput("i_sep",
                           "Separator",
                           c("Comma (CSV)" = ',',
                             "Tab (TSV)" = '\t')),
               htmlOutput("warning"),
               hr(),
               h4("Cite:"),
               p("If you used this tool in your work, please cite it as:"),
               p("{citation}")
             )),

      # center panel for displaying+customizing plot
      column(6,

             # plot can take a while to load, so let the user know it's working
             withLoader(plotOutput("mutationPlot"), type="html", loader="dnaspin"),
             br(), br(), br(),

             # for borders, enclose the following block in a wellPanel()
             tabsetPanel(type="tabs",

                         # one tab for customizing labels
                         tabPanel(
                           h5("Customize Labels"),
                           fluidRow(
                             column(4,
                                    textInput("x_label",
                                              h5("Label of x-axis"),
                                              "x-axis")
                                    ),
                             column(4,
                                    textInput("y_label",
                                              h5("Label of y-axis"),
                                              "y-axis")
                                    ),
                             column(4,
                                    textInput("title",
                                              h5("Plot Title"),
                                              "title")
                                    )
                             ),
                           fluidRow(
                             column(5,
                                    sliderInput("title_height",
                                                h5("Title Position"),
                                                min=0,
                                                max=1,
                                                value=0.95)
                                    ),
                             column(3,
                                    h5("Shrink Plot"),
                                    checkboxInput("shrink",
                                                  "(useful for large number of drugs)",
                                                  value=FALSE)
                                    ),
                             column(4,
                                    downloadButton("download_plot",
                                                   "Download Plot",
                                                   style="display: flex;
                                                          align-items: center;
                                                          justify-content: center;")
                                      )
                             )
                           ),

                         # another one for customizing colors
                         tabPanel(
                           h5("Customize Colors"),

                           # handled server-side
                           uiOutput("drug_color_selector")
                         )
             )),

      # right panel for help
      column(3,
             wellPanel(
                       h4("Input File Format"),
                       p("Mutation List needs to be in CSV/TSV format, with columns",
                         " containing the mutation name and the corresponding",
                         " drug resistance respectively. If a mutation provides",
                         " resistance against a drug combination, add the name of",
                         " each drug separated by '+'. The file should look like: ",
                         style="text-align: justify"),
                       tableOutput("exampleList"),
                       p("Different types of mutations are illustrated in the plot",
                         "using different symbols, shown as:"),
                       tableOutput("mutationMap"),

                       # edit/remove accordingly
                       h4("Help:"),
                       p("For any query, please refer to:",
                         tags$a(href="https://github.com",
                                "{github-link}",
                                target="_blank")))),
      ),
)

# Server-side algorithm
server = function(input, output) {

    # house-keeping: delete any temporary files older than 24 hours
    for (f in list.files(pattern="\\.pdf$")) {

      # we name the files so that their first part is a timestamp
      timestamp = as.numeric(strsplit(f, split='_')[[1]][1])
      if (!is.na(timestamp)) {

        # check if that timestamp is > 24 hours old
        if ((as.numeric(Sys.time()) - as.integer(timestamp)) > 86400) {
          file.remove(f)
        }
      }
    }

    # a new row is added dynamically for each domain
    output$pos_inputs = renderUI({
      n = as.integer(input$num_domains)
      lapply(1:n, function(i) {
        fluidRow(
          column(3,
                 numericInput(paste0("dom_start", i),
                              h6("Start"),
                              value=1),
                 style="padding-right: 0px;"),
          column(3,
                 numericInput(paste0("dom_width", i),
                              h6("Width"),
                              value=200),
                 style="padding: 0px;"),
          column(6,
                 textInput(paste0("dom_name", i),
                           h6("Domain Name"), paste0("Domain ", i)))
          )
      })
    })

    # a warning will be displayed in case the list contains more than 12
    # unique drug names, since the resulting palette colors will not be
    # easily distinguishable
    test = reactive({
      req(input$mutation_list)
      check = length(unique(read.csv(input$mutation_list$datapath)$Resistance))
      ifelse(check > 12,1,0)
    })

    # the warning that will be displayed
    output$warning = renderText({
      if(test()==1){
        paste0("<b style='color: red;'>WARNING:</b><br>",
               "Your list contains more than 12 unique drug names. The default ",
               "colors may thus not be distinguishable from each other. Please ",
               "try using a smaller dataset or manually setting the colors.")
      }
    })

    # if no file is uploaded, a default plot will be created
    get_default = reactive({
      if (is.null(input$mutation_list)) {
        read.csv("example_mutations.csv")
      }
      else {
        df = read.csv(input$mutation_list$datapath, sep=input$i_sep)

        # return it after some clean-up: remove NAs, duplicates, and whitespaces
        df = na.omit(df) %>% mutate(across(where(is.character), str_trim))
        df[!duplicated(df), ]
      }
    })

    # generate default drug colors
    drug_colors = reactive({
      drug_list = unique(get_default()$Resistance)

      # RColorBrewer only supports up to 12 colors for making palette, so we
      # manually interpolate the color values to get a larger palette
      if (length(drug_list) > 12) {
        drug_col = data.frame(drug=drug_list,
                              color=colorRampPalette(brewer.pal(12, "Paired")
                                                     )(length(drug_list)))
      } else if (length(drug_list) < 3) {
        drug_col = data.frame(drug=drug_list,
                              color=brewer.pal(3, "Paired")[1:length(drug_list)])
      } else {
        drug_col = data.frame(drug=drug_list,
                              color=brewer.pal(length(drug_list), "Paired"))
      }
      drug_col
    })

    # unless the user wants to change the colors
    output$drug_color_selector = renderUI({

      # split by odd-even
      drug_list = split(unique(get_default()$Resistance),
                        rep(1:2,
                            length=length(unique(get_default()$Resistance))))
      colors = drug_colors()

      # split into 2 columns for aesthetics
      fluidRow(
        column(6,
               lapply(drug_list[[1]], function(x){
                 colorPickr(inputId=paste0("color_", x),
                            x,
                            width="100%",
                            selected=colors$color[colors$drug == x])
               })
               ),
        column(6,
               lapply(drug_list[[2]], function(x){
                 colorPickr(inputId=paste0("color_", x),
                            x,
                            width="100%",
                            selected=colors$color[colors$drug == x])
                 })
        )
      )
    })

    # make sure colors get loaded as soon as page is opened (needs to be forced)
    outputOptions(output, "drug_color_selector", suspendWhenHidden=FALSE)

    # to be accessible by different functions
    temp_fname = reactiveValues(a='')

    # create a random filename for a temporary file (static function)
    get_filename = function() {
        f = get_default()

        # so that we don't keep saving example plot with different names
        if (identical(f, read.csv("example_mutations.csv"))) {
          temp_fname$a = "default_plot.pdf"
        } else {

          if (temp_fname$a != "default_plot.pdf") {

            # check if it's been at least 2 seconds since we last saved a file
            prev_tstamp = as.integer(strsplit(temp_fname$a, split='_')[[1]][1])
            if ((as.numeric(Sys.time()) - as.integer(prev_tstamp)) < 2) {
              temp_fname$a = temp_fname$a
            }
          } else {

            # keep a few components in the name, for housekeeping
            temp_fname$a = paste0(as.integer(Sys.time()), # timestamp
                                  '_',
                                  input$title,            # plot title
                                  '_',
                                  do.call(paste0,         # random tail
                                          replicate(10,
                                                    sample(c(0:9,
                                                             letters),
                                                           1,
                                                           TRUE),
                                                    FALSE)),
                                  ".pdf")
          }
        }
    }

    # main code
    mutationPlot = reactive({

      # get input file
      mutation_list = get_default()

      # merge all domain starts, ends, and names into data frame
      n = as.integer(input$num_domains)
      domain_locs = unlist(lapply(1:n,
                                  function(i) {
                                    input[[paste0("dom_start", i)]]
                                    }))
      domain_widths = unlist(lapply(1:n,
                                    function(i) {
                                      input[[paste0("dom_width", i)]]
                                      }))
      domain_names = lapply(1:n,
                            function(i) {
                              input[[paste0("dom_name", i)]]
                              })

      # split mutation name into amino acids and locations; here's how it works:
      # first we look for the index where behind us is a letter and ahead of us
      # is a number; this is the location for first split (unless we hit an
      # in-frame insertion; there we leave first column as NA and move on), then
      # we look for the index such that behind us is a number and ahead of us
      # is a letter (or a '*' in case of a non-sense mutation); this is the
      # location for second split
      mutation_names = separate(mutation_list,
                                "Mutation",
                                into=c("orig", "rest"),
                                sep=c("(?<=[A-Za-z])(?=[0-9])"), fill="left")
      mutation_names = separate(mutation_names,
                                "rest",
                                into=c("loc", "new"),
                                sep=("(?<=[0-9])(?=[A-Za-z*])"))

      # only get the drug names for now
      drug_names = unique(mutation_list$Resistance)

      # fetch the color palette for drugs
      get_color_pal = reactive({data.frame(drug=drug_names,
                                      color=unlist(lapply(drug_names,
                                                          function(x) {
                                                            input[[paste0("color_", x)]]
                                                          }))
        )
      })

      color_pal = get_color_pal()

      # map these colors back to mutations
      mut_colors = merge(color_pal,
                         mutation_names,
                         by.x="drug",
                         by.y="Resistance",
                         all=TRUE)

      # now assign mutation type to each mutation
      mut_colors$shape = rep("circle", length(mut_colors$loc))
      mut_colors[mut_colors$new == "fs", "shape"] = "diamond"
      mut_colors[mut_colors$new == '*', "shape"] = "square"
      mut_colors[is.na(mut_colors$orig), "shape"] = "triangle_point_down"

      # in-frame insertions need more work
      if (length(mut_colors[is.na(mut_colors$orig), "loc"])) {
        mut_colors[is.na(mut_colors$orig), "loc"] = strsplit(
                                                             mut_colors[is.na(mut_colors$orig),
                                                                        "loc"],
                                                             split='_')[[1]][1]
        mut_colors[is.na(mut_colors$orig), "orig"] = ''
      }

      # mutation order is now changed, now use the updated lists
      original_aa = mut_colors$orig
      PM = mut_colors$loc
      replaced_aa = mut_colors$new

      # create separate color palette for domains
      if (n < 3) {
        dom_color = brewer.pal(3, "Dark2")[1:n]
      } else if (n > 8) {
        dom_color = colorRampPalette(brewer.pal(8, "Dark2"))(n)
      } else {
        dom_color = brewer.pal(n, "Dark2")
      }

      # change plot height if user tells us to
      if (input$shrink) {
        max_h = 2
      }
      else {
        max_h = 10
      }

      # create the domains
      # (and add a dummy domain as big as the whole sequence; I couldn't find
      # another way to keep the whole sequence visible in the plot)
      features = GRanges("chr1",
                          IRanges(append(1, c(domain_locs)),
                                  width=append(input$seq_length, c(domain_widths)),
                                  names=append('', domain_names)),
                          fill=append("#FFFFFF", dom_color),
                          seqlengths=seqlengths(Seqinfo("chr1", c(input$seq_length))))

      # and put mutations on top of it
      mutations = GRanges("chr1",
                          IRanges(as.integer(PM),
                                  width=1,
                                  names=paste0(c(original_aa),
                                               PM,
                                               c(replaced_aa))),
                                  color=mut_colors$color,
                                  shape=mut_colors$shape,
                                  score=runif(length(PM))*max_h)

      # legend for the drugs used here (domain legend is made automatically)
      legend = list(labels=color_pal$drug, fill=color_pal$color)

      # create axis with ticks uniformly distributed to sequence length
      # (rounded to the nearest hundred)
      xaxis = round(seq(1,
                        round(input$seq_length,
                              digits=-(length(input$seq_length)+1)),
                              length.out=5)
                    )

      # if the last tick is too close to the axis length, skip the last one
      if ((xaxis[length(xaxis)] - input$seq_length) < 20) {
        xaxis = append(xaxis[1:length(xaxis)-1], input$seq_length)
      } else {

        # otherwise, just add a last tick for the total sequence length
        xaxis = append(xaxis, input$seq_length)
      }

      # make the plot and add labels
      lolliplot(mutations, features, yaxis=FALSE, rescale=FALSE,
                xaxis=xaxis, legend=legend, ylab=input$y_label)
      grid.text(input$x_label, x=.5, y=0, just="bottom")
      grid.text(input$title, x=.5, y=input$title_height, just="top",
                gp=gpar(cex=1.5, fontface="bold"))
    })

    # send the plotting function to renderer
    output$mutationPlot = renderPlot({

      # make a temporary filename
      get_filename()

      # make the plot
      mutationPlot()

      # also save a temporary copy of the plot as soon as it is made
      dev.copy2pdf(file=temp_fname$a, out.type="pdf")
      })

    # the example file used for default plot
    output$exampleList = renderTable({read.csv("example_mutations.csv")})

    # and map for mutation types
    output$mutationMap = renderTable({read.csv("mutation_map.csv")})

    # download the plot
    output$download_plot = downloadHandler(
      filename = function() {
        paste(input$title, ".pdf", sep='')
      },
      content = function(file) {

        # now just copy the temporary file made earlier
        file.copy(temp_fname$a, file, overwrite=TRUE)
      }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
