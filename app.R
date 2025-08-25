library(shiny)
library(data.table)
library(highcharter)

quickbinomtest <- function(x, n, conf.level = 0.95, mult = 100, minX = 0, minN = 50) {
  alpha <- (1 - conf.level)/2
  if(x < minX || n < minN) list(value = NA_real_, lci = NA_real_, uci = NA_real_) else
    list(value = x/n*mult,
         lci = if (x == 0) 0 * mult else qbeta(alpha, x, n - x + 1) * mult,
         uci = if (x == n) 1 * mult else qbeta(1 - alpha, x + 1, n - x) * mult)
}

roundwzero <- function(num, dig = 1) format(round(num, digits = dig), nsmall = dig)

formatwci <- function(value, lci, uci, cond, uom, dig = 1) {
  paste0(roundwzero(value, dig), uom,
         if(cond) paste0(" (95% CI: ", roundwzero(lci, dig), uom,
                         " - ", roundwzero(uci, dig), uom, ")") else "")
}

hcoptslang <- getOption("highcharter.lang")
hcoptslang$contextButtonTitle <- "Helyi menü"
hcoptslang$exitFullscreen <- "Kilépés a teljes képernyős módból"
hcoptslang$hideData <- "Adatok elrejtése"
hcoptslang$loading <- "Betöltés..."
hcoptslang$mainBreadcrumb <- "Fő ábra"
hcoptslang$noData <- "Nincs megjeleníthető adat"
hcoptslang$printChart <- "Ábra nyomtatása"
hcoptslang$viewData <- "Adatok megtekintése"
hcoptslang$viewFullscreen <- "Teljes képernyős nézet"
hcoptslang$months <- c("január", "február", "március", "április", "május","június", "július",
                       "augusztus", "szeptember", "október", "november", "december")
hcoptslang$shortMonths <- c("jan", "febr", "márc", "ápr", "máj", "jún", "júl", "aug", "szept",
                            "okt", "nov", "dec")
hcoptslang$weekdays <- c("vasárnap", "hétfő", "kedd", "szerda", "csütörtök", "péntek",
                         "szombat")
hcoptslang$shortWeekdays <- c("Vas", "Hét", "Ked", "Sze", "Csü", "Pén", "Szo", "Vas")
hcoptslang$exportButtonTitle <- "Exportál"
hcoptslang$printButtonTitle <- "Importál"
hcoptslang$rangeSelectorFrom <- "ettől"
hcoptslang$rangeSelectorTo <- "eddig"
hcoptslang$rangeSelectorZoom <- "mutat:"
hcoptslang$downloadPNG <- "Letöltés PNG képként"
hcoptslang$downloadJPEG <- "Letöltés JPEG képként"
hcoptslang$downloadPDF <- "Letöltés PDF dokumentumként"
hcoptslang$downloadSVG <- "Letöltés SVG formátumban"
hcoptslang$downloadCSV <- "Letöltés CSV formátumú táblázatként"
hcoptslang$downloadXLS <- "Letöltés XLS formátumú táblázatként"
hcoptslang$resetZoom <- "Nagyítás alaphelyzetbe állítása"
hcoptslang$resetZoomTitle <- "Nagyítás alaphelyzetbe állítása"
hcoptslang$thousandsSep <- " "
hcoptslang$decimalPoint <- ","
hcoptslang$numericSymbols <- NA
options(highcharter.lang = hcoptslang)
options(highcharter.download_map_data = FALSE)

mapdata <- readLines("hu-all-custom.geojson", warn = FALSE, encoding = "UTF-8")
mapdata[1] <- gsub(".* = ", "", mapdata[1])
mapdata <- paste(mapdata, collapse = "\n")
mapdata <- stringr::str_remove(mapdata, ";$")
mapdata <- jsonlite::fromJSON(mapdata, simplifyVector = FALSE)

AgeTable <- data.table(AGEcode = c("Y_LT5", paste0("Y", seq(5, 85, 5), "-", seq(9, 89, 5)),
                                   "Y_GE90"),
                       AGE = c("<5", paste0(seq(5, 85, 5), "-", seq(9, 89, 5)), ">=90"))
AgeTable$AGE <- factor(AgeTable$AGE, levels = AgeTable$AGE)

NameTable <- data.table(
  variable = c("AMI", "STEMI", "NSTEMI", "MYOCARDIALIS", "SZIVELEGTELENSEG", "COPD",
               "HYPERTONIA", "STROKE", "DIABETES", "PERIFERIALIS_ERBETEGSEG", "HYPERLIPIDAEMIA",
               "ESEMENY_SZINTU_PCI", "MORT30DAY", "MORT1YEAR", "DATEFORMATTED", "MEGYE", "NEM",
               "AGE", "KORHAZI_DIAGNOZIS",
               "varname", "valueformatted"),
  varname = c("Heveny szívinfarktus", "STEMI", "NSTEMI", "Kórelőzményben szereplő infarktus",
              "Szívelégtelenség", "Krónikus obstruktív tüdőbetegség", "Magasvérnyomás-betegség",
              "Kórelőzményben szereplő stroke", "Cukorbetegség", "Perifériás érbetegség",
              "Magas vérzsír-szint", "Katéteres érmegnyitás (PCI)", "30 napos halálozás",
              "1 éves halálozás", "Időszak", "Megye", "Nem", "Életkor", "Diagnózis",
              "varname", "valueformatted"),
  type = c("incidence", "incidence", "incidence", "anamnestic", "anamnestic", "anamnestic",
           "comorb", "anamnestic", "comorb", "comorb", "comorb", "treatment", "survival",
           "survival", "header", "header", "header", "header", "header",
           "header", "header")
)

leirasurl <- paste0("https://github.com/tamas-ferenci/NSZR-Lekerdezo?tab=readme-ov-file#",
                    "a-nemzeti-sz%C3%ADvinfarktusregiszter-nszr-",
                    "interakt%C3%ADv-lek%C3%A9rdez%C5%91-fel%C3%BClete")

RawData <- readRDS("RawData.rds")
datadate <- RawData$datadate
deathcutoffdate <- RawData$deathcutoffdate
RawData <- RawData$RawData

PopData <- readRDS("PopData.rds")

ESP2013 <- data.table(
  AGEcode = c("Y_LT5", "Y5-9", "Y10-14", "Y15-19", "Y20-24", "Y25-29", "Y30-34", "Y35-39",
              "Y40-44", "Y45-49", "Y50-54", "Y55-59", "Y60-64", "Y65-69", "Y70-74", "Y75-79",
              "Y80-84", "Y85-89", "Y_GE90"),
  STDPOP = c(5000L, 5500L, 5500L, 5500L, 6000L, 6000L, 6500L, 7000L, 7000L, 7000L, 7000L, 6500L,
             6000L, 5500L, 5000L, 4000L, 2500L, 1500L, 1000L)
)

StdPopData <- merge(ESP2013, AgeTable, by = "AGEcode", sort = FALSE)[, .(AGE, STDPOP)]

pickeropts <- shinyWidgets::pickerOptions(
  actionsBox = TRUE,
  liveSearch = TRUE,
  noneSelectedText = "Válasszon!",
  noneResultsText = "Nincs találat {0}",
  countSelectedText = "{0} elem kiválasztva",
  maxOptionsText ="Legfeljebb {n} elem választható",
  selectAllText = "Mindegyik",
  deselectAllText = "Egyik sem",
  multipleSeparator = ", ",
  style = "btn btn-outline-dark"
)

pickeroptsWOSearch <- shinyWidgets::pickerOptions(
  actionsBox = TRUE,
  liveSearch = FALSE,
  noneSelectedText = "Válasszon!",
  noneResultsText = "Nincs találat {0}",
  countSelectedText = "{0} elem kiválasztva",
  maxOptionsText ="Legfeljebb {n} elem választható",
  selectAllText = "Mindegyik",
  deselectAllText = "Egyik sem",
  multipleSeparator = ", ",
  style = "btn btn-outline-dark"
)

pickeroptsWOSearchWOACtion <- shinyWidgets::pickerOptions(
  actionsBox = FALSE,
  liveSearch = FALSE,
  noneSelectedText = "Válasszon!",
  noneResultsText = "Nincs találat {0}",
  countSelectedText = "{0} elem kiválasztva",
  maxOptionsText ="Legfeljebb {n} elem választható",
  selectAllText = "Mindegyik",
  deselectAllText = "Egyik sem",
  multipleSeparator = ", ",
  style = "btn btn-outline-dark"
)

ownpanel <- function(id, metrics, primary, primarylabel, primaryvalues,
                     timestrat, comorbSelEnable, pciSelEnable, spaceLimitYear) {
  c(
    list(
      radioButtons(paste0(id, "Subject"), "Vizsgálat tárgya",
                   c("Időbeli alakulás" = "time", "Területi alakulás" = "space")),
      selectInput(paste0(id, "Metric"), "Mutató", metrics),
      conditionalPanel(
        paste0("input.", id, "Subject == 'time'"), 
        shinyWidgets::pickerInput(
          paste0(id, "Time", primary),
          div(primarylabel, bslib::tooltip(
            bsicons::bs_icon("question-circle"),
            paste0("Egy idejűleg több elem is kiválasztható, de ez esetben nem érhető el lebontás"),
            placement = "right")),
          primaryvalues, unlist(primaryvalues)[1], multiple = TRUE, options = if(length(primaryvalues) <= 3) pickeroptsWOSearchWOACtion else pickeroptsWOSearch),
        conditionalPanel(paste0("input.", id, "Time", primary, ".length == 1 & (input.", id, "Metric == 'absolute' | input.", id, "Metric == 'cruderate' | input.", id, "Metric == 'crude')"),
                         selectInput(paste0(id, "TimeStrat"), "Lebontás", timestrat)),
        conditionalPanel(paste0("input.", id, "Time", primary, ".length == 1 & input.", id, "Metric == 'adjrate'"),
                         selectInput(paste0(id, "TimeStratAdjrate"), "Lebontás", timestrat[timestrat != "AGE"])),
        conditionalPanel(paste0("input.", id, "Time", primary, ".length == 1 & input.", id, "Metric == 'agesexadj'"),
                         selectInput(paste0(id, "TimeStratAgesexadj"), "Lebontás", timestrat[!timestrat %in% c("AGE", "NEM")])),
        conditionalPanel(paste0("input.", id, "Time", primary, ".length == 1 & input.", id, "Metric == 'agesexcomorbadj'"),
                         selectInput(paste0(id, "TimeStratAgesexcomorbadj"), "Lebontás", timestrat[!timestrat %in% c("AGE", "NEM", NameTable[type %in% c("comorb", "anamnestic")]$variable)])),
        conditionalPanel(
          paste0("(input.", id, "Time", primary, ".length == 1 & input.", id, "Metric != 'adjrate' & ",
                 "input.", id, "Metric != 'agesexadj' & input.", id, "Metric != 'agesexcomorbadj' & ",
                 "(input.", id, "TimeStrat == 'AGE' | input.", id, "TimeStrat == 'MEGYE')) | ",
                 "(input.", id, "Time", primary, ".length == 1 & input.", id, "Metric == 'adjrate' & input.", id, "TimeStratAdjrate == 'MEGYE')"),
          actionButton(paste0(id, "ShowButton"), "Minden görbe megjelenítése"),
          actionButton(paste0(id, "HideButton"), "Minden görbe elrejtése")),
        checkboxInput(paste0(id, "TimeSelEnable"), "Szűkítés"),
        conditionalPanel(paste0("input.", id, "TimeSelEnable == 1"),
                         conditionalPanel(paste0("['crude', 'absolute', 'cruderate', 'adjrate'].includes(input.", id, "Metric)"),
                                          fluidRow(
                                            column(width = 4, checkboxInput(paste0(id, "TimeSexSelEnable"), "Nem szerint")),
                                            column(width = 8, conditionalPanel(paste0("input.", id, "TimeSexSelEnable == 1"),
                                                                               selectInput(paste0(id, "TimeSexSel"), NULL,
                                                                                           unique(RawData$NEM)))))),
                         conditionalPanel(paste0("['crude', 'absolute', 'cruderate'].includes(input.", id, "Metric)"),
                                          fluidRow(
                                            column(width = 4, checkboxInput(paste0(id, "TimeAgeSelEnable"), "Életkor szerint")),
                                            column(width = 8, conditionalPanel(paste0("input.", id, "TimeAgeSelEnable == 1"),
                                                                               selectInput(paste0(id, "TimeAgeSel"), NULL,
                                                                                           sort(unique(RawData$AGE))))))),
                         fluidRow(
                           column(width = 4, checkboxInput(paste0(id, "TimeCountySelEnable"), "Megye szerint")),
                           column(width = 8, conditionalPanel(paste0("input.", id, "TimeCountySelEnable == 1"),
                                                              selectInput(paste0(id, "TimeCountySel"), NULL,
                                                                          sort(unique(RawData$MEGYE)))))),
                         conditionalPanel(paste0("['crude', 'agesexadj', 'agesexcomorbadj'].includes(input.", id, "Metric)"),
                                          fluidRow(
                                            column(width = 4, checkboxInput(paste0(id, "TimeDgSelEnable"), "Diagnózis szerint")),
                                            column(width = 8, conditionalPanel(paste0("input.", id, "TimeDgSelEnable == 1"),
                                                                               selectInput(paste0(id, "TimeDgSel"), NULL,
                                                                                           c("STEMI", "NSTEMI")))))),
                         if(pciSelEnable) {
                           fluidRow(
                             column(width = 4, checkboxInput(paste0(id, "TimePciSelEnable"), "Katéteres érmegnyitás (PCI) szerint")),
                             column(width = 8, conditionalPanel(paste0("input.", id, "TimePciSelEnable == 1"),
                                                                selectInput(paste0(id, "TimePciSel"), NULL,
                                                                            unique(RawData$ESEMENY_SZINTU_PCI)))))
                         },
                         if(comorbSelEnable) {
                           conditionalPanel(paste0("input.", id, "Metric != 'agesexcomorbadj'"),
                                            lapply(1:nrow(NameTable[type%in%c("comorb", "anamnestic")]), function(i) {
                                              fluidRow(
                                                column(width = 5, checkboxInput(paste0(id, "Time", NameTable[type%in%c("comorb", "anamnestic")]$variable[i], "SelEnable"),
                                                                                paste0(NameTable[type%in%c("comorb", "anamnestic")]$varname[i], " szerint"))),
                                                column(width = 7, conditionalPanel(paste0("input.", id, "Time", NameTable[type%in%c("comorb", "anamnestic")]$variable[i], "SelEnable == 1"),
                                                                                   selectInput(paste0(id, "Time", NameTable[type%in%c("comorb", "anamnestic")]$variable[i], "Sel"), NULL,
                                                                                               c("Igen", "Nem")))))
                                            }))
                         }),
        selectInput(paste0(id, "TimeFreq"), "Időbeli sűrűség",
                    c("Évi" = "year", "Havi" = "month")),
        conditionalPanel(paste0("input.", id, "Metric != 'absolute'"),
                         checkboxInput(paste0(id, "TimeCI"), "Konfidenciaintervallum feltüntetése")),
        checkboxInput(paste0(id, "TimeIncludeZero"), "Függőleges tengely 0-nál kezdődik"),
        checkboxInput(paste0(id, "TimeDisplayallvalues"), "Minden érték megjelenítése")
      ),
      conditionalPanel(
        paste0("input.", id, "Subject == 'space'"),
        shinyWidgets::pickerInput(paste0(id, "Space", primary), primarylabel,
                                  primaryvalues, primaryvalues[1], multiple = FALSE,
                                  options = pickeroptsWOSearch),
        radioButtons(paste0(id, "SpacePlotType"), "Ábrázolás módja",
                     c("Térkép" = "map", "Oszlopdiagram" = "barchart")),
        conditionalPanel(paste0("input.", id, "SpacePlotType == 'barchart'"),
                         checkboxInput(paste0(id, "SpacePlotBarOrder"), "Nagyság szerint sorbarendezve"),
                         checkboxInput(paste0(id, "SpacePlotBarHorizontal"), "Vízszintes diagram")),
        conditionalPanel(paste0("input.", id, "Metric != 'absolute'"),
                         checkboxInput(paste0(id, "SpaceCI"), "Konfidenciaintervallum feltüntetése")),
        checkboxInput(paste0(id, "SpaceDisplayallvalues"), "Minden érték megjelenítése"),
        checkboxInput(paste0(id, "SpaceFixedcoloraxis"), "Rögzített színskála"),
        conditionalPanel(paste0("input.", id, "SpaceFixedcoloraxis == 1"),
                         numericInput(paste0(id, "SpaceFixedcoloraxisMin"), "Minimum", 0, 0),
                         numericInput(paste0(id, "SpaceFixedcoloraxisMax"), "Maximum", 100, 0)),
        tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
        sliderInput(paste0(id, "SpaceYearSel"), div("Vizsgált év(tartomány)", bslib::tooltip(
          bsicons::bs_icon("question-circle"),
          paste0("Amennyiben a csúszka két vége nem esik egybe, úgy a megjelenített ",
                 "adat a tartomány összesített eredménye. A két végpont egy évre is ",
                 "összehúzható, ez esetben a kérdéses év adata fog látszódni."),
          placement = "right")), min(lubridate::year(RawData$DATE)),
          if(!spaceLimitYear) max(lubridate::year(RawData$DATE)) else lubridate::year(datadate) - 1,
          c(2019, 2023), 1, sep = "", width = "100%"),
        checkboxInput(paste0(id, "SpaceSelEnable"), "Szűkítés"),
        conditionalPanel(paste0("input.", id, "SpaceSelEnable == 1"),
                         conditionalPanel(paste0("['crude', 'absolute', 'cruderate', 'adjrate'].includes(input.", id, "Metric)"),
                                          fluidRow(
                                            column(width = 4, checkboxInput(paste0(id, "SpaceSexSelEnable"), "Nem szerint")),
                                            column(width = 8, conditionalPanel(paste0("input.", id, "SpaceSexSelEnable == 1"),
                                                                               selectInput(paste0(id, "SpaceSexSel"), NULL,
                                                                                           unique(RawData$NEM)))))),
                         conditionalPanel(paste0("['crude', 'absolute', 'cruderate'].includes(input.", id, "Metric)"),
                                          fluidRow(
                                            column(width = 4, checkboxInput(paste0(id, "SpaceAgeSelEnable"), "Életkor szerint")),
                                            column(width = 8, conditionalPanel(paste0("input.", id, "SpaceAgeSelEnable == 1"),
                                                                               selectInput(paste0(id, "SpaceAgeSel"), NULL,
                                                                                           sort(unique(RawData$AGE))))))),
                         conditionalPanel(paste0("['crude', 'agesexadj', 'agesexcomorbadj'].includes(input.", id, "Metric)"),
                                          fluidRow(
                                            column(width = 5, checkboxInput(paste0(id, "SpaceDgSelEnable"), "Diagnózis szerint")),
                                            column(width = 7, conditionalPanel(paste0("input.", id, "SpaceDgSelEnable == 1"),
                                                                               selectInput(paste0(id, "SpaceDgSel"), NULL,
                                                                                           c("STEMI", "NSTEMI")))))),
                         if(pciSelEnable) {
                           fluidRow(
                             column(width = 4, checkboxInput(paste0(id, "SpacePciSelEnable"), "Katéteres érmegnyitás (PCI) szerint")),
                             column(width = 8, conditionalPanel(paste0("input.", id, "SpacePciSelEnable == 1"),
                                                                selectInput(paste0(id, "SpacePciSel"), NULL,
                                                                            unique(RawData$ESEMENY_SZINTU_PCI)))))
                         },
                         if(comorbSelEnable) {
                           conditionalPanel(paste0("input.", id, "Metric != 'agesexcomorbadj'"),
                                            lapply(1:nrow(NameTable[type%in%c("comorb", "anamnestic")]), function(i) {
                                              fluidRow(
                                                column(width = 5, checkboxInput(paste0(id, "Space", NameTable[type%in%c("comorb", "anamnestic")]$variable[i], "SelEnable"),
                                                                                paste0(NameTable[type%in%c("comorb", "anamnestic")]$varname[i], " szerint"))),
                                                column(width = 7, conditionalPanel(paste0("input.", id, "Space", NameTable[type%in%c("comorb", "anamnestic")]$variable[i], "SelEnable == 1"),
                                                                                   selectInput(paste0(id, "Space", NameTable[type%in%c("comorb", "anamnestic")]$variable[i], "Sel"), NULL,
                                                                                               c("Igen", "Nem")))))
                                            }))
                         }
        )
      )
    )
  )
}

ui <- navbarPage(
  theme = bslib::bs_theme(bootswatch = "default"),
  title = "Nemzeti Szívinfarktus Regiszter",
  header = list(
    p("A program használatát részletesen bemutató súgó, valamint a technikai részletek",
      a("itt", href = leirasurl, target = "_blank" ), "olvashatóak el.")
  ),
  
  footer = list(
    hr(),
    p(paste0("Az adatok lezárásának dátuma: ", format(datadate, "%Y. %m. %d"), ". A halálozási ",
             "adatok teljeskörűek a következő dátumig: ",
             format(deathcutoffdate, "%Y. %m. %d.")))
  ),
  
  tabPanel(
    shinyjs::useShinyjs(),
    title = "Előfordulás",
    sidebarLayout(
      sidebarPanel(
        ownpanel("incidence", c("Abszolút szám" = "absolute", "Nyers ráta" = "cruderate",
                                "Standardizált ráta" = "adjrate"),
                 "Incidence", "Betegség",
                 c("Heveny szívinfarktus" = "AMI", "STEMI" = "STEMI",
                   "NSTEMI" = "NSTEMI"),
                 c("Nincs" = "None", "Nem" = "NEM", "Életkor" = "AGE", "Megye" = "MEGYE"),
                 FALSE, FALSE, TRUE)
      ),
      
      mainPanel(
        bslib::navset_tab(
          bslib::nav_panel("Ábra", shinycssloaders::withSpinner(highchartOutput("incidencePlot"))),
          bslib::nav_panel("Táblázat", DT::DTOutput("incidenceTab")),
          bslib::nav_panel("Magyarázat", withMathJax(includeMarkdown("incidenceExplanation.md")))
        )
      )
    )
  ),
  
  tabPanel(
    title = "Társbetegségek",
    sidebarLayout(
      sidebarPanel(
        ownpanel("comorb", c("Nyers arány" = "crude"), "Comorb", "Kórelőzmény, társbetegség",
                 list(`Kórelőzményben szereplő betegség` = setNames(NameTable[type == "anamnestic"]$variable, NameTable[type == "anamnestic"]$varname),
                      `Társbetegség` = setNames(NameTable[type == "comorb"]$variable, NameTable[type == "comorb"]$varname)),
                 c("Nincs" = "None", "Nem" = "NEM", "Életkor" = "AGE", "Megye" = "MEGYE", "Diagnózis" = "KORHAZI_DIAGNOZIS"),
                 FALSE, FALSE, FALSE)
      ),
      
      mainPanel(
        bslib::navset_tab(
          bslib::nav_panel("Ábra", shinycssloaders::withSpinner(highchartOutput("comorbPlot"))),
          bslib::nav_panel("Táblázat", DT::DTOutput("comorbTab")),
          bslib::nav_panel("Magyarázat", withMathJax(includeMarkdown("comorbExplanation.md")))
        )
      )
    )
  ),
  
  tabPanel(
    title = "Ellátás",
    sidebarLayout(
      sidebarPanel(
        ownpanel("treatment",
                 c("Nyers arány" = "crude", "Életkorra és nemre korrigált arány" = "agesexadj",
                   "Életkorra, nemre és társbetegségekre korrigált arány" = "agesexcomorbadj"),
                 "Treatment", "Ellátási jellemző",
                 c("Katéteres érmegnyitás (PCI)" = "ESEMENY_SZINTU_PCI"),
                 c("Nincs" = "None", "Nem" = "NEM", "Életkor" = "AGE", "Megye" = "MEGYE", "Diagnózis" = "KORHAZI_DIAGNOZIS",
                   setNames(NameTable[type %in%c("anamnestic", "comorb")]$variable, NameTable[type %in%c("anamnestic", "comorb")]$varname)),
                 TRUE, FALSE, FALSE)
      ),
      
      mainPanel(
        bslib::navset_tab(
          bslib::nav_panel("Ábra", shinycssloaders::withSpinner(highchartOutput("treatmentPlot"))),
          bslib::nav_panel("Táblázat", DT::DTOutput("treatmentTab")),
          bslib::nav_panel("Magyarázat", withMathJax(includeMarkdown("treatmentExplanation.md")))
        )
      )
    )
  ),
  
  tabPanel(
    title = "Túlélés",
    sidebarLayout(
      sidebarPanel(
        ownpanel("survival",
                 c("Nyers arány" = "crude", "Életkorra és nemre korrigált arány" = "agesexadj",
                   "Életkorra, nemre és társbetegségekre korrigált arány" = "agesexcomorbadj"),
                 "Survival", "Túlélési mutató",
                 c("30 napos halálozás" = "MORT30DAY", "1 éves halálozás" = "MORT1YEAR"),
                 c("Nincs" = "None", "Nem" = "NEM", "Életkor" = "AGE", "Megye" = "MEGYE",
                   "Diagnózis" = "KORHAZI_DIAGNOZIS", "Katéteres érmegnyitás (PCI)" = "ESEMENY_SZINTU_PCI", 
                   setNames(NameTable[type %in%c("anamnestic", "comorb")]$variable, NameTable[type %in%c("anamnestic", "comorb")]$varname)),
                 TRUE, TRUE, FALSE)
      ),
      
      mainPanel(
        bslib::navset_tab(
          bslib::nav_panel("Ábra", shinycssloaders::withSpinner(highchartOutput("survivalPlot"))),
          bslib::nav_panel("Táblázat", DT::DTOutput("survivalTab")),
          bslib::nav_panel("Magyarázat", withMathJax(includeMarkdown("survivalExplanation.md")))
        )
      )
    )
  )
)

server <- function(input, output) {
  
  observeEvent("", {
    showModal(modalDialog(
      h1("Nemzeti Szívinfarktus Regiszter - Interaktív Lekérdező"),
      p(paste0("Ez a weboldal lehetőséget biztosít a Nemzeti Szívinfarktus Regiszter aggregált ",
               "adatainak interaktív lekérdezésére, elemzésére és vizualizálására. ",
               "Az adatok lehetővé teszik a szívinfarktus hazai helyzetének, időbeli és ",
               "térbeli alakulásának vizsgálatát többféle, népegészségügyileg fontos ",
               "vetületben (előfordulás, társbetegségek, az ellátás jellemzői, túlélés).")),
      p(paste0("Az adatok értelmezése, valamint a különféle összehasonlítások esetén ",
               "érdemes bármilyen következtetés ",
               "levonása előtt tanulmányozni az oldalhoz kapcsolódó "),
        a(href = leirasurl, target = "_blank", "leírást", .noWS = "outside"),
        ", mely igyekszik közérthetően összefoglalni a legfontosabb szempontokat. ",
        "Ezen kívül a weboldal maga is biztosít magyarázatot minden elemzési vetülethez; ",
        "ezen szempontok ismerete fontos a helyes következtetések levonásához."),
      p(paste0("Ugyanezen a linken elérhető részletes technikai magyarázat is a weboldal ",
               "működéséhez, valamint a transzparencia jegyében letölthető a weboldal ",
               " és a számításokat végző program teljes forráskódja.")),
      p(paste0("Az alapadatok a Nemzeti Szívinfarktus Regiszterből származnak (Gottsegen György ",
               "Országos Kardiovaszkuláris Intézet, szakmai vezető: ",
               "Prof. Dr. Jánosi András), a weboldalt készítette és a statisztikai ",
               "számításokat végezte Ferenci Tamás (Óbudai Egyetem).")),
      easyClose = TRUE,
      footer = tagList(
        actionButton(inputId = "opening", label = "Bezárás", icon = icon("info-circle"))
      )
    ))
  })
  
  observeEvent(input$opening, {
    removeModal()
  })
  
  observeEvent(input$incidenceMetric, {
    updateNumericInput(inputId = "incidenceSpaceFixedcoloraxisMax",
                       value = if(input$incidenceMetric == "absolute") 20000 else 200)
  })
  
  observeEvent(input$incidenceShowButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.show();});")
  })
  observeEvent(input$incidenceHideButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.hide();});")
  })
  observeEvent(input$comorbShowButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.show();});")
  })
  observeEvent(input$comorbHideButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.hide();});")
  })
  observeEvent(input$treatmentShowButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.show();});")
  })
  observeEvent(input$treatmentHideButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.hide();});")
  })
  observeEvent(input$survivalShowButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.show();});")
  })
  observeEvent(input$survivalHideButton, {
    shinyjs::runjs("$.each(Highcharts.charts.slice(-1)[0].series, function(i, ser) {ser.hide();});")
  })
  
  dataInput <- function(metric, primary, subject, strat, freq, ci,
                        yearSel, sexSel, dgSel, ageSel, comorbSel, pciSel, countySel,
                        uom, titleText) {
    dat <- RawData
    
    if(!is.null(yearSel)) dat <- dat[lubridate::year(DATE) >= yearSel[1] & 
                                       lubridate::year(DATE) <= yearSel[2]]
    if(!is.null(sexSel)) dat <- dat[NEM == sexSel]
    if(!is.null(dgSel)) dat <- dat[KORHAZI_DIAGNOZIS == dgSel]
    if(!is.null(ageSel)) dat <- dat[AGE == ageSel]
    if(!is.null(comorbSel))
      for(i in 1:nrow(NameTable[type %in% c("comorb", "anamnestic")]))
        if(!is.na(comorbSel[i])) dat <- dat[ dat[[ NameTable[type %in% c("comorb", "anamnestic")]$variable[i] ]] == comorbSel[i] ]
    if(!is.null(pciSel)) dat <- dat[ESEMENY_SZINTU_PCI == pciSel]
    if(!is.null(countySel)) dat <- dat[MEGYE == countySel]
    
    dat$DATE <- switch(freq, "year" = dat$YEAR, "month" = dat$YEARMON)
    
    if(strat == "None") strat <- NULL
    stratvar <- c(subject, strat)
    
    dat <- rbindlist(lapply(primary, function(prim) {
      if(metric == "crude") {
        temp <- cbind(dat[!is.na(get(prim)), quickbinomtest(sum(get(prim) == "Igen"), .N),
                          stratvar])[order(get(subject))]
        if(subject == "MEGYE")
          temp <- merge(CJ(MEGYE = unique(RawData$MEGYE)), temp, by = "MEGYE", all.x = TRUE)
        temp$variable <- prim
        # if(length(primary) == 1 && strat != "None") temp <- temp[order(get(strat))]
        temp
      } else if(metric %in% c("agesexadj", "agesexcomorbadj")) {
        as.data.table(emmeans::emmeans(
          glm(as.formula(switch(
            metric,
            "agesexadj" = paste0(prim, " == 'Igen' ~ as.factor(", paste0(stratvar, collapse = ") + as.factor("), ") + NEM + AGEcont"),
            "agesexcomorbadj" = paste0(prim, " == 'Igen' ~ as.factor(", paste0(stratvar, collapse = ") + as.factor("), ") + NEM + ",
                                       "AGEcont + MYOCARDIALIS + SZIVELEGTELENSEG + ",
                                       "HYPERTONIA + STROKE + DIABETES"))),
            data = dat, family = binomial(link = "logit")),
          stratvar, type = "response", weights = "proportional", rg.limit = 15000))[
            , c(.SD, .(value = if(nrow(dat) < 50) NA_real_ else response * 100,
                       lci = if(nrow(dat) < 50) NA_real_ else asymp.LCL * 100,
                       uci = if(nrow(dat) < 50) NA_real_ else asymp.UCL * 100, variable = prim))]
      } else {
        temp <- dat
        if(prim != "AMI") temp <- temp[KORHAZI_DIAGNOZIS == prim]
        
        pd <- PopData
        if(!is.null(sexSel)) pd <- pd[NEM == sexSel]
        if(!is.null(ageSel)) pd <- pd[AGE == ageSel]
        if(!is.null(countySel)) pd <- pd[MEGYE == countySel]
        
        temp <- merge(CJ(AGE = if(is.null(ageSel)) AgeTable$AGE else ageSel,
                         DATE = unique(switch(freq,
                                              "year" = if(is.null(yearSel)) RawData$YEAR else RawData[lubridate::year(DATE) >= yearSel[1] & lubridate::year(DATE) <= yearSel[2]]$YEAR,
                                              "month" = if(is.null(yearSel)) RawData$YEARMON else RawData[lubridate::year(DATE) >= yearSel[1] & lubridate::year(DATE) <= yearSel[2]]$YEARMON)),
                         NEM = if(is.null(sexSel)) unique(RawData$NEM) else sexSel,
                         MEGYE = if(is.null(countySel)) unique(RawData$MEGYE) else countySel),
                      temp[, .N, .(AGE, DATE, NEM, MEGYE)],
                      by = c("AGE", "DATE", "NEM", "MEGYE"), all.x = TRUE)[
                        , .(N = sum(N, na.rm = TRUE)), c(union(stratvar, c("DATE", "AGE")))]
        
        # if(subject == "MEGYE")
        #   temp <- merge(CJ(MEGYE = unique(RawData$MEGYE)), temp, by = "MEGYE", all.x = TRUE)
        
        temp <- merge(temp, pd[FREQ == freq, .(POP = sum(POP)), c(union(stratvar, c("DATE", "AGE")))],
                      by = union(stratvar, c("DATE", "AGE")))
        
        temp <- switch(metric,
                       "absolute" = temp[, .(value = if(sum(N) < 50) NA_integer_ else sum(N), lci = NA_integer_, uci = NA_integer_), stratvar],
                       "cruderate" = temp[, quickbinomtest(sum(N), sum(POP), mult = 1e5, minX = 50), stratvar],
                       "adjrate" = merge(temp, StdPopData, by = "AGE")[
                         , if(sum(N) < 50) list(value = NA_real_, lci = NA_real_, uci = NA_real_) else
                           with(as.list(epitools::ageadjust.direct(N, POP, stdpop = STDPOP)),
                                list(value = adj.rate * 1e5, lci = lci * 1e5, uci = uci * 1e5)),
                         stratvar])
        
        temp$variable <- prim
        
        temp
      }
    }))
    
    if(subject == "DATE") {
      dat <- dat[order(DATE)]
      dat$DATEFORMATTED <- format(dat$DATE, switch(freq, "year" = "%Y", "month" = "%Y. %m."))
    }
    
    if(!is.null(yearSel) && subject == "MEGYE")
      dat$DATEFORMATTED <- if(diff(yearSel) == 0) yearSel[1] else paste0(yearSel, collapse = " - ")
    if(!is.null(sexSel)) dat$NEM <- sexSel
    if(!is.null(dgSel)) dat$KORHAZI_DIAGNOZIS <- dgSel
    if(!is.null(ageSel)) dat$AGE <- ageSel
    if(!is.null(comorbSel))
      for(i in 1:nrow(NameTable[type %in% c("comorb", "anamnestic")]))
        if(!is.na(comorbSel[i])) dat[[ NameTable[type %in% c("comorb", "anamnestic")]$variable[i] ]] <- comorbSel[i]
    if(!is.null(pciSel)) dat$ESEMENY_SZINTU_PCI <- pciSel
    
    dat$valueformatted <- ifelse(is.na(dat$value), NA,
                                 formatwci(dat$value, dat$lci, dat$uci, ci && metric != "absolute",
                                           uom, if(metric == "absolute") 0 else 1))
    
    list(data = dat,
         title = paste0("<b>", if(subject == "DATE" && length(primary) > 1) "Adott" else NameTable[variable %in% primary]$varname, " ", titleText,
                        switch(metric,
                               "absolute" = " (esetszám)",
                               "cruderate" = " (nyers ráta)",
                               "adjrate" = " (standardizált ráta)",
                               "crude" = " (nyers arány)",
                               "agesexadj" = " (életkorra és nemre korrigált arány)",
                               "agesexcomorbadj" = " (életkorra, nemre és társbetegségekre korrigált arány)"),
                        if(!is.null(yearSel)) paste0(", ", if(yearSel[1] == yearSel[2]) yearSel[1] else paste0(yearSel, collapse = " - ")),
                        "</b><br>",
                        if(!is.null(sexSel)) paste0("Nem: ", sexSel),
                        if(!is.null(ageSel)) paste0(" Életkor: ", ageSel),
                        if(!is.null(countySel)) paste0(" Megye: ", countySel),
                        if(!is.null(dgSel)) paste0(" Diagnózis: ", dgSel),
                        if(!is.null(pciSel)) paste0(" Katéteres érmegnyitás: ", pciSel),
                        if(!is.null(comorbSel)) do.call(paste0, lapply(1:nrow(NameTable[type %in% c("comorb", "anamnestic")]), function(i)
                          if(!is.na(comorbSel[i])) paste0(" ", NameTable[type %in% c("comorb", "anamnestic")]$varname[i], ": ", comorbSel[i]))))
    )
  }
  
  owntab <- function(dat, nametext, valuetext, titletext) {
    dat <- merge(dat, NameTable, by = "variable")
    dat <- dat[, colnames(dat) %in% NameTable$variable, with = FALSE]
    setcolorder(dat, c("varname", if("DATEFORMATTED" %in% colnames(dat)) "DATEFORMATTED",
                       if("MEGYE" %in% colnames(dat)) "MEGYE",
                       setdiff(colnames(dat),
                               c("varname", "DATEFORMATTED", "MEGYE", "valueformatted")),
                       "valueformatted"))
    dat <- na.omit(dat)
    colnames(dat) <- NameTable$varname[match(colnames(dat), NameTable$variable)]
    colnames(dat)[colnames(dat) == "varname"] <- nametext
    colnames(dat)[colnames(dat) == "valueformatted"] <- valuetext
    
    DT::datatable(
      dat, rownames = FALSE, extensions = "Buttons",
      caption = paste0("50-nél kevesebb betegen alapuló adatok nem ",
                       "jelennek meg a táblázatban."),
      options = list(
        language = list(url = "https://cdn.datatables.net/plug-ins/2.0.8/i18n/hu.json"),
        dom = paste0("<'row'<'col-sm-4'l><'col-sm-8'f><'col-sm-12'tr>><'row'<'col-sm-4'i>",
                     "<'col-sm-8'p>>B"),
        buttons = list(list(extend = "copy", title = titletext), "spacer",
                       list(extend = "excel", title = titletext), "spacer",
                       list(extend = "csv", bom = TRUE, fieldSeparator = ";"), "spacer",
                       list(extend = "pdf", title = titletext))),
    )
  }
  
  ownplot <- function(dat, type, timeIncludeZero, timeCI, freq,
                      spacePlotType, spacePlotBarOrder, spacePlotBarHorizontal, spaceCI,
                      ytitle, uom, dig, displayallvalues, fixedcoloraxis,
                      fixedcoloraxisMin, fixedcoloraxisMax, titleText) {
    dat <- dat[!is.na(strat)]
    if(type == "space") dat <- dat[!is.na(MEGYE)]
    dat$value2 <- ifelse(is.na(dat$value), 0, dat$value)
    
    switch(type,
           "time" = {
             p <- highchart() |>
               hc_add_series(dat, type = "line",
                             hcaes(x = datetime_to_timestamp(DATE), y = value, group = strat),
                             id = unique(dat$strat)) |>
               hc_xAxis(type = "datetime",
                        dateTimeLabelFormats = list(year = "%Y", month = "%Y. %m.")) |>
               hc_yAxis(softMin = if(timeIncludeZero) 0 else NULL,
                        title = list(text = ytitle)) |>
               hc_tooltip(dateTimeLabelFormats = list(year = "%Y", month = "%Y. %B"),
                          pointFormat = paste0('<span style="color:{point.color}">\u25CF</span> {series.name}: <b>{point.value:.', dig, "f}", uom, "</b>",
                                               if(timeCI) paste0(" (95% CI: {point.lci:.", dig, "f}", uom, " - {point.uci:.", dig, "f}", uom, ")"), "<br>")) |>
               hc_title(text = titleText)
             if(displayallvalues) p <- p |>
                 hc_plotOptions(line = list(dataLabels = list(enabled = TRUE, allowOverlap = TRUE, format = paste0(
                   "{point.y:.", dig, "f}", if(timeCI) paste0(" (95% CI: {point.lci:.", dig, "f}", " - {point.uci:.", dig, "f}", ")") else "")), enableMouseTracking = FALSE))
             
             colorlist <- getOption("highcharter.color_palette")[
               rep(1:length(getOption("highcharter.color_palette")),
                   length(unique(dat$strat)))[1:length(unique(dat$strat))]]
             
             if(timeCI) {
               if(freq == "year") p <- p |>
                   hc_add_series(dat, type = "errorbar",
                                 hcaes(x = datetime_to_timestamp(DATE), y = value, group = strat,
                                       low = lci, high = uci),
                                 linkedTo = unique(dat$strat),
                                 color = colorlist,
                                 tooltip = list(enabled = FALSE))
               if(freq == "month") p <- p |>
                   hc_add_series(dat, type = "arearange",
                                 hcaes(x = datetime_to_timestamp(DATE), y = value, group = strat,
                                       low = lci, high = uci), fillOpacity = 0.1, lineWidth = 0,
                                 marker = list(enabled = FALSE,
                                               states = list(hover = list(enabled = FALSE))),
                                 linkedTo = unique(dat$strat), name = "95% CI",
                                 color = colorlist,
                                 tooltip = list(enabled = FALSE))
             }
           },
           "space" = {
             switch(spacePlotType,
                    "map" = {
                      p <- highchart(type = "map") |>
                        hc_add_series(mapData = mapdata, data = dat,
                                      joinBy = c("woe-name", "MEGYE"),
                                      name = unique(dat$strat),
                                      nullInteraction = TRUE) |>
                        hc_colorAxis(auxpar = NULL)
                      if(fixedcoloraxis) p <- p |>
                          hc_colorAxis(min = fixedcoloraxisMin, max = fixedcoloraxisMax)
                    },
                    "barchart" = {
                      p <- hchart(dat[order(if(spacePlotBarOrder) value2 else MEGYE)],
                                  type = if(spacePlotBarHorizontal) "bar" else "column",
                                  hcaes(x = MEGYE, y = value2), name = unique(dat$strat)) |>
                        hc_xAxis(title = list(text = "")) |>
                        hc_yAxis(title = list(text = ytitle))
                      if(spaceCI) {
                        p <- p |>
                          hc_add_series(dat, type = "errorbar",
                                        hcaes(x = MEGYE, y = value, low = lci, high = uci))
                      }
                    })
             dataLabelFormatString <- paste0("{point.value:.", dig, "f}", if(spaceCI)
               paste0(" (95% CI: {point.lci:.", dig, "f}", " - {point.uci:.", dig, "f}", ")") else "")
             p <- p |>
               hc_title(text = titleText) |>
               hc_plotOptions(
                 bar = list(minPointLength = 1, dataLabels = list(enabled = displayallvalues, format = dataLabelFormatString), enableMouseTracking = !displayallvalues),
                 column = list(minPointLength = 1, dataLabels = list(enabled = displayallvalues, format = dataLabelFormatString), enableMouseTracking = !displayallvalues),
                 map = list(dataLabels = list(enabled = displayallvalues, format = dataLabelFormatString, allowOverlap = TRUE), enableMouseTracking = !displayallvalues)
               ) |>
               # hc_tooltip(pointFormat = paste0("{point.name}: <b>{point.value:.", dig, "f}", uom, "</b>",
               #                                 if(spaceCI) paste0(" (95% CI: {point.lci:.", dig, "f}", uom, " - {point.uci:.", dig, "f}", uom, ")"), "<br>"),
               #            nullFormat = "{point.name}: Hiányzó adat, vagy 50-nél kevesebb betegen alapuló és ezért nem megjelenített adat")
               hc_tooltip(formatter = JS(paste0(
                 "function() {if (this.point.value === null) return(this.point.name + ': 50-nél kevesebb betegen alapuló és ezért nem megjelenített adat');",
                 "else return('<span style=\"color:' + this.point.color + '\">' + '\u25CF' + '</span>' + this.series.name + '<br>'",
                 "+ this.point.name + ': <b>' + this.point.value.toLocaleString(undefined, {minimumFractionDigits: ", dig, ", maximumFractionDigits: ", dig, "}) + '", uom, "</b>'",
                 if(spaceCI) paste0("+ ' (95% CI: ' + this.point.lci.toLocaleString(undefined, {minimumFractionDigits: ", dig, ", maximumFractionDigits: ", dig, "}) +'", uom, " - ' + this.point.uci.toLocaleString(undefined, {minimumFractionDigits: ", dig, ", maximumFractionDigits: ", dig, "}) +'", uom, ")'") else "", ")}")))
             # hc_tooltip(pointFormatter = JS(paste0("function() {console.log(this); if (this.missing) return('a'); else return(this.name + ': <b>' + this.value.toFixed(", dig, ") + '", uom, "</b>');}")))
             # hc_tooltip(formatter = JS(paste0("function() {console.log(this);}")))
           })
    
    if(sum(is.na(dat$value)) > 0) p <- p |>
      hc_caption(text = paste0("50-nél kevesebb betegen alapuló adatok nem ",
                               "jelennek meg az ábrán."),
                 align = "right")
    
    p <- p |>
      hc_subtitle(text = "Nemzeti Szívinfarktus Regiszter<br>https://nszr.gokvi.hu/",
                  align = "left", verticalAlign = "bottom") |>
      hc_credits(enabled = TRUE) |>
      hc_add_theme(hc_theme(chart = list(backgroundColor = "white"))) |>
      hc_exporting(enabled = TRUE, chartOptions = list(legend = TRUE),
                   sourceWidth = 1600/2, sourceHeight = 900/2,
                   pdfFont = list(
                     normal = "https://nszr.gokvi.hu/_upload/fonts/NotoSans-Regular.ttf",
                     bold = "https://nszr.gokvi.hu/_upload/fonts/NotoSans-Bold.ttf",
                     bolditalic = "https://nszr.gokvi.hu/_upload/fonts/NotoSans-BoldItalic.ttf",
                     italic = "https://nszr.gokvi.hu/_upload/fonts/NotoSans-Italic.ttf"
                   ),
                   buttons = list(contextButton = list(menuItems = list(
                     "viewFullscreen", "printChart", "separator", "downloadPNG", "downloadJPEG",
                     "downloadPDF", "downloadSVG"))))
    
    p
  }
  
  dataInputTimeIncidence <- reactive(dataInput(
    metric = input$incidenceMetric,
    primary = input$incidenceTimeIncidence,
    subject = "DATE",
    strat = if(length(input$incidenceTimeIncidence) > 1) "None" else (if(input$incidenceMetric != "adjrate") input$incidenceTimeStrat else input$incidenceTimeStratAdjrate),
    freq = input$incidenceTimeFreq,
    ci = input$incidenceTimeCI,
    yearSel = c(lubridate::year(min(RawData$YEAR)), lubridate::year(datadate) - 1),
    sexSel = if(input$incidenceTimeSelEnable && input$incidenceTimeSexSelEnable) input$incidenceTimeSexSel else NULL,
    dgSel = NULL,
    ageSel = if(input$incidenceTimeSelEnable && input$incidenceTimeAgeSelEnable && input$incidenceMetric != "adjrate") input$incidenceTimeAgeSel else NULL,
    comorbSel = NULL,
    pciSel = NULL,
    countySel = if(input$incidenceTimeSelEnable && input$incidenceTimeCountySelEnable) input$incidenceTimeCountySel else NULL,
    uom = "",
    titleText = "betegség incidenciája"))
  dataInputSpaceIncidence <- reactive(dataInput(
    metric = input$incidenceMetric,
    primary = input$incidenceSpaceIncidence,
    subject = "MEGYE",
    strat = "None",
    freq = "year",
    ci = input$incidenceSpaceCI,
    yearSel = input$incidenceSpaceYearSel,
    sexSel = if(input$incidenceSpaceSelEnable && input$incidenceSpaceSexSelEnable) input$incidenceSpaceSexSel else NULL,
    dgSel = NULL,
    ageSel = if(input$incidenceSpaceSelEnable && input$incidenceSpaceAgeSelEnable && input$incidenceMetric != "adjrate") input$incidenceSpaceAgeSel else NULL,
    comorbSel = NULL,
    pciSel = NULL,
    countySel = NULL,
    uom = "",
    titleText = "betegség incidenciája"))
  
  dataInputIncidence <- reactive(switch(input$incidenceSubject,
                                        "time" = dataInputTimeIncidence(),
                                        "space" = dataInputSpaceIncidence()))
  
  output$incidencePlot <- renderHighchart({
    if(input$incidenceSubject == "time" && length(input$incidenceTimeIncidence) == 0) return(NULL)
    
    di <- dataInputIncidence()
    dat <- merge(di$data, NameTable[type == "incidence"], by = "variable")
    strat <- if(input$incidenceMetric != "adjrate") input$incidenceTimeStrat else input$incidenceTimeStratAdjrate
    dat$strat <- if(input$incidenceSubject == "time" && length(input$incidenceTimeIncidence) == 1 &&
                    strat != "None") dat[[strat]] else dat$varname
    
    ownplot(dat, input$incidenceSubject, input$incidenceTimeIncludeZero, input$incidenceTimeCI && input$incidenceMetric != "absolute",
            input$incidenceTimeFreq, input$incidenceSpacePlotType,
            input$incidenceSpacePlotBarOrder, input$incidenceSpacePlotBarHorizontal,
            input$incidenceSpaceCI && input$incidenceMetric != "absolute",
            switch(input$incidenceMetric, "absolute" = "Esetszám", "cruderate" = "Nyers ráta [/100 ezer fő/év]", "adjrate" = "Standardizált ráta [/100 ezer fő/év]"),
            switch(input$incidenceMetric, "absolute" = "", "cruderate" = "/100 ezer fő/év", "adjrate" = "/100 ezer fő/év"),
            if(input$incidenceMetric == "absolute") 0 else 1,
            switch(input$incidenceSubject, "time" = input$incidenceTimeDisplayallvalues, "space" = input$incidenceSpaceDisplayallvalues),
            input$incidenceSpaceFixedcoloraxis, input$incidenceSpaceFixedcoloraxisMin, input$incidenceSpaceFixedcoloraxisMax, di$title)
  })
  
  output$incidenceTab <- DT::renderDT(
    server = FALSE,
    {
      if(input$incidenceSubject == "time" && length(input$incidenceTimeIncidence) == 0) return(NULL)
      owntab(dataInputIncidence()$data,
             "Betegség", switch(input$incidenceMetric, "absolute" = "Esetszám [fő]",
                                "cruderate" = "Nyers ráta [/100 ezer fő/év]",
                                "adjrate" = "Standardizált ráta [/100 ezer fő/év]"),
             "Betegség előfordulása, NSZR")
    })
  
  dataInputTimeComorb <- reactive(dataInput(
    metric = input$comorbMetric,
    primary = input$comorbTimeComorb,
    subject = "DATE",
    strat = if(length(input$comorbTimeComorb) > 1) "None" else input$comorbTimeStrat,
    freq = input$comorbTimeFreq,
    ci = input$comorbTimeCI,
    yearSel = NULL,
    sexSel = if(input$comorbTimeSelEnable && input$comorbTimeSexSelEnable) input$comorbTimeSexSel else NULL,
    dgSel = if(input$comorbTimeSelEnable && input$comorbTimeDgSelEnable) input$comorbTimeDgSel else NULL,
    ageSel = if(input$comorbTimeSelEnable && input$comorbTimeAgeSelEnable) input$comorbTimeAgeSel else NULL,
    comorbSel = NULL,
    pciSel = NULL,
    countySel = if(input$comorbTimeSelEnable && input$comorbTimeCountySelEnable) input$comorbTimeCountySel else NULL,
    uom = "%",
    titleText = "kórelőzmény, társbetegség aránya"))
  dataInputSpaceComorb <- reactive(dataInput(
    metric = input$comorbMetric,
    primary = input$comorbSpaceComorb,
    subject = "MEGYE",
    strat = "None",
    freq = NA,
    ci = input$comorbSpaceCI,
    yearSel = input$comorbSpaceYearSel,
    sexSel = if(input$comorbSpaceSelEnable && input$comorbSpaceSexSelEnable) input$comorbSpaceSexSel else NULL,
    dgSel = if(input$comorbSpaceSelEnable && input$comorbSpaceDgSelEnable) input$comorbSpaceDgSel else NULL,
    ageSel = if(input$comorbSpaceSelEnable && input$comorbSpaceAgeSelEnable) input$comorbSpaceAgeSel else NULL,
    comorbSel = NULL,
    pciSel = NULL,
    countySel = NULL,
    uom = "%",
    titleText = "kórelőzmény, társbetegség aránya"))
  
  dataInputComorb <- reactive(switch(input$comorbSubject,
                                     "time" = dataInputTimeComorb(),
                                     "space" = dataInputSpaceComorb()))
  
  output$comorbPlot <- renderHighchart({
    if(input$comorbSubject == "time" && length(input$comorbTimeComorb) == 0) return(NULL)
    
    di <- dataInputComorb()
    dat <- merge(di$data, NameTable[type %in%c("anamnestic", "comorb")], by = "variable")
    dat$strat <- if(input$comorbSubject == "time" && length(input$comorbTimeComorb) == 1 &&
                    input$comorbTimeStrat != "None") dat[[input$comorbTimeStrat]] else dat$varname
    
    ownplot(dat, input$comorbSubject, input$comorbTimeIncludeZero, input$comorbTimeCI,
            input$comorbTimeFreq, input$comorbSpacePlotType,
            input$comorbSpacePlotBarOrder, input$comorbSpacePlotBarHorizontal,
            input$comorbSpaceCI && input$comorbMetric != "absolute",
            "Arány [%]", "%", 1,
            switch(input$comorbSubject, "time" = input$comorbTimeDisplayallvalues, "space" = input$comorbSpaceDisplayallvalues),
            input$comorbSpaceFixedcoloraxis, input$comorbSpaceFixedcoloraxisMin, input$comorbSpaceFixedcoloraxisMax, di$title)
  })
  
  output$comorbTab <- DT::renderDT(
    server = FALSE,
    {
      owntab(dataInputComorb()$data, "Társbetegség", "Nyers arány [%]",
             "Adott kórelőzménnyel, társbetegséggel rendelkezők aránya, NSZR")
    })
  
  dataInputTimeTreatment <- reactive(dataInput(
    metric = input$treatmentMetric,
    primary = input$treatmentTimeTreatment,
    subject = "DATE",
    strat = if(length(input$treatmentTimeTreatment) > 1) "None" else switch(
      input$treatmentMetric,
      "crude" = input$treatmentTimeStrat,
      "agesexadj" = input$treatmentTimeStratAgesexadj,
      "agesexcomorbadj" = input$treatmentTimeStratAgesexcomorbadj),
    freq = input$treatmentTimeFreq,
    ci = input$treatmentTimeCI,
    yearSel = NULL,
    sexSel = if(input$treatmentTimeSelEnable && input$treatmentTimeSexSelEnable && input$treatmentMetric == "crude") input$treatmentTimeSexSel else NULL,
    dgSel = if(input$treatmentTimeSelEnable && input$treatmentTimeDgSelEnable) input$treatmentTimeDgSel else NULL,
    ageSel = if(input$treatmentTimeSelEnable && input$treatmentTimeAgeSelEnable && input$treatmentMetric == "crude") input$treatmentTimeAgeSel else NULL,
    comorbSel = sapply(NameTable[type %in%c("anamnestic", "comorb")]$variable, function(com)
      if(input$treatmentTimeSelEnable && input$treatmentMetric != "agesexcomorbadj" && input[[paste0("treatmentTime", com, "SelEnable")]])
        input[[paste0("treatmentTime", com, "Sel")]] else NA),
    pciSel = NULL,
    countySel = if(input$treatmentTimeSelEnable && input$treatmentTimeCountySelEnable) input$treatmentTimeCountySel else NULL,
    uom = "%",
    titleText = "ellátásban részesülők aránya"))
  dataInputSpaceTreatment <- reactive(dataInput(
    metric = input$treatmentMetric,
    primary = input$treatmentSpaceTreatment,
    subject = "MEGYE",
    strat = "None",
    freq = NA,
    ci = input$treatmentSpaceCI,
    yearSel = input$treatmentSpaceYearSel,
    sexSel = if(input$treatmentSpaceSelEnable && input$treatmentMetric == "crude" && input$treatmentSpaceSexSelEnable) input$treatmentSpaceSexSel else NULL,
    dgSel = if(input$treatmentSpaceSelEnable && input$treatmentSpaceDgSelEnable) input$treatmentSpaceDgSel else NULL,
    ageSel = if(input$treatmentSpaceSelEnable && input$treatmentMetric == "crude" && input$treatmentSpaceAgeSelEnable) input$treatmentSpaceAgeSel else NULL,
    comorbSel = sapply(NameTable[type %in%c("anamnestic", "comorb")]$variable, function(com)
      if(input$treatmentSpaceSelEnable && input$treatmentMetric != "agesexcomorbadj" && input[[paste0("treatmentSpace", com, "SelEnable")]]) input[[paste0("treatmentSpace", com, "Sel")]] else NA),
    pciSel = NULL,
    countySel = NULL,
    uom = "%",
    titleText = "ellátásban részesülők aránya"))
  
  dataInputTreatment <- reactive(switch(input$treatmentSubject,
                                        "time" = dataInputTimeTreatment(),
                                        "space" = dataInputSpaceTreatment()))
  
  output$treatmentPlot <- renderHighchart({
    if(input$treatmentSubject == "time" && length(input$treatmentTimeTreatment) == 0) return(NULL)
    
    di <- dataInputTreatment()
    dat <- merge(di$data, NameTable[type == "treatment"], by = "variable")
    strat <- switch(input$treatmentMetric,
                    "crude" = input$treatmentTimeStrat,
                    "agesexadj" = input$treatmentTimeStratAgesexadj,
                    "agesexcomorbadj" = input$treatmentTimeStratAgesexcomorbadj)
    dat$strat <- if(input$treatmentSubject == "time" && length(input$treatmentTimeTreatment) == 1 &&
                    strat != "None") dat[[strat]] else dat$varname
    
    ownplot(
      dat, input$treatmentSubject, input$treatmentTimeIncludeZero, input$treatmentTimeCI,
      input$treatmentTimeFreq, input$treatmentSpacePlotType,
      input$treatmentSpacePlotBarOrder, input$treatmentSpacePlotBarHorizontal,
      input$treatmentSpaceCI && input$treatmentMetric != "absolute",
      "Arány [%]", "%", 1,
      switch(input$treatmentSubject, "time" = input$treatmentTimeDisplayallvalues, "space" = input$treatmentSpaceDisplayallvalues),
      input$treatmentSpaceFixedcoloraxis, input$treatmentSpaceFixedcoloraxisMin, input$treatmentSpaceFixedcoloraxisMax, di$title)
  })
  
  output$treatmentTab <- DT::renderDT(
    server = FALSE,
    {
      owntab(dataInputTreatment()$data, "Ellátás", paste0("Arány ", switch(
        input$treatmentMetric, "crude" = "(nyers)",
        "agesexadj" = "(életkorra és nemre korrigált)",
        "agesexcomorbadj" = "(életkorra, nemre és társbetegségekre korrigált)", " [%]")),
        "Adott ellátásban részesülők aránya, NSZR")
    })
  
  dataInputTimeSurvival <- reactive(dataInput(
    metric = input$survivalMetric,
    primary = input$survivalTimeSurvival,
    subject = "DATE",
    strat = if(length(input$survivalTimeSurvival) > 1) "None" else switch(input$survivalMetric,
                                                                          "crude" = input$survivalTimeStrat,
                                                                          "agesexadj" = input$survivalTimeStratAgesexadj,
                                                                          "agesexcomorbadj" = input$survivalTimeStratAgesexcomorbadj),
    freq = input$survivalTimeFreq,
    ci = input$survivalTimeCI,
    yearSel = NULL,
    sexSel = if(input$survivalTimeSelEnable && input$survivalTimeSexSelEnable && input$survivalMetric == "crude") input$survivalTimeSexSel else NULL,
    dgSel = if(input$survivalTimeSelEnable && input$survivalTimeDgSelEnable) input$survivalTimeDgSel else NULL,
    ageSel = if(input$survivalTimeSelEnable && input$survivalTimeAgeSelEnable && input$survivalMetric == "crude") input$survivalTimeAgeSel else NULL,
    comorbSel = sapply(NameTable[type %in%c("anamnestic", "comorb")]$variable, function(com)
      if(input$survivalTimeSelEnable && input$survivalMetric != "agesexcomorbadj" && input[[paste0("survivalTime", com, "SelEnable")]])
        input[[paste0("survivalTime", com, "Sel")]] else NA),
    pciSel = if(input$survivalTimeSelEnable && input$survivalTimePciSelEnable) input$survivalTimePciSel else NULL,
    countySel = if(input$survivalTimeSelEnable && input$survivalTimeCountySelEnable) input$survivalTimeCountySel else NULL,
    uom = "%",
    titleText = "kimenet aránya"))
  dataInputSpaceSurvival <- reactive(dataInput(
    metric = input$survivalMetric,
    primary = input$survivalSpaceSurvival,
    subject = "MEGYE",
    strat = "None",
    freq = NA,
    ci = input$survivalSpaceCI,
    yearSel = input$survivalSpaceYearSel,
    sexSel = if(input$survivalSpaceSelEnable && input$survivalMetric == "crude" && input$survivalSpaceSexSelEnable) input$survivalSpaceSexSel else NULL,
    dgSel = if(input$survivalSpaceSelEnable && input$survivalSpaceDgSelEnable) input$survivalSpaceDgSel else NULL,
    ageSel = if(input$survivalSpaceSelEnable && input$survivalMetric == "crude" && input$survivalSpaceAgeSelEnable) input$survivalSpaceAgeSel else NULL,
    comorbSel = sapply(NameTable[type %in%c("anamnestic", "comorb")]$variable, function(com)
      if(input$survivalSpaceSelEnable && input$survivalMetric != "agesexcomorbadj" && input[[paste0("survivalSpace", com, "SelEnable")]]) input[[paste0("survivalSpace", com, "Sel")]] else NA),
    pciSel = if(input$survivalSpaceSelEnable && input$survivalSpacePciSelEnable) input$survivalSpacePciSel else NULL,
    countySel = NULL,
    uom = "%",
    titleText = "kimenet aránya"))
  
  dataInputSurvival <- reactive(switch(input$survivalSubject,
                                       "time" = dataInputTimeSurvival(),
                                       "space" = dataInputSpaceSurvival()))
  
  output$survivalPlot <- renderHighchart({
    if(input$survivalSubject == "time" && length(input$survivalTimeSurvival) == 0) return(NULL)
    
    di <- dataInputSurvival()
    dat <- merge(di$data, NameTable[type == "survival"], by = "variable")
    strat <- switch(input$survivalMetric,
                    "crude" = input$survivalTimeStrat,
                    "agesexadj" = input$survivalTimeStratAgesexadj,
                    "agesexcomorbadj" = input$survivalTimeStratAgesexcomorbadj)
    dat$strat <- if(input$survivalSubject == "time" && length(input$survivalTimeSurvival) == 1 &&
                    strat != "None") dat[[strat]] else dat$varname
    
    ownplot(
      dat, input$survivalSubject, input$survivalTimeIncludeZero, input$survivalTimeCI,
      input$survivalTimeFreq, input$survivalSpacePlotType,
      input$survivalSpacePlotBarOrder, input$survivalSpacePlotBarHorizontal,
      input$survivalSpaceCI && input$survivalMetric != "absolute",
      "Arány [%]", "%", 1,
      switch(input$survivalSubject, "time" = input$survivalTimeDisplayallvalues, "space" = input$survivalSpaceDisplayallvalues),
      input$survivalSpaceFixedcoloraxis, input$survivalSpaceFixedcoloraxisMin, input$survivalSpaceFixedcoloraxisMax, di$title)
  })
  
  output$survivalTab <- DT::renderDT(
    server = FALSE,
    {
      owntab(dataInputSurvival()$data, "Túlélési mutató", paste0("Arány ", switch(
        input$survivalMetric, "crude" = "(nyers)",
        "agesexadj" = "(életkorra és nemre korrigált)",
        "agesexcomorbadj" = "(életkorra, nemre és társbetegségekre korrigált)"), " [%]"),
        "Halálozási arány, NSZR")
    })
}

shinyApp(ui = ui, server = server)
