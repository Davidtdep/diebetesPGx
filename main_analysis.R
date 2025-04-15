################################################################################
# This script reproduces the original workflow.  Each numbered
# section is self‑contained so collaborators can jump to what they need.
################################################################################

################################################################################
# 1.  Load Required Libraries ---------------------------------------------------
################################################################################
# Uncomment to install any that are missing
# pkgs <- c("dplyr", "ggplot2", "ggExtra", "ggsci", "reshape2", "circlize", 
#           "openxlsx", "rentrez", "httr", "jsonlite", "xml2", "rcrossref", 
#           "fuzzyjoin")
# install.packages(setdiff(pkgs, rownames(installed.packages())))

library(dplyr)
library(ggplot2)
library(ggExtra)
library(ggsci)
library(reshape2)
library(circlize)
library(openxlsx)
library(rentrez)
library(httr)
library(jsonlite)
library(xml2)
library(rcrossref)
library(fuzzyjoin)
library(stringr)

################################################################################
# 2.  Import Primary Dataset ----------------------------------------------------
################################################################################
# Pharmacogenomic associations for diabetes (change the path as needed)
raw_path <- "~/data/data.tsv"
stopifnot(file.exists(raw_path))

data <- read.delim(raw_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

################################################################################
# 3.  Data Cleaning & Feature Engineering --------------------------------------
################################################################################
# 3.1 Extract simple alleles / homozygous genotypes -----------------------------
extract_allele <- function(txt) {
  if (grepl("^Allele [A-Z] ", txt)) {
    sub("^Allele ([A-Z]) .*", "\\1", txt)
  } else if (grepl("^Genotype [A-Z]{2} ", txt)) {
    geno <- sub("^Genotype ([A-Z]{2}).*", "\\1", txt)
    if (substr(geno, 1, 1) == substr(geno, 2, 2)) geno else NA
  } else {
    NA
  }
}

# 3.2 Extract effect direction (increased / decreased) --------------------------
extract_effect <- function(txt) {
  m <- regmatches(txt, regexpr("associated with (\\w+)", txt))
  ifelse(length(m) > 0, sub("associated with ", "", m), NA)
}

# 3.3 Apply transformations -----------------------------------------------------

data <- data %>%
  mutate(Allele = sapply(Association, extract_allele),
         Effect = sapply(Association, extract_effect))

# 3.4 Keep significant rows with single‑letter alleles -------------------------

data_significant <- data %>%
  filter(Significance == "yes", !is.na(Allele), nchar(Allele) == 1)

# (Optional) Persist the cleaned subset
# write.csv(data_significant, "~/Desktop/endocongress2025/data/diabetes_significant.csv", row.names = FALSE)

################################################################################
# 4.  Correlation Analysis (Numeric Variables) ----------------------------------
################################################################################
num_vars <- c("ATQCES", "ATQPGC", "CLM", "PLQ", "CHG")

# Import data_significant with allele frequencies (change the path as needed)
data_significant <- read.csv("~/data/data_significant.csv", header = TRUE)

# Compute pairwise Pearson correlations (excluding missing pairs)
cor_matrix <- matrix(NA, nrow = length(num_vars), ncol = length(num_vars))
rownames(cor_matrix) <- num_vars
colnames(cor_matrix) <- num_vars

for (i in 1:length(num_vars)) {
  for (j in 1:length(num_vars)) {
    if (i != j) {  # Evitar la diagonal principal (que es siempre 1)
      cor_result <- cor.test(data_significant[[num_vars[i]]], 
                             data_significant[[num_vars[j]]], 
                             use = "complete.obs", 
                             method = "pearson")
      cor_matrix[i, j] <- cor_result$estimate  # Guardar el coeficiente de correlación
    } else {
      cor_matrix[i, j] <- 1  # La diagonal principal siempre es 1
    }
  }
}

# Heat‑map (lower triangle)
cor_matrix <- as.data.frame(cor_matrix)
cor_long <- cor_matrix %>%
  tibble::rownames_to_column("Var1") %>%
  reshape2::melt(id.vars = "Var1", variable.name = "Var2", value.name = "r") %>%
  filter(as.numeric(factor(Var1, levels = num_vars)) >=
           as.numeric(factor(Var2, levels = num_vars)))

cor_long$Var1 <- factor(cor_long$Var1, levels = rownames(cor_matrix))
cor_long$Var2 <- factor(cor_long$Var2, levels = rownames(cor_matrix))

heatmap_plot <- ggplot(cor_long, aes(Var1, Var2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
print(heatmap_plot)


################################################################################
# 5.  Scatter Plot by Ancestry --------------------------------------------------
################################################################################
plot_df <- data_significant %>%
  mutate(African_Freq  = rowMeans(select(., PLQ, CHG), na.rm = TRUE),
         European_Freq = rowMeans(select(., ATQCES, CLM, ATQPGC), na.rm = TRUE),
         Phenotype     = paste(Phenotype.Categories, Effect)) %>%
  filter(!is.na(African_Freq), !is.na(European_Freq))

scatter_base <- ggplot(plot_df, aes(African_Freq, European_Freq, colour = Drugs)) +
  geom_point(size = 3) +
  scale_colour_npg() +
  theme_classic(base_size = 13) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "African Ancestry", y = "European Ancestry")

ggMarginal(scatter_base, type = "histogram", groupColour = TRUE, groupFill = TRUE)

################################################################################
# 6.  Chord Diagram: Cases & Controls by Ancestry ------------------------------
################################################################################
# 6.1 Recode ancestry labels

data <- data %>%
  mutate(Ancestry = case_when(
    Biogeographical.Groups %in% c("African American/Afro-Caribbean")                     ~ "African",
    Biogeographical.Groups %in% c("American", "Latino")                                 ~ "American",
    Biogeographical.Groups %in% c("Central/South Asian", "East Asian", "Near Eastern") ~ "Asian",
    Biogeographical.Groups == "European"                                                 ~ "European",
    Biogeographical.Groups == "Multiple groups"                                          ~ "Mixed",
    Biogeographical.Groups == "Unknown"                                                 ~ "Unknown",
    TRUE                                                                                 ~ "Other"))

# 6.2 Summaries
ancestry_totals <- data %>%
  group_by(Ancestry) %>%
  summarise(Cases    = sum(X..of.Cases,     na.rm = TRUE),
            Controls = sum(X..of.Controls,  na.rm = TRUE), .groups = "drop")

links <- tidyr::pivot_longer(ancestry_totals, Cases:Controls, names_to = "to", values_to = "value") %>%
  rename(from = Ancestry) %>%
  filter(value > 0)

order_from <- rev(c("African", "Mixed", "American", "Asian", "Unknown", "European"))
order_to   <- c("Cases", "Controls")

links$from <- factor(links$from, levels = order_from)
links$to   <- factor(links$to,   levels = order_to)

grid.col <- c("European" = "#5fa8d3", "Unknown" = "#cae9ff", "Asian" = "#1b4965",
              "American"  = "#62b6cb", "Mixed"  = "#bee9e8", "African" = "#e63946",
              "Cases"     = "black",   "Controls" = "darkgrey")
transparency <- c("European" = 1, "Unknown" = 1, "Asian" = 1,
                  "American" = 1, "Mixed" = 1, "African" = 0)

circos.clear()
chordDiagram(links, grid.col = grid.col,
             transparency = transparency[as.character(links$from)],
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.01))

################################################################################
# 7.  Descriptive Bar Plots -----------------------------------------------------
################################################################################
bar_palette <- c("Efficacy decreased" = "#fdca40", "Efficacy increased" = "#df2935",
                 "Toxicity decreased" = "#3772ff", "Toxicity increased" = "#540d6e")

ggplot(plot_df, aes(reorder(Drugs, table(Drugs)[Drugs]), fill = Phenotype)) +
  geom_bar() +
  scale_fill_manual(values = bar_palette) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggplot(plot_df, aes(reorder(Genes, table(Genes)[Genes]), fill = Phenotype)) +
  geom_bar() +
  scale_fill_manual(values = bar_palette) +
  coord_flip() +
  labs(x = "Genes") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

################################################################################
# 8.  Bibliometric Metadata -----------------------------------------------------
################################################################################
# 8.1 Separate PubMed & PMC identifiers ----------------------------------------
ids    <- unique(data$Literature)
pmids  <- sub("PMID:",  "", ids[grepl("^PMID:",  ids)])
pmcids <- sub("PMCID:", "", ids[grepl("^PMCID:", ids)])

# 8.2 Helper functions ---------------------------------------------------------
get_info_pmids <- function(pmid_vec) {
  purrr::map_dfr(pmid_vec, function(pmid) {
    # Basic summary
    s <- tryCatch(entrez_summary(db = "pubmed", id = pmid), error = function(e) NULL)
    if (is.null(s)) return(data.frame())
    doi <- if (!is.null(s$elocationid)) s$elocationid else NA
    
    # Author & country
    xml_raw <- tryCatch(entrez_fetch(db = "pubmed", id = pmid, rettype = "xml"), error = function(e) NULL)
    if (is.null(xml_raw)) return(data.frame())
    xml <- read_xml(xml_raw)
    author_node <- xml_find_first(xml, ".//AuthorList/Author[1]")
    first_author <- paste(xml_text(xml_find_first(author_node, "./ForeName")),
                          xml_text(xml_find_first(author_node, "./LastName")))
    aff <- xml_text(xml_find_first(author_node, ".//Affiliation"))
    country <- if (!is.na(aff)) trimws(tail(strsplit(strsplit(aff, ";")[[1]][1], ",")[[1]], 1)) else NA
    
    tibble::tibble(id = paste0("PMID:", pmid),
                   title = s$title,
                   journal = s$fulljournalname,
                   pubdate = s$pubdate,
                   first_author = first_author,
                   first_author_country = country,
                   doi = doi)
  })
}

get_info_pmcids <- function(pmc_vec) {
  pmid_lookup <- purrr::map_chr(pmc_vec, function(pmcid) {
    links <- tryCatch(entrez_link(dbfrom = "pmc", id = sub("PMC", "", pmcid), db = "pubmed"),
                      error = function(e) NULL)
    if (!is.null(links) && length(links$links$pmc_pubmed) > 0) links$links$pmc_pubmed[1] else NA
  })
  info <- get_info_pmids(na.omit(pmid_lookup))
  info$id <- paste0("PMCID:", pmc_vec[match(info$id, paste0("PMID:", pmid_lookup))])
  info
}

# 8.3 Retrieve metadata --------------------------------------------------------
articles <- bind_rows(get_info_pmids(pmids), get_info_pmcids(pmcids))

# 8.4 Open Access & citations ---------------------------------------------------
articles$doi <- sub("^doi:\\s*", "", articles$doi, ignore.case = TRUE)

get_open_access <- function(doi) {
  if (is.na(doi) || doi == "" || !grepl("^10\\.", doi)) return(NA)
  
  url <- paste0("https://api.unpaywall.org/v2/", doi, "?email=david.tucorreo@ejemplo.edu.co")
  
  res <- tryCatch(GET(url), error = function(e) return(NULL))
  
  if (!is.null(res) && status_code(res) == 200) {
    data <- fromJSON(content(res, as = "text", encoding = "UTF-8"))
    return(data$is_oa)
  } else {
    return(NA)
  }
}

get_citations <- function(doi) {
  if (is.na(doi) || doi == "") return(NA)
  tryCatch(cr_citation_count(doi)$count, error = function(e) NA)
}

articles$doi[articles$id == "PMID:21857094"] <- "10.1016/s1734-1140(11)70595-3"
articles$doi[articles$id == "PMCID:PMC2352037"] <- "10.1136/bmj.313.7057.591"
articles$doi[articles$id == "PMID:8562390"] <- "10.1046/j.1365-2141.1996.275810.x"
articles$doi[articles$id == "PMID:27271184"] <- "10.2337/dc15-2464"
articles$doi[articles$id == "PMID:9551410"] <- "10.1111/j.1523-1755.1998.00847.x"
articles$doi[articles$id == "PMID:11007831"] <- "10.1093/ndt/15.10.1617"
articles$doi[articles$id == "PMID:15608561"] <- "10.1097/00008571-200412000-00005"

articles <- articles %>%
  mutate(open_access = sapply(doi, get_open_access),
         citations   = sapply(doi, get_citations))

# 8.5 Journal metrics (SJR) -----------------------------------------------------
sjr_path <- "~/Desktop/endocongress2025/data/scimagojr 2024.csv"
stopifnot(file.exists(sjr_path))

sjr <- read.csv(sjr_path, sep = ";", stringsAsFactors = FALSE)

setdiff(articles$journal, sjr$Title)

articles$journal <- gsub("Diabetes care", "Diabetes Care",
                          articles$journal)
articles$journal <- gsub("Transplantation proceedings", "Transplantation Proceedings",
                          articles$journal)
articles$journal <- gsub("Pharmacogenetics and genomics", "Pharmacogenetics and Genomics",
                          articles$journal)
articles$journal <- gsub("American journal of hypertension", "American Journal of Hypertension",
                          articles$journal)
articles$journal <- gsub("Pediatric nephrology \\(Berlin, Germany\\)", "Pediatric Nephrology",
                         articles$journal)
articles$journal <- gsub("Pharmacogenetics", "Pharmacogenetics and Genomics",
                         articles$journal)
articles$journal <- gsub("Pharmacogenetics and Genomics and Genomics", "Pharmacogenetics and Genomics",
                         articles$journal)
articles$journal <- gsub("Clinical transplantation", "Clinical Transplantation",
                         articles$journal)
articles$journal <- gsub("Nephrology, dialysis, transplantation : official publication of the European Dialysis and Transplant Association - European Renal Association", "Nephrology Dialysis Transplantation",
                         articles$journal)
articles$journal <- gsub("Pharmacological reports : PR", "Pharmacological Reports",
                         articles$journal)
articles$journal <- gsub("British journal of haematology", "British Journal of Haematology",
                         articles$journal)
articles$journal <- gsub("The world journal of biological psychiatry : the official journal of the World Federation of Societies of Biological Psychiatry", "World Journal of Biological Psychiatry",
                         articles$journal)
articles$journal <- gsub("Kidney international", "Kidney International",
                         articles$journal)
articles$journal <- gsub("PloS one", "PLoS ONE",
                         articles$journal)
articles$journal <- gsub("Clinical pharmacology and therapeutics", "Clinical Pharmacology and Therapeutics",
                         articles$journal)
articles$journal <- gsub("Indian heart journal", "Indian Heart Journal",
                         articles$journal)
articles$journal <- gsub("Journal of the American Society of Nephrology : JASN", "Journal of the American Society of Nephrology",
                         articles$journal)
articles$journal <- gsub("Frontiers in pharmacology", "Frontiers in Pharmacology",
                         articles$journal)
articles$journal <- gsub("BMJ \\(Clinical research ed.\\)", "BMJ",
                         articles$journal)
articles$journal <- gsub("AMIA Joint Summits on Translational Science proceedings. AMIA Joint Summits on Translational Science", "AMIA Joint Summits on Translational Science proceedings",
                         articles$journal)
articles$journal <- gsub("Nature genetics", "Nature Genetics",
                         articles$journal)


articles <- merge(articles,
                  sjr %>% select(Title, H.index, SJR.Best.Quartile),
                  by.x = "journal", by.y = "Title", all.x = TRUE)




# 8.6 Data depuration -----------------------------------------------------------

# Extract only the year (4 consecutive digits) from the 'pubdate' column
articles$pubdate <- gsub(".*(\\d{4}).*", "\\1", articles$pubdate)



# Step 1: Extract the first "word-like" country portion before any punctuation or email
articles$first_author_country <- str_trim(str_extract(articles$first_author_country, "^[^,\\.@]+"))

# Step 2: Manual standardization (you can expand this list as needed)
standardize_country <- function(country) {
  country <- str_trim(country)  # Remove extra spaces
  
  country <- case_when(
    country %in% c("USA", "US", "U.S.", "U.S.A") ~ "United States",
    country %in% c("UK", "U.K.") ~ "United Kingdom",
    country %in% c("The Netherlands", "the Netherlands", "Netherlands") ~ "Netherlands",
    country %in% c("People's Republic of China", "China") ~ "China",
    country %in% c("Belgique", "Belgium") ~ "Belgium",
    country %in% c("Korea", "South Korea") ~ "South Korea",
    country %in% c("Okinawa Japan", "Japan") ~ "Japan",
    TRUE ~ country  # Leave it as-is if no match
  )
  
  return(country)
}

# Apply the function to standardize the country names
articles$first_author_country <- sapply(articles$first_author_country, standardize_country)

# Change some country names
articles$first_author_country <- gsub("MN", "United States",
                         articles$first_author_country)

articles$first_author_country <- gsub("ivan", "Slovakia",
                         articles$first_author_country)

library(dplyr)

# Add WHO region column
articles <- articles %>%
  mutate(who_region = case_when(
    first_author_country %in% c("France", "Germany", "Italy", "Netherlands", "Spain", 
                                "United Kingdom", "Sweden", "Denmark", "Belgium", 
                                "Poland", "Slovenia", "Slovakia", "Scotland", "Switzerland") ~ "Europe",
    
    first_author_country %in% c("United States", "Canada") ~ "Americas",
    
    first_author_country %in% c("Mexico") ~ "Americas",
    
    first_author_country %in% c("China", "Japan", "South Korea") ~ "Western Pacific",
    
    first_author_country %in% c("India") ~ "South-East Asia",
    
    first_author_country %in% c("Australia") ~ "Western Pacific",
    
    first_author_country %in% c("Tunisia") ~ "Eastern Mediterranean",
    
    # Puedes agregar más países si aparecen en tu dataset
    TRUE ~ "Other/Unknown"
  ))



################################################################################
# 9.  Descriptive bibliometric analysis
################################################################################

# 9.1 Number of articles by year ------------------------------------------------
# Ensure 'pubdate' is numeric
articles$pubdate <- as.numeric(articles$pubdate)

# Create a 5-year interval column
articles <- articles %>%
  mutate(year_group = cut(pubdate,
                          breaks = seq(floor(min(pubdate, na.rm = TRUE) / 5) * 5,
                                       ceiling(max(pubdate, na.rm = TRUE) / 5) * 5,
                                       by = 5),
                          include.lowest = TRUE,
                          right = FALSE,
                          labels = paste(seq(floor(min(pubdate, na.rm = TRUE) / 5) * 5,
                                             ceiling(max(pubdate, na.rm = TRUE) / 5) * 5 - 5,
                                             by = 5),
                                         seq(floor(min(pubdate, na.rm = TRUE) / 5) * 5 + 4,
                                             ceiling(max(pubdate, na.rm = TRUE) / 5) * 5 - 1,
                                             by = 5),
                                         sep = "-")))


# Set the stacking order: most frequent at the bottom
articles$who_region <- factor(articles$who_region,
                              levels = rev(c("Europe", "Western Pacific", "Americas", 
                                         "South-East Asia", "Eastern Mediterranean")))


# Define custom colors for each WHO region
region_colors <- c(
  "Europe" =                "#ff595e",               # Blue
  "Western Pacific" =       "#ffca3a",      # Purple
  "Americas" =              "#8ac926",             # Orange
  "South-East Asia" =       "#1982c4",      # Green
  "Eastern Mediterranean" = "#6a4c93" # Red
)

# Plot
articles %>%
  group_by(year_group, who_region) %>%
  summarise(n_articles = n(), .groups = "drop") %>%
  ggplot(aes(x = year_group, y = n_articles, fill = who_region)) +
  
  # Add horizontal dashed lines at y = 5, 10, 15
  geom_hline(yintercept = c(5, 10, 15), linetype = "dashed", color = "gray60", size = 0.2) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  # Add solid baseline at y = 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) +
  
  # Custom fill colors
  scale_fill_manual(values = region_colors) +
  
  # Labels
  labs(x = "5-Year Interval",
       y = NULL,
       fill = "WHO Region") +
  
  # Theme adjustments
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove all gridlines
    axis.text.x = element_text(angle = 45, hjust = 1)
    ,legend.position = "none",
  )




# 9.2 Number of articles by income group ----------------------------------------

# Assign income groups manually
articles$income_group <- dplyr::case_when(
  articles$first_author_country %in% c("United States", "United Kingdom", "Netherlands", 
                                       "Sweden", "Switzerland", "Denmark", "Japan", 
                                       "Australia", "Belgium", "Italy", "South Korea", 
                                       "France", "Spain", "Poland", "Slovenia", 
                                       "Slovakia", "Scotland") ~ "High income",
  
  articles$first_author_country == "China" ~ "Upper-middle income",
  articles$first_author_country == "Mexico" ~ "Upper-middle income",
  articles$first_author_country == "India" ~ "Lower-middle income",
  articles$first_author_country == "Tunisia" ~ "Lower-middle income",
  
  TRUE ~ "Unknown"
)

articles$income_group <- factor(articles$income_group,
                              levels = rev(c("High income", "Upper-middle income", "Lower-middle income", 
                                             "Unknown")))

# Plot
articles %>%
  group_by(year_group, income_group) %>%
  summarise(n_articles = n(), .groups = "drop") %>%
  ggplot(aes(x = year_group, y = n_articles, fill = income_group)) +
  geom_hline(yintercept = c(5, 10, 15), linetype = "dashed", color = "gray60", size = 0.3) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) +
  scale_fill_manual(values = c(
    "High income" = "#1d3557",
    "Upper-middle income" = "#457b9d",
    "Lower-middle income" = "#a8dadc",
    "Unknown" = "#cccccc"
  )) +
  theme_minimal() +
  labs(x = "5-Year Interval",
       y = NULL,
       fill = "Income Group") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
    ,legend.position = "none",
  )



# 9.3 Number of articles by journal quartile ------------------------------------


# Plot
articles %>%
  group_by(year_group, SJR.Best.Quartile) %>%
  summarise(n_articles = n(), .groups = "drop") %>%
  ggplot(aes(x = year_group, y = n_articles, fill = SJR.Best.Quartile)) +
  geom_hline(yintercept = c(5, 10, 15), linetype = "dashed", color = "gray60", size = 0.3) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) +
  scale_fill_manual(values = c(
    "Q1" = "#1d3557",
    "Q2" = "#457b9d",
    "Q3" = "#a8dadc",
    "NA" = "#cccccc"
  )) +
  theme_minimal() +
  labs(x = "5-Year Interval",
       y = NULL,
       fill = "Journal Quartil") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
    ,legend.position = "none",
  )




# 9.4 Number of articles by open access -----------------------------------------

# Plot
articles %>%
  group_by(year_group, open_access) %>%
  summarise(n_articles = n(), .groups = "drop") %>%
  ggplot(aes(x = year_group, y = n_articles, fill = open_access)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.3) +
  scale_fill_manual(values = c(
    "TRUE" = "#e63946",
    "FALSE" = "#1d3557",
    "NA" = "#cccccc"
  )) +
  theme_minimal() +
  labs(x = "5-Year Interval",
       y = NULL,
       fill = "Open Access Yes No") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
    #,legend.position = "none",
  )


################################################################################
# End of Script ----------------------------------------------------------------
################################################################################
