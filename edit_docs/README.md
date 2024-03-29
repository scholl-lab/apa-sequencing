# APA-Sequencing Documentation

The repository subfolder for the APA-Sequencing documentation.

To build the documentation execute following commands:

```
## load libraries
library(bookdown)
library(config)

project_topic <- "development"
project_name <- "apa-sequencing"
script_path <- "/edit_docs/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

# render the book
bookdown::render_book("index.Rmd", "all")
```

The CSL file is [American Psychological Association 7th edition, from Zotero](https://www.zotero.org/styles/apa).

## required packages
Rendering HTML widgets for PDF requires webshot and phantomjs [FROM:](https://bookdown.org/yihui/bookdown/html-widgets.html).
```
install.packages("webshot")
webshot::install_phantomjs()
```
