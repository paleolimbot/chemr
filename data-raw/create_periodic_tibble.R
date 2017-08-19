
# create the periodic tibble
library(dplyr)
library(stringr)

# we need a list of elements, which can be found at
# https://en.wikipedia.org/wiki/List_of_chemical_elements
elements <- htmltab::htmltab("https://en.wikipedia.org/wiki/List_of_chemical_elements", which = 2)

# clean names
names(elements) <- str_replace(names(elements), "^List of chemical elements >> ", "")
names(elements) <- str_replace(tolower(names(elements)), "\\(.*\\)", "")

# select columns, fix types
pt <- elements %>%
  filter(!str_detect(z, "^Notes")) %>%
  select(z, symbol, name = element, group, period, mass = `atomic weight`) %>%
  mutate(z = as.integer(z), group = as.integer(group), period = as.integer(period)) %>%
  # fix references in atomic weight
  mutate(mass = str_replace_all(mass, "[\\[\\(](.*?)[\\]\\)]", "\\1") %>% as.double()) %>%
  mutate(symbol = factor(symbol, levels = symbol)) %>%
  tibble::as.tibble()

# add to package
devtools::use_data(pt, overwrite = TRUE)


