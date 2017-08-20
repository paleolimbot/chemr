
# create the periodic tibble
library(dplyr)
library(stringr)
library(tidyr)

# we need a list of elements, which can be found at
# https://en.wikipedia.org/wiki/List_of_chemical_elements
elements <- htmltab::htmltab("https://en.wikipedia.org/wiki/List_of_chemical_elements", which = 2)

# clean names
names(elements) <- str_replace(names(elements), "^List of chemical elements >> ", "")
names(elements) <- str_replace(tolower(names(elements)), "\\(.*\\)", "")

# also, a list of oxidation states
oxidation <- htmltab::htmltab("https://en.wikipedia.org/wiki/List_of_oxidation_states_of_the_elements",
                              which = 1, rm_nodata_cols = FALSE)
# clean names
names(oxidation) <- str_replace(names(oxidation), "^Positive oxidation states >> ", "")
names(oxidation) <- str_replace(names(oxidation), "^Negative oxidation states >> ", "")

# select columns, types from oxidation states table
oxidation <- oxidation %>%
  select(-Element, -Group, -Notes) %>%
  select(z = Z, symbol = `0`, everything()) %>%
  mutate(`0` = "0") %>%
  gather(-z, -symbol, key = "valence", value = "value") %>%
  filter(!is.na(value)) %>%
  # get rid of en dashes in valence column
  mutate(valence = str_replace(valence, "^[^+]([0-9])", "-\\1")) %>%
  mutate(valence = as.integer(valence)) %>%
  mutate(z = as.integer(z)) %>%
  select(-value) %>%
  group_by(z, symbol) %>%
  summarise(valence = list(sort(valence)))

# select columns, fix types
pt <- elements %>%
  filter(!str_detect(z, "^Notes")) %>%
  select(z, symbol, name = element, group, period, mass = `atomic weight`) %>%
  mutate(z = as.integer(z), group = as.integer(group), period = as.integer(period)) %>%
  # fix references in atomic weight
  mutate(mass = str_replace_all(mass, "[\\[\\(](.*?)[\\]\\)]", "\\1") %>% as.double()) %>%
  left_join(oxidation, by = c("z", "symbol")) %>%
  tibble::as.tibble()

class(pt) <- c("periodic_tibble", class(pt))


# add to package
devtools::use_data(pt, overwrite = TRUE)


