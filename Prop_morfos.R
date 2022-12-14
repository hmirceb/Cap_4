a<-data.frame(level = c('Order', 'Family', 'Genus', 'Species'),
           N = c(PERDIVER %>%
                   filter(genus == order) %>%
                   pull(Species) %>%
                   unique(.) %>%
                   length(),
                 PERDIVER %>%
                   filter(genus == family) %>%
                   pull(Species) %>%
                   unique(.) %>%
                   length(),
                 PERDIVER %>%
                   filter(genus != family & genus != order) %>%
                   filter(grepl('_sp\\.', Species)) %>%
                   pull(Species) %>%
                   unique(.) %>%
                   length(),
                 PERDIVER %>%
                   filter(genus != family & genus != order) %>%
                   filter(!grepl('_sp\\.', Species)) %>%
                   pull(Species) %>%
                   unique(.) %>%
                   length()))
a$prop <- round(100*a$N/sum(a$N), 1)
a

PERDIVER %>%
  filter(genus == order) %>%
  dplyr::select(c(Species, order)) %>%
  distinct() %>%
  group_by(order) %>%
  summarize(n = n())

PERDIVER %>%
  filter(genus == family) %>%
  dplyr::select(c(Species, order)) %>%
  distinct() %>%
  group_by(order) %>%
  summarize(n = n())

PERDIVER %>%
  filter(genus != family & genus != order) %>%
  filter(grepl('_sp\\.', Species)) %>%
  dplyr::select(c(Species, order)) %>%
  distinct() %>%
  group_by(order) %>%
  summarize(n = n())

PERDIVER %>%
  filter(genus != family & genus != order) %>%
  filter(!grepl('_sp\\.', Species)) %>%
  dplyr::select(c(Species, order)) %>%
  distinct() %>%
  group_by(order) %>%
  summarize(n = n())
