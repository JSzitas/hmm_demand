

intermittent_examples <- tsibbledata::PBS %>%
  as.data.frame() %>%
  dplyr::group_by(Concession, Type, ATC1, ATC1_desc, ATC2, ATC2_desc) %>%
  dplyr::summarise( how_many_zero = sum(Scripts == 0), 
                    series_length = length(Scripts),
                    fraction_zero = how_many_zero / series_length
                    ) %>%
  dplyr::filter(fraction_zero > 0.3 & fraction_zero < 0.6)
