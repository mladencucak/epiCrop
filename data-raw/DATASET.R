## code to prepare `weather` dataset goes here
weather <-
  read.csv(here::here("data-raw", "weather_OP_2016.csv"),
           row.names= NULL, stringsAsFactors = FALSE)

usethis::use_data("weather", overwrite = TRUE)
