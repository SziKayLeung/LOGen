
# 
read_dea_files <- function(dir){
  f = excel_sheets(dir) %>% purrr::set_names() %>% purrr::map(read_excel, path = dir)
  return(f)
}
