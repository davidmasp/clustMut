# check output status

list_of_files = readLines("required_files_after_tests.list")

files_exist = purrr::map_lgl(list_of_files,fs::file_exists)

if (any(!files_exist)){
  message(glue::glue("Files required don't exist"))
  for (i in list_of_files[!files_exist]){
    warning(glue::glue("file {i} not present after the tests."))
  }

  stop("Quitting")
} else {
  message("All required files found")
  checksums = purrr::map_chr(list_of_files,tools::md5sum)
  df = data.frame(
    files = list_of_files,
    checksum = checksums
  )

  checksums_file = glue::glue(".chcksum_{lubridate::today()}.txt")
  message(glue::glue("Writting checksums in {checksums_file} file."))
  readr::write_tsv(df,path = checksums_file,col_names = FALSE)
}
