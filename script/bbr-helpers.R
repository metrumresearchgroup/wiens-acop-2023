library(bbr)

edit_model <- function(.mod){
  get_model_path(.mod) %>% 
    file.edit()
}

view_output_file <- function(.mod, .extension){
  ## grab file with .extension
  tmpdir <- get_output_dir(.mod)
  tmpfile <- list.files(tmpdir, pattern=regex(.extension, ignore_case = TRUE))
  
  ## stop if file doesn't exist
  if(length(tmpfile) == 0){
    stop("File does not exist")
  }
  
  ## open file
  file.edit(file.path(tmpdir,tmpfile))
}

view_lst <- function(.mod){
  view_output_file(.mod, .extension="lst")
}

view_prderr <- function(.mod){
  view_output_file(.mod, .extension="PRDERR")
}

update_run_number <- function(
    .mod,
    .suffixes = c(
      '_its.ext',
      '_saem.ext',
      '_imp.ext',
      '.MSF',
      '.ext',
      '.tab',
      '.chn',
      'par.tab'
    )
){
  runno <- get_model_id(.mod)
  modelfile <- get_model_path(.mod)
  based_on <- .mod$based_on[1]
  message(glue::glue("replacing {based_on} with {runno} in {modelfile}"))
  
  txt <- readLines(modelfile) 
  ## edit text of new model file
  for (.s in .suffixes) {
    txt <- gsub(
      paste0(based_on, .s), 
      paste0(runno, .s), 
      txt, 
      ignore.case = T
    )
  }
  ## write updated model file
  cat(txt, file=modelfile, sep='\n')
  ## return model to make this a pipeable function
  return(.mod)
}

copy_mod <- function(...){
  copy_model_from(...,
                  .inherit_tags = TRUE, 
                  .update_model_file=FALSE) %>% 
    update_run_number()
}


mod_diff <- function(.x,.y){
  model_diff(
    file.path(MODEL_DIR, .x) %>% read_model,
    file.path(MODEL_DIR, .y) %>% read_model
  )
}

#' Arrange tags on a model object
#'
#' Arrange the tags on a model object so that they appear in the same order as
#' in the corresponding tags YAML file.
#'
#' @param .mod The `bbi_{.model_type}_model` object to modify
#' @param glossary Tags glossary list (e.g., value from
#'   `yaml::read_yaml("tags.yaml")`)
arrange_tags <- function(.mod, glossary = yaml::read_yaml("tags.yaml")) {
  x <- .mod[["tags"]]
  y <- unlist(glossary)
  return(replace_all_tags(.mod, .tags = x[order(match(x, y))]))
}

arrange_dir_tags <- function(.dir, .tags){
  walk(pull(run_log(.dir, .recurse=FALSE), "run"), ~ {
    read_model(file.path(.dir, .x)) %>%
      arrange_tags(glossary=.tags)
  })
}

add_aic <- function(.runlog) mutate(.runlog, aic=2*param_count + ofv)
