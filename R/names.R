# follow Census processing rules at <https://www2.census.gov/topics/genealogy/2010surnames/surnames.pdf>
# misspellings list not implemented
# double-barrel names not removed yet
proc_names = function(x) {
    x = stringr::str_to_upper(x)
    x = if_else(stringr::str_starts(x, "[A-Z] [A-Z]$"), NA_character_, x)
    x = if_else(x %in% c("JUNIOR", "SENIOR", "THIRD", "CRUTHIRD"), x,
                stringr::str_remove(x, "(JUNIOR|SENIOR|THIRD|JR|III| II| J R| S R)$"))
    x = stringr::str_remove(x, "^(JR|III|II |J R |S R )")
    x = if_else(stringr::str_length(x) >= 7, stringr::str_remove(x, "SR$"), x)
    x = stringr::str_remove_all(x, "[.,\\'\"!@#$%^&*/?~`]")
    x = stringr::str_replace_all(x, "^(MC|MAC|O) ", "\\1")
    x = stringr::str_replace_all(x, "-", " ")
    x = stringr::str_squish(x)
    x
}

is_double_name = function(x) {
    stringr::str_detect(x, "^\\S+ \\S+$")
}
