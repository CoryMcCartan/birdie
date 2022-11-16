# follow Census processing rules at <https://www2.census.gov/topics/genealogy/2010surnames/surnames.pdf>
# misspellings list not implemented
# double-barrel names not removed yet
proc_names = function(x) {
    x = stringr::str_to_upper(x)
    x = stringi::stri_trans_general(x, "Latin-ASCII")
    x = if_else(str_starts(x, "[A-Z] [A-Z]$"), NA_character_, x)
    x = if_else(str_starts(x, "([A-Z] ){2,}[A-Z]$"), str_remove_all(x, " "), x)
    x = if_else(x %in% c("JUNIOR", "SENIOR", "THIRD", "CRUTHIRD"), x,
                str_remove(x, "(JUNIOR|SENIOR|THIRD|JR|III| II| IV| J R| S R)$"))
    x = str_remove(x, "^(JR|III|II |J R |S R )")
    x = if_else(stringr::str_length(x) >= 7, str_remove(x, "SR$"), x)
    x = str_remove_all(x, "[.,\\'\"!@#$%^&*/?~`]")
    x = str_replace_all(x, "^(MC|MAC|O) ", "\\1")
    x = str_replace_all(x, "-", " ")
    x = stringr::str_squish(x)
    x = if_else(stringr::str_length(x) == 2 & x != "NG" &
                    !str_detect(x, "[AEIOUY]"), NA_character_, x)
    x = if_else(x %in% c("DECLINE TO STATE", "DONT KNOW", "DECLINE",
                         "NO NAME", "NO NOMBRE", "SAME AS ABOVE"),
                NA_character_, x)
    x
}

is_double_name = function(x) {
    str_detect(x, "^\\S+ \\S+$")
}
