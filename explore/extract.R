library(readr)

dict = read_csv("nc-dict.csv", skip=15, n_max=771)

cols_pick = c(state_id="Voters_StateVoterID",
              first_name="Voters_FirstName",
              middle_name="Voters_MiddleName",
              last_name="Voters_LastName",
              party="Parties_Description",
              turnout="Voters_VotingPerformanceEvenYearGeneral",
              age="Voters_Age",
              ethn="Ethnic_Description",
              gender="Voters_Gender",
              zip="Residence_Addresses_Zip",
              county="County",
              tract="Residence_Addresses_CensusTract",
              bgroup="Residence_Addresses_CensusBlockGroup",
              block="Residence_Addresses_CensusBlock")

d = read_tsv("nc_ex.tab", lazy=TRUE, num_threads=4, col_select=all_of(cols_pick),
             col_types=cols(.default="c"))
if (nrow(d) > 10e6) stop("Too many voters.")

write_rds(d, "~/docs/nc_voters.rds", compress="gz")
