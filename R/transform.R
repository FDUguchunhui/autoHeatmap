#'@import magrittr dplyr tibble tidyr xlsx

data_demo <- 'data/trans_demo.xlsx'

#'@title transform a perth mfi excel file
#'@description transform a perth mfi excel file
#'@param data a excel file
#'@export
perth_table_trans <- function(data){

  # read-in data, with only Subject as char, others as numeric
  dat <- read.xlsx(data, sheetIndex = 1,
                         colClasses = (c('character',
                                          rep('numeric', 20)))) %>%
    as_tibble() %>%
    select(-contains('NA')) %>% # delete NA column which is used as delim
    drop_na() %>% # drop row without information
    # select useful column
    select(-c('X.Pacific.Blue.A.', 'X.Alexa.Fluor.700.A.', 'X.PE.Cy7.A.')) %>%
    select(Subject = `Sample..Set.1.blood.PMN.old`, everything())

  # need to create new var panel first
  dat$panel = regmatches(dat$Subject,
                         regexpr(dat$Subject,
                                 pattern = '(?U)((?<=FL_)[^_]+(?=_.+.fcs))',
                                 perl = T))

  # extract meaningful name to muate original Subject var
  dat$Subject = regmatches(dat$Subject,
                           regexpr(dat$Subject,
                                   pattern = '(?U)(?<=_)[^_]+((?=__WB)|(?=_WB))',
                                   perl = T))

  # transform the tibble into a long format, new columns to indicate FITC or APC
  dat.long <- dat %>% pivot_longer(cols = c('X.FITC.A.', 'X.APC.A.'),
                                   names_to = 'FITC or APC',
                                   values_to = 'intensity',
  )

  # re-transform the tibble into a wide format, use `FITC or APC`and panel
  #togather to create multiple columns
  dat.wide <-  dat.long %>% pivot_wider(names_from = c(`FITC or APC`, panel),
                                        values_from = intensity
  )

  # summarize(average) all the columns grouped by Subject
  dat.tranformed <- dat.wide %>%
    group_by(Subject) %>%
    summarise_all(mean, na.rm = T) %>%
    mutate(
           HNE = matches('X.FITC.A._9PR', ignore.case = T),
           CD91 = matches('X.FITC.A._10Tcel', ignore.case = T),
           `PD-1` = X.FITC.A._11PD1,
           CD47 = rowSums(.[grep('X.FITC.A._12Ap', ignore.case = T),], na.rm = T),
           EXO = X.FITC.A._13EV,
           ###################
           CD163 = X.APC.A._9PR,
           CDD172a = X.APC.A._10Tcel,
           CD36 = X.APC.A._11PD1,
           FLICA = X.APC.A._12Ap,
           `TIMP-1` = X.APC.A._13EV,
           )
  # %>%
  #   select(Subject, CD45,	CD63,	CD16,	CD66b,	CD33,	HNE,	CD91,
  #          `PD-1`,	CD47,	EXO,	CD163,	CDD172a,	CD36,	FLICA,	`TIMP-1`,
  #          everything()
  #          )

  # for FITC
  # 9PR HNE
  # 10Tcel CD91
  # 11PD1  PD-1
  # 12Ap CD47
  # 13EV EXO

  ## For APC
  # 9PR CD163
  # 10Tcel CDD172a
  # 11PD1  CD36
  # 12Ap FLICA
  # 13EV TIMP-1



  return(dat.tranformed)
}



