library(stringr)
context('String length')

test_that('Str_length is number of characters',
          expect_equal(str_length('a'), 1)
          )
