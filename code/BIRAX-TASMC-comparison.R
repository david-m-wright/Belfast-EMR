# Inferential comparison of Belfast and Tel Aviv analysis cohorts
# Using summary statistics to generate the test statistics


# Function to calculate the two sample t-test from summary data
# From https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha/30450#30450
t.test.from.summary.data <- function(mean1, sd1, n1, mean2, sd2, n2, ...) {
  data1 <- scale(1:n1)*sd1 + mean1
  data2 <- scale(1:n2)*sd2 + mean2
  t.test(data1, data2, ...)
}

