returnCritValKDE <- function(x, alpha) {
  kdeObj <- kde1d::kde1d(x)
  kde1d::qkde1d(1 - alpha, kdeObj)
}
