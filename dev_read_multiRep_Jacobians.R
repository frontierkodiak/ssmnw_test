

require('R.matlab')

data.dir <- 'C:/Users/icurtsd/Documents/Project Atlanta/Project Structure and stability of mutualistic networks/data/main_test/testresults_2018Jul11T1636/'  
setwd(data.dir)

file1  <- dir(pattern = 'allJacobian*')[1]

foo <- readMat(file1)
str(foo)

m3d <- foo$repMat
str(m3d)

m1 <- m3d[ , , 1 ]

for (r in 1:dim(m3d)[3]){
  m <- m3d[ , , r ]
  print(r)
  str(m)
  print(sum(m))
}