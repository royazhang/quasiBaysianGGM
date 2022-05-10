library(Matrix)

p = 100
zeroP2 = as(matrix(0,p^2,p^2), "sparseMatrix")

c1_coef_1 = as(diag(1, p^2), "sparseMatrix")
c1_coef_2 = as(diag(-1, p^2), "sparseMatrix")
c1_coef_3 = as(matrix(0, p^2, 1), "sparseMatrix")
c1_coef = cbind(c1_coef_1, c1_coef_2, c1_coef_3)

c2_coef_1 = as(diag(1, p^2), "sparseMatrix")
c2_coef_2 = as(diag(1, p^2), "sparseMatrix")
c2_coef = cbind(c2_coef_1, c2_coef_2, c1_coef_3)

c3_coef_2 = as(diag(1, p^2), "sparseMatrix")
c3_coef = cbind(zeroP2, c3_coef_2, c1_coef_3)

c4_coef_1 = as(matrix(0, (p^2-p)/2, p^2), "sparseMatrix")
count = 0
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    count = count + 1
    c4_coef_1[count,c((j-1)*p+i,(i-1)*p+j)] = c(1,-1)
  }
}
c4_coef = cbind(c4_coef_1, as(matrix(0, (p^2-p)/2, p^2+1), "sparseMatrix"))

c5_coef_2 = as(matrix(0, p, p^2), "sparseMatrix")
for (i in 1:p) {
  c5_coef_2[i, ((i-1)*p+1):(i*p)] = 1
}
c5_coef = cbind(as(matrix(0, p, p^2), "sparseMatrix"), 
                c5_coef_2, 
                as(matrix(-1, p, 1), "sparseMatrix"))


constraint_mat = rbind(c1_coef, c2_coef, c3_coef, c4_coef, c5_coef)
constraint_num = c(dim(c1_coef)[1], dim(c2_coef)[1], dim(c3_coef)[1], dim(c4_coef)[1], dim(c5_coef)[1])
constraint_dir = rep(c("<=",">=", ">=", "=", "<="),
                     constraint_num)

library(lpSolveAPI)
Sym.LP_fun = function(W) {
  constraint_rhs = c(as.vector(W)+1, as.vector(W)+1,
                     rep(0, sum(constraint_num[3:5])))
  
  lprec = make.lp(nrow=0, ncol=ncol(constraint_mat))
  set.objfn(lprec, obj=c(rep(0,2*p^2),1))
  #Add Constraints
  #Note Direction and RHS is included along with Constraint Value
  for(i in 1:nrow(constraint_mat)){
    add.constraint(lprec,
                   xt= constraint_mat[i,], 
                   type = constraint_dir[i], 
                   rhs = constraint_rhs[i])
    #print(i)
  }
  set.type(lprec, c(1:ncol(constraint_mat)), type = c("real"))
  solve(lprec)
  #get.total.iter(lprec)
  B = matrix(get.variables(lprec)[1:p^2], p, p)
  return (list(Wsym = B-1, loss = get.objective(lprec)))
}
