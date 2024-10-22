
##----------------------------General utility functions---------------------------
neighbor_coords <- function(j, i, l){
  if (2^(j-1)+1 <= l){
    return('not enough elements to build parent sets')
  }
  candidates = seq(max(floor((i-2*l)/2)*2,0), min(ceiling((i+2*l)/2)*2,2^j), 2)
  dists = abs(candidates-i)
  id_pa = candidates[sort(dists,index.return=TRUE)$ix[1:(l+1)]]
  x_pa = id_pa / 2^j
  return(x_pa)
}

neighbor_coords_x <- function(j, x, l){
  if (2^(j-1)+1 <= l){
    return('not enough elements to build parent sets')
  }
  i = round(x*2^j)
  candidates = seq(max(floor((i-2*l)/2)*2,0), min(ceiling((i+2*l)/2)*2,2^j), 2)
  dists = abs(candidates/2^j-x)
  id_pa = candidates[sort(dists,index.return=TRUE)$ix[1:(l+1)]]
  x_pa = id_pa / 2^j
  return(x_pa)
}

#' Convert a R DAG to Rcpp DAG
#'
#' The indices in R start from 1, while indices in c++ start from 0. This function
#' simply minus 1 for all the parent sets indices in a dag.
#'
#' @param dag_ord A list, where each element of the list is a numeric vector containing indices of parents
#' @return A list similar to input, except that all indices of parents get minus 1.
#' @export
Rdag_to_Cppdag <- function(dag_ord){
  n = length(dag_ord)
  Cppdag = vector('list',length=n)
  for (i in 1:n){
    if (length(dag_ord[[i]])<1){
      Cppdag[[i]] = vector('numeric',0)
    } else{
      Cppdag[[i]] = dag_ord[[i]] - 1
    }
  }
  return(Cppdag)
}

##-----------------------------d=1--------------------------------
#' Norming DAG on grid, 1-d
#'
#' This function constructs norming DAG on a unidimesional grid of cardinality 2^J+1. For grid of other cardinality, please call
#' the function DAGgrid
#'
#' @param J An integer, number of layers. J determins the cardinality of the grid.
#' @param alpha A numeric value, the smoothness of the Matern process, determining the cardinality of parent sets.
#' @param sort_in_ord A numeric vector representing the ordering of elements in Xcoords. One can simply input the "sort_in_ord" variable from the output of the function DAGgrid_per_hd
#' @return A list of three entries:
#'        X_ord: coordinates of the grid in DAG ordering;
#'        dag_ord: a list of numeric vectors, each entry contains the indices of its parent sets in DAG ordering;
#'        sort_in_ord: a numeric vector, the ordering of DAG with respect to the coordinate ordering.
#' @export
DAGgrid_per_1d <- function(J, alpha){
  n = 2^J+1
  nprime = 2^J
  l = ceiling(alpha)-1   ## order of polynomials
  m = choose(l+1, l)     ## cardinality of parent set
  X = seq(0,2^J,1)/2^J
  dag = vector('list',n)
  ord_2_id = vector('numeric',n)
  ord = 1
  for (x in 0:1){
    id = round(x*nprime)+1
    ord_2_id[[ord]] = id
    ord = ord + 1
    dag[[id]] = vector('numeric',0)
  }
  for (j in 1:J){
    if (2^(j-1)+1 <= l){
      build = FALSE
    } else{
      build = TRUE
    }
    for (i in 0:(2^j)){
      if (i %% 2 == 0){
        next
      }
      x = i/2^j
      id = round(x*nprime)+1
      ord_2_id[[ord]] = id
      ord = ord + 1
      if (build == TRUE){
        dag[[id]] = vector('numeric',m)
        x_coords = neighbor_coords(j,i,l)
        i_nei = 1
        for (ix in 0:l){
          dag[[id]][i_nei] = round(x_coords[ix+1]*nprime)+1
          i_nei = i_nei + 1
        }
      } else{
        dag[[id]] = vector('numeric',2^(j-1)+1)
        i_nei = 1
        for (ix in 0:(2^(j-1))){
          dag[[id]][i_nei] = round(ix/2^(j-1)*nprime)+1
          i_nei = i_nei + 1
        }
      }
    }
  }
  id_2_ord = sort(ord_2_id,index.return=TRUE)$ix
  X_ord = X[ord_2_id]
  dag_ord = vector('list',n)
  for (i in 1:n){
    pa = dag[[i]]
    ord = id_2_ord[i]
    dag_ord[[ord]] = pa
    for (ipa in 1:length(pa)){
      dag_ord[[ord]][ipa] = id_2_ord[pa[ipa]]
    }
  }
  for (i in 1:2){
    dag_ord[[i]] = vector('numeric',0)
  }
  return(list('X_ord'=X_ord, 'dag_ord'=dag_ord, 'sort_in_ord'=ord_2_id))
}

#' Norming DAG for test data, 1-d grid algorithm
#'
#' This function constructs norming DAG on a test set using 1-d grid algorithm. Only
#' the training set is required to be on a grid, while test set can be arbitrary locations.
#'
#' @param Xcoords List of numeric vectors, each element of the list contains the coordinate values of training data
#' @param Xtest A numeric matrix, each row corresponds to the coordinates of a test location
#' @param alpha A numeric value, the smoothness of the Matern process, determining the cardinality of parent sets.
#' @return A list of numeric vectors representing the test DAG. Each element of the list contains the indices of parents in DAG ordering.
#' @export
DAGgrid_test_1d <- function(X, Xtest, alpha){
  l = ceiling(alpha)-1   ## degree of polynomials
  m = choose(l+1, l)     ## cardinality of parent set
  dag_test_mat = FNN::get.knnx(X, Xtest, m)$nn.index
  dag_test = vector('list',length(Xtest))
  for (i in 1:length(Xtest)){
    dag_test[[i]] = dag_test_mat[i,]
  }
  return(dag_test)
}




##------------------------------d=2--------------------------------
## convert x values to id, only used in perfect DAG cases.
x_coords_to_id <- function(J, x1, x2){
  nprime = 2^J
  return((nprime+1) * round(x1*nprime) + round(x2*nprime) + 1)
}

## convert ids in each coordinate to the overall id
coord_ids_to_id <- function(n1, n2, id1, id2){
  return( n2*(id1-1)+id2 )
}

## generate j1 * j2 grid
gridX <- function(n1, n2=NULL){
  ## create coordinates for a two dimensional (2^J+1)^2 grid; require J\ge 1
  if (is.null(n2)){
    n2 = n1
  }
  sep = 1/(max(n1,n2)-1)
  n = n1*n2
  X = matrix(0,nrow=n,ncol=2)
  id = 1
  for (i1 in 0:(n1-1)){
    for (i2 in 0:(n2-1)){
      X[id,] = c(i1*sep, i2*sep)
      id = id + 1
    }
  }
  return(X)
}

#' Norming DAG on grid, 2-d
#'
#' This function constructs norming DAG on a 2-dimesional grid of cardinality (2^J+1)^2. For grid of other cardinality, please call
#' the function DAGgrid
#'
#' @param J An integer, number of layers. J determins the cardinality of the grid.
#' @param alpha A numeric value, the smoothness of the Matern process, determining the cardinality of parent sets.
#' @param sort_in_ord A numeric vector representing the ordering of elements in Xcoords. One can simply input the "sort_in_ord" variable from the output of the function DAGgrid_per_hd
#' @return A list of three entries:
#'        X_ord: coordinates of the grid in DAG ordering;
#'        dag_ord: a list of numeric vectors, each entry contains the indices of its parent sets in DAG ordering;
#'        sort_in_ord: a numeric vector, the ordering of DAG with respect to the coordinate ordering.
#' @export
DAGgrid_per_hd <- function(J, alpha){
  n = (2^J+1)^2
  l = ceiling(alpha)-1   ## order of polynomials
  m = choose(l+2, l)     ## cardinality of parent set
  X = gridX(2^J+1)
  dag = vector('list',n)
  ord_2_id = vector('numeric',n)
  ord = 1
  for (x1 in 0:1){
    for (x2 in 0:1){
      id = x_coords_to_id(J,x1,x2)
      ord_2_id[[ord]] = id
      ord = ord + 1
      dag[[id]] = vector('numeric',0)
    }
  }
  for (j in 1:J){
    if (2^(j-1)+1 <= l){
      build = FALSE
    } else{
      build = TRUE
    }
    for (i1 in 0:(2^j)){
      for (i2 in 0:(2^j)){
        if ((i1 %% 2 == 0) & (i2 %% 2 == 0)){
          next
        }
        x = c(i1/2^j, i2/2^j)
        id = x_coords_to_id(J,x[1],x[2])
        ord_2_id[[ord]] = id
        ord = ord + 1
        if (build == TRUE){
          dag[[id]] = vector('numeric',m)
          x1_coords = neighbor_coords(j,i1,l)
          x2_coords = neighbor_coords(j,i2,l)
          i_nei = 1
          for (ix1 in 0:l){
            for (ix2 in 0:(l-ix1)){
              dag[[id]][i_nei] = x_coords_to_id(J,x1_coords[ix1+1],x2_coords[ix2+1])
              i_nei = i_nei + 1
            }
          }
        } else{
          dag[[id]] = vector('numeric',(2^(j-1)+1)^2)
          i_nei = 1
          for (ix1 in 0:(2^(j-1))){
            for (ix2 in 0:(2^(j-1))){
              dag[[id]][i_nei] = x_coords_to_id(J,ix1/2^(j-1),ix2/2^(j-1))
              i_nei = i_nei + 1
            }
          }
        }
      }
    }
  }
  id_2_ord = sort(ord_2_id,index.return=TRUE)$ix
  X_ord = X[ord_2_id,]
  dag_ord = vector('list',n)
  for (i in 1:n){
    pa = dag[[i]]
    ord = id_2_ord[i]
    dag_ord[[ord]] = pa
    for (ipa in 1:length(pa)){
      dag_ord[[ord]][ipa] = id_2_ord[pa[ipa]]
    }
  }
  for (i in 1:4){
    dag_ord[[i]] = vector('numeric',0)
  }
  return(list('X_ord'=X_ord, 'dag_ord'=dag_ord, 'sort_in_order'=ord_2_id))
}


#' @export
DAGgrid_Not_Norm <- function(J){
  n = (2^J+1)^2
  l = 2                  ## order of polynomials
  m = choose(l+2, l)     ## cardinality of parent set
  X = gridX(2^J+1)
  dag = vector('list',n)
  ord_2_id = vector('numeric',n)
  ord = 1
  for (x1 in 0:1){
    for (x2 in 0:1){
      id = x_coords_to_id(J,x1,x2)
      ord_2_id[[ord]] = id
      ord = ord + 1
      dag[[id]] = vector('numeric',0)
    }
  }
  for (j in 1:J){
    if (2^(j-1)+1 <= l){
      build = FALSE
    } else{
      build = TRUE
    }
    for (i1 in 0:(2^j)){
      for (i2 in 0:(2^j)){
        if ((i1 %% 2 == 0) & (i2 %% 2 == 0)){
          next
        }
        x = c(i1/2^j, i2/2^j)
        id = x_coords_to_id(J,x[1],x[2])
        ord_2_id[[ord]] = id
        ord = ord + 1
        if (build == TRUE){
          dag[[id]] = vector('numeric',m)
          x1_coords = neighbor_coords(j,i1,l)
          x2_coords = neighbor_coords(j,i2,l)
          i_nei = 1
          for (ix1 in 0:l){
            for (ix2 in 0:(l-1)){
              dag[[id]][i_nei] = x_coords_to_id(J,x1_coords[ix1+1],x2_coords[ix2+1])
              i_nei = i_nei + 1
            }
          }
        } else{
          dag[[id]] = vector('numeric',(2^(j-1)+1)^2)
          i_nei = 1
          for (ix1 in 0:(2^(j-1))){
            for (ix2 in 0:(2^(j-1))){
              dag[[id]][i_nei] = x_coords_to_id(J,ix1/2^(j-1),ix2/2^(j-1))
              i_nei = i_nei + 1
            }
          }
        }
      }
    }
  }
  id_2_ord = sort(ord_2_id,index.return=TRUE)$ix
  X_ord = X[ord_2_id,]
  dag_ord = vector('list',n)
  for (i in 1:n){
    pa = dag[[i]]
    ord = id_2_ord[i]
    dag_ord[[ord]] = pa
    for (ipa in 1:length(pa)){
      dag_ord[[ord]][ipa] = id_2_ord[pa[ipa]]
    }
  }
  for (i in 1:4){
    dag_ord[[i]] = vector('numeric',0)
  }
  return(list('X_ord'=X_ord, 'dag_ord'=dag_ord))
}

#' Norming DAG for test data, 2-d grid algorithm
#'
#' This function constructs norming DAG on a test set using 2-d grid algorithm. Only
#' the training set is required to be on a grid, while test set can be arbitrary locations.
#'
#' @param Xcoords List of numeric vectors, each element of the list contains the coordinate values of training data
#' @param Xtest A numeric matrix, each row corresponds to the coordinates of a test location
#' @param alpha A numeric value, the smoothness of the Matern process, determining the cardinality of parent sets.
#' @param sort_in_ord A numeric vector representing the ordering of elements in Xcoords. One can simply input the "sort_in_ord" variable from the output of the function DAGgrid_per_hd
#' @return A list of numeric vectors representing the test DAG. Each element of the list contains the indices of parents in DAG ordering.
#' @export
DAGgrid_test_hd <- function(Xcoords, Xtest, alpha, sort_in_ord){
  sort_in_id = sort(sort_in_ord,index.return=TRUE)$ix
  ntrain_coord = length(Xcoords[[1]])
  ntest = nrow(Xtest)
  d = ncol(Xtest)
  l = ceiling(alpha)-1   ## order of polynomials
  m = choose(l+d, l)     ## cardinality of parent set
  coords_nei = vector('list',d)
  for (j in 1:d){
    coords_nei[[j]] = FNN::get.knnx(Xcoords[[j]], Xtest[,j], l+1)$nn.index
  }
  dag_test = vector('list',ntest)
  for (i in 1:ntest){
    nei_lst_1 = coords_nei[[1]][i,]
    nei_lst_2 = coords_nei[[2]][i,]
    loc = 1
    dag_test[[i]] = rep(0,m)
    for (i1 in 0:l){
      for (i2 in 0:(l-i1)){
        now_id = (nei_lst_1[i1+1]-1)*ntrain_coord + nei_lst_2[i2+1]
        dag_test[[i]][loc] = sort_in_id[now_id]
        loc = loc + 1
      }
    }
  }
  return(dag_test)
}



##-------------------------------output function-------------------------------

#' Norming DAG on grid
#'
#' This function computes DAG for grid data. You do not need to provide coordinates
#' of grid, only the observed response variable on the grid. The size of the grid is
#' automatically extracted from the shape of the response Y.
#'
#' @param Y A numeric vector or matrix, the observed response variable on the grid.
#' @param minsep A numeric value, Minimal separation distance of the grid.
#' @param alpha A numeric value, the smoothness of the Matern process, determining the cardinality of parent sets.
#' @return A list of four entries:
#'        X_ord: coordinates of the grid in DAG ordering;
#'        dag_ord: a list of numeric vectors, each entry contains the indices of its parent sets in DAG ordering;
#'        Y_ord: a numeric vector, the Y values in DAG ordering
#'        sort_in_ord: a numeric vector, the ordering of DAG with respect to the coordinate ordering.
#' @export
DAGgrid <- function(Y, minsep, alpha=3/2){
  if (is.vector(Y)){
    d = 1
    n = length(Y)
  } else {
    dds = dim(Y)
    if (length(dds)==2 & min(dds)==1){
      Y = as.vector(Y)
      d = 1
      n = length(Y)
    } else if (length(dds)==2 & min(dds)>1){
      d = 2
      ns = c(nrow(Y),ncol(Y))
    }
  }

  if (d==1){
    scaling = minsep * (n-1)
    X = seq(0,1,1/(n-1))
    J = floor(log(n-1,base=2))
    n_per = 2^J + 1
    n_res = n - n_per
    dag_per_obj = DAGgrid_per_1d(J, alpha)
    if (n_res>0){
      gap = floor((n+1)/(n_res+1))-1
      id_res = seq(1,n_res)*(gap+1)
      X_res = X[id_res]
      X_per = X[-id_res]
      dag_res = DAGgrid_test_1d(X_per,X_res,alpha)
      X_ord = c(X_per[dag_per_obj$sort_in_ord], X_res)
      dag_ord = c(dag_per_obj$dag_ord, dag_res)
      sort_in_ord = c(seq(1,n,1)[-id_res][dag_per_obj$sort_in_ord], id_res)
    } else{
      X_ord = dag_per_obj$X_ord
      dag_ord = dag_per_obj$dag_ord
      sort_in_ord = dag_per_obj$sort_in_ord
    }
    # sort_in_id = sort(X_ord,index.return=TRUE)$ix
    # sort_in_ord = sort(sort_in_id,index.return=TRUE)$ix
    Y_ord = Y[sort_in_ord]
  } else if (d==2){
    scaling = minsep * (max(ns)-1)
    n = prod(ns)
    X_all_coords = vector('list',d)
    Js = rep(0,d)
    for (j in 1:d){
      Js[j] = floor(log(ns[j]-1,base=2))
      X_all_coords[[j]] = seq(0,ns[j]-1,1)/(max(ns)-1)
    }
    J = min(Js)
    n_per = 2^J + 1
    id_pers = vector('list',d)
    X_per_coords = vector('list',d)
    dag_per_obj = DAGgrid_per_hd(J,alpha)
    if (n_per^d == n){
      X_ord = dag_per_obj$X_ord
      dag_ord = dag_per_obj$dag_ord
      sort_in_ord = dag_per_obj$sort_in_ord
      Y_ord = as.vector(t(Y))[sort_in_ord]
    } else{
      for (j in 1:d){
        if (n_per < ns[j]/2){
          has_res = TRUE
          gap = floor((ns[j]-1)/(n_per-1))
          id_pers[[j]] = c(seq(1,by=gap,length.out=n_per-1),ns[j])
          X_per_coords[[j]] = X_all_coords[[j]][id_pers[[j]]]
        } else{
          gap = floor((ns[j]+1)/(ns[j]-n_per+1))-1
          id_pers[[j]] = seq(1,ns[j])[-seq(1,ns[j]-n_per)*(gap+1)]
          X_per_coords[[j]] = X_all_coords[[j]][id_pers[[j]]]
        }
      }
      label_res = rep(1, n)
      for (i1 in id_pers[[1]]){
        for (i2 in id_pers[[2]]){
          label_res[ns[2]*(i1-1)+i2] = 0
        }
      }
      X_res = gridX(ns[1],ns[2])[which(label_res==1),]
      dag_res = DAGgrid_test_hd(X_per_coords, X_res, alpha, dag_per_obj$sort_in_ord)
      dag_ord = c(dag_per_obj$dag_ord, dag_res)
      sort_in_ord = c(which(label_res==0)[dag_per_obj$sort_in_ord], which(label_res==1))
      X_ord = gridX(ns[1],ns[2])[sort_in_ord,]
      Y_ord = as.vector(t(Y))[sort_in_ord]
    }
  }
  return(list(X_ord=X_ord*scaling, dag_ord=dag_ord, Y_ord=Y_ord, sort_in_ord=sort_in_ord))
}







##---------------------------dump functions--------------------------------
# DAGgrid_test_1d <- function(J, alpha, Xtest, X, dag){
#   ntrain = length(X)
#   ntest = length(Xtest)
#   X_input = data.frame(X=X,input_ord=seq(1,ntrain))
#   X_input$coord_ord = sort(X,index.return=TRUE)$
#
#   ## calling X_input[coords_trans] will provide X sorted in coords
#   coords_trans = sort(X_input$coord_ord,index.return=TRUE)$ix
#   X_coord = X_input[coords_trans]
#
#   ## totalsort maps the input order of a test location to the input order of the closest training location to its left
#   totalsort = sort(c(X,Xtest),index.return=TRUE)$ix
#   Xtest_id_left = rep(0,ntest)
#   for (i in 1:(ntrain+ntest)){
#     if (totalsort[i]>ntrain){
#       for (j in seq(i-1,1,-1)){
#         if (totalsort[j]<=ntrain){
#           Xtest_id_left[totalsort[i]-ntrain] = totalsort[j]
#         }
#       }
#     }
#   }
#
#   l = ceiling(alpha)-1   ## degree of polynomials
#   m = choose(l+1, l)     ## cardinality of parent set
#   dag_all = vector('list',ntrain+ntest)
#   dag_all[1:ntrain] = dag
#   for (i in 1:ntest){
#     dag[[i+ntrain]] = vector('numeric',l+1)
#     x = Xtest[i]
#     x_left_input_ord = Xtest_id_left[i]
#     x_left_coord_ord = X_input$coord_ord[i]
#     nei_candidates = vector('numeric',0)
#     if (x_left_coord_ord>0){
#       left_coords = seq(max(1,x_left_coord_ord-m), x_left_coord_ord, 1)
#     } else{
#       left_coords = vector('numeric',0)
#     }
#     if (x_left_coord_ord<ntrain){
#       right_coords = seq(x_left_coord_ord+1, min(n,x_left_coord_ord+1+m), 1)
#     } else{
#       right_coords = vector('numeric',0)
#     }
#     both_coords = c(left_inputs, right_inputs)
#
#     ## convert coord order into input order
#     nei_candidates = X_coord$input_ord[both_coords]
#
#     ## continue writing
#     # dag[[i+ntrain]] =
#     #
#     # i_nei = 1
#     # for (ix in 0:l){
#     #   dag[[i+ntrain]][i_nei] = ord_2_id[round(x_coords[ix+1]*nprime)+1]
#     #   i_nei = i_nei + 1
#     # }
#
#   }
#   return(list('X'=c(X,Xtest), 'dag'=dag))
# }











