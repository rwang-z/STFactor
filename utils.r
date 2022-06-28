# functions
library(stringr)

# calculate RMSE of complete tensor
cal_rmse = function(x, y){
    se = (x-y)^2
    rmse = sqrt(mean(se))
    
    return(rmse)
}

res_processing = function(res, filter_size = 1){
    # params = res$params
    A = res$A$mu
    X = res$X$mom1
    num_comp = nrow(X)
    PIP = round(res$X$gamma)   # membership matrix
    X_PIP = X * PIP   # loadings
    
    # removing empty components
    member_count = apply(X_PIP, 1, function(x) length(which(x != 0)))
    removed_components = which(member_count < filter_size)
    if(length(removed_components) == num_comp){
        stop('No factor identified, please try to set the parameter num_components to a larger value')
    }else if(length(removed_components) > 0){
        group_member_loadings = X_PIP[-removed_components,]
        loc_loadings = A[, -removed_components]
    }else{
        group_member_loadings = X_PIP
        loc_loadings = A
    }
    if(is.matrix(group_member_loadings)){
        num_group = nrow(group_member_loadings)
    }else{
        num_group = 1
    }
    return(list(group_loadings = group_member_loadings, loc_loadings = loc_loadings, num_group = num_group))
}

cal_mat_inv = function(mat){
    mat_inv = tryCatch(chol2inv(chol(mat)),
                        error = function(e){
                            # print('Using solve function instead')
                            tryCatch(return(solve(mat)),
                                    error = function(e){
                                        # add random noise to cov_mat
                                        print('Failed in using solve, adding noise to the matrix')
                                        mat_noise = mat + diag(runif(dim(mat)[1], -1e-6, 0))
                                        return(solve(mat_noise))
                                    })
                        })
    return(mat_inv)
}

cal_cov_inv_mat = function(r, dist_mat){
    cov_mat = exp(-0.5 * r * dist_mat)
    cov_inv = cal_mat_inv(cov_mat)
    return(cov_inv)
}

cal_cov_inv_list = function(r, dist_mat){
    num_comp = dim(r)[1]
    cov_list = list()
    for(c in c(1:num_comp)){
        cov_list[[c]] = cal_cov_inv_mat(r[c], dist_mat)
    }
    return(cov_list)
}

cal_precision_log_det = function(prec_mat){
    det_log = tryCatch(determinant(prec_mat)$modulus[1],
                        error = function(e){
                            return(2 * determinant(chol(prec_mat))$modulus[1])
                        })
    return(det_log)
}

get_prec_mat_diag = function(prec_list, dimension, num_comp){
    prec_diag_mat = matrix(0, dimension, num_comp)
    for(c in c(1:num_comp)){
        prec_diag_mat[, c] = diag(prec_list[[c]])
    }
    return(prec_diag_mat)
}

convert_str_2_loc = function(location_str, loc_sep = 'x'){
    locations = c()
    i = 1
    for(loc_str in location_str){
        splitted = strsplit(loc_str, loc_sep)
        x_pos = str_trim(splitted[[1]][1], 'both')
        y_pos = str_trim(splitted[[1]][2], 'both')
        loc_name = paste0('loc_', i)
        new_loc_str = str_replace(loc_str, loc_sep, 'x')
        locations = rbind(locations, c(loc_name, new_loc_str, x_pos, y_pos))
        i = i + 1
    }
    locations = as.data.frame(locations)
    colnames(locations) = c('loc', 'raw_name', 'x_pos', 'y_pos')
    rownames(locations) = locations$raw_name
    return(locations)
}

cal_distance_mat = function(locations, dimension = '2D'){
    # locations: dataframe including columns 'x_pos' and 'y_pos' to calcualte the distances
    num_loc = dim(locations)[1]
    if(dimension == '2D'){
        mat = locations[, c('x_pos', 'y_pos')]
    }else if(dimension == '3D'){
        print('calculating distance for 3D data')
        mat = locations[, c('x_pos', 'y_pos', 'z_pos')]
    }
    dist_mat = as.matrix(dist(mat, method = "euclidean"))
    dist_mat = dist_mat^2
    return(dist_mat)
}

smooth_extreme_values = function(loading_vec){
    sorted_loadings = sort(loading_vec, decreasing = TRUE)
    extreme_list = c(5:1)
    num_loc = length(loading_vec)
    for(extreme_ind in extreme_list){
        if(sorted_loadings[extreme_ind] - sorted_loadings[extreme_ind + 1] > 2){
            loading_vec[which(loading_vec == sorted_loadings[extreme_ind])] = sorted_loadings[extreme_ind + 1] + 0.5
        }
        if(abs(sorted_loadings[num_loc - extreme_ind + 1] - sorted_loadings[num_loc - extreme_ind]) > 2){
            loading_vec[which(loading_vec == sorted_loadings[num_loc - extreme_ind + 1])] = sorted_loadings[num_loc - extreme_ind] - 0.5
        }
    }
    return(loading_vec)
}

