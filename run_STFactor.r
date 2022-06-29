# Run STFactor on spatially resolved transcriptomics data

# Format of count matrix file: 
#       - a matrix of spots (rows) * genes (columns)
#       - rownames: coordinates of the spots, x and y position are connected with a separator (e.g., '10x10')
#       - colnames: genes

### Parameters:
###     - file_name: file of the count matrix
###     - num_components: number of underlying factors
###     - output_flag: a string used in the names of the output files
###     - file_sep: separator used to read the count matrix file
###     - top_hvg: number of highly variable genes to select
###     - gene_filtering: used to remove genes expressed in less than gene_filtering locations
###     - loc_sep: the separator between x and y position of the spot names

####### Examples
# source('run_STFactor.r')
# STFactor('data/Layer2_BC_count_matrix-1.tsv', 10, 'bc2')
# STFactor('data/Rep11_MOB_count_matrix-1.tsv', 10, 'mob11')


source('utils.r')
source('data_preprocessing.r')
source('matrix_factorization_spatial.r')
library(ggplot2)

STFactor = function(file_name, num_components, output_flag, file_sep = '\t', top_hvg = 2000, gene_filtering = 0.1, loc_sep = 'x'){
    # generate data from count matrix in the file file_name
    data = data_preprocessing(file_name, top_hvg = top_hvg, file_sep = file_sep, loc_sep = loc_sep, gene_filtering = gene_filtering)
    exp_mat = as.matrix(data$exp_mat)
    dist_mat = data$dist_mat
    num_locs = nrow(exp_mat)
    num_genes = ncol(exp_mat)
    print(paste0('Number of locations: ', num_locs))
    print(paste0('Number of genes: ', num_genes))
    print(paste0('Number of factors: ', num_components))

    # check data
    if(any(is.na(exp_mat))){
    stop('ERROR: the data contains missing data')
    }
    if(any(is.infinite(as.matrix(exp_mat)))){
    stop('ERROR: the data contains infinite data')
    }
    if(!is.numeric(exp_mat)){
    stop('ERROR: the data contains non-numeric data')
    }

    # run factorization
    params=list(N=num_locs, L=num_genes, C=num_components, num_locs = num_locs, 
                top_hvg = top_hvg, gene_filtering = gene_filtering,
                a=1e-6, b=1e6, c=1e-6, d=1e6, e=1e-6, f=1e6, g=0, h=0, u=1e-6, v=1e6, r=1, z=1,
                check_FE = FALSE, decrease_stop = TRUE)
    res <- matrix_factorization(params, exp_mat, dist_mat)
    print("Estimation finished!")
    print(paste(res$maximumiteration,' Iterations were carried out.'))

    # post processing
    post_res = res_processing(res)
    num_group = post_res$num_group
    print(paste0('Number of underlying factors identified: ', num_group))
    group_member_loadings = post_res$group_loadings  # factor by gene
    loc_loadings = post_res$loc_loadings   # loc by factor
    gene_list = data$gene_list
    loc_list = data$loc_list
    factor_list = sapply(c(1:num_group), function(x) paste0('factor_', x))

    # save associated genes of the underlying factors
    factor_genes = c()
    for(g in c(1:num_group)){
        group_loadings = group_member_loadings[g,]
        member = which(group_loadings != 0)
        member_info = data.frame(factor = array(factor_list[g], c(length(member))), gene = gene_list[member], activity = group_loadings[member])
        factor_genes = rbind(factor_genes, member_info)
    }
    write.table(factor_genes, file = paste0('results/', output_flag, '_factor_associated_genes.txt'), sep = '\t', row.names = FALSE)

    # save spatial patterns of the underlying factors
    pattern_loadings = apply(loc_loadings, 2, scale) # scale the activities across the locations of each pattern
    rownames(pattern_loadings) = loc_list
    colnames(pattern_loadings) = factor_list
    write.table(pattern_loadings, file = paste0('results/', output_flag, '_factor_spatial_patterns.txt'), sep = '\t')

    # visulaization of the spatial patterns, output to a pdf file
    spot_locations = data$loc_info
    spot_locations$y_pos = as.numeric(spot_locations$y_pos)
    spot_locations$x_pos = as.numeric(spot_locations$x_pos)
    y_max = max(spot_locations$y_pos) + 1
    y_min = min(spot_locations$y_pos) - 1
    x_max = max(spot_locations$x_pos) + 1
    x_min = min(spot_locations$x_pos) - 1

    print('Plotting spatial patterns')
    pdf(paste0('results/', output_flag, '_factor_spatial_patterns_visualization.pdf'))
    for(g in c(1:num_group)){
        loading_vec = pattern_loadings[,g]
        loading_vec = smooth_extreme_values(loading_vec)
        pattern_df = spot_locations
        pattern_df$Activity = loading_vec
        scatter <- ggplot(pattern_df, aes(x = x_pos, y = y_pos, color = Activity)) +
                    geom_point(size = 6) +
                    coord_fixed() + 
                    xlim(x_min, x_max) + 
                    ylim(y_max, y_min) +
                    scale_color_gradient2(low = '#0f388a', high = '#a11212') +
                    ggtitle(paste0('Spatial pattern of underlying factor ', g)) +
                    theme_void() +
                    theme(legend.position = "right")
        plot(scatter)
    }
    dev.off()
}


