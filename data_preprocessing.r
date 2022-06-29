### spatial transcriptomics data preprocessing
### Generate the data for matrix factorization

library(Seurat)
source('utils.r')

data_preprocessing = function(file_name, top_hvg = 2000, file_sep = '\t', loc_sep = 'x', gene_filtering = 0.1){
    # load count matrix
    count_mat = read.table(file_name, header = TRUE, sep = file_sep, row.names = 1)
    print(paste0('Raw count mat: ', nrow(count_mat), ' spots, ', ncol(count_mat), ' genes'))

    # generate locations
    location_str = rownames(count_mat)
    loc_data = convert_str_2_loc(location_str, loc_sep)
    # loc_data: dataframe with columns 'loc', 'raw_name', 'x_pos', 'y_pos'

    # remove locations with no expression (zero count in total)
    loc_sum = rowSums(count_mat)
    remove_locs = which(loc_sum == 0)
    if(length(remove_locs) > 0){
        count_mat = count_mat[-remove_locs,]
        print(paste0(length(remove_locs), ' spots removed due to zero expression'))
        print('The following locations are removed')
        print(location_str[remove_locs])
    }

    # remove genes expressed in less than gene_filtering locations
    if(gene_filtering > 0){
        print(paste0('Removing genes expressed in less than ', gene_filtering * 100, '% locations'))
        gene_loc_exp = apply(count_mat, 2, function(x) length(which(x > 0))) # number of locations that each gene expressed
        remove_genes = which(gene_loc_exp < gene_filtering * nrow(count_mat))
        if(length(remove_genes) > 0){
            count_mat = count_mat[, -remove_genes]
            removed_gene_names = colnames(count_mat)[remove_genes]
        }
        print(paste0(length(remove_genes), ' genes removed'))
        print(paste0(ncol(count_mat), ' genes remained'))
    }

    # normalization and HVG selection
    print('Normalizing filtered count matrix')
    dataobj <- CreateSeuratObject(counts = t(count_mat))
    dataobj <- NormalizeData(dataobj, normalization.method = "RC", scale.factor = 10000) # default setting
    dataobj <- FindVariableFeatures(dataobj, selection.method = 'vst', nfeatures = top_hvg)
    hvgs <- VariableFeatures(dataobj)
    transformed_mat = GetAssayData(object = dataobj, slot = 'data')  # transformed data of HVGs
    transformed_mat = t(as.matrix(transformed_mat[hvgs,]))  # loc by gene
    
    # calculate distance matrix of locations
    locations = loc_data[rownames(transformed_mat),]
    dist_mat = cal_distance_mat(locations, '2D')

    # result data
    data = list(exp_mat = transformed_mat, 
                dist_mat = dist_mat, 
                loc_info = locations,
                gene_list = colnames(transformed_mat), 
                loc_list = rownames(transformed_mat), 
                top_hvg = top_hvg, 
                gene_filtering = gene_filtering)

    return(data)
}


