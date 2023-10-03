#created by Tina Del Carpio

#these functions were used to summarize recombination rate data
#that was infered with pyrho (https://github.com/popgenmethods/pyrho)


########################
# calculated a weighted value of r for a chromosome

calc_weighted_r <- function(df) { #df = your rmap from pyrho
        windows <- nrow(df)
        window_len <- vector(length = windows)
        r_vals <- vector(length = windows)
        
        for(ii in 1:(windows)){
                window_len[ii] <- (df[ii,2] - df[ii,1])
        }
        
        for(ii in 1:(windows)){
                r_vals[ii] <- (df[(ii),3])
        }
        
        weighted_r <- weighted.mean(r_vals, window_len)
        
        return(weighted_r)
}

#############################################
#comparing recombination rates between 2 populations in specific windows of the genome


#modified version of function above to work with our comparisons
pyrho_weighted_mean_r <- function(pyrho_chr_data){ 
        #pyrho_chr_data is your dataframe of data from a chromosomes rmap 
        #it needs column names: 1 = start, col 2 = end, col 3 = r 
        windows <- nrow(pyrho_chr_data)
        window_len <- vector(length = windows)
        chr_r <- vector(length = windows)
        
        for(jj in 1:windows){
                window_len[jj] <- (pyrho_chr_data[jj, "end"] - pyrho_chr_data[jj, "start"] + 1 ) #subtract start from end to get length
        }
        
        for(kk in 1:windows){
                chr_r[kk] <- (pyrho_chr_data[(kk), "r"]) #make a vector with all the values of r
        }
        
        chr_weighted_mean_r <- weighted.mean(chr_r, window_len) #calc avg of r weighted by window length
        return(chr_weighted_mean_r)
        
}

#function that partitions your data into the same genomic windows
#and calculates the weighted value of r for each population in those windows

binning <- function(pop1_df, pop2_df, bin_size ){
        #pop1_df and pop2_df are the dataframes of your rmap data for the same
        #chromosome, 1 for each of the populations you will compare 
       
         #find the min and max for the chromosome for both populations
        
        start_val <- max(pop1_df[1,1], pop2_df[1,1])
        
        end_val <- min(pop1_df[nrow(pop1_df),2], pop2_df[nrow(pop2_df),2])
        
        #define the windows
        chrm_start <-(signif(start_val, (floor(log10(start_val)) - 1)) + 100) 
        chrm_end <-(signif(end_val, (floor(log10(end_val)) - 1)) - 100) 
        
        win_size = bin_size #in bp
        
        win_starts <-vector()
        win_ends <- vector()
        nwindows <- round(chrm_end/win_size)
        
        for (ii in 1:nwindows){
                temp_start <- (chrm_start + (win_size * (ii-1)))
                temp_end <- (temp_start + win_size - 1)
                if (temp_end > chrm_end){
                        break
                }
                else {
                        win_starts[ii] <- temp_start
                        win_ends[ii] <- temp_end
                }
        }
        
        #make empty vector to store the avg r of each window
        
        binned_df <- as.data.frame(matrix(nrow = length(win_starts), ncol = 4))
        colnames(binned_df) <- c("win_start", "win_end", "pop1_win_r", "pop2_win_r")
        
        binned_df$win_start <- win_starts
        binned_df$win_end <- win_ends
        
        #subset the lines for window one into new df
        for (ii in 1:length(win_starts)){
                temp_pop1_df <- pop1_df[pop1_df$start < win_ends[ii] & pop1_df$end > win_starts[ii], ]
                temp_pop2_df <- pop2_df[pop2_df$start < win_ends[ii] & pop2_df$end > win_starts[ii], ]    
                
                #change the first start locatin
                
                temp_pop1_df[1, "start"] <- win_starts[ii]
                temp_pop2_df[1, "start"] <- win_starts[ii]
                
                #change the last end location
                temp_pop1_df[nrow(temp_pop1_df), "end"] <- win_ends[ii]
                temp_pop2_df[nrow(temp_pop2_df), "end"] <- win_ends[ii]
                
                # print(temp_pop1_df)
                # print(temp_pop2_df)
                
                #calc weighted r for that window
                #store in the vector of window values
                
                binned_df[ii, "pop1_win_r"] <- pyrho_weighted_mean_r(temp_pop1_df)
                binned_df[ii, "pop2_win_r"] <- pyrho_weighted_mean_r(temp_pop2_df)
                #when last loop return the df of weighted rs by window
                
                if (ii == length(win_starts)) {
                        return(binned_df)
                }
                
        }
}

#function to conver the dataframe made with the "binning" function to calculate ranks
#the ranks are generated within a given population so both populations
#will have windows ranked from 1 to n number of windows

rank_df <- function(df){
        
        rank_matrix <- as.data.frame(matrix(nrow = nrow(df), ncol = ncol(df)))
        colnames(rank_matrix) <- colnames(df)
        rank_matrix$win_start <- df$win_start
        rank_matrix$win_end <- df$win_end
        rank_matrix$pop1_win_r <- rank(df$pop1_win_r)
        rank_matrix$pop2_win_r <- rank(df$pop2_win_r)
        
        return(rank_matrix)
}










